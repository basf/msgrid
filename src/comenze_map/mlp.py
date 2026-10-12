"""MLP - Molecular Lipophilic Potential grid (pure Python/NumPy).

Reimplementation in modern Python of pyMLP by Julien Lefeuvre
(https://github.com/plantigrade/pyMLP, BSD 3-Clause, see LICENSES/pyMLP.txt)
    Broto P et al. Eur J Med Chem 1984, 19, 71-78
    Laguerre M et al. Pharm Sci 1997, 3, 217-222
    Ghose AK et al. J Phys Chem A 1998, 102, 3762-3772
Amended by Dmitry Suplatov (BASF SE): Crippen table, ``Readin`` option.

Each non-hydrogen ATOM of the structure contributes ``100 * fi * f(d)`` to
every grid point, where ``fi`` is the atomic lipophilicity and ``d`` the
distance. Available ``f``:

    dubost      1 / (1 + d)
    fauchere    exp(-d)            (default)
    brasseur    exp(-d / 3.1)
    buckingham  1 / d**3
    type5       exp(-sqrt(d))

``fauchere`` is evaluated within a per-atom cutoff (``DEFAULT_CUTOFF``); the
neglected tail is at most ``100 * 2 * exp(-25) ~ 3e-9`` per atom.
"""

from __future__ import annotations

import argparse
import ast
import logging
import sys
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np

from .cliutil import exit_on_error
from .grid import GridMap, accumulate, require_positive
from .mlp_tables import TABLES
from .structure import Structure

log = logging.getLogger(__name__)

DEFAULT_CUTOFF = 25.0
DEFAULT_SPACING = 1.0
DEFAULT_MARGIN = 10.0
# Buckingham diverges at d = 0; distances are floored to keep the grid finite.
BUCKINGHAM_MIN_DIST = 0.1

METHODS = {
    "dubost": lambda d: 1.0 / (1.0 + d),
    "fauchere": lambda d: np.exp(-d),
    "brasseur": lambda d: np.exp(-d / 3.1),
    "buckingham": lambda d: np.maximum(d, BUCKINGHAM_MIN_DIST) ** -3.0,
    "type5": lambda d: np.exp(-np.sqrt(d)),
}

# Residue-name variants (force-field specific) mapped onto the table names.
RESNAME_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS",
    "HSP": "HIS", "HISD": "HIS", "HISE": "HIS", "HISH": "HIS",
    "ASH": "ASP", "ASPH": "ASP", "GLH": "GLU", "GLUH": "GLU",
    "CYX": "CYS", "CYM": "CYS", "LYN": "LYS", "ARN": "ARG",
}  # fmt: skip


FiTable = Dict[str, Dict[str, float]]


def load_logp_file(path: Union[str, Path]) -> FiTable:
    """Read a ``Crippen = {...}`` dictionary from a Python file.

    The file is parsed, not executed: only the literal assigned to ``Crippen``
    is evaluated.
    """
    tree = ast.parse(Path(path).read_text(), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "Crippen" for t in node.targets
        ):
            return ast.literal_eval(node.value)
    raise ValueError(f"No 'Crippen' dictionary found in {path}")


def get_table(name: str = "Crippen", logp_file: Optional[Union[str, Path]] = None) -> FiTable:
    if name == "Readin":
        if not logp_file:
            raise ValueError("Table 'Readin' requires a logP file")
        return load_logp_file(logp_file)
    try:
        return TABLES[name]
    except KeyError:
        raise ValueError(
            f"Unknown reference table {name!r}; choose from {[*TABLES, 'Readin']}"
        ) from None


def _lookup_names(resname: str):
    yield resname
    if resname in RESNAME_ALIASES:
        yield RESNAME_ALIASES[resname]
    if len(resname) == 4 and resname[0] in "NC":  # terminal-residue prefix, e.g. NALA
        yield resname[1:]
        if resname[1:] in RESNAME_ALIASES:
            yield RESNAME_ALIASES[resname[1:]]


def assign_fi(structure: Structure, table: FiTable):
    """Return (xyz, fi) for the non-hydrogen ATOM records of ``structure``."""
    sel = np.flatnonzero((structure.record == "ATOM") & ~structure.is_hydrogen)
    if sel.size == 0:
        raise ValueError(f"No non-hydrogen ATOM records in {structure.source or 'structure'}")
    cache: Dict[tuple, float] = {}
    missing = []
    fi = np.empty(sel.size)
    for n, i in enumerate(sel):
        key = (str(structure.resname[i]), str(structure.name[i]))
        if key not in cache:
            for res in _lookup_names(key[0]):
                if res in table and key[1] in table[res]:
                    cache[key] = float(table[res][key[1]])
                    break
            else:
                cache[key] = np.nan
                missing.append(key)
        fi[n] = cache[key]
    if missing:
        shown = ", ".join(f"{r}:{a}" for r, a in missing[:10])
        raise KeyError(f"No fi value for {len(missing)} residue/atom pair(s): {shown}")
    return structure.xyz[sel], fi


def grid_geometry(xyz: np.ndarray, spacing: float, margin: float = DEFAULT_MARGIN):
    """Grid origin and shape: bounding box + margin, padded to whole spacings."""
    hi = xyz.max(axis=0) + margin
    lo = xyz.min(axis=0) - margin
    adj = ((spacing - (hi - lo)) % spacing) / 2
    hi, lo = hi + adj, lo - adj
    n = np.rint((hi - lo) / spacing).astype(int) + 1
    return lo, tuple(int(v) for v in n)


def compute_mlp(
    structure: Structure,
    table: str = "Crippen",
    method: str = "fauchere",
    spacing: float = DEFAULT_SPACING,
    logp_file: Optional[Union[str, Path]] = None,
    cutoff: Optional[float] = DEFAULT_CUTOFF,
    n_jobs: Optional[int] = 1,
) -> GridMap:
    """Compute the MLP grid of ``structure``.

    ``cutoff`` (Angstrom) only applies to ``fauchere``; other methods are
    evaluated over the whole grid.
    """
    if method not in METHODS:
        raise ValueError(f"Unknown method {method!r}; choose from {list(METHODS)}")
    require_positive(spacing=spacing, cutoff=cutoff)
    xyz, fi = assign_fi(structure, get_table(table, logp_file))
    origin, shape = grid_geometry(xyz, spacing)
    log.info(
        "MLP: %d atoms, table %s, method %s -> grid %dx%dx%d (%s points, spacing %s A)",
        len(fi), table, method, *shape, f"{np.prod(shape):,}", spacing,
    )
    data = accumulate(
        xyz,
        100.0 * fi,
        METHODS[method],
        origin,
        spacing,
        shape,
        cutoff=cutoff if method == "fauchere" else None,
        dtype=np.float64,
        n_jobs=n_jobs,
    )
    return GridMap(data, origin, float(spacing))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="MLP",
        description="Molecular Lipophilic Potential evaluator",
        epilog="Reimplementation of pyMLP by Julien Lefeuvre (https://github.com/plantigrade/pyMLP)",
    )
    p.add_argument("-f", metavar="input file", required=True,
                   help="PDB (or PQR) coordinate file, e.g. input.pdb")
    p.add_argument("-o", metavar="output file", required=True,
                   help="output file in dx format, e.g. output.dx")
    p.add_argument("-t", metavar="reference table", default="Crippen",
                   choices=[*TABLES, "Readin"],
                   help="Laguerre; Ghose; Crippen (default); Readin")
    p.add_argument("-l", metavar="precalculated logp data", default="",
                   help="logP data file for use with -t Readin")
    p.add_argument("-m", metavar="calculation method", default="fauchere", choices=list(METHODS),
                   help="dubost [1/(1+d)]; fauchere (default) [exp(-d)]; brasseur [exp(-d/3.1)]; "
                        "buckingham [1/d**3]; type5 [exp(-sqrt(d))]")
    p.add_argument("-s", metavar="grid spacing", type=float, default=DEFAULT_SPACING,
                   help="grid spacing in Angstrom (default: 1.0)")
    p.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF,
                   help=f"per-atom cutoff in Angstrom for fauchere (default: {DEFAULT_CUTOFF})")
    p.add_argument("--threads", type=int, default=1,
                   help="number of threads (0 = all available; default: 1)")
    return p


@exit_on_error
def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    structure = Structure.from_file(args.f)
    result = compute_mlp(
        structure, table=args.t, method=args.m, spacing=args.s,
        logp_file=args.l or None, cutoff=args.cutoff, n_jobs=args.threads,
    )
    result.write_dx(args.o)
    log.info("MLP written to %s", args.o)
    return 0


if __name__ == "__main__":
    sys.exit(main())
