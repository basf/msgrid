"""MEP - Molecular Electrostatic Potential grid from a PQR file (pure Python/NumPy).

Screened Coulomb (Yukawa / Debye-Hueckel) potential of the heavy-atom group
charges (heavy atom + bonded hydrogens within 1.3 A, same residue):

    V(r_i) = scale * sum_k q_k * exp(-kappa * r_eff) / (eps * r_eff)

    r_eff     = max(r_ik, softening)
    eps       protein interior dielectric            (default 4.0)
    kappa     = 1 / kappa_inv, kappa_inv = 8.0 A     (~150 mM NaCl)
    softening minimum distance for 1/r               (default 0.5 A)
    scale     conversion to kT/e units               (default 1200.0)

Each atom contributes within ``cutoff_factor * kappa_inv`` Angstrom only
(default 7 x 8 = 56 A, exp(-7) ~ 9e-4), making the cost linear in protein size.

Validated values (kT/e at self-position, softening=0.5 A, scale=1200):
    Asp OD1  (q_grp = -0.80 e) ->  -452 kT/e
    Arg NH2  (q_grp = +0.03 e) ->   +19 kT/e
    Arg NE   (q_grp = -0.18 e) ->  -104 kT/e
    bb N-H   (q_grp = -0.07 e) ->   -41 kT/e
    bb C=O   (q_grp = -0.59 e) ->  -332 kT/e
    Lys NZ   (q_grp = +0.63 e) ->  +358 kT/e

The output is OpenDX, readable by ``gridData.Grid`` like an APBS map.
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

import numpy as np

from .cliutil import exit_on_error
from .grid import GridMap, accumulate, require_positive
from .structure import Structure

log = logging.getLogger(__name__)

DEFAULT_EPS = 4.0
DEFAULT_KAPPA_INV = 8.0
DEFAULT_SOFTENING = 0.5
DEFAULT_SCALE = 1200.0
DEFAULT_SPACING = 1.0
DEFAULT_PADDING = 5.0
DEFAULT_CUTOFF_FACTOR = 7.0


def compute_mep(
    structure: Structure,
    spacing: float = DEFAULT_SPACING,
    padding: float = DEFAULT_PADDING,
    eps: float = DEFAULT_EPS,
    kappa_inv: float = DEFAULT_KAPPA_INV,
    softening: float = DEFAULT_SOFTENING,
    scale: float = DEFAULT_SCALE,
    cutoff_factor: float = DEFAULT_CUTOFF_FACTOR,
    n_jobs: Optional[int] = 1,
) -> GridMap:
    """Compute the MEP grid (float32, kT/e) of a structure with charges (PQR)."""
    require_positive(
        spacing=spacing, padding=padding, eps=eps, kappa_inv=kappa_inv,
        softening=softening, cutoff_factor=cutoff_factor,
    )
    positions, group_q = structure.heavy_group_charges()

    kappa = 1.0 / kappa_inv
    cutoff = cutoff_factor * kappa_inv
    lo = positions.min(axis=0) - padding
    hi = positions.max(axis=0) + padding
    shape = tuple(int(v) for v in np.ceil((hi - lo) / spacing).astype(int) + 1)
    log.info(
        "MEP: %d heavy atoms -> grid %dx%dx%d (%s points, spacing %s A, cutoff %.1f A)",
        len(positions), *shape, f"{np.prod(shape):,}", spacing, cutoff,
    )

    def kernel(r: np.ndarray) -> np.ndarray:
        r_eff = np.maximum(r, softening)
        return np.exp(-kappa * r_eff) / (eps * r_eff)

    data = accumulate(
        positions, group_q, kernel, lo, spacing, shape,
        cutoff=cutoff, dtype=np.float32, n_jobs=n_jobs,
    )
    data *= np.float32(scale)
    return GridMap(data, lo, float(spacing))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="MEP",
        description="Screened Coulomb molecular electrostatic potential grid from a PQR file.",
    )
    p.add_argument("pqr", help="input PQR file (PDB2PQR AMBER force-field output)")
    p.add_argument("output", help="output .dx file (OpenDX format)")
    p.add_argument("--spacing", type=float, default=DEFAULT_SPACING,
                   help=f"grid spacing in Angstrom (default: {DEFAULT_SPACING})")
    p.add_argument("--padding", type=float, default=DEFAULT_PADDING,
                   help=f"padding around atom bounding box in Angstrom (default: {DEFAULT_PADDING})")
    p.add_argument("--eps", type=float, default=DEFAULT_EPS,
                   help=f"dielectric constant (default: {DEFAULT_EPS})")
    p.add_argument("--kappa-inv", type=float, default=DEFAULT_KAPPA_INV,
                   help=f"Debye screening length in Angstrom (default: {DEFAULT_KAPPA_INV})")
    p.add_argument("--softening", type=float, default=DEFAULT_SOFTENING,
                   help=f"minimum effective distance for 1/r in Angstrom (default: {DEFAULT_SOFTENING})")
    p.add_argument("--scale", type=float, default=DEFAULT_SCALE,
                   help=f"scale factor to kT/e units (default: {DEFAULT_SCALE})")
    p.add_argument("--cutoff-factor", type=float, default=DEFAULT_CUTOFF_FACTOR,
                   help="per-atom contribution radius in units of kappa_inv "
                        f"(default: {DEFAULT_CUTOFF_FACTOR}; larger = more accurate but slower)")
    p.add_argument("--threads", type=int, default=1,
                   help="number of threads (0 = all available; default: 1)")
    return p


@exit_on_error
def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    structure = Structure.from_file(args.pqr)
    result = compute_mep(
        structure, spacing=args.spacing, padding=args.padding, eps=args.eps,
        kappa_inv=args.kappa_inv, softening=args.softening, scale=args.scale,
        cutoff_factor=args.cutoff_factor, n_jobs=args.threads,
    )
    result.write_dx(args.output)
    log.info(
        "MEP written to %s (range %.1f to %.1f kT/e)",
        args.output, result.data.min(), result.data.max(),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
