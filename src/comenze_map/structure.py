"""Minimal PDB / PQR reader for a single molecular object (first model only)."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Union

import numpy as np

# Covalent X-H bond cutoff (Angstrom) used to assign hydrogens to heavy atoms.
H_BOND_CUTOFF = 1.3

_PQR_ATOM_RE = re.compile(
    r"^(ATOM|HETATM)\s+"
    r"(\d+)\s+"  # serial
    r"(\S+)\s+"  # atom name
    r"(\S+?)\s+"  # residue name
    r"(?:([A-Za-z])\s+)?"  # optional chain id
    r"(-?\d+)([A-Za-z]?)\s+"  # residue number + optional insertion code
    r"(\S+)\s+(\S+)\s+(\S+)\s+"  # x y z
    r"(\S+)\s+(\S+)"  # charge radius
)


def _is_hydrogen_name(name: str) -> bool:
    """Name-based hydrogen test (``H...`` or PDB-style ``1H..``)."""
    return name[:1] == "H" or (len(name) > 1 and name[0].isdigit() and name[1] == "H")


@dataclass
class Structure:
    """Atoms of a single molecule.

    ``charge`` and ``radius`` are NaN when read from a plain PDB file.
    """

    record: np.ndarray  # 'ATOM' / 'HETATM'
    name: np.ndarray
    resname: np.ndarray
    chain: np.ndarray
    resnum: np.ndarray
    icode: np.ndarray
    xyz: np.ndarray  # (N, 3) float64
    charge: np.ndarray
    radius: np.ndarray
    is_hydrogen: np.ndarray  # bool
    source: str = ""

    def __len__(self) -> int:
        return len(self.name)

    @property
    def has_charges(self) -> bool:
        return bool(len(self)) and not np.isnan(self.charge).any()

    @classmethod
    def from_file(cls, path: Union[str, Path], fmt: str = "auto") -> "Structure":
        """Read ``path`` as ``pdb`` or ``pqr`` (``auto`` uses the file extension)."""
        path = Path(path)
        if fmt == "auto":
            fmt = "pqr" if path.suffix.lower() == ".pqr" else "pdb"
        if fmt not in ("pdb", "pqr"):
            raise ValueError(f"Unknown format {fmt!r}; expected 'pdb', 'pqr' or 'auto'")
        with open(path) as fh:
            rows = _parse_pqr(fh, str(path)) if fmt == "pqr" else _parse_pdb(fh, str(path))
        if not rows:
            raise ValueError(f"No ATOM/HETATM records found in {path}")
        return cls._from_rows(rows, str(path))

    @classmethod
    def _from_rows(cls, rows: List[tuple], source: str) -> "Structure":
        rec, name, resname, chain, resnum, icode, x, y, z, q, r, h = zip(*rows)
        return cls(
            record=np.array(rec),
            name=np.array(name),
            resname=np.array(resname),
            chain=np.array(chain),
            resnum=np.array(resnum, dtype=np.int64),
            icode=np.array(icode),
            xyz=np.array([x, y, z], dtype=np.float64).T,
            charge=np.array(q, dtype=np.float64),
            radius=np.array(r, dtype=np.float64),
            is_hydrogen=np.array(h, dtype=bool),
            source=source,
        )

    def heavy_group_charges(self) -> Tuple[np.ndarray, np.ndarray]:
        """Heavy-atom positions and group charges (atom + bonded H charges).

        A hydrogen belongs to a heavy atom of the same residue when closer
        than ``H_BOND_CUTOFF``.
        """
        if not self.has_charges:
            raise ValueError(
                "Atomic charges are required (use a PQR file); "
                f"none found in {self.source or 'structure'}"
            )
        heavy = np.flatnonzero(~self.is_hydrogen)
        if heavy.size == 0:
            raise ValueError(f"No heavy atoms found in {self.source or 'structure'}")
        hyd = np.flatnonzero(self.is_hydrogen)

        positions = self.xyz[heavy]
        group_q = self.charge[heavy].copy()

        if hyd.size:
            keys = [
                (c, n, i)
                for c, n, i in zip(self.chain, self.resnum, self.icode)
            ]
            h_by_res: dict = {}
            for idx in hyd:
                h_by_res.setdefault(keys[idx], []).append(idx)
            for pos_i, atom_i in enumerate(heavy):
                h_idx = h_by_res.get(keys[atom_i])
                if not h_idx:
                    continue
                diff = self.xyz[h_idx] - self.xyz[atom_i]
                close = (diff * diff).sum(axis=1) < H_BOND_CUTOFF**2
                group_q[pos_i] += self.charge[h_idx][close].sum()
        return positions, group_q


def _parse_pdb(lines: Iterable[str], source: str = "") -> List[tuple]:
    rows = []
    altloc: dict = {}  # atom key -> alternate-location id kept (first one seen)
    for lineno, line in enumerate(lines, 1):
        rec = line[:6].strip()
        if rec == "ENDMDL":
            break
        if rec not in ("ATOM", "HETATM"):
            continue
        try:
            name = line[12:16].strip()
            chain = line[21:22].strip()
            resnum = int(line[22:26])
            icode = line[26:27].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            raise ValueError(f"{source}:{lineno}: cannot parse PDB record") from None
        alt = line[16:17].strip()
        if alt and altloc.setdefault((chain, resnum, icode, name), alt) != alt:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        is_h = element == "H" if element else _is_hydrogen_name(name)
        rows.append(
            (rec, name, line[17:20].strip(), chain, resnum, icode, x, y, z, np.nan, np.nan, is_h)
        )
    return rows


def _is_pqr_hydrogen(record: str, name: str, resname: str) -> bool:
    # Single-atom HETATM residues named like their element (HG, HO, ...) are not hydrogens.
    if record == "HETATM" and name.upper() == resname.upper():
        return False
    return _is_hydrogen_name(name)


def _parse_pqr_fixed(line: str) -> Optional[tuple]:
    """Fixed-column fallback for PQR lines whose fields run together."""
    try:
        return (
            line[:6].strip(),
            line[12:16].strip(),
            line[17:20].strip(),
            line[21:22].strip(),
            int(line[22:26]),
            line[26:27].strip(),
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
            float(line[54:62]),
            float(line[62:70]),
        )
    except ValueError:
        return None


def _looks_like_pdb(line: str) -> bool:
    """True for PDB records with occupancy, B-factor and element columns."""
    def is_float(text: str) -> bool:
        try:
            float(text)
        except ValueError:
            return False
        return True

    return (
        len(line) >= 78
        and line[76:78].strip().isalpha()
        and not line[66:76].strip()
        and is_float(line[54:60])
        and is_float(line[60:66])
    )


def _parse_pqr(lines: Iterable[str], source: str) -> List[tuple]:
    rows = []
    for lineno, line in enumerate(lines, 1):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if not rows and _looks_like_pdb(line):
            raise ValueError(
                f"{source}: looks like a PDB file (occupancy/B-factor columns), not a PQR, "
                "so it has no atomic charges; convert it with PDB2PQR first"
            )
        m = _PQR_ATOM_RE.match(line)
        if m:
            row = (
                m.group(1), m.group(3), m.group(4), m.group(5) or "", int(m.group(6)), m.group(7),
                float(m.group(8)), float(m.group(9)), float(m.group(10)),
                float(m.group(11)), float(m.group(12)),
            )
        else:
            row = _parse_pqr_fixed(line)
            if row is None:
                raise ValueError(f"{source}:{lineno}: cannot parse PQR record")
        rows.append((*row, _is_pqr_hydrogen(row[0], row[1], row[2])))
    return rows
