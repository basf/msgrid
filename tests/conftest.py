from pathlib import Path

import pytest

# name, resname, resnum, xyz, charge ; hydrogens are named H*
_ATOMS = [
    ("N", "ALA", 1, (0.000, 0.000, 0.000), -0.4157),
    ("H", "ALA", 1, (0.000, 1.000, 0.000), 0.2719),
    ("CA", "ALA", 1, (1.458, -0.300, 0.000), 0.0337),
    ("CB", "ALA", 1, (1.900, -1.100, 1.200), -0.1825),
    ("HB1", "ALA", 1, (2.900, -1.100, 1.200), 0.0603),
    ("C", "ALA", 1, (2.100, 1.050, 0.000), 0.5973),
    ("O", "ALA", 1, (1.500, 2.100, 0.300), -0.5679),
    ("N", "ASP", 2, (3.400, 1.000, -0.500), -0.5163),
    ("H", "ASP", 2, (3.900, 1.900, -0.500), 0.2936),
    ("CA", "ASP", 2, (4.200, -0.200, -1.000), 0.0381),
    ("CB", "ASP", 2, (5.600, 0.100, -1.500), -0.0303),
    ("CG", "ASP", 2, (6.500, -1.000, -2.000), 0.7994),
    ("OD1", "ASP", 2, (6.300, -2.100, -1.500), -0.8014),
    ("OD2", "ASP", 2, (7.400, -0.700, -2.800), -0.8014),
    ("C", "ASP", 2, (4.300, -1.200, 0.200), 0.5366),
    ("O", "ASP", 2, (4.800, -2.300, 0.100), -0.5819),
    ("OXT", "ASP", 2, (3.800, -0.900, 1.400), -0.5819),
]


def pqr_text() -> str:
    lines = []
    for i, (name, res, num, (x, y, z), q) in enumerate(_ATOMS, 1):
        lines.append(
            f"ATOM  {i:5d} {name:<4s} {res:>3s} A{num:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f} {q:7.4f} 1.8000"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def pdb_text() -> str:
    lines = []
    for i, (name, res, num, (x, y, z), _q) in enumerate(_ATOMS, 1):
        elem = "H" if name.startswith("H") else name[0]
        lines.append(
            f"ATOM  {i:5d}  {name:<3s} {res:>3s} A{num:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


@pytest.fixture
def pqr_file(tmp_path: Path) -> Path:
    p = tmp_path / "mol.pqr"
    p.write_text(pqr_text())
    return p


@pytest.fixture
def pdb_file(tmp_path: Path) -> Path:
    p = tmp_path / "mol.pdb"
    p.write_text(pdb_text())
    return p
