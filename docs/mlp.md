# MLP: Molecular Lipophilic Potential

## Theory

Each non-hydrogen atom $k$ of the molecule has an atomic lipophilicity contribution $f_k$ (the `fi`
value of a reference table; positive = lipophilic, negative = hydrophilic). The potential at a grid
point $i$ at distance $d_{ik}$ from atom $k$ is

$$\mathrm{MLP}_i = 100 \sum_k f_k \, g(d_{ik})$$

where $g$ is a distance-decay function chosen with `-m`. The result is a relative, unit-less lipophilicity
field; only comparisons between grids computed with the same table and method are meaningful.

### Distance functions

| `-m` | $g(d)$ | Notes |
|---|---|---|
| `fauchere` (default) | $e^{-d}$ | Evaluated only within `--cutoff` of each atom |
| `dubost` | $1/(1+d)$ | Long-range; evaluated on the whole grid |
| `brasseur` | $e^{-d/3.1}$ | Evaluated on the whole grid |
| `buckingham` | $d^{-3}$ | Distances below 0.1 Å are set to 0.1 Å to keep the grid finite; evaluated on the whole grid |
| `type5` | $e^{-\sqrt{d}}$ | Evaluated on the whole grid |

For `fauchere` the cutoff keeps the cost linear in the number of atoms. With the default 25 Å the neglected
tail is at most $100 \cdot 2 e^{-25} \approx 3\times10^{-9}$ per atom. The other methods have no cutoff,
so their cost grows with (atoms x grid points).

## Atom selection

- Only `ATOM` records that are not hydrogens are used (`HETATM` records and hydrogens are ignored).
- Only the first model of the file is read.
- For atoms with alternate locations only the first one found in the file is kept.
- Hydrogens are recognised from the element column, or from the atom name when the element is missing.
- Atoms with $|100 f_k| < 10^{-9}$ are skipped.
- Residue-name variants are mapped to the standard residue when the exact name is not in the table:
  `HID/HIE/HIP/HSD/HSE/HSP/HISD/HISE/HISH` -> `HIS`, `ASH/ASPH` -> `ASP`, `GLH/GLUH` -> `GLU`,
  `CYX/CYM` -> `CYS`, `LYN` -> `LYS`, `ARN` -> `ARG`; a leading `N`/`C` on a four-letter name
  (terminal residues such as `NALA`) is also removed.
- If an atom has no value in the selected table, the calculation stops with an error listing the
  residue/atom pairs.

## Reference tables (`-t`)

| Table | Description |
|---|---|
| `Crippen` (default) | Wildman-Crippen atomic contributions calculated with RDKit; table introduced by Dmitry Suplatov (BASF SE) |
| `Laguerre` | Table of Laguerre et al., *Pharm Sci* 1997, 3, 217-222 |
| `Ghose` | Table based on Ghose et al., *J Phys Chem A* 1998, 102, 3762-3772 |
| `Readin` | User-supplied table read from the file given with `-l` |

### Custom table (`-t Readin -l logp.py`)

The file is a Python file containing a dictionary named `Crippen`, mapping residue name -> atom name -> value:

```python
Crippen = {
    "ALA": {"N": -0.7096, "CA": -0.2051, "C": -0.2783, "O": -0.1526, "CB": 0.1441, "OXT": -0.2893},
    # ...
}
```

The file is parsed, not executed: only the literal assigned to `Crippen` is evaluated, so other code in the
file is ignored.

## Grid

- The grid is the bounding box of the atoms plus a fixed 10 Å margin on every side, adjusted so that its
  size is a whole number of spacings.
- `-s` (default 1.0 Å) is the spacing. Halving it multiplies the number of points by 8.
- Values are stored as 64-bit floats in an OpenDX file (readable by VMD, PyMOL, Chimera, APBS tools and
  `gridData.Grid`).

## Options

| Option | Default | Meaning |
|---|---|---|
| `-f` | - | Input PDB (or PQR) file |
| `-o` | - | Output `.dx` file |
| `-t` | `Crippen` | Reference table: `Crippen`, `Laguerre`, `Ghose`, `Readin` |
| `-l` | - | logP data file for `-t Readin` |
| `-m` | `fauchere` | Distance function (see above) |
| `-s` | 1.0 | Grid spacing in Å |
| `--cutoff` | 25.0 | Per-atom cutoff in Å (`fauchere` only) |
| `--threads` | 1 | Number of threads (0 = all available) |

Within `msgrid` the same options are `--mlp-table`, `--mlp-logp`, `--mlp-method`, `--mlp-spacing`
and `--mlp-cutoff`.

## References

See [the README](../README.md#credits-and-references). The Crippen table is described in
[the README](../README.md#crippen-table).
