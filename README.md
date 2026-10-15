# msgrid

`msgrid` (molecular surface grid) computes physicochemical property maps of a single molecule, typically a
protein, on a regular 3D grid surrounding it. Two properties are available:

- **MLP** – Molecular Lipophilic Potential: where the surface is lipophilic (hydrophobic) or hydrophilic,
  from per-atom lipophilicity contributions.
- **MEP** – Molecular Electrostatic Potential (screened Coulomb): where the surface is positively or negatively
  charged, from the atomic partial charges.

The input is a PDB or PQR file. Each map is written as an OpenDX (`.dx`) file that can be loaded in molecular
viewers such as PyMOL, VMD or Chimera, to colour the molecular surface by the property or to display isosurfaces.
The maps are meant to describe and compare surface patches, for example to analyse hydrophobic or charged
regions of a protein. A single command computes both maps from the same parsed structure, and the package is
pure Python (NumPy only), so no compiler is needed to install it.

**More details are available in the [docs](docs/) folder:** the theory, all options and implementation notes
are described in [docs/mlp.md](docs/mlp.md) (lipophilic potential) and [docs/mep.md](docs/mep.md)
(electrostatic potential).

## Installation

```bash
pip install .            # or: pip install -e ".[test]" for development
```

Requires Python >= 3.9 and `numpy`. `gridDataFormats` is only used by the tests; output `.dx` files are
readable with `gridData.Grid`.

## Usage

### Wrapper: both maps on the same object

```bash
msgrid protein.pqr                       # -> protein_MLP.dx, protein_MEP.dx
msgrid protein.pqr -d out/ -p run1 --threads 8
msgrid protein.pdb --no-mep              # MLP only (a PDB has no charges)
```

MEP needs atomic charges and radii, so its input must be a PQR file. MLP only uses
coordinates and atom/residue names, so a PQR works for both. Run `msgrid -h` for all
`--mlp-*` / `--mep-*` options.

#### Starting from a PDB: PDB2PQR conversion

`--pdb2pqr` runs an external `pdb2pqr` on the PDB, writes `<prefix>.pqr`, and uses that PQR (the same
protonated, charged structure) for both maps. PDB2PQR is not a hard dependency; install it with
`pip install pdb2pqr` (or `pip install ".[pdb2pqr]"`), or the APBS-PDB2PQR distribution, so that `pdb2pqr`
is on `PATH`.

```bash
msgrid protein.pdb --pdb2pqr                          # AMBER, pH 7.0, PROPKA
msgrid protein.pdb --pdb2pqr --ff CHARMM --ph 6.5 --no-propka
```

```python
run_maps("protein.pdb", prefix="out/protein", pdb2pqr_options={"forcefield": "AMBER", "ph": 7.0})
```

Without `--pdb2pqr`, a PDB input only works with `--no-mep`. You can also run `pdb2pqr` yourself and pass
the PQR; the charges (and therefore the MEP) depend on the force field and protonation state chosen.
The conversion can add hydrogens and rename residues (e.g. `HIS` -> `HIE`), which MLP handles.

```python
from msgrid import run_maps, Structure, compute_mlp, compute_mep

run_maps("protein.pqr", prefix="out/protein")          # writes both maps
s = Structure.from_file("protein.pqr")
mlp = compute_mlp(s, method="fauchere", spacing=1.0)   # GridMap(data, origin, spacing)
mep = compute_mep(s, spacing=1.0)
mep.write_dx("mep.dx")
```

### Individual tools

```bash
MLP -f input.pdb -o output.dx [-t Crippen|Laguerre|Ghose|Readin] [-l logp.py] [-m fauchere|dubost|brasseur|buckingham|type5] [-s 1.0]
MEP input.pqr output.dx [--spacing 1.0] [--padding 5.0] [--eps 4.0] [--kappa-inv 8.0] [--softening 0.5] [--scale 1200] [--cutoff-factor 7]
```

Both tools accept `--threads N` (0 = all CPUs).

Theory, options and implementation details: [docs/mlp.md](docs/mlp.md) and [docs/mep.md](docs/mep.md).

## Visualizing the maps in PyMOL

The `.dx` files use the same coordinates as the input structure, so load them together:

```
load protein.pqr, protein
load protein_MLP.dx, mlp
load protein_MEP.dx, mep
```

Isosurfaces (the levels depend on the system; adjust them while looking at the map):

```
# MLP: positive = lipophilic, negative = hydrophilic
isosurface mlp_lipo,  mlp,  5
isosurface mlp_hydro, mlp, -5
color orange, mlp_lipo
color cyan,   mlp_hydro

# MEP in kT/e: blue = positive, red = negative
isosurface mep_pos, mep,  1
isosurface mep_neg, mep, -1
color blue, mep_pos
color red,  mep_neg
```

To colour the molecular surface by the MEP instead:

```
ramp_new mep_ramp, mep, [-5, 0, 5], [red, white, blue]
show surface, protein
set surface_color, mep_ramp, protein
```

Useful settings: `set transparency, 0.5, mep_pos` (and likewise for the other objects) to see the atoms through a
surface, and `isomesh` instead of `isosurface` for a mesh. The standalone `MEP` tool logs the value range of the
map when it finishes, which helps to choose the levels.

## Notes

- Built-in lightweight PDB/PQR reader. MLP uses non-hydrogen `ATOM` records (not `HETATM`); only the first
  model is read. Hydrogens are detected from the element column, or from the atom name when it is missing.
- PQR files with or without a chain ID are accepted.
- MLP residue-name variants (`HID/HIE/HIP`, `ASH`, `GLH`, `CYX/CYM`, `LYN`, ...) fall back to the standard residue.
- `-t Readin` parses the `Crippen = {...}` literal from the logP file instead of executing it.
- An unknown MLP method or reference table is an error; missing atoms in the table raise an error listing
  the residue/atom pairs.
- Both maps are written as OpenDX files by the same writer.

## Authors

Miguel Angel Bayo Valero, Dmitry Suplatov, Valerio Ferrario (BASF SE).

## Acknowledgements

This work was carried out within the COMENZE project, funded by the European Union under
Grant Agreement No. 101169327.

## MLP and pyMLP

The MLP part of this package is a reimplementation, in modern Python (Python 3, vectorised NumPy), of
[pyMLP](https://github.com/plantigrade/pyMLP) by Julien Lefeuvre, itself a platform-independent Python
equivalent of the MLPP program of Laguerre et al. It keeps the same lipophilicity-potential methods and
extends them:

- installable package with command-line tools and a Python API, instead of a single script for Python 2.4/2.5;
- vectorised NumPy calculation with an optional multi-threaded mode (`--threads`) and a per-atom distance
  cutoff;
- additional atomic lipophilicity tables (Crippen) and a `Readin` option for user-supplied logP data;
- support for PQR input, residue-name variants (`HID/HIE/HIP`, `ASH`, ...) and clear errors for unknown atoms;
- the same input structure can also be used for the MEP map in a single run.

If you use the MLP map, please cite the references below. pyMLP is distributed under the BSD 3-Clause
license; its copyright notice and license text are reproduced in [LICENSES/pyMLP.txt](LICENSES/pyMLP.txt).

## Crippen table

The built-in `Crippen` lipophilicity table (the default for `-t` / `--mlp-table`) was introduced by Dmitry
Suplatov (BASF SE) and calculated with [RDKit](https://www.rdkit.org), using its implementation of the
Wildman-Crippen atomic contributions (Wildman SA, Crippen GM. *J Chem Inf Comput Sci* 1999, 39, 868-873).
RDKit is not needed to use the package.

## Credits and references

- pyMLP: Julien Lefeuvre, https://github.com/plantigrade/pyMLP (BSD 3-Clause).
- Broto P, Moreau G, Vandycke C. Molecular structures: perception, autocorrelation descriptor and SAR studies.
  System of atomic contributions for the calculation of the n-octanol/water partition coefficients.
  *Eur J Med Chem* 1984, 19(1), 71-78.
- Laguerre M, Saux M, Dubost JP, Carpy A. MLPP: a program for the calculation of molecular lipophilicity
  potential in proteins. *Pharm Sci* 1997, 3(5-6), 217-222.
- Ghose AK, Viswanadhan VN, Wendoloski JJ. Prediction of hydrophobic (lipophilic) properties of small organic
  molecules using fragmental methods: an analysis of ALOGP and CLOGP methods. *J Phys Chem A* 1998, 102,
  3762-3772.
- Crippen table (introduced in the MLP calculation) and cutoff amendments: Dmitry Suplatov (BASF SE).

## License

BSD 3-Clause, see [LICENSE](LICENSE). The MLP part is a reimplementation of pyMLP, also BSD 3-Clause;
its notice is in [LICENSES/pyMLP.txt](LICENSES/pyMLP.txt).

## Tests

```bash
pip install -e ".[test]" && pytest
```
