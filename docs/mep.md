# MEP: Molecular Electrostatic Potential

## Theory

The MEP is a screened Coulomb (Yukawa / Debye-Hückel) potential of the heavy-atom group charges. At a grid
point $i$:

$$V_i = s \sum_k \frac{q_k \, e^{-r_{\mathrm{eff}}/\lambda}}{\varepsilon \, r_{\mathrm{eff}}},
\qquad r_{\mathrm{eff}} = \max(r_{ik},\ \sigma)$$

| Symbol | Option | Default | Meaning |
|---|---|---|---|
| $\varepsilon$ | `--eps` | 4.0 | Dielectric constant (protein interior) |
| $\lambda = \kappa^{-1}$ | `--kappa-inv` | 8.0 Å | Debye screening length (about 150 mM NaCl) |
| $\sigma$ | `--softening` | 0.5 Å | Minimum distance in the $1/r$ term |
| $s$ | `--scale` | 1200 | Conversion factor to kT/e units |
| - | `--cutoff-factor` | 7 | Per-atom contribution radius, in units of $\lambda$ |
| - | `--padding` | 5.0 Å | Margin around the atoms' bounding box |
| - | `--spacing` | 1.0 Å | Grid spacing |

The result is in kT/e. It is an analytical screened potential: it does not solve the Poisson-Boltzmann
equation and does not model the solvent boundary.

### Effect of each parameter

- **`eps`**: larger values weaken the whole potential (it enters as $1/\varepsilon$).
- **`kappa-inv`**: the screening length $\lambda$. Smaller values screen more strongly (higher ionic
  strength) and give a shorter-ranged potential; larger values give a longer-ranged one. It also scales the
  cutoff radius.
- **`softening`**: prevents the divergence of $1/r$ at the atom positions. Grid points closer than $\sigma$
  to an atom see the same value as at distance $\sigma$.
- **`scale`**: a pure multiplicative factor; it only changes the units.
- **`cutoff-factor`**: each atom contributes only within `cutoff-factor x kappa-inv` (56 Å with the defaults).
  The neglected tail is of the order $e^{-7}\approx 9\times10^{-4}$ of the full term. Larger values are more
  accurate but slower; the cost is linear in the number of atoms for a fixed cutoff.
- **`padding`**: sets the box size. The potential decays with distance, so a larger padding adds points that
  are farther from the molecule and have small values.
- **`spacing`**: halving it multiplies the number of grid points by 8.

## Group charges

The potential is computed from the **heavy atoms only**. All heavy atoms of the PQR are used, including
`HETATM` records (ligands, ions, waters), unlike MLP, which uses `ATOM` records only. Remove atoms you do not
want in the potential from the PQR first (PDB2PQR can drop waters with `--drop-water` in a manual run).

Each heavy atom carries its own charge plus the charges of the hydrogens bonded to it: a hydrogen belongs to
a heavy atom of the same residue (same chain, residue number and insertion code) when it is closer than 1.3 Å.
Hydrogens are recognised from the element column or, if it is missing, from the atom name. A structure
without hydrogens is processed with the heavy atom charges as they are.

Heavy atoms with a group charge of absolute value below $10^{-9}$ are skipped.

## Input: why PQR

The charges come from the input file, so a plain PDB (no charges or radii) cannot be used. A PQR file holds a
partial charge and a radius for every atom. The charge values depend on the force field and on the
protonation state used to create the PQR, and so does the MEP.

You can either:

- create the PQR yourself with PDB2PQR and pass it to `MEP` or `comenze-map`, or
- let `comenze-map` do it with `--pdb2pqr` (needs `pdb2pqr` installed; options `--ff`, `--ph`, `--no-propka`).
  The generated `<prefix>.pqr` is kept and used for both MLP and MEP. The default force field for the
  conversion is `AMBER`.

The atomic radii in the PQR file are not used by the calculation.

## Grid and output

- The grid is the bounding box of the heavy atoms plus `--padding` on every side. Its origin is the lower
  corner of that box.
- Values are computed and stored as 32-bit floats and written to an OpenDX file (readable by VMD, PyMOL,
  Chimera, APBS tools and `gridData.Grid`).

## Reference values

With the default parameters ($\sigma=0.5$ Å, $s=1200$), the potential at the position of the atom itself is
approximately:

| Group | Group charge | MEP at its own position |
|---|---|---|
| Asp OD1 | -0.80 e | -452 kT/e |
| Arg NH2 | +0.03 e | +19 kT/e |
| Arg NE | -0.18 e | -104 kT/e |

These values are used as a sanity check.

## Options

| Option | Default | Meaning |
|---|---|---|
| `pqr` | - | Input PQR file |
| `output` | - | Output `.dx` file |
| `--spacing` | 1.0 | Grid spacing in Å |
| `--padding` | 5.0 | Padding around the bounding box in Å |
| `--eps` | 4.0 | Dielectric constant |
| `--kappa-inv` | 8.0 | Debye screening length in Å |
| `--softening` | 0.5 | Minimum distance in $1/r$ in Å |
| `--scale` | 1200 | Scale factor to kT/e |
| `--cutoff-factor` | 7 | Cutoff radius in units of `kappa-inv` |
| `--threads` | 1 | Number of threads (0 = all available) |

Within `comenze-map` the same options are `--mep-spacing`, `--mep-padding`, `--mep-eps`, `--mep-kappa-inv`,
`--mep-softening`, `--mep-scale` and `--mep-cutoff-factor`.
