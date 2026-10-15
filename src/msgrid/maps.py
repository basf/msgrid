"""Run MLP and MEP on the same single structure."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Union

from .mep import compute_mep
from .mlp import compute_mlp
from .pqr import pdb_to_pqr
from .structure import Structure

log = logging.getLogger(__name__)


def run_maps(
    structure: Union[str, Path, Structure],
    prefix: Optional[Union[str, Path]] = None,
    mlp: bool = True,
    mep: bool = True,
    mlp_options: Optional[dict] = None,
    mep_options: Optional[dict] = None,
    n_jobs: Optional[int] = 1,
    pdb2pqr_options: Optional[dict] = None,
) -> Dict[str, Path]:
    """Compute MLP and/or MEP for one structure and write them as OpenDX files.

    The structure is parsed once and shared by both calculations. MEP needs
    atomic charges, so for ``mep=True`` the input must be a PQR file, unless
    ``pdb2pqr_options`` (a dict of :func:`pdb_to_pqr` keyword arguments, may be
    empty) is given: then a PDB path is first converted to ``<prefix>.pqr`` with
    PDB2PQR and that PQR is used for both maps. MLP alone only needs coordinates
    and atom/residue names (PDB or PQR).

    ``prefix`` defaults to the input path without extension; outputs are
    ``<prefix>_MLP.dx`` and ``<prefix>_MEP.dx``. Returns ``{'MLP': path, 'MEP': path}``
    for the maps written, plus ``'PQR'`` when a conversion was run.
    """
    if not (mlp or mep):
        raise ValueError("Nothing to do: both MLP and MEP are disabled")
    outputs: Dict[str, Path] = {}
    if not isinstance(structure, Structure):
        path = Path(structure)
        if prefix is None:
            prefix = path.with_suffix("")
        if mep and pdb2pqr_options is not None and path.suffix.lower() != ".pqr":
            pqr_path = Path(str(prefix) + ".pqr")
            if pqr_path.exists():
                log.warning("Overwriting existing %s", pqr_path)
            outputs["PQR"] = pdb_to_pqr(path, pqr_path, **pdb2pqr_options)
            path = pqr_path
        structure = Structure.from_file(path)
    if prefix is None:
        raise ValueError("'prefix' is required when passing a Structure object")
    prefix = Path(prefix)

    if mep and not structure.has_charges:
        raise ValueError(
            f"MEP requires atomic charges (PQR input); {structure.source or 'structure'} has none. "
            "Provide a PQR file, enable the PDB2PQR conversion, or disable MEP."
        )

    if mlp:
        out = prefix.with_name(prefix.name + "_MLP.dx")
        compute_mlp(structure, n_jobs=n_jobs, **(mlp_options or {})).write_dx(out)
        outputs["MLP"] = out
        log.info("MLP written to %s", out)
    if mep:
        out = prefix.with_name(prefix.name + "_MEP.dx")
        compute_mep(structure, n_jobs=n_jobs, **(mep_options or {})).write_dx(out)
        outputs["MEP"] = out
        log.info("MEP written to %s", out)
    return outputs
