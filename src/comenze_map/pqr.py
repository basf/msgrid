"""Optional PDB -> PQR conversion through an external PDB2PQR installation."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Sequence, Union

log = logging.getLogger(__name__)

DEFAULT_FORCEFIELD = "AMBER"
DEFAULT_PH = 7.0


def pdb_to_pqr(
    pdb: Union[str, Path],
    pqr: Union[str, Path],
    forcefield: str = DEFAULT_FORCEFIELD,
    ph: Optional[float] = DEFAULT_PH,
    propka: bool = True,
    extra_args: Sequence[str] = (),
    executable: Optional[str] = None,
) -> Path:
    """Run ``pdb2pqr`` on ``pdb`` and write the result to ``pqr``.

    ``ph=None`` skips the pH-dependent protonation (``propka`` is then ignored).
    """
    exe = executable or shutil.which("pdb2pqr") or shutil.which("pdb2pqr30")
    if not exe:
        raise RuntimeError(
            "pdb2pqr not found on PATH; install it (pip install pdb2pqr) "
            "or provide a PQR file instead of a PDB"
        )
    pdb, pqr = Path(pdb), Path(pqr)
    cmd = [exe, f"--ff={forcefield}", "--keep-chain"]
    if ph is not None:
        cmd.append(f"--with-ph={ph}")
        if propka:
            cmd.append("--titration-state-method=propka")
    cmd += [*extra_args, str(pdb), str(pqr)]

    log.info("Running: %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"pdb2pqr failed (exit {res.returncode}):\n{res.stderr.strip()}")
    if not pqr.exists() or pqr.stat().st_size == 0:
        raise RuntimeError(f"pdb2pqr produced no output file {pqr}")
    return pqr
