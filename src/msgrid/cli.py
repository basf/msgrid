"""Command line wrapper: MLP + MEP on a single structure."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import mep as _mep
from . import mlp as _mlp
from . import pqr as _pqr
from .cliutil import exit_on_error
from .maps import run_maps
from .mlp_tables import TABLES


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="msgrid",
        description="Compute the MLP (lipophilic) and MEP (electrostatic) maps of a single "
                    "structure. MEP needs charges, so use a PQR input or add --pdb2pqr to convert "
                    "a PDB (MLP works on PDB or PQR).",
    )
    p.add_argument("input", help="input structure (.pqr for both maps; .pdb for MLP only, "
                                 "or for both with --pdb2pqr)")
    p.add_argument("-d", "--outdir", help="output directory (default: next to the input)")
    p.add_argument("-p", "--prefix", help="output name prefix (default: input file stem)")
    p.add_argument("--no-mlp", action="store_true", help="skip the MLP map")
    p.add_argument("--no-mep", action="store_true", help="skip the MEP map")
    p.add_argument("--threads", type=int, default=1,
                   help="number of threads (0 = all available; default: 1)")

    g = p.add_argument_group("PDB2PQR options (requires an external pdb2pqr installation)")
    g.add_argument("--pdb2pqr", action="store_true",
                   help="convert the PDB input to <prefix>.pqr with pdb2pqr and use it for both maps")
    g.add_argument("--ff", default=_pqr.DEFAULT_FORCEFIELD, help="force field (default: %(default)s)")
    g.add_argument("--ph", type=float, default=_pqr.DEFAULT_PH,
                   help="pH for protonation (default: %(default)s)")
    g.add_argument("--no-propka", action="store_true", help="do not run PROPKA for pKa estimation")

    g = p.add_argument_group("MLP options")
    g.add_argument("--mlp-table", default="Crippen", choices=[*TABLES, "Readin"])
    g.add_argument("--mlp-logp", help="logP data file for --mlp-table Readin")
    g.add_argument("--mlp-method", default="fauchere", choices=list(_mlp.METHODS))
    g.add_argument("--mlp-spacing", type=float, default=_mlp.DEFAULT_SPACING)
    g.add_argument("--mlp-cutoff", type=float, default=_mlp.DEFAULT_CUTOFF)

    g = p.add_argument_group("MEP options")
    g.add_argument("--mep-spacing", type=float, default=_mep.DEFAULT_SPACING)
    g.add_argument("--mep-padding", type=float, default=_mep.DEFAULT_PADDING)
    g.add_argument("--mep-eps", type=float, default=_mep.DEFAULT_EPS)
    g.add_argument("--mep-kappa-inv", type=float, default=_mep.DEFAULT_KAPPA_INV)
    g.add_argument("--mep-softening", type=float, default=_mep.DEFAULT_SOFTENING)
    g.add_argument("--mep-scale", type=float, default=_mep.DEFAULT_SCALE)
    g.add_argument("--mep-cutoff-factor", type=float, default=_mep.DEFAULT_CUTOFF_FACTOR)
    return p


@exit_on_error
def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    src = Path(args.input)
    prefix = Path(args.prefix) if args.prefix else Path(src.stem)
    outdir = Path(args.outdir) if args.outdir else src.parent
    outdir.mkdir(parents=True, exist_ok=True)

    run_maps(
        src,
        prefix=outdir / prefix,
        mlp=not args.no_mlp,
        mep=not args.no_mep,
        mlp_options=dict(
            table=args.mlp_table, logp_file=args.mlp_logp, method=args.mlp_method,
            spacing=args.mlp_spacing, cutoff=args.mlp_cutoff,
        ),
        mep_options=dict(
            spacing=args.mep_spacing, padding=args.mep_padding, eps=args.mep_eps,
            kappa_inv=args.mep_kappa_inv, softening=args.mep_softening,
            scale=args.mep_scale, cutoff_factor=args.mep_cutoff_factor,
        ),
        n_jobs=args.threads,
        pdb2pqr_options=(
            dict(forcefield=args.ff, ph=args.ph, propka=not args.no_propka)
            if args.pdb2pqr else None
        ),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
