"""Grid accumulation (pure NumPy) and OpenDX output shared by MLP and MEP."""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Tuple, Union

import numpy as np

# Atoms whose weight is below this are skipped (no contribution).
MIN_WEIGHT = 1e-9


@dataclass
class GridMap:
    """A regular 3-D scalar field."""

    data: np.ndarray  # (NX, NY, NZ)
    origin: np.ndarray  # (3,) Angstrom
    spacing: float  # Angstrom

    def write_dx(self, path: Union[str, Path]) -> None:
        write_dx(self.data, self.origin, self.spacing, path)


def require_positive(**values: Optional[float]) -> None:
    """Raise ``ValueError`` for any given option that is not a finite number > 0."""
    for name, value in values.items():
        if value is not None and not (np.isfinite(value) and value > 0):
            raise ValueError(f"{name} must be a positive number (got {value})")


def resolve_threads(n_jobs: Optional[int]) -> int:
    """``None``/``0`` means all CPUs available to this process."""
    if n_jobs:
        return max(1, int(n_jobs))
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def accumulate(
    positions: np.ndarray,
    weights: np.ndarray,
    kernel: Callable[[np.ndarray], np.ndarray],
    origin: np.ndarray,
    spacing: float,
    shape: Tuple[int, int, int],
    cutoff: Optional[float] = None,
    dtype=np.float64,
    n_jobs: Optional[int] = 1,
) -> np.ndarray:
    """Sum ``weights[k] * kernel(r_k)`` over all atoms ``k`` on a regular grid.

    ``kernel`` maps distances (Angstrom) to values. If ``cutoff`` is given,
    each atom only contributes within that radius (applied on the atom's
    local grid box, so cost is linear in the number of atoms); otherwise the
    whole grid is evaluated for every atom.

    The x axis is split into slabs processed by ``n_jobs`` threads; slabs are
    disjoint, so no locking is needed.
    """
    shape_a = np.array(shape, dtype=np.int64)
    origin = np.asarray(origin, dtype=np.float64)
    grid = np.zeros(shape, dtype=dtype)

    keep = np.abs(weights) >= MIN_WEIGHT
    pos = np.asarray(positions, dtype=np.float64)[keep]
    w = np.asarray(weights, dtype=np.float64)[keep]
    if len(w) == 0:
        return grid

    n_cut = int(shape_a.max()) if cutoff is None else int(np.ceil(cutoff / spacing))
    idx = np.rint((pos - origin) / spacing).astype(np.int64)
    lo = np.clip(idx - n_cut, 0, shape_a)
    hi = np.clip(idx + n_cut + 1, 0, shape_a)
    nonempty = np.all(lo < hi, axis=1)
    pos, w, lo, hi = pos[nonempty], w[nonempty], lo[nonempty], hi[nonempty]

    axes = [(origin[a] + np.arange(shape[a]) * spacing).astype(dtype) for a in range(3)]
    pos_t = pos.astype(dtype)
    r_max = np.inf if cutoff is None else float(cutoff)
    NX = shape[0]

    # Slab thickness so that one atom-slab temporary holds ~2M elements.
    box_yz = min(shape[1], 2 * n_cut + 1) * min(shape[2], 2 * n_cut + 1)
    slab = max(1, int(2_000_000 // box_yz))
    threads = resolve_threads(n_jobs)
    if threads > 1:
        slab = min(slab, max(1, -(-NX // (threads * 2))))

    def process(x0: int) -> None:
        x1 = min(x0 + slab, NX)
        for k in np.flatnonzero((lo[:, 0] < x1) & (hi[:, 0] > x0)):
            xa, xb = max(lo[k, 0], x0), min(hi[k, 0], x1)
            ya, yb = lo[k, 1], hi[k, 1]
            za, zb = lo[k, 2], hi[k, 2]
            dx = (axes[0][xa:xb] - pos_t[k, 0]).reshape(-1, 1, 1)
            dy = (axes[1][ya:yb] - pos_t[k, 1]).reshape(1, -1, 1)
            dz = (axes[2][za:zb] - pos_t[k, 2]).reshape(1, 1, -1)
            r = np.sqrt(dx * dx + dy * dy + dz * dz)
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                val = kernel(r)
                val *= float(w[k])
            if cutoff is not None:
                val[r > r_max] = 0.0
            grid[xa:xb, ya:yb, za:zb] += val

    starts = range(0, NX, slab)
    if threads == 1:
        for x0 in starts:
            process(x0)
    else:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            list(pool.map(process, starts))
    return grid


_DX_FOOTER = (
    'attribute "dep" string "positions"\n'
    'object "regular positions regular connections" class field\n'
    'component "positions" value 1\n'
    'component "connections" value 2\n'
    'component "data" value 3\n'
)


def write_dx(grid: np.ndarray, origin, spacing: float, path: Union[str, Path]) -> None:
    """Write a grid in OpenDX format (APBS-compatible, readable by ``gridData.Grid``)."""
    NX, NY, NZ = grid.shape
    ox, oy, oz = (float(v) for v in origin)
    n_total = NX * NY * NZ
    header = (
        f"object 1 class gridpositions counts {NX} {NY} {NZ}\n"
        f"origin {ox:.6f} {oy:.6f} {oz:.6f}\n"
        f"delta {spacing:.6f} 0.000000 0.000000\n"
        f"delta 0.000000 {spacing:.6f} 0.000000\n"
        f"delta 0.000000 0.000000 {spacing:.6f}\n"
        f"object 2 class gridconnections counts {NX} {NY} {NZ}\n"
        f"object 3 class array type double rank 0 items {n_total} data follows\n"
    )
    flat = grid.ravel(order="C")
    full = (n_total // 3) * 3
    chunk = 3 * 200_000
    with open(path, "w") as fh:
        fh.write(header)
        for s in range(0, full, chunk):
            block = flat[s : min(s + chunk, full)].reshape(-1, 3)
            np.savetxt(fh, block, fmt="%.6e", delimiter="  ")
        if full < n_total:
            fh.write("  ".join(f"{float(v):.6e}" for v in flat[full:]) + "\n")
        fh.write(_DX_FOOTER)
