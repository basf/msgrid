"""Helpers shared by the command-line entry points."""

from __future__ import annotations

import functools
import sys
from typing import Callable


def exit_on_error(main: Callable[..., int]) -> Callable[..., int]:
    """Turn expected user errors into a one-line message and exit status 1."""

    @functools.wraps(main)
    def wrapper(*args, **kwargs) -> int:
        try:
            return main(*args, **kwargs)
        except (ValueError, KeyError, RuntimeError, OSError) as exc:
            msg = exc.args[0] if isinstance(exc, KeyError) and exc.args else exc
            print(f"error: {msg}", file=sys.stderr)
            return 1

    return wrapper
