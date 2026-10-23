"""ab1kit entry point: the ``ab1kit`` command launches the GUI."""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    from ui.main_window import run

    return run(argv if argv is not None else sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main())
