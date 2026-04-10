#!/usr/bin/env python3
from __future__ import annotations

from figure_g_panel import parse_args, run_panel


def main() -> None:
    args = parse_args(panel_required=False, show_panel=False)
    run_panel("D", args)


if __name__ == "__main__":
    main()
