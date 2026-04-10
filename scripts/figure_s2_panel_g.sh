#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${SCRNA_PYTHON:-python}"

exec "$PYTHON_BIN" "$SCRIPT_DIR/figure_s2_panel_g.py" "$@"
