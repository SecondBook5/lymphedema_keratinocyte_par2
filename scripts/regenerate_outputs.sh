#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: regenerate_outputs.sh [--cohort COHORT] [--annotated PATH] [--python PYTHON]

Regenerate tracked keratinocyte/PAR2 figures, tables, and notes from an annotated h5ad file.
Defaults:
  --cohort    manuscript_cohort
  --annotated work/<cohort>/results/adata_annotated.h5ad
  --python    $SCRNA_PYTHON or python
EOF
}

COHORT="manuscript_cohort"
PYTHON_BIN="${SCRNA_PYTHON:-python}"
ANNOTATED=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cohort)
            COHORT="$2"
            shift 2
            ;;
        --annotated)
            ANNOTATED="$2"
            shift 2
            ;;
        --python)
            PYTHON_BIN="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -z "$ANNOTATED" ]]; then
    ANNOTATED="$REPO_ROOT/work/$COHORT/results/adata_annotated.h5ad"
fi

echo "Regenerating curated keratinocyte outputs from: $ANNOTATED"

"$PYTHON_BIN" "$SCRIPT_DIR/keratinocyte_f2rl1_analysis.py" --adata "$ANNOTATED"
"$PYTHON_BIN" "$SCRIPT_DIR/keratinocyte_state_analysis.py" --adata "$ANNOTATED"
