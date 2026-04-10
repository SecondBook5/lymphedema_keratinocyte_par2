#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_core_rerun.sh [--cohort COHORT] [--matrix-root PATH] [--python PYTHON]

Build runtime assets for a cohort rerun and execute the paired_sc core workflow.
Defaults:
  --cohort      manuscript_cohort
  --matrix-root $SCRNA_MATRIX_ROOT
  --python      $SCRNA_PYTHON or python
EOF
}

COHORT="manuscript_cohort"
MATRIX_ROOT="${SCRNA_MATRIX_ROOT:-}"
PYTHON_BIN="${SCRNA_PYTHON:-python}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cohort)
            COHORT="$2"
            shift 2
            ;;
        --matrix-root)
            MATRIX_ROOT="$2"
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
WORKDIR="$REPO_ROOT/work/$COHORT"

build_args=(
    "$SCRIPT_DIR/build_runtime_assets.py"
    --cohort "$COHORT"
    --workdir "$WORKDIR"
    --validate
)

if [[ -n "$MATRIX_ROOT" ]]; then
    build_args+=(--matrix-root "$MATRIX_ROOT")
fi

"$PYTHON_BIN" "${build_args[@]}"

PROJECT_PATH="$WORKDIR/runtime/project.yaml"
MANIFEST_PATH="$WORKDIR/runtime/manifest.csv"

echo "Running paired_sc core workflow..."
echo "  project: $PROJECT_PATH"
echo "  manifest: $MANIFEST_PATH"
echo "  workdir: $WORKDIR"

"$PYTHON_BIN" -m paired_sc.cli.app run core \
    --project "$PROJECT_PATH" \
    --manifest "$MANIFEST_PATH" \
    --workdir "$WORKDIR"

