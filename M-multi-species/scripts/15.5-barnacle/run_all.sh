#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

INPUT_DIR="$PROJECT_ROOT/M-multi-species/output/14-pca-orthologs"
OUT_BASE="$PROJECT_ROOT/M-multi-species/output/15.5-barnacle"

mkdir -p "$OUT_BASE"

ranks=(3 8 10 12 16 20)
for r in "${ranks[@]}"; do
  OUT_DIR="$OUT_BASE/rank-$r"
  mkdir -p "$OUT_DIR"
  echo "[15-barnacle] Running rank=$r -> $OUT_DIR"
  uv run python "$PROJECT_ROOT/M-multi-species/scripts/15.5-barnacle/build_tensor_and_run.py" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUT_DIR" \
    --rank "$r" --lambda-gene 0.1 --lambda-sample 0.1 --lambda-time 0.05 \
    --max-iter 1000 --tol 1e-5 --seed 42 | tee "$OUT_DIR/run.log"
done

echo "[15.5-barnacle] Generating summary"
uv run python "$PROJECT_ROOT/M-multi-species/scripts/15.5-barnacle/summarize_runs.py" \
  --base-dir "$OUT_BASE"

echo "[15.5-barnacle] Done. Summary at $OUT_BASE/SUMMARY.md"

