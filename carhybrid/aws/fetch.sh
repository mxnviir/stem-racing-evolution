#!/bin/bash
# Copy a run's results from S3 into carhybrid/results/<batch>/ (the same files the GitHub workflow's collect step
# commits), plus the instance's log and status under carhybrid/results/<batch>/_aws/<run_id>/. On your computer,
# from the repository root; then review and commit:
#   bash carhybrid/aws/fetch.sh <run_id>
#   git add carhybrid/results/<batch> && git commit -m "carhybrid <batch>: results (AWS run <run_id>)" && git push
set -euo pipefail
CONF=${CONF:-$HOME/.config/carhybrid/aws.env}; . "$CONF"; export AWS_DEFAULT_REGION=$REGION
RUN=${1:?usage: fetch.sh <run_id>}; P=s3://$BUCKET/runs/$RUN
cd "$(git rev-parse --show-toplevel)"
BATCH=$(aws s3 cp --quiet "$P/manifest.json" - | python3 -c "import json,sys; print(json.load(sys.stdin)['batch'])")
D=carhybrid/results/$BATCH
mkdir -p "$D/_aws/$RUN"
aws s3 sync --only-show-errors "$P/out/" "$D/" --exclude "*" --include "*.json" --include "*.log" --include "*.out" \
    --include "*.md" --include "*.png" --include "*.npz"
for f in manifest.json status.json instance.log boot.log DONE; do aws s3 cp --quiet "$P/$f" "$D/_aws/$RUN/$f" 2>/dev/null || true; done
# GitHub refuses files over 100 MB; the workflow keeps npz up to 90 MB
find "$D" -name "*.npz" -size +90M -print -delete | sed 's/^/too large for git, left in S3: /'
echo "fetched into $D:"; ls "$D"
python3 - "$D" <<'PY'
import json, sys
from pathlib import Path
for f in sorted(Path(sys.argv[1]).glob("*/fine_adjoint.json")):
    j = json.load(open(f)); r = j.get("resources") or {}
    print(f"{f.parent.name}: cells {j.get('cells')}, peak {r.get('peak_GB_whole_machine')} GB, F1 {j.get('F1_pass')} "
          f"({j.get('F1_dafoam_vs_pipeline_drag')}), F2 {j.get('F2_pass')} ({j.get('F2_span_mN')} mN), F3 {j.get('F3_pass')}, "
          f"grad {j.get('grad_json')}")
PY
