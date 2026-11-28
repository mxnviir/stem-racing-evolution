#!/bin/bash
# Where a launched run is (carhybrid/aws/README.md). On your computer:
#   bash carhybrid/aws/status.sh [run_id]        (default: the latest run launched from this computer)
#   TAIL=200 bash carhybrid/aws/status.sh <run_id>   (more of the running case's log)
set -uo pipefail
CONF=${CONF:-$HOME/.config/carhybrid/aws.env}; . "$CONF"; export AWS_DEFAULT_REGION=$REGION
RUN=${1:-$(ls -t "$(dirname "$CONF")/runs" 2>/dev/null | head -1)}
[ -n "$RUN" ] || { echo "no run id"; exit 1; }
P=s3://$BUCKET/runs/$RUN
IID=$(aws s3 cp --quiet "$P/instance.json" - 2>/dev/null | python3 -c "import json,sys; print(json.load(sys.stdin)['instance_id'])" 2>/dev/null)
echo "run $RUN  instance ${IID:-?}: $( [ -n "$IID" ] && aws ec2 describe-instances --instance-ids "$IID" \
     --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null || echo unknown)"
if aws s3 cp --quiet "$P/DONE" - 2>/dev/null; then echo "(batch ended)"; fi
echo "=== status"; aws s3 cp --quiet "$P/status.json" - 2>/dev/null || echo "(no status yet: the machine is still booting / installing)"
echo; echo "=== instance log (tail)"; aws s3 cp --quiet "$P/instance.log" - 2>/dev/null | tail -15
echo "=== running case (tail)"; aws s3 cp --quiet "$P/current_job.log" - 2>/dev/null | tail -"${TAIL:-25}"
aws s3 cp --quiet "$P/boot.log" - 2>/dev/null | tail -20 | sed 's/^/boot: /'
echo "=== results in S3"; aws s3 ls "$P/out/" 2>/dev/null || echo "(none yet)"
