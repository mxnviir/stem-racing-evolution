#!/bin/bash
# Launch one batch of fine-mesh DAFoam cases on a fresh EC2 instance with the AWS CLI (carhybrid/aws/README.md).
# Run on YOUR computer, from a clone of stem-racing-evolution on branch claude/wizardly-brahmagupta-1edalc, after
# aws_setup.sh. The instance runs the batch's cases one after another, copies each case's results to S3 as soon as it
# ends, and terminates itself when the batch is done (or after MAX_HOURS, whatever happens).
#
#   bash carhybrid/aws/launch.sh <batch>            # a key of carhybrid/request.json "batches", e.g. aws_smoke, f0, p1f
#
# Environment (all optional):
#   TYPE        instance type (default r7i.16xlarge: 64 vCPU, 512 GiB; aws_smoke defaults to r7i.4xlarge, 128 GiB)
#   SPOT=1      spot instead of on-demand (cheaper; AWS can reclaim it mid-run: finished cases are already in S3)
#   MAX_HOURS   hard stop: the instance powers off (and terminates) after this many hours (default: the batch's
#               "aws_max_hours" note, else 8 h per case, at most 48)
#   DISK_GB     root volume (default 500 GB gp3)
#   ADJ_REPO    local clone of mxnviir/stem-racing-adjoint (default ../stem-racing-adjoint)
#   STOP_ON_FAIL=0  run every case even after one fails a gate (default 1: a failed case or gate F1/F3 stops the batch)
#   DRY_RUN=1   package and print what would be sent, call no AWS API
# What is shipped: this repository at HEAD (uncommitted changes are NOT included: commit first), the pinned builder
# (EVO_SHA) and adjoint code (ADJ_SHA) as pinned in .github/workflows/carhybrid.yml, so the instance runs exactly what
# the GitHub workflow would.
set -euo pipefail
BATCH=${1:?usage: launch.sh <batch in carhybrid/request.json>}
ROOT=$(git rev-parse --show-toplevel); cd "$ROOT"
CONF=${CONF:-$HOME/.config/carhybrid/aws.env}
if [ "${DRY_RUN:-0}" != 1 ]; then
  [ -f "$CONF" ] || { echo "no $CONF: run carhybrid/aws/aws_setup.sh first"; exit 1; }
fi
[ -f "$CONF" ] && . "$CONF"
REGION=${REGION:-us-east-1}; BUCKET=${BUCKET:-dry-run-bucket}; PROFILE_NAME=${PROFILE_NAME:-carhybrid-dafoam}
export AWS_DEFAULT_REGION=$REGION
WF=.github/workflows/carhybrid.yml
EVO_SHA=$(sed -n 's/^ *EVO_SHA: *//p' $WF | head -1); ADJ_SHA=$(sed -n 's/^ *ADJ_SHA: *//p' $WF | head -1)
PY_DEPS=$(sed -n 's/^ *PY_DEPS: *//p' $WF | head -1)
ADJ_REPO=${ADJ_REPO:-$ROOT/../stem-racing-adjoint}
RUN_ID=${BATCH}_$(date -u +%Y%m%dT%H%M%SZ)
PREFIX=s3://$BUCKET/runs/$RUN_ID
say() { echo "[launch] $*"; }

# ---------------------------------------------------------------- the batch and its manifest
WORKD=$(mktemp -d); trap 'rm -rf "$WORKD"' EXIT
python3 - "$BATCH" "$RUN_ID" "$WORKD/manifest.json" "$EVO_SHA" "$ADJ_SHA" "$(git rev-parse HEAD)" "$PY_DEPS" \
         "${STOP_ON_FAIL:-1}" <<'PY'
import json, sys
batch, run_id, out, evo, adj, head, deps, stop = sys.argv[1:9]
R = json.load(open("carhybrid/request.json"))
if batch not in R["batches"]:
    sys.exit(f"batch {batch!r} is not in carhybrid/request.json (have: {', '.join(R['batches'])})")
cases = R["batches"][batch]
bad = [c.get("name") for c in cases if c.get("type") != "fine_adjoint"]
if bad:
    sys.exit(f"launch.sh runs fine_adjoint cases only; not fine_adjoint: {bad}")
names = [c["name"] for c in cases]
if len(set(names)) != len(names):
    sys.exit("duplicate case names")
json.dump({"run_id": run_id, "batch": batch, "cases": cases, "domain_m": R["domain_m"], "bg_div": R["bg_div"],
           "evo_sha": evo, "adj_sha": adj, "repo_sha": head, "py_deps": deps, "stop_on_fail": stop == "1",
           "note": R.get(f"_{batch}_note", ""), "aws_max_hours": R.get(f"_{batch}_aws_max_hours")},
          open(out, "w"), indent=1)
print(f"{len(cases)} cases: " + ", ".join(names))
PY
NCASES=$(python3 -c "import json; print(len(json.load(open('$WORKD/manifest.json'))['cases']))")
MH=$(python3 -c "import json; print(json.load(open('$WORKD/manifest.json')).get('aws_max_hours') or '')")
MAX_HOURS=${MAX_HOURS:-${MH:-$(( NCASES * 8 > 48 ? 48 : NCASES * 8 ))}}
DEFAULT_TYPE=r7i.16xlarge; [ "$BATCH" = aws_smoke ] && DEFAULT_TYPE=r7i.4xlarge
TYPE=${TYPE:-$DEFAULT_TYPE}

# ---------------------------------------------------------------- the code (three pinned trees in one archive)
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  say "WARNING: uncommitted changes in $ROOT are not shipped (the archive is HEAD $(git rev-parse --short HEAD))"
fi
git cat-file -e "$EVO_SHA^{commit}" 2>/dev/null || git fetch -q origin "$EVO_SHA" || git fetch -q origin claude/loving-keller-q72mvv
git cat-file -e "$EVO_SHA^{commit}" || { echo "builder commit $EVO_SHA not found (git fetch origin claude/loving-keller-q72mvv)"; exit 1; }
[ -d "$ADJ_REPO/.git" ] || { echo "no adjoint clone at $ADJ_REPO (set ADJ_REPO)"; exit 1; }
git -C "$ADJ_REPO" cat-file -e "$ADJ_SHA^{commit}" 2>/dev/null || git -C "$ADJ_REPO" fetch -q origin "$ADJ_SHA" \
  || git -C "$ADJ_REPO" fetch -q origin
git -C "$ADJ_REPO" cat-file -e "$ADJ_SHA^{commit}" || { echo "adjoint commit $ADJ_SHA not in $ADJ_REPO"; exit 1; }
mkdir -p "$WORKD/repo/evo" "$WORKD/repo/adj"
# results folders are not needed on the instance (it writes its own out/); everything else as the workflow checks out
git archive HEAD | tar -x -C "$WORKD/repo" --exclude="carhybrid/results" --exclude="audit/results"
git archive "$EVO_SHA" | tar -x -C "$WORKD/repo/evo"
git -C "$ADJ_REPO" archive "$ADJ_SHA" | tar -x -C "$WORKD/repo/adj"
tar -czf "$WORKD/code.tgz" -C "$WORKD" repo
say "code: repo $(git rev-parse --short HEAD), builder ${EVO_SHA:0:7}, adjoint ${ADJ_SHA:0:7}; $(du -h "$WORKD/code.tgz" | cut -f1)"

# ---------------------------------------------------------------- user data: fetch the code, start the runner
cat > "$WORKD/user_data.sh" <<EOF
#!/bin/bash
# carhybrid fine-mesh DAFoam instance, run $RUN_ID (carhybrid/aws/launch.sh)
exec > /var/log/carhybrid-boot.log 2>&1
set -x
shutdown -h +$(( MAX_HOURS * 60 )) "carhybrid MAX_HOURS=$MAX_HOURS reached"   # hard stop; the instance terminates
export BUCKET=$BUCKET RUN_ID=$RUN_ID AWS_DEFAULT_REGION=$REGION MAX_HOURS=$MAX_HOURS PATH=\$PATH:/snap/bin
for i in 1 2 3 4 5; do snap install aws-cli --classic && break; sleep 20; done
command -v aws || { apt-get update -q && apt-get install -y -q unzip curl && curl -sfL -o /tmp/awscli.zip \
  https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip && unzip -q /tmp/awscli.zip -d /tmp && /tmp/aws/install; }
mkdir -p /opt/ch && cd /opt/ch
if ! { aws s3 cp $PREFIX/code.tgz . && aws s3 cp $PREFIX/manifest.json . && tar -xzf code.tgz; }; then
  aws s3 cp /var/log/carhybrid-boot.log $PREFIX/boot.log; shutdown -h now; exit 1
fi
nohup bash /opt/ch/repo/carhybrid/aws/instance_run.sh > /opt/ch/instance.log 2>&1 &
EOF

AMI_PARAM=/aws/service/canonical/ubuntu/server/24.04/stable/current/amd64/hvm/ebs-gp3/ami-id
MARKET=(); [ "${SPOT:-0}" = 1 ] && MARKET=(--instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}')
NET=(--security-group-ids "${SG_ID:-sg-dryrun}")
[ -n "${SUBNET_ID:-}" ] && NET=(--network-interfaces "DeviceIndex=0,SubnetId=$SUBNET_ID,Groups=$SG_ID,AssociatePublicIpAddress=true")
KEY=(); [ -n "${KEY_NAME:-}" ] && KEY=(--key-name "$KEY_NAME")

say "batch $BATCH -> run $RUN_ID: $NCASES cases on $TYPE ($([ "${SPOT:-0}" = 1 ] && echo spot || echo on-demand)), ${DISK_GB:-500} GB," \
    "hard stop after $MAX_HOURS h; results to $PREFIX/out/"
if [ "${DRY_RUN:-0}" = 1 ]; then
  say "DRY_RUN: manifest"; cat "$WORKD/manifest.json"; say "DRY_RUN: user data"; cat "$WORKD/user_data.sh"
  say "DRY_RUN: archive contents (top)"; tar -tzf "$WORKD/code.tgz" | awk -F/ 'NF<=3' | head -30
  [ -n "${KEEP_DRY:-}" ] && cp "$WORKD/code.tgz" "$WORKD/manifest.json" "$WORKD/user_data.sh" "$KEEP_DRY/"
  exit 0
fi

aws s3 cp --only-show-errors "$WORKD/code.tgz" "$PREFIX/code.tgz"
aws s3 cp --only-show-errors "$WORKD/manifest.json" "$PREFIX/manifest.json"
AMI=$(aws ssm get-parameter --name "$AMI_PARAM" --query Parameter.Value --output text)
IID=$(aws ec2 run-instances --image-id "$AMI" --instance-type "$TYPE" --count 1 \
  --iam-instance-profile Name="$PROFILE_NAME" "${NET[@]}" "${KEY[@]}" "${MARKET[@]}" \
  --instance-initiated-shutdown-behavior terminate \
  --metadata-options HttpTokens=required,HttpEndpoint=enabled \
  --block-device-mappings "DeviceName=/dev/sda1,Ebs={VolumeSize=${DISK_GB:-500},VolumeType=gp3,DeleteOnTermination=true}" \
  --user-data "file://$WORKD/user_data.sh" \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=carhybrid-$RUN_ID},{Key=carhybrid-run,Value=$RUN_ID}]" \
  --query 'Instances[0].InstanceId' --output text)
echo "{\"run_id\": \"$RUN_ID\", \"instance_id\": \"$IID\", \"type\": \"$TYPE\", \"region\": \"$REGION\", \"ami\": \"$AMI\"}" \
  | aws s3 cp --only-show-errors - "$PREFIX/instance.json"
mkdir -p "$(dirname "$CONF")/runs"; echo "$IID" > "$(dirname "$CONF")/runs/$RUN_ID"
say "instance $IID started"
echo
echo "  status:   bash carhybrid/aws/status.sh $RUN_ID"
echo "  results:  bash carhybrid/aws/fetch.sh $RUN_ID        (into carhybrid/results/$BATCH/)"
echo "  stop now: aws ec2 terminate-instances --instance-ids $IID"
