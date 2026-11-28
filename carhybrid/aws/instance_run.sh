#!/bin/bash
# On the EC2 instance (started by launch.sh's user data, as root): set the machine up exactly as the GitHub workflow's
# case job does, run every case of the manifest with carhybrid/fine_adjoint.sh one after another, copy each case's
# results to S3 the moment it ends, and power off (= terminate) at the end.
#
# Layout: /opt/ch/manifest.json, /opt/ch/repo (this repository at the launched commit) with repo/evo (builder at
# EVO_SHA) and repo/adj (adjoint code at ADJ_SHA) inside it, as the workflow checks them out.
# S3 (s3://$BUCKET/runs/$RUN_ID/): instance.log and status.json (every minute), out/<case>/ (each case's outputs
# as the workflow's collect step keeps them: json, logs, npz, png, md), DONE when the batch has ended.
#   SKIP_SETUP=1  (tests) skip package installs; NO_POWEROFF=1 keep the machine up at the end; AWS_CMD (tests) the
#   aws binary.
set -uo pipefail
W=${W:-/opt/ch}; R=$W/repo; M=$W/manifest.json
PREFIX=s3://${BUCKET:?}/runs/${RUN_ID:?}
AWS=${AWS_CMD:-aws}
export PATH=$PATH:/snap/bin
say() { echo "[instance $(date -u +%FT%TZ)] $*"; }
status() {   # status.json: phase, current case, finished cases with their exit codes and gates
  python3 - "$W/status.json" "$1" "${2:-}" <<'PY'
import json, sys, time
p, phase, case = sys.argv[1:4]
try: s = json.load(open(p))
except Exception: s = {"cases": {}}
s.update({"phase": phase, "current": case or None, "updated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
json.dump(s, open(p, "w"), indent=1)
PY
}
push_status() { $AWS s3 cp --only-show-errors "$W/instance.log" "$PREFIX/instance.log" 2>/dev/null
                $AWS s3 cp --only-show-errors "$W/status.json" "$PREFIX/status.json" 2>/dev/null
                [ -f "$W/current_job.log" ] && $AWS s3 cp --only-show-errors "$W/current_job.log" "$PREFIX/current_job.log" 2>/dev/null; true; }
finish() {
  status "${1:-done}"; say "end: ${1:-done}"; push_status
  echo "${1:-done}" | $AWS s3 cp --only-show-errors - "$PREFIX/DONE"
  [ "${NO_POWEROFF:-0}" = 1 ] || { sync; shutdown -h now; }
  exit 0
}
status setup
( while :; do sleep 60; push_status; done ) & PUSHER=$!
trap 'kill $PUSHER 2>/dev/null' EXIT
say "run $RUN_ID on $(nproc) cores, $(awk '/MemTotal/{printf "%.0f", $2/1048576}' /proc/meminfo) GB RAM;" \
    "$(df -h / | awk 'NR==2{print $4}') disk free; hard stop after ${MAX_HOURS:-?} h"

# ---------------------------------------------------------------- 1. the machine (as the workflow's setup step)
cd "$R" || finish "failed: no code at $R"
PY_DEPS=$(python3 -c "import json; print(json.load(open('$M'))['py_deps'])")
if [ "${SKIP_SETUP:-0}" != 1 ]; then
  export DEBIAN_FRONTEND=noninteractive
  for i in 1 2 3 4 5 6; do timeout 300 apt-get update -q && break; sleep 20; done
  for i in 1 2 3 4; do apt-get install -y -q git curl ca-certificates docker.io libopengl0 libglu1-mesa libegl1 libgl1 \
       libgomp1 libxkbcommon0 && break; sleep 30; done
  systemctl enable --now docker
  # Python 3.11, as actions/setup-python gives the workflow (the pinned deps are tested there): uv's standalone build
  curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh || finish "failed: uv install"
  uv python install 3.11 && uv venv -q -p 3.11 "$W/venv" || finish "failed: python 3.11"
  . "$W/venv/bin/activate"
  # shellcheck disable=SC2086
  uv pip install -q $PY_DEPS || finish "failed: python deps"
  OF_DEBS=$W/ofdebs bash adj/.github/scripts/install_openfoam.sh || finish "failed: OpenFOAM install"
  for i in 1 2 3; do docker pull -q dafoam/opt-packages:v5.1.1 && break; sleep 30; done
  docker image inspect dafoam/opt-packages:v5.1.1 > /dev/null || finish "failed: DAFoam image"
else
  [ -f "$W/venv/bin/activate" ] && . "$W/venv/bin/activate"
fi
python seedhybrid/patch_adj.py adj || finish "failed: patch_adj"
python evo_patches/apply.py evo || finish "failed: evo_patches"
# OpenFOAM's bashrc reads unset variables: source it without -u (as the workflow, whose steps have no -u)
set +u; [ -f /usr/lib/openfoam/openfoam2412/etc/bashrc ] && . /usr/lib/openfoam/openfoam2412/etc/bashrc; set -u
export OMPI_MCA_rmaps_base_oversubscribe=1 OMPI_MCA_hwloc_base_use_hwthreads_as_cpus=1 \
       OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
RUNNER_NAME="ec2-$(cat /sys/devices/virtual/dmi/id/board_asset_tag 2>/dev/null || hostname)"; export RUNNER_NAME RUNNER_ENVIRONMENT=aws-ec2
say "setup done: $(python --version 2>&1), OpenFOAM $([ -f /usr/lib/openfoam/openfoam2412/etc/bashrc ] && echo v2412 || echo MISSING)," \
    "docker $(docker --version 2>/dev/null | cut -d, -f1)"
push_status

# ---------------------------------------------------------------- 2. the cases, one at a time
# each case's environment = the workflow's fine_adjoint step env (same defaults)
python3 - "$M" "$W/cases" <<'PY'
import json, shlex, sys
from pathlib import Path
m = json.load(open(sys.argv[1])); d = Path(sys.argv[2]); d.mkdir(exist_ok=True)
dom = json.dumps(m["domain_m"])
for i, c in enumerate(m["cases"]):
    g = lambda k, v="": "" if c.get(k) is None else str(c.get(k)).lower() if isinstance(c.get(k), bool) else str(c.get(k))
    u = c.get("u", "zero"); u = u if isinstance(u, str) else json.dumps(u)
    lay = c.get("layout", "team"); lay = lay if isinstance(lay, str) else json.dumps(lay)
    env = {"NAME": c["name"], "LAYOUT": lay, "U": u or "zero", "K": g("k") or "0",
           "BG_DIV": g("bg_div") or str(m["bg_div"]), "GRID": g("grid"), "NP": g("np"), "SCORE_NP": g("score_np") or "4",
           "GATE_DOMAIN_JSON": dom, "DAF_PRIMAL_ONLY": "1" if c.get("primal_only") else "",
           "DAF_FUNCS": g("funcs") or "D,L", "DAF_PRIMAL_ITERS": g("daf_primal_iters") or "300",
           "DAF_GMRES_TOL": g("gmres_tol") or "1e-6", "DAF_GMRES_ITERS": g("gmres_iters") or "6000",
           "DAF_GMRES_RESTART": g("gmres_restart") or "200", "DAF_PC_FILL": g("pc_fill") or "0",
           "DAF_PC_LEVEL": g("pc_level") or "1"}
    (d / f"{i:02d}_{c['name']}.env").write_text("".join(f"export {k}={shlex.quote(v)}\n" for k, v in env.items()))
PY
STOP_ON_FAIL=$(python3 -c "import json; print(1 if json.load(open('$M')).get('stop_on_fail', True) else 0)")
mkdir -p "$R/out" "$W/tmp"
for E in "$W"/cases/*.env; do
  ( . "$E"; echo "$NAME" ) > "$W/current_name"; NAME=$(cat "$W/current_name")
  status running "$NAME"; say "case $NAME: start"; push_status
  ( . "$E"; [ -z "$NP" ] && unset NP
    export CAR=$W/car_$NAME RUNNER_TEMP=$W/tmp ADJ=adj EVO=evo OUT=out/$NAME
    rm -rf "$CAR"; bash carhybrid/fine_adjoint.sh ) > "$W/current_job.log" 2>&1
  RC=$?
  mkdir -p "out/$NAME"; cp "$W/current_job.log" "out/$NAME/job.log"; rm -f "$W/current_job.log"
  # the results as the workflow's collect step keeps them (json, logs, npz up to 90 MB, renders), then to S3
  $AWS s3 sync --only-show-errors "out/$NAME" "$PREFIX/out/$NAME" --exclude "*" --include "*.json" --include "*.log" \
      --include "*.out" --include "*.md" --include "*.png" --include "*.npz"
  GATES=$(python3 - "out/$NAME/fine_adjoint.json" <<'PY'
import json, sys
try: f = json.load(open(sys.argv[1]))
except Exception: f = {}
print(json.dumps({k: f.get(k) for k in ("F1_pass", "F2_pass", "F3_pass", "F1_dafoam_vs_pipeline_drag", "F2_span_mN",
                                        "pipeline_drag_mN", "dafoam_D_mN", "cells")}
                 | {"peak_GB": (f.get("resources") or {}).get("peak_GB_whole_machine"),
                    "wall_s": (f.get("resources") or {}).get("job_wall_s")}))
PY
)
  python3 - "$W/status.json" "$NAME" "$RC" "$GATES" <<'PY'
import json, sys
p, n, rc, g = sys.argv[1:5]
s = json.load(open(p)); s.setdefault("cases", {})[n] = {"exit": int(rc), **json.loads(g)}
json.dump(s, open(p, "w"), indent=1)
PY
  say "case $NAME: exit $RC $GATES"; push_status
  rm -rf "$W/car_$NAME" "$W/tmp/fine_adjoint"
  # a gate known to be missed (False) or a failed job stops the batch: nothing downstream runs on a failed gate
  if [ "$STOP_ON_FAIL" = 1 ]; then
    if [ "$RC" != 0 ] || echo "$GATES" | grep -qE '"F[13]_pass": false'; then finish "stopped after $NAME (exit $RC)"; fi
  fi
done
finish "done"
