#!/bin/bash
# The fine-mesh adjoint job body (case type "fine_adjoint", carhybrid/FINE_ADJOINT.md): build -> the audit-protocol
# score of frame k on the production mesh, its converged state kept -> steadiness (gate F2) -> DAFoam on that same
# state (mesh, schemes, BCs; dafoam/make_case.py) with whole-machine memory sampled every 10 s -> the gradient mapped
# to the free-form numbers (carhybrid/grad.py). Every key number and the tail of DAFoam's output go to stdout (the
# job log is what we can read; artifacts are not).
#
#   bash carhybrid/fine_adjoint.sh          (from the repository root, inputs from the environment)
#
# Inputs: NAME, LAYOUT (team | genes json | Lxxxx), U ("zero" | json list), K (frame, default 0), BG_DIV (default
#   12.889), GRID (optional, passed to make_car.py and to grad.py if it takes --grid), SMOOTH (grad.py --smooth mm;
#   default 10 on any grid but 5,2,3,4), NP (DAFoam MPI ranks, default all cores),
#   SCORE_NP (score's ranks, default 4 = the audit protocol, as every hosted score), GATE_DOMAIN_JSON (meshing frame),
#   ADJ (adjoint checkout, default adj), EVO (evolution checkout, default evo), CAR (built car, default /tmp/car; built
#   here if CAR/rec.json is missing), OUT (default out/$NAME), WORK (scratch, default $RUNNER_TEMP/fine_adjoint),
#   DAF_PRIMAL_ONLY=1 (DAFoam primal only, no adjoint), DAF_FUNCS (default D,L), DAF_PRIMAL_ITERS (default 300),
#   DAF_GMRES_TOL / DAF_GMRES_ITERS / DAF_GMRES_RESTART / DAF_PC_FILL / DAF_PC_LEVEL (adjoint settings).
# Outputs in OUT: result.json (score), steady.json (F2), dafoam_sens.npz/.json (or dafoam_primal.json), dafoam.log
#   (DAFoam's output, condensed), resources.json (cells, ranks, peak memory, wall times), grad.json, fine_adjoint.json
#   (the gates F1 F2 F3 and the key numbers in one place).
# Exit: 0 when the requested product exists (grad.json; with DAF_PRIMAL_ONLY=1 dafoam_primal.json), else 1.
set -uo pipefail
ROOT=$(pwd)
NAME=${NAME:?NAME}; LAYOUT=${LAYOUT:-team}; U=${U:-zero}; K=${K:-0}; BG_DIV=${BG_DIV:-12.889}; GRID=${GRID:-}
ADJ=$(realpath "${ADJ:-adj}"); EVO=$(realpath "${EVO:-evo}"); CAR=${CAR:-/tmp/car}
OUT=$(mkdir -p "${OUT:-out/$NAME}" && realpath "${OUT:-out/$NAME}")
WORK=${WORK:-${RUNNER_TEMP:-/tmp}/fine_adjoint}
NP=${NP:-$(nproc)}; SCORE_NP=${SCORE_NP:-4}
IMAGE=${DAFOAM_IMAGE:-dafoam/opt-packages:v5.1.1}
PY=${PY:-python}
rm -rf "$WORK"; mkdir -p "$WORK"
T0=$(date +%s)
say() { echo "[fine_adjoint $(date -u +%T)] $*"; }
phase() { echo "$1 $(date +%s)" >> "$WORK/phases.txt"; say "phase $1"; }

# whole-machine memory in use (RAM + swap, MB) every 10 s, for the whole job; peaks per phase at the end
( while :; do awk -v t="$(date +%s)" '/MemTotal/{m=$2}/MemAvailable/{a=$2}/SwapTotal/{s=$2}/SwapFree/{f=$2}
    END{printf "%d %d\n", t, (m-a+s-f)/1024}' /proc/meminfo >> "$WORK/mem.tsv"; sleep 10; done ) &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null' EXIT

say "case $NAME: layout ${LAYOUT:0:60} u ${U:0:60} k $K bg_div $BG_DIV grid ${GRID:-default} ranks $NP (score $SCORE_NP)" \
    "host $(nproc) cores, $(awk '/MemTotal/{printf "%.0f", $2/1048576}' /proc/meminfo) GB RAM," \
    "$(awk '/SwapTotal/{printf "%.0f", $2/1048576}' /proc/meminfo) GB swap, runner ${RUNNER_NAME:-?} (${RUNNER_ENVIRONMENT:-?})"

# ---------------------------------------------------------------- 0. the car (the workflow's build step normally did it)
phase build
if [ ! -f "$CAR/rec.json" ] || [ ! -d "$CAR/cfd" ]; then
  rm -rf "$CAR"
  $PY -W ignore carhybrid/make_car.py "$EVO" "$LAYOUT" "$U" "$CAR" ${GRID:+--grid "$GRID"} > "$OUT/make_car.log" 2>&1 || { tail -c 2500 "$OUT/make_car.log"; say "BUILD FAILED"; }
fi
[ -d "$CAR/cfd" ] || { say "no car folder $CAR/cfd: stop"; exit 1; }

# ---------------------------------------------------------------- 1. the score (audit protocol), state kept
phase score
$PY - "$WORK/case.json" <<PY
import json, sys
json.dump({"name": "$NAME", "design": "folder:$CAR/cfd", "k": int("$K"), "bg_div": float("$BG_DIV"),
           "np": int("$SCORE_NP"), "keep_state": "$WORK/state", "timeout_s": int("${SCORE_TIMEOUT_S:-60000}")},
          open(sys.argv[1], "w"))
PY
for TRY in 1 2 3; do
  $PY -W ignore audit/run_case.py "$ADJ" "$WORK/case.json" "$OUT" | tail -1
  # decomposePar of a coarse layered mesh dies with SIGFPE in boundBox now and then, on a decomposition that changes
  # from run to run (seen at bg_div 4 on the same car: passed, then failed); it fails in the first minutes, so retry
  if grep -q "decomposePar -force' died with <Signals.SIGFPE" "$OUT/prodcheck.log" 2>/dev/null && [ $TRY -lt 3 ]; then
    say "decomposePar SIGFPE (try $TRY): meshing again"; continue
  fi
  break
done
echo "=== result"; head -c 1500 "$OUT/result.json"; echo
if ! $PY -c "import json,sys; r=json.load(open('$OUT/result.json')); sys.exit(0 if r.get('converged') and 'path' in (r.get('kept_state') or {}) else 1)"; then
  echo "=== prodcheck tail"; tail -c 3000 "$OUT/prodcheck.log" 2>/dev/null
  # the pipeline keeps a failed run folder: print what its meshing and decomposition said (artifacts are not readable)
  RUN=$(grep -o "run dir KEPT for diagnosis: [^ ]*" "$OUT/prodcheck.log" 2>/dev/null | tail -1 | awk '{print $NF}')
  if [ -n "$RUN" ] && [ -d "$RUN/logs" ]; then
    mkdir -p "$OUT/failed_logs"; cp "$RUN"/logs/*.log "$OUT/failed_logs/" 2>/dev/null
    echo "=== failed run $RUN: snappy layers"; grep -E "Extruding|Layer mesh|cells:" "$RUN/logs/snappyHexMesh.log" 2>/dev/null | tail -8
    grep -A8 "patch  faces        layers" "$RUN/logs/snappyHexMesh.log" 2>/dev/null | tail -9
    echo "=== snappy tail"; tail -25 "$RUN/logs/snappyHexMesh.log" 2>/dev/null
    echo "=== checkMesh"; grep -E "cells:|points:|faces:|\*\*\*|Mesh OK|Failed|^ +[A-Za-z]+ +[0-9]+ +[0-9]+ " "$RUN/logs/checkMesh.log" 2>/dev/null | head -30
    echo "=== decomposePar tail"; tail -40 "$RUN/logs/decomposePar.log" 2>/dev/null
  fi
  say "SCORE FAILED or no kept state: stop"; STOP=1
fi

# ---------------------------------------------------------------- 2. steadiness (gate F2 = mark A1)
phase steady
[ -z "${STOP:-}" ] && $PY - "$OUT" <<'PY'
import json, sys
from pathlib import Path
o = Path(sys.argv[1]); r = json.loads((o / "result.json").read_text())
d, l = r.get("steady_drag") or {}, r.get("steady_lift") or {}
s = {"gate": "F2 (mark A1): 100-iteration block means span <= 0.10 mN over the last 500 iterations (half car)",
     "drag": d, "lift": l, "span_mN": d.get("span_mN"), "lift_span_mN": l.get("span_mN"),
     "pass": d.get("span_mN") is not None and d["span_mN"] <= 0.10}
(o / "steady.json").write_text(json.dumps(s, indent=1))
print("=== F2 steadiness: drag block means", [round(x, 3) for x in d.get("block_means_mN", [])], "span",
      round(d.get("span_mN") or float("nan"), 3), "mN ->", "PASS" if s["pass"] else "FAIL",
      "| lift span", round(l.get("span_mN") or float("nan"), 3), "mN")
PY

# ---------------------------------------------------------------- 3. DAFoam on the kept state
phase make_case
if [ -z "${STOP:-}" ]; then
  $PY carhybrid/patch_daf.py "$ADJ" || { say "patch_daf failed"; STOP=1; }
fi
if [ -z "${STOP:-}" ]; then
  ( cd "$ADJ/research/sidepod_opt" && $PY dafoam/make_case.py --src "$WORK/state" --out "$WORK/daf" --iters "${DAF_PRIMAL_ITERS:-300}" ) \
    > "$WORK/make_case.log" 2>&1 || { cat "$WORK/make_case.log"; say "make_case failed"; STOP=1; }
  head -1 "$WORK/make_case.log"
  rm -rf "$WORK/state"                                    # the case holds its own copy
fi
phase dafoam
if [ -z "${STOP:-}" ]; then
  say "DAFoam: $NP ranks, functions ${DAF_FUNCS:-D,L}, primal only ${DAF_PRIMAL_ONLY:-0}"
  env DAF_NP="$NP" DAF_FUNCS="${DAF_FUNCS:-D,L}" DAF_ADJ="${DAF_ADJ:-Krylov}" DAF_PC_FILL="${DAF_PC_FILL:-0}" \
      DAF_PC_LEVEL="${DAF_PC_LEVEL:-1}" DAF_GMRES_ITERS="${DAF_GMRES_ITERS:-2000}" DAF_GMRES_RESTART="${DAF_GMRES_RESTART:-200}" \
      DAF_GMRES_TOL="${DAF_GMRES_TOL:-1e-6}" DAF_PRIMAL_ONLY="${DAF_PRIMAL_ONLY:-}" \
      "$ADJ/research/sidepod_opt/dafoam/run_docker.sh" "$WORK/daf" "$NP" "$IMAGE" > "$WORK/dafoam.out" 2>&1
  DRC=$?
  say "DAFoam exit code $DRC"
  # the container writes as its own user: make the case removable by the runner user (self-hosted workspaces persist)
  docker run --rm -u root -v "$WORK/daf:/c" --entrypoint chmod "$IMAGE" -R a+rwX /c 2>/dev/null || true
  F="$WORK/dafoam.out"
  { echo "=== head"; head -40 "$F"
    echo "=== key lines"; grep -nE "DAFOAM|Killed|killed|Error|ERROR|rror:|nCells|Global Adjoint States|AllNonZeros|Ncolors|Calculating dRdW Coloring\.\.\. Completed|Main iteration|Completed\*\*|PetscConvergedReason|Residual tolerance|signal|Maximum resident|Total Residual Norm2|ExecutionTime = " "$F" | awk 'NR<=4000'
    echo "=== tail"; tail -150 "$F"; } > "$OUT/dafoam.log"
  echo "=== DAFoam key lines"
  grep -E "DAFOAM|Killed|killed|FATAL|Error|Global Adjoint States|Ncolors|Coloring\.\.\. Completed|Completed\*\*|PetscConvergedReason|Total Residual Norm2" "$F" | tail -40
  echo "=== DAFoam GMRES (every 500th line)"; grep -E "Main iteration" "$F" | awk 'NR%5==1' | tail -40
  echo "=== DAFoam tail"; tail -40 "$F"
  cp "$WORK/daf/dafoam_sens.npz" "$WORK/daf/dafoam_sens.json" "$WORK/daf/dafoam_primal.json" "$OUT/" 2>/dev/null
fi
phase grad

# ---------------------------------------------------------------- 4. GMRES from the log (always), gates F1 F3, sizes
$PY - "$OUT" "$WORK" "${DAF_FUNCS:-D,L}" <<'PY'
import json, re, sys
from pathlib import Path
o, w, funcs = Path(sys.argv[1]), Path(sys.argv[2]), [f for f in sys.argv[3].split(",") if f]
log = (w / "dafoam.out").read_text(errors="replace") if (w / "dafoam.out").exists() else ""
# one block per adjoint solve: "Main iteration N KSP Residual norm R t s." lines closed by "**Completed**! ..."
solves, cur = [], []
for line in log.splitlines():
    m = re.match(r"\s*Main iteration\s+(\d+)\s+KSP Residual norm\s+([0-9.eE+-]+)", line)
    if m:
        cur.append((int(m.group(1)), float(m.group(2))))
    m = re.search(r"\*\*Completed\*\*! Total iterations: (\d+)\. PetscConvergedReason: (-?\d+)", line)
    if m:
        solves.append({"iterations": int(m.group(1)), "reason": int(m.group(2)), "converged": int(m.group(2)) > 0,
                       "residual_first": cur[0][1] if cur else None, "residual_last": cur[-1][1] if cur else None,
                       "relative_drop": (cur[-1][1] / cur[0][1]) if cur and cur[0][1] else None})
        cur = []
gl = {f: s for f, s in zip(funcs, solves)}
sj = o / "dafoam_sens.json"
if sj.exists():
    d = json.loads(sj.read_text()); d["gmres_log"] = gl; sj.write_text(json.dumps(d, indent=1))
m = re.search(r"Global Adjoint States:\s*(\d+)", log)
res = {"adjoint_states": int(m.group(1)) if m else None, "gmres_log": gl}
json.dump(res, open(w / "daf_log.json", "w"))
PY

if [ -z "${STOP:-}" ] && [ -f "$OUT/dafoam_sens.npz" ]; then
  read D20 L20 < <($PY -c "import json; r=json.load(open('$OUT/result.json')); print(2*r['drag_mN']/1e3, 2*(r.get('lift_mN') or 0)/1e3)")
  GA=(--body "$CAR/evo" --d20 "$D20")
  HELP=$($PY carhybrid/grad.py --help 2>&1 | head -c 4000 || true)
  if echo "$HELP" | grep -q -- "--l20" || grep -q -- '"--l20"' carhybrid/grad.py; then GA+=(--l20 "$L20"); fi
  if [ -n "$GRID" ] && { echo "$HELP" | grep -q -- "--grid" || grep -q -- '"--grid"' carhybrid/grad.py; }; then GA+=(--grid "$GRID"); fi
  # the 720-number grid needs smoothing of the surface sensitivities (grad.py: --smooth 10 mm or more there)
  if [ -z "${SMOOTH:-}" ] && [ -n "$GRID" ] && [ "$GRID" != "5,2,3,4" ]; then SMOOTH=10; fi
  if [ -n "${SMOOTH:-}" ] && echo "$HELP" | grep -q -- "--smooth"; then GA+=(--smooth "$SMOOTH"); fi
  say "grad.py ${GA[*]}"
  $PY -W ignore carhybrid/grad.py "$EVO" "$CAR/rec.json" "$OUT/dafoam_sens.npz" "$U" "$OUT/grad.json" "${GA[@]}" 2>&1 | tail -20 \
    || say "grad mapping failed"
  echo "=== grad.json"; head -c 2000 "$OUT/grad.json" 2>/dev/null; echo
fi
phase end

# ---------------------------------------------------------------- 5. resources.json and the one-place summary
$PY - "$OUT" "$WORK" "$NP" "$SCORE_NP" "$T0" "${DAF_PRIMAL_ONLY:-0}" <<'PY'
import json, os, sys
from pathlib import Path
o, w, npr, snp, t0, primal_only = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6] == "1"
def J(p):
    try: return json.loads(Path(p).read_text())
    except Exception: return {}
r, st, sens, prim, lg, gr = J(o / "result.json"), J(o / "steady.json"), J(o / "dafoam_sens.json"), J(o / "dafoam_primal.json"), J(w / "daf_log.json"), J(o / "grad.json")
mem = [tuple(map(int, l.split())) for l in (w / "mem.tsv").read_text().split("\n") if len(l.split()) == 2] if (w / "mem.tsv").exists() else []
ph = [l.split() for l in (w / "phases.txt").read_text().split("\n") if l.strip()] if (w / "phases.txt").exists() else []
ph = [(n, int(t)) for n, t in ph]
spans = {n: (t, ph[i + 1][1] if i + 1 < len(ph) else 10**12) for i, (n, t) in enumerate(ph)}
def peak(a, b):
    v = [m for t, m in mem if a <= t <= b + 10]
    return round(max(v) / 1024, 2) if v else None
cells = r.get("cells")
meminfo = {l.split(":")[0]: int(l.split()[1]) for l in open("/proc/meminfo") if l.split()[0].rstrip(":") in ("MemTotal", "SwapTotal")}
pk = peak(0, 10**12)
res = {"cells": cells, "ranks": npr, "score_ranks": snp, "adjoint_states": lg.get("adjoint_states"),
       "peak_GB_whole_machine": pk, "peak_GB_score": peak(*spans["score"]) if "score" in spans else None,
       "peak_GB_dafoam": peak(*spans["dafoam"]) if "dafoam" in spans else None,
       "baseline_GB": round(mem[0][1] / 1024, 2) if mem else None,
       "kB_per_cell_dafoam_peak": round((peak(*spans["dafoam"]) or 0) * 1048576 / cells, 1) if cells and "dafoam" in spans and peak(*spans["dafoam"]) else None,
       "estimate_GB_at_150kB_per_cell": round(cells * 150e3 / 1e9, 1) if cells else None,
       "wall_s": {n: (b - a if b < 10**12 else None) for n, (a, b) in spans.items() if n != "end"},
       "score_wall_s": r.get("wall_s"), "dafoam_primal_s": (sens or prim).get("primal_s"), "dafoam_adjoint_s": sens.get("adjoint_s"),
       "dafoam_adjoint_s_per_function": sens.get("adjoint_s_per_function"), "job_wall_s": int(ph[-1][1]) - t0 if ph else None,
       "host_cores": os.cpu_count(), "host_mem_GB": round(meminfo.get("MemTotal", 0) / 1048576, 1),
       "host_swap_GB": round(meminfo.get("SwapTotal", 0) / 1048576, 1),
       "runner": os.environ.get("RUNNER_NAME"), "runner_environment": os.environ.get("RUNNER_ENVIRONMENT"),
       "memory_samples": len(mem), "sample_interval_s": 10}
(o / "resources.json").write_text(json.dumps(res, indent=1))
d_pipe = r.get("drag_mN"); l_pipe = r.get("lift_mN")
src = sens or prim
d_daf = src.get("D_half_N"); l_daf = src.get("L_half_N")
f1 = (d_daf * 1e3 / d_pipe - 1) if d_daf and d_pipe else None
gm = sens.get("gmres") or {}
f3 = {f: bool(v.get("converged")) for f, v in gm.items()} if gm else {f: v.get("converged") for f, v in (lg.get("gmres_log") or {}).items()}
summ = {"name": (r.get("case") or {}).get("name"), "cells": cells, "y_plus_max": r.get("y_plus_max"),
        "pipeline_drag_mN": d_pipe, "pipeline_drag_se_mN": r.get("drag_block_se_mN"), "pipeline_lift_mN": l_pipe,
        "dafoam_D_mN": d_daf * 1e3 if d_daf else None, "dafoam_L_mN": l_daf * 1e3 if l_daf else None,
        "F1_dafoam_vs_pipeline_drag": f1, "F1_pass": abs(f1) <= 0.02 if f1 is not None else None,
        "lift_dafoam_vs_pipeline": (l_daf * 1e3 / l_pipe - 1) if l_daf and l_pipe else None,
        "F2_span_mN": st.get("span_mN"), "F2_pass": st.get("pass"),
        "F3_gmres_converged": f3, "F3_pass": (all(f3.values()) if f3 else None),
        "gmres": gm or lg.get("gmres_log"), "gradients": sens.get("gradients"), "primal_only": primal_only,
        "grad_json": bool(gr), "norm_dT_du_s": gr.get("norm_dT_du_s"), "resources": res}
(o / "fine_adjoint.json").write_text(json.dumps(summ, indent=1, default=str))
print("=== fine_adjoint summary")
for k, v in summ.items():
    if k != "resources":
        print(f"  {k}: {json.dumps(v, default=str)[:600]}")
print("=== resources"); print(json.dumps(res, indent=1)[:3000])
PY

if [ "${DAF_PRIMAL_ONLY:-0}" = 1 ]; then [ -f "$OUT/dafoam_primal.json" ]; else [ -f "$OUT/grad.json" ]; fi
