#!/bin/bash
# The fine adjoint job inside the container (AWS Batch, or `docker run` anywhere).
#   in:  S3_IN  = s3://.../in.tgz   (batch.py package: job.json + evolution/ + adjoint/ code snapshots [+ car/])
#        or IN_TGZ = a local path (docker run with a bind mount; the GitHub end-to-end test)
#   out: S3_OUT = s3://.../<job id>/  (out.tgz, result.json, and progress/ every PROGRESS_S seconds)
#        or OUT_DIR = a local folder
#   NP: MPI ranks (default: every CPU the container sees); KEEP_STATE=1 also returns the steady flow (warm starts)
set -uo pipefail
NP=${NP:-$(nproc)}
W=/scratch/job
mkdir -p $W/in $W/run
say() { echo "[entrypoint $(date -u +%H:%M:%S)] $*"; }

if [ -n "${S3_IN:-}" ]; then
  say "fetch $S3_IN"; python /usr/local/bin/s3io.py get "$S3_IN" $W/in.tgz || exit 2
else
  cp "${IN_TGZ:?S3_IN or IN_TGZ}" $W/in.tgz || exit 2
fi
tar xzf $W/in.tgz -C $W/in || exit 2
say "job $(python -c "import json;print(json.load(open('$W/in/job.json')).get('id'))"), NP $NP, $(nproc) CPUs, $(free -g | awk '/Mem/{print $2}') GB"

pack() {   # the results (small) and, with KEEP_STATE=1, the steady flow
  cd $W/run
  files=$(ls -d result.json fine_adjoint.log gradient.json production.json dafoam.json steady/steady.json \
          dafoam_case/dafoam_sens.npz dafoam_case/dafoam_sens.json dafoam_case/log.* \
          steady/steady_*/log.* production/prodcheck.log car/layout_car.json car/run_car.log 2>/dev/null)
  tar czf $W/out.tgz $files
  if [ "${KEEP_STATE:-0}" = 1 ]; then
    st=$(python -c "import json;print(json.load(open('steady/steady.json'))['state'])" 2>/dev/null) \
      && tar czf $W/state.tgz -C "$(dirname "$st")" "$(basename "$st")"
  fi
}
progress() {   # partial results while the job runs (a long job can be watched and a lost one diagnosed)
  while sleep "${PROGRESS_S:-600}"; do
    pack >/dev/null 2>&1
    [ -n "${S3_OUT:-}" ] && python /usr/local/bin/s3io.py put $W/out.tgz "${S3_OUT%/}/progress/out.tgz" 2>/dev/null
  done
}
progress & PROG=$!

cd $W/in/evolution
python -W ignore search/adjoint/fine_adjoint.py $W/in/job.json $W/run --adjoint-root $W/in/adjoint --np "$NP" 2>&1 | tee $W/run/fine_adjoint.log
kill $PROG 2>/dev/null
pack
if [ -n "${S3_OUT:-}" ]; then
  python /usr/local/bin/s3io.py put $W/out.tgz "${S3_OUT%/}/out.tgz"
  python /usr/local/bin/s3io.py put $W/run/result.json "${S3_OUT%/}/result.json"
  python /usr/local/bin/s3io.py put $W/run/fine_adjoint.log "${S3_OUT%/}/fine_adjoint.log"
  [ -f $W/state.tgz ] && python /usr/local/bin/s3io.py put $W/state.tgz "${S3_OUT%/}/state.tgz"
else
  mkdir -p "${OUT_DIR:-/out}" && cp $W/out.tgz $W/run/result.json $W/run/fine_adjoint.log "${OUT_DIR:-/out}/" && \
    { [ ! -f $W/state.tgz ] || cp $W/state.tgz "${OUT_DIR:-/out}/"; }
fi
python -c "import json,sys; r=json.load(open('$W/run/result.json')); print('USABLE' if r.get('usable') else 'NOT USABLE', r.get('error',''))"
python -c "import json,sys; sys.exit(0 if json.load(open('$W/run/result.json')).get('usable') else 3)"
