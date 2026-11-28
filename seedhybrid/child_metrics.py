"""Per-child geometry metrics for the analysis (PROTOCOL.md E1), run inside the adjoint checkout's research/sidepod_opt
with a lineage's base surfaces at gate/original(_rm) and STEP_REGION=1.

    python <this file> <in.json> <out.json>
in.json: {child: {"field": ..., "amp": ..., "grad": <npz path rel. to sidepod_opt or null>}}
out.json: {child: {"rms_mm": RMS normal move over the moved region, "peak_mm": ..., "pred_mN": first-order drag change
           predicted by the gradient (dafoam/judge.py's mapping, as esi/predict_candidates.py does) or null}}

The move is rebuilt with esi/predict_candidates.move (make_step's own construction: smoothing, region taper,
volume-neutral, 95th-percentile normalisation, clip), so it is deterministic from the stored field and amplitude.
"""
import json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "esi"); sys.path.insert(0, ".")
import predict_candidates as PC          # noqa: E402
from dafoam.judge import displacement_at  # noqa: E402

spec = json.loads(Path(sys.argv[1]).read_text())
out, grads = {}, {}
moved = PC.taper > 0.5
for name, c in spec.items():
    w = PC.move(c["field"], c["amp"])
    rms = float(np.sqrt((PC.A[moved] * w[moved] ** 2).sum() / PC.A[moved].sum()) * 1e3)
    pred = None
    if c.get("grad"):
        if c["grad"] not in grads:
            grads[c["grad"]] = np.load(c["grad"])
        S = grads[c["grad"]]
        pred = float((S["dD_dXs"] * displacement_at(S["Xs"], PC.V, PC.N * w[:, None])).sum()) * 1e3
    out[name] = {"rms_mm": rms, "peak_mm": float(np.abs(w).max() * 1e3), "pred_mN": pred}
Path(sys.argv[2]).write_text(json.dumps(out, indent=1))
print(f"{len(out)} children measured", flush=True)
