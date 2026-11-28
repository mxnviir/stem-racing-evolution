"""Build one whole car with the pinned evolution builder and adapt it for OpenFOAM (carhybrid/DESIGN.md P0).

    python carhybrid/make_car.py <evolution checkout> <layout: "team" | genes json | Lxxxx (search/layout/results)>
                                 <stage-2 numbers u: "zero" | json list> <out folder>

Writes <out>/evo/ (the evolution car folder: summary.json, layout_car.json, body_half.stl, parts/, body_full.stl,
wheels), <out>/rec.json (the build record: genes, decoded layout incl. the free-form seed, legal, mass, cfd report) and
<out>/cfd/ (the OpenFOAM half-car folder, carhybrid/adapt.py). Exit code 0 only if the car is legal and adapted.
The stage-2 numbers move the free-form body and nose exactly as stage 2 does (simscale_shape.ShapeCampaign.job:
freeform.from_unit(seed, u[:45], AMP_MM, *COARSE)); the seed is the stage-1 car's own organic free-form surfaces.
"""
import json, subprocess, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def genes_of(evo, layout):
    import layout_genes as LG
    if layout == "team":
        return LG.x_default()
    if layout.strip().startswith("["):
        return np.asarray(json.loads(layout), float)
    for f in (evo / "search" / "layout" / "results" / "designs.jsonl",):
        for line in f.read_text().splitlines():
            r = json.loads(line)
            if r.get("id") == layout and r.get("genes") is not None:
                return np.asarray(r["genes"], float)
    raise SystemExit(f"layout {layout} not found")


def main(evo, layout, u_s, out):
    evo, out = Path(evo).resolve(), Path(out).resolve()
    sys.path.insert(0, str(evo)); import repo_paths; repo_paths.add()
    import layout_car as LC, layout_genes as LG, freeform as fr, simscale_shape as SS
    x = genes_of(evo, layout)
    d = LG.decode(x)
    freeform = None
    if u_s != "zero":
        u = np.asarray(json.loads(u_s), float)
        seed = fr.Freeform.from_dict(d["freeform"])
        n = fr.n_genes(*SS.COARSE)
        freeform = fr.from_unit(seed, u[:n], SS.AMP_MM, *SS.COARSE).as_dict()
    car = out / "evo"
    rec = LC.build(x, car, freeform=freeform)
    rec["u"] = None if u_s == "zero" else json.loads(u_s)
    rec["layout_name"] = layout if not layout.strip().startswith("[") else "genes"
    if rec.get("legal"):
        try:
            rec.update(LC.export_cfd(car, rec["decoded"]))
        except Exception as e:                       # noqa: BLE001  (refused: not one closed piece, etc.)
            rec["cfd_error"] = repr(e)[-400:]
    (out / "rec.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps({k: rec.get(k) for k in ("legal", "failed", "why", "mass", "build_s")}, default=float)[:1500])
    if not rec.get("legal") or "cfd" not in rec:
        return 2
    r = subprocess.run([sys.executable, str(HERE / "adapt.py"), str(car), str(out / "cfd")])
    return r.returncode


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
