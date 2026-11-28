"""Legalise a morphed design (PROTOCOL.md E4): a half-car surface the loop produced -> the same legal builder every
seed went through (rules, 3-axis 6.35 mm ball-end mill, 3 mm walls, one piece, cargo inside, ground clearance;
everything outside the sidepod region reset to the team's CAD) -> a legal car folder.

    python seedhybrid/legalise.py <adjoint checkout> <car_half.stl (car frame, m, y >= 0)> <out folder> [label]

Writes <out folder>/car/ (Part 2 car folder) and <out folder>/legalise.json: whether the builder accepted it, the
volume the rules changed inside the region, the volume differing outside the region before the reset, pod volume, mass.
"""
import json, sys
from pathlib import Path

import numpy as np


def main(adj, stl, out, label="legal"):
    sp = Path(adj).resolve() / "research" / "sidepod_opt"
    stl, out = str(Path(stl).resolve()), Path(out).resolve()
    sys.path.insert(0, str(sp)); sys.path.insert(0, str(sp.parents[1] / "part2-simulation"))
    import os
    os.chdir(sp)
    import trimesh
    import manifold3d as mf
    import cadio as CI
    import shape_loop as SL
    from seedevo import ingest
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    car = SL.load_car()
    m = trimesh.load(stl, process=True)
    # car frame (m, x from the nose, mirrored) -> STEP frame (mm); the inverse of geometry.Car.to_car_frame
    T = np.diag([-1e3, 1e3, 1e3, 1.0]); T[0, 3] = car.s.x_nose_mm
    m.apply_transform(T)
    if m.volume < 0:
        m.invert()
    other = m.copy(); other.apply_transform(np.diag([1.0, -1.0, 1.0, 1.0]))
    if other.volume < 0:
        other.invert()
    full = CI.from_manifold(CI.to_manifold(m) + CI.to_manifold(other))       # the whole body, watertight
    f, info = ingest.field_from_mesh(car, full)
    b = SL.build_design(car, f)
    rec = {"source": str(stl), "accepted": bool(b["check"]["ok"]), "check": b["check"],
           "rules_changed_in_region_mm3": float(b["rules_changed_mm3"]),
           "outside_region_differs_mm3_before_reset": info["outside_region_changed_mm3"],
           "pod_cm3": float((b["voxels"] & car.region).sum() * car.h ** 3 / 1e3), "mass": b.get("mass"),
           "full_body_volume_cm3": float(full.volume / 1e3)}
    if rec["accepted"]:
        SL.export_design(car, b, out, label)          # writes <out>/car/
    (out / "legalise.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps({k: rec[k] for k in ("accepted", "rules_changed_in_region_mm3", "outside_region_differs_mm3_before_reset",
                                          "pod_cm3")}, default=float))
    return 0 if rec["accepted"] else 2


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
