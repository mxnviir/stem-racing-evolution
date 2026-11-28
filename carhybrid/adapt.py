"""An evolution car folder's CFD export -> an OpenFOAM half-car folder for stem-racing-adjoint's pipeline.

    python carhybrid/adapt.py <evolution car folder> <out folder>

Input (written by search/layout/layout_car.export_cfd): body_full.stl (the whole car as ONE closed solid with its
supports, wings, nose, halo and cartridge; metres, car frame: x from the nose tip, y left, z up), wheelF.stl /
wheelR.stl (the +y wheels) and wheels.json (their centres).
Output (Part 2's car folder): car_half.stl (solid "car": the y >= 0 half, closed exactly at y = 0 by a Boolean),
wheelF.stl, wheelR.stl, assembly.json (rotating wheels: axis +y, omega = -U/r), adapt.json (bounds, volumes, checks).
"""
import json, sys
from pathlib import Path

import numpy as np

U_MPS = 20.0
SINK_M = 0.3e-3                 # Part 2's wheel sink (stem-racing-adjoint geometry.WHEEL_SINK_MM)


def main(src, out):
    import trimesh
    import manifold3d as mf
    src, out = Path(src), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    m = trimesh.load(src / "body_full.stl", force="mesh")
    m.apply_scale(1e3)                                             # mm for the Boolean (manifold is float32)
    M = mf.Manifold(mf.Mesh(vert_properties=np.asarray(m.vertices, np.float32), tri_verts=np.asarray(m.faces, np.uint32)))
    half = M ^ mf.Manifold.cube([2e3, 1e3, 1e3]).translate([-1e3, 0.0, -5e2])
    hm = half.to_mesh()
    # manifold's output is closed and 2-manifold by construction: keep its exact topology (process=False). trimesh's
    # own vertex merge (process=True, tolerance 1e-8 in metres, near float32 spacing at 0.27 m) could fuse distinct
    # vertices and open edges: P1 run 46132124252 lost three legal cars to "half_watertight": false that way.
    V, F = np.asarray(hm.vert_properties)[:, :3].astype(np.float64), np.asarray(hm.tri_verts)
    h = trimesh.Trimesh(V * 1e-3, F, process=False)
    how = "manifold topology"
    if not h.is_watertight:                                        # fallback: merge + repair, recorded
        h = trimesh.Trimesh(V * 1e-3, F, process=True)
        trimesh.repair.fill_holes(h); trimesh.repair.fix_normals(h)
        how = "trimesh merge + fill_holes (manifold topology was not watertight)"
    if h.volume < 0:
        h.invert()
    rep = {"full_volume_mm3": float(m.volume), "half_volume_mm3": float(h.volume * 1e9),
           "half_watertight": bool(h.is_watertight), "pieces": int(len(half.decompose())), "half_mesh": how}
    rep["half_ratio"] = rep["half_volume_mm3"] / rep["full_volume_mm3"]
    (out / "car_half.stl").write_text(_ascii(h, "car"))           # ASCII, solid named "car" (Part 2 reads patch names)
    wheels = json.loads((src / "wheels.json").read_text())
    surfaces, bounds = [], [h.bounds]
    for name in ("wheelF", "wheelR"):
        w = trimesh.load(src / f"{name}.stl", force="mesh")
        # the evolution route floats the wheels 0.15 mm over the road (SimScale); Part 2's sinks them 0.3 mm into it
        # (a contact patch snappyHexMesh can mesh): moved down so the lowest point is 0.3 mm under the road
        dz = -SINK_M - float(w.bounds[0][2])
        w.apply_translation([0.0, 0.0, dz])
        (out / f"{name}.stl").write_text(_ascii(w, name))
        bounds.append(w.bounds)
        c = np.array(wheels[name], float); c[2] += dz
        r = 0.5 * float(w.bounds[1][2] - w.bounds[0][2])
        surfaces.append({"name": name, "stl": f"{name}.stl",
                         "rotating": {"origin": c.tolist(), "axis": [0, 1, 0], "omega": -U_MPS / r}})
        rep[f"{name}_moved_down_mm"] = -dz * 1e3
    (out / "assembly.json").write_text(json.dumps({"extra_surfaces": surfaces, "label": src.name}, indent=1))
    b = np.array(bounds)
    rep["bounds_m"] = [b[:, 0].min(0).tolist(), b[:, 1].max(0).tolist()]
    rep["ok"] = bool(rep["half_watertight"] and abs(rep["half_ratio"] - 0.5) < 0.01 and rep["pieces"] == 1)
    (out / "adapt.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep))
    return 0 if rep["ok"] else 2


def _ascii(m, name):
    v, f, n = m.vertices, m.faces, m.face_normals
    lines = [f"solid {name}"]
    for (a, b, c), nn in zip(f, n):
        lines += [f"facet normal {nn[0]:.9g} {nn[1]:.9g} {nn[2]:.9g}", "outer loop",
                  f"vertex {v[a][0]:.9g} {v[a][1]:.9g} {v[a][2]:.9g}", f"vertex {v[b][0]:.9g} {v[b][1]:.9g} {v[b][2]:.9g}",
                  f"vertex {v[c][0]:.9g} {v[c][1]:.9g} {v[c][2]:.9g}", "endloop", "endfacet"]
    lines.append(f"endsolid {name}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
