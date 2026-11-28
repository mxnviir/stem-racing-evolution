"""Squareness and hard-coded-feature audit of a built car's body (owner 2026-11-28: no extreme squareness, nothing
extremely hard-coded in the cars the program designs).

    python carhybrid/shape_audit.py <evolution car folder> <out json>

Measured on body_half.stl (the milled body, mm or m):
- section rectangularity: at stations every 5 mm from Ref A + 5 to the rear face - 25 (clear of the block's flat rear
  face and the bore), the cross-section's area / its bounding box's area. An ellipse is 0.785, the old loft's
  superellipse (p = 2.6) about 0.85, a rectangle 1.0. Flag: any station > 0.93 ("box section").
- flat patches: faces whose normal is within 2 deg of an axis, grouped by the plane they lie on (0.05 mm bins);
  patches over 20 mm2 listed with orientation, position and area. Known rule planes are named (the centre plane y = 0,
  the block's rear face, the halo-pocket floor z = 24, the T4.4/sidepod-top cap z = 27.175); any other large flat patch
  is flagged.
- creases: total length of edges whose faces meet at more than 60 deg, excluding the closing face at y = 0.
"""
import json, sys
from pathlib import Path

import numpy as np

CAP_Z_MM = 27.175
POCKET_FLOOR_Z_MM = 24.0


def _near_halo(pts, halo, reach_mm=15.0):
    """Most of a flat patch's area within reach_mm of the halo seat in x (the ring the seat's rim makes)."""
    x = pts[:, 0]
    return np.mean((x > halo[0] - reach_mm) & (x < halo[1] + 2.0)) > 0.5


def main(folder, out_json):
    import trimesh
    folder = Path(folder)
    rec = json.loads((folder / "layout_car.json").read_text()) if (folder / "layout_car.json").exists() else {}
    d = rec.get("decoded") or {}
    m = trimesh.load(folder / "body_half.stl", force="mesh")
    if m.bounds[1][0] < 1.0:
        m.apply_scale(1e3)
    ref_a = d.get("x_front", m.bounds[0][0] + 20.0) - 16.0
    rear = d.get("x_front", 0) + d.get("W", 0) + d.get("rear", 0) if d else m.bounds[1][0]
    # the halo's seat (halo_seat.py): the body meets the pocket's rim (24 + 3.175 = 27.175 mm) flush, by the owner's
    # 2026-11-23 rule (no trench), and T4.4.2 keeps it at or under the rim ahead of / beside the halo: a flat ring there
    halo = (ref_a + d["d_halo"], ref_a + d["d_halo"] + 50.0) if "d_halo" in d else None
    # the nose's foam/print joint (evo_patches/nose_material.py): the glued face at x_cut and the tenon's pocket
    x_cut = None
    for src in (folder / "nose_split.json", folder / "summary.json"):
        if src.exists():
            js = json.loads(src.read_text())
            x_cut = js.get("x_cut_mm") or ((js.get("mass") or {}).get("nose_split") or {}).get("x_cut_mm")
            if x_cut:
                break
    axles = (d["x_front"], d["x_front"] + d["W"]) if "x_front" in d and "W" in d else ()
    # sections
    rect, stations = [], np.arange(ref_a + 5.0, rear - 25.0, 5.0)
    for x in stations:
        s = m.section(plane_origin=[x, 0, 0], plane_normal=[1, 0, 0])
        if s is None:
            continue
        try:
            p2, _ = s.to_2D(normal=[1, 0, 0]) if hasattr(s, 'to_2D') else s.to_planar(normal=[1, 0, 0])
        except Exception:
            continue
        A = sum(abs(poly.area) for poly in p2.polygons_full)
        if A <= 1.0:
            continue
        (y0, z0), (y1, z1) = p2.bounds
        rect.append((float(x), float(A / max((y1 - y0) * (z1 - z0), 1e-9))))
    r = np.array([v for _, v in rect]) if rect else np.array([np.nan])
    # flat patches
    n, c, A = m.face_normals, m.triangles_center, m.area_faces
    patches = []
    for ax, name in ((0, "x"), (1, "y"), (2, "z")):
        al = np.abs(n[:, ax]) > np.cos(np.radians(2.0))
        if not al.any():
            continue
        key = np.round(c[al, ax] / 0.05).astype(int)
        sgn = np.sign(n[al, ax]).astype(int)
        for (k, sg) in set(zip(key.tolist(), sgn.tolist())):
            sel = (key == k) & (sgn == sg)
            area = float(A[al][sel].sum())
            if area < 20.0:
                continue
            pos = k * 0.05
            label = None
            cc = c[al][sel]
            (cc_x0, cc_y0, cc_z0), (cc_x1, cc_y1, cc_z1) = cc.min(0), np.abs(cc).max(0) * [0, 1, 1] + cc.max(0) * [1, 0, 0]
            if name == "y" and abs(pos) < 0.1:
                label = "centre plane (closing face)"
            elif name == "x" and abs(pos - rear) < 0.6:
                label = "block rear face"
            elif name == "z" and abs(pos - POCKET_FLOOR_Z_MM) < 0.2:
                label = "halo pocket floor"
            elif name == "z" and abs(pos - CAP_Z_MM) < 0.2 and halo is not None and _near_halo(c[al][sel], halo):
                label = "halo rim, body flush with the seat (T4.4.2, owner 2026-11-23)"
            elif name == "z" and abs(pos - CAP_Z_MM) < 0.2:
                label = "27.175 mm: the team sidepod cap (dropped 2026-11-28 by sidepod_rules S1; flat here = an older builder)"
            elif name == "z" and abs(pos - (POCKET_FLOOR_Z_MM - 6.35)) < 0.2:
                label = "ballast container floor (App. ix, internal)"
            elif name == "x" and abs(pos - (rear - 52.0)) < 2.5:
                label = "cartridge chamber end (52 mm bore, internal)"
            elif name == "x" and sg < 0 and x_cut and abs(pos - x_cut) < 0.3:
                label = "nose glued face at x_cut (foam meets print; covered by the print)"
            elif x_cut and cc_x0 >= x_cut - 0.5 and cc_x1 <= x_cut + 20.0 and (
                    (name == "y" and abs(abs(pos) - 4.0) < 0.15) or (name == "z" and sg < 0)):
                label = "tenon pocket (nose joint, internal)"
            elif (name == "y" and abs(abs(pos) - 3.0) < 0.15) or (name == "z" and sg < 0 and abs(pos - 8.0) < 0.3):
                label = "official block underside slot (6 x 6 mm, kept unmachined)"
            elif name == "x" and sg < 0 and cc_y1 < 4.5 and cc_z1 < 12.5:
                label = "official block underside slot end"
            elif name == "x" and sg < 0 and halo is not None and abs(pos - halo[1]) < 0.6:
                label = "halo rear face seat (body backs the halo, team rule)"
            elif name == "z" and sg > 0 and halo is not None and halo[1] <= cc_x0 <= halo[1] + 15.0:
                label = "halo-canister loft deck (required body behind the halo)"
            elif name == "y" and abs(abs(pos) - 6.35) < 0.15:
                label = "ballast container side wall (App. ix, internal)"
            elif name == "z" and abs(pos - 14.12) < 0.35 and cc_x1 - cc_x0 < 6.0:
                label = "axle-height feature (fixed hardware)"
            elif name == "y" and axles and any(abs(0.5 * (cc_x0 + cc_x1) - a) < 15.0 for a in axles):
                label = "wheel zone face (T7.9 clearance round a wheel)"
            elif name == "z" and pos < 2.0:
                label = "floor near the track"
            patches.append({"axis": name, "normal_sign": sg, "at_mm": round(pos, 2), "area_mm2": round(area, 1),
                            "label": label or "UNEXPLAINED",
                            "bbox_mm": [np.round(cc.min(0), 1).tolist(), np.round(cc.max(0), 1).tolist()]})
    patches.sort(key=lambda p: -p["area_mm2"])
    # creases
    ang = m.face_adjacency_angles
    e = m.face_adjacency_edges
    L = np.linalg.norm(m.vertices[e[:, 0]] - m.vertices[e[:, 1]], axis=1)
    mid_y = 0.5 * (m.vertices[e[:, 0], 1] + m.vertices[e[:, 1], 1])
    sharp = (ang > np.radians(60)) & (np.abs(mid_y) > 0.05)
    out = {"stations": len(rect), "rectangularity_median": float(np.nanmedian(r)), "rectangularity_max": float(np.nanmax(r)),
           "rectangularity_by_x": [[round(x, 1), round(v, 3)] for x, v in rect],
           "box_sections": [[round(x, 1), round(v, 3)] for x, v in rect if v > 0.93],
           "flat_patches": patches[:20],
           "unexplained_flat_mm2": round(sum(p["area_mm2"] for p in patches if p["label"] == "UNEXPLAINED"), 1),
           "cap_area_mm2": round(sum(p["area_mm2"] for p in patches if p["label"] and p["label"].startswith("sidepod-top cap")), 1),
           "halo_rim_area_mm2": round(sum(p["area_mm2"] for p in patches if p["label"] and p["label"].startswith("halo rim")), 1),
           "crease_length_mm": round(float(L[sharp].sum()), 1),
           "body_area_mm2": round(float(m.area), 1)}
    out["flags"] = ([f"box section at x {x} (rectangularity {v})" for x, v in out["box_sections"]] +
                    ([f"unexplained flat patches {out['unexplained_flat_mm2']} mm2"] if out["unexplained_flat_mm2"] > 50 else []) +
                    ([f"hard-coded sidepod-top cap {out['cap_area_mm2']} mm2"] if out["cap_area_mm2"] > 20 else []))
    Path(out_json).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("stations", "rectangularity_median", "rectangularity_max", "unexplained_flat_mm2",
                                          "cap_area_mm2", "crease_length_mm", "flags")}))


if __name__ == "__main__":
    main(*sys.argv[1:3])
