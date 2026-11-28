"""One car of the layout search: genes -> Parts 1-5 build -> legality -> mass -> SimScale folder.

The whole car is made by the repo's own pipeline (car/build/run_car.py, --gate-only):
Part 1's parametric body (param_body: the halo-canister loft through eight stations from Ref A
to the rear face, sidepods, the cartridge bore, the halo pocket and ballast container, milled
top and bottom with the 3.175 mm ball, 3 mm walls), Part 4's nose cone, front and rear wings,
the official halo and the team's wheels, and Part 5's referee. Wheel supports and tether guides
are left off (user, 2026-11-22: added later); the tether guides and the T4.2 virtual cargo are on,
the cargo where the genes put it (front to back, wide end forward or rearward).

The body is built at a fixed skin (no mass sizing loop): what it weighs is what it weighs, and
below the 48.2 g floor the rest is lead in the T1.22 container under the halo (as in the
sidepod search), so a layout that saves foam is not hidden by the sizing.

    python layout_car.py <out_dir> [--x-json '[...]']      # default: the team-like start car
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
P5 = ROOT / "car/build"
import layout_genes as LG  # noqa: E402

FLOOR_G = 48.2                  # T3.6 48.0 + 0.2 margin, as the sidepod search
CARTRIDGE_SHELL_G = 21.15       # Part 2's car_weight_kg includes the empty cartridge (Appendix vi)
CARTRIDGE_OUT_MM = 12.5         # shortest cartridge 64.5 - 52 mm bore: how far it sticks out
CHAMBER_Z_MM = 31.0             # Part 1's chamber axis (geometry_contract.CHAMBER_AXIS_Z_MM)
# Not counted against a search car: T3.6 is met by the lead below; T5.1 the chamber is the block's own
# drilled bore (18.10-18.15 mm), which this 1 mm build grid draws 0.6 mm small -- the finalists are rebuilt
# on the fine grid, where it is measured properly. The wheel supports are on every car (support_engine.py)
# and every one of their rules, loads and machining checks counts.
WAIVED_PREFIXES = ("T3.6_mass", "T5.1_diameter")
LEAD_CAP_G = None               # read from Part 1's ballast.capacity_kg("lead") at run time


def run_car_args(d: dict, out: Path, lattice: dict | None = None) -> list[str]:
    bp = LG.body_params(d)
    extra = ["--lattice-json", json.dumps(lattice)] if lattice else []
    extra += ["--support", "engine", "--support-json", json.dumps(d["support"])]
    if d.get("freeform"):                     # the free-form body and nose (car/body/freeform.py)
        extra += ["--freeform-json", json.dumps(d["freeform"])]
    return [sys.executable, "-W", "ignore", str(P5 / "run_car.py"), "--out", str(out),
            "--W", str(d["W"]), "--x-front", str(d["x_front"]), "--d-halo", str(d["d_halo"]),
            "--rear-mm", str(d["rear"]), "--no-joints", "--skin-mm", str(d["skin_mm"]),
            "--cargo-json", json.dumps(d["cargo"]), "--tether-json", json.dumps(list(d["tethers"])),
            "--body-json", json.dumps(bp), "--nose-json", json.dumps(d["nose"]),
            "--fwing-json", json.dumps(d["fwing"]), "--rwing-json", json.dumps(d["rwing"])] + extra


def build(x, out: Path, timeout_s=3600, lattice: dict | None = None, support: dict | None = None,
          fwing: dict | None = None, freeform: dict | None = None, rwing: dict | None = None) -> dict:
    """Genes -> a built, checked car in `out`. Returns a record (legal, mass, files, why).
    lattice: the shape search's free-form change of the body (part1 param_body.build lattice=).
    support: the shape search's wheel-support design (support_engine.SupportDesign numbers), else the
    catalogue's (layout_genes.support_design)."""
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)     # absolute: run_car runs in car/build
    d = LG.decode(x)
    if support:
        d["support"] = dict(d["support"], **support)
    if fwing:
        d["fwing"] = dict(d["fwing"], **fwing)
    if rwing:                                 # the searches' rear wing and its struts (rear_wing_engine.py)
        d["rwing"] = rwing
    if freeform:
        d["freeform"] = freeform
    rec = {"genes": [float(v) for v in np.clip(x, 0, 1)], "layout": {k: d[k] for k in ("W", "rear", "halo_rear", "x_front", "d_halo")},
           "decoded": d, "legal": False, "lattice": lattice}
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT / "car/parts"), os.environ.get("PYTHONPATH", "")]))
    t0 = time.time()
    with open(out / "run_car.log", "w") as log:
        try:
            p = subprocess.run(run_car_args(d, out, lattice), cwd=P5, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout_s)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    rec["build_s"] = round(time.time() - t0)
    return judge(out, rec, d, rc)


def judge(out: Path, rec: dict, d: dict, rc=0) -> dict:
    """The referee's verdict, the mass with lead, and the CFD solid, from a finished run_car folder."""
    sp = out / "summary.json"
    if rc != 0 or not sp.exists():
        tail = (out / "run_car.log").read_text()[-600:] if (out / "run_car.log").exists() else ""
        rec["why"] = f"build failed ({rc}): {tail}"
        return _save(out, rec)
    S = json.loads(sp.read_text())
    failed = S["legality"]["summary"]["failed"]
    rec["failed_all"] = failed
    rec["failed"] = [f for f in failed if not f.startswith(WAIVED_PREFIXES) and not _grid_bore(f, S)]
    rec["mass"] = mass(out, d)
    if rec["mass"]["lead_g"] > rec["mass"]["lead_cap_g"] + 1e-9:
        rec["failed"].append("T1.22_ballast_capacity (lead to reach 48.2 g)")
    asm = json.loads((out / "parts" / "assembly.json").read_text())
    se_info = (asm.get("info") or {}).get("support_engine") or {}
    rec["support"] = {"design": se_info.get("design"),
                      **{tag: {k: se_info[tag].get(k) for k in ("t_mm", "saddle_h_mm", "chord_used_mm", "net_g_each",
                                                                 "arm_span_mm", "sag_mm", "sizing_iterations")}
                         for tag in ("front", "rear") if tag in se_info}}
    rec["legal"] = not rec["failed"]
    try:                                                    # every car's supports, to inspect by eye
        import support_render
        support_render.render(out, out / "supports.png")
    except Exception as e:                                  # noqa: BLE001  (a picture never fails a car)
        rec["support_render_error"] = repr(e)[:300]
    if rec["legal"]:
        rec.update(export_cfd(out, d))
        if not (rec["cfd"]["one_piece"] and rec["cfd"]["watertight"]):
            rec["legal"] = False
            rec["failed"].append(f"CFD car is not one closed solid: floating {rec['cfd']['floating_mm3']} mm3, "
                                 f"watertight {rec['cfd']['watertight']}")
    return _save(out, rec)


BLOCK_SPECK_MM3 = 1.0                    # outside the block by less than this in one place: a grid speck, not a body
BORE_BOX = ((-9.2, 9.2), (21.8, 40.2))   # |y|, z of the 18.15 mm bore round the 31 mm chamber axis, 0.1 mm over


def _grid_bore(key, S) -> bool:
    """The 1 mm build grid draws the bore 0.6 mm small (T5.1, waived): its foam inside the bore also fails
    T5.5_chamber_clear (only the body in the chamber) and block_body_inside_official_block (everything
    outside the block inside the bore's own envelope). Waived only when that is all it is."""
    import re
    c = S["legality"]["checks"].get(key) or {}
    note = c.get("note", "")
    if key == "T5.5_chamber_clear":
        return bool(note) and all(part.strip().startswith("machined_body") for part in note.split(","))
    if key == "block_body_inside_official_block":
        # each "<v> mm3 in [[x, y, z], [x, y, z]]": the bore's own envelope is waived, and so are specks under
        # SPECK_MM3 anywhere (a free-form surface grazing the block's floor by a cell edge: 0.2 mm3, 2026-11-28)
        boxes = re.findall(r"([\d.]+) mm3 in \[\[([-\d.]+), ([-\d.]+), ([-\d.]+)\], \[([-\d.]+), ([-\d.]+), ([-\d.]+)\]\]", note)
        (ylo, yhi), (zlo, zhi) = BORE_BOX
        return bool(boxes) and all(float(b[0]) < BLOCK_SPECK_MM3 or
                                   (ylo <= float(b[2]) and float(b[5]) <= yhi and zlo <= float(b[3]) and float(b[6]) <= zhi)
                                   for b in boxes)
    return False


def _save(out, rec):
    (Path(out) / "layout_car.json").write_text(json.dumps(rec, indent=1, default=float))
    return rec


FOAM_G_CM3 = 0.163


def mass(out: Path, d) -> dict:
    """Race mass and COM (car frame: x from the nose tip, z from the track), from the parts as built:
    the milled body (its mesh's volume x 0.163 g/cm3), the printed nose, wings and guides, the halo,
    the wheels with the wheel supports (support_engine.py: each part less the foam its slot removes, at
    the axle; Part 4 puts them in the wheels' mass), lead in the T1.22 container up to 48.2 g, and the
    cartridge."""
    import trimesh
    out = Path(out)
    asm = json.loads((out / "parts" / "assembly.json").read_text())
    f = asm["fixed_hardware_kwargs"]
    body = trimesh.load(str(out / "body_half.stl"), force="mesh")
    foam = 2 * abs(body.volume) * 1e6 * FOAM_G_CM3                     # m3 -> cm3
    bc = body.center_mass * 1e3
    xf, xr = d["x_front"], d["x_front"] + d["W"]
    items = [("foam", foam, bc[0], bc[2]),
             ("printed", f["rear_wing_mass_kg"] * 1e3, f["rear_wing_com"][0] * 1e3, f["rear_wing_com"][2] * 1e3),
             ("wheels_front", f["wheels_front_mass_kg"] * 1e3, xf, f["wheels_front_com"][2] * 1e3),
             ("wheels_rear", f["wheels_rear_mass_kg"] * 1e3, xr, f["wheels_rear_com"][2] * 1e3),
             ("halo", f["halo_mass_kg"] * 1e3, f["halo_com"][0] * 1e3, f["halo_com"][2] * 1e3)]
    car_g = sum(i[1] for i in items)
    sys.path.insert(0, str(ROOT / "car/body"))
    import ballast
    cap = ballast.capacity_kg("lead") * 1e3
    lead = max(0.0, FLOOR_G - car_g)
    # lead sits in the container under the halo pocket (centre 34 mm behind the pocket front, 6.35 deep)
    items.append(("lead", lead, d["x_front"] - 16.0 + d["d_halo"] + 34.0, 24.0 - 6.35 / 2))
    rear_face = d["x_front"] + d["W"] + d["rear"]
    items.append(("cartridge", CARTRIDGE_SHELL_G, rear_face + CARTRIDGE_OUT_MM - 32.5, CHAMBER_Z_MM))
    tot = sum(i[1] for i in items)
    return dict(car_g=car_g, foam_g=foam, lead_g=lead, lead_cap_g=cap, race_mass_g=tot,
                com_x_car_mm=sum(i[1] * i[2] for i in items) / tot, com_z_mm=sum(i[1] * i[3] for i in items) / tot,
                items={i[0]: [round(float(v), 3) for v in i[1:]] for i in items})


CARTRIDGE_R_MM = 9.0             # Appendix vi: 18 mm cartridge; the bore is 18.1-18.15
CARTRIDGE_LEN_MM = 64.5         # the shortest legal cartridge
POCKET_MM = (50.0, 25.0, 3.175)  # Appendix ix: halo pocket, length x width x depth; floor at z = 24
FLOATING_MAX_MM3 = 2.0          # loose specks tolerated in the merged car (dropped); more is a refusal
WRAP_PITCH_MM = 0.3             # grid the CFD surface is rebuilt on


CARTRIDGE_STL = ROOT / "car" / "body" / "hardware_cad" / "canister_steel.stl"   # half (y <= 0.95), neck at x = 0
CARTRIDGE_FULL_LEN_MM = 65.0
BORE_DEPTH_MM = 52.0            # the block's bore, from the rear face (rules_constants.BLOCK_BORE_DEPTH_MM)


def cartridge(rear_face_mm: float):
    """The official cartridge (whole, mm, manifold3d) seated in a car whose rear face is at x = rear_face_mm:
    its base at the bottom of the 52 mm bore, its neck out behind the rear face, axis at the chamber height."""
    import trimesh
    sys.path.insert(0, str(HERE))
    import cadio as CI
    full = None
    for f in (CARTRIDGE_STL, CARTRIDGE_STL.with_name("canister_co2_charge.stl")):   # the shell, filled with its charge
        h = CI.to_manifold(trimesh.load(str(f), force="mesh"))
        h = h + h.mirror([0.0, 1.0, 0.0]) if f == CARTRIDGE_STL else h + h.mirror([0.0, 1.0, 0.0]).translate([0.0, 0.01, 0.0])
        full = h if full is None else full + h
    full = max(full.decompose(), key=lambda q: q.volume())          # coincident shell/charge faces leave specks
    # cartridge x (neck 0 -> base 65) to car x: base at rear_face - bore depth, neck towards +x
    return full.mirror([1.0, 0.0, 0.0]).translate([rear_face_mm - BORE_DEPTH_MM + CARTRIDGE_FULL_LEN_MM, 0.0, CHAMBER_Z_MM])


SUPPORT_CUT_MM = 0.3            # the support is cut this far short of the wheel's inner face ...
SUPPORT_INTO_WHEEL_MM = 0.8     # ... and its section there runs this far into the (solid) CFD wheel


def cfd_supports(S, asm):
    """The supports as the CFD sees them (both halves, mm, manifold3d): each axle's part cut SUPPORT_CUT_MM
    short of its wheels' inner faces, its section there extruded SUPPORT_INTO_WHEEL_MM into the wheel.
    The disc sits just inside the wheel's face (0.05 mm proud of it at its inner face); kept, it would
    leave a 0.05 mm ring step on the wheel no mesh can resolve, and the wheels are solid in the CFD."""
    import manifold3d as mf
    import trimesh
    out = None
    for s in asm["extra_surfaces"]:
        if s["name"] not in ("wheelF", "wheelR"):
            continue
        p = Path(s["stl"])
        w = trimesh.load(str(p if p.exists() else Path(asm.get("_dir", ".")) / Path(s["stl"]).name), force="mesh")
        b = w.bounds * 1e3
        xc, y_in, R = 0.5 * (b[0, 0] + b[1, 0]), float(b[0, 1]), 0.5 * (b[1, 0] - b[0, 0])
        yc = y_in - SUPPORT_CUT_MM
        box = mf.Manifold.cube([2 * R + 2, 2 * yc, 80.0]).translate([xc - R - 1, -yc, -1.0])
        piece = S ^ box
        # the section just inside the cut, extruded into the wheel (each side)
        half = S ^ mf.Manifold.cube([2 * R + 2, 60.0, 80.0]).translate([xc - R - 1, 0.0, -1.0])
        cs = half.warp_batch(lambda V: np.c_[V[:, 2], V[:, 0], V[:, 1]]).slice(yc - 0.05)   # (z, x) at y
        if cs.area() > 0.01:
            L = SUPPORT_CUT_MM + 0.05 + SUPPORT_INTO_WHEEL_MM
            plug = mf.Manifold.extrude(cs, L).translate([0.0, 0.0, yc - 0.05])
            plug = plug.warp_batch(lambda V: np.c_[V[:, 1], V[:, 2], V[:, 0]])                # back to (x, y, z)
            piece = piece + plug + plug.mirror([0.0, 1.0, 0.0])
        out = piece if out is None else out + piece
    return out if out is not None else mf.Manifold()


def export_cfd(out: Path, d, pitch_mm: float | None = None) -> dict:
    """ONE closed solid of the whole car for the CFD, and the wheels.

    The parts as built only touch (the halo sits in its pocket, the rear wing's pylon stands on the
    deck, a guide hangs under the floor): cut from the air as they are they would leave slits, a
    pocket of trapped air under the halo and separate solids. Here they are fused the way the glued
    car is: the halo pocket (a box in Part 1, the halo's own outline when milled) and the ballast
    container under it are filled, the cartridge is in the chamber (a plain cylinder, 12.5 mm out),
    the rear wing and the guides are sunk / stalked into the body, everything is unioned, enclosed
    cavities are removed, and the result must be one closed piece or the car is refused.
    Writes body_full.stl (m, car frame; body_wrap_mm.stl is the same in mm), wheelF/wheelR.stl
    (+y wheels), wheels.json."""
    import manifold3d as mf
    import trimesh
    sys.path.insert(0, str(HERE))
    import cadio as CI
    out = Path(out)
    asm = json.loads((out / "parts" / "assembly.json").read_text())
    asm["_dir"] = str(out / "parts")

    def full(p):
        m = trimesh.load(str(p), force="mesh"); m.apply_scale(1e3)          # m -> mm (manifold is float32)
        # 20 um across the centre plane, so the two halves overlap: halves that only share the
        # plane's face leave zero-thickness leftovers there (edges shared by four faces)
        try:                                                    # as written; cleaning can break a good one
            h = CI.to_manifold(m)                               # (dropping 2 zero-area faces did, car L0005)
        except ValueError:
            h = CI.to_manifold(CI.clean(m))
        h = h.translate([0.0, -0.02, 0.0])
        return h + h.mirror([0.0, 1.0, 0.0])

    def box(x0, x1, y, z0, z1):
        return mf.Manifold.cube([x1 - x0, 2 * y, z1 - z0]).translate([x0, -y, z0])

    rep = {}
    body = full(out / "body_half.stl")
    bp = sorted(body.decompose(), key=lambda q: -q.volume())
    rep["body_specks_mm3"] = round(sum(q.volume() for q in bp[1:]), 3)
    body = bp[0]
    parts = {n: full(out / "parts" / f"{n}.stl") for n in ("nose", "halo", "fwing", "rwing", "tethers")
             if (out / "parts" / f"{n}.stl").exists()}
    whole = body
    # halo pocket + ballast container: solid from under the container to the pocket's rim
    px0 = d["x_front"] - 16.0 + d["d_halo"]
    # 0.05 mm into the foam all round and 0.03 under the rim: faces exactly on the pocket's own walls
    # leave zero-area leftovers in the union
    fill = box(px0 - 0.05, px0 + POCKET_MM[0] + 0.05, POCKET_MM[1] / 2 + 0.05, 24.0 - 6.35 - 0.6, 24.0 + POCKET_MM[2] - 0.03)
    rep["pocket_fill_mm3"] = round((fill - body - parts["halo"]).volume(), 1)
    whole += fill + parts["halo"]
    # the slot between the halo's back and the body rising behind it (the pocket is cut as a box, ~1.5 mm
    # longer than the halo's foot): closed up to the body's own height just behind, never above the halo's back
    xb = px0 + POCKET_MM[0]
    behind = box(xb + 2.0, xb + 4.0, POCKET_MM[1] / 2, 24.0, 60.0) ^ body
    if behind.volume() > 0:
        zt = min(33.9, behind.bounding_box()[5])
        slot = box(xb - 0.05, xb + 2.5, POCKET_MM[1] / 2 + 0.05, 24.0, zt)
        rep["halo_back_slot_mm3"] = round((slot - body - parts["halo"]).volume(), 1)
        whole += slot
    # the cartridge in the chamber
    rear_face = d["x_front"] + d["W"] + d["rear"]
    def cyl(x0, x1, r):
        return mf.Manifold.cylinder(x1 - x0, r, r, 64).rotate([0.0, 90.0, 0.0]).translate([x0, 0.0, CHAMBER_Z_MM])

    # inside the chamber 0.5 mm into the wall (no 0.07 mm crevice round the cartridge), 18 mm where it sticks out
    # the chamber filled 0.5 mm into its wall (no crevice round the cartridge), then the real cartridge in it:
    # the official shell (Appendix vi; car/body/hardware_cad/canister_steel.stl), base first, its neck out
    # behind the rear face -- an 18 mm body, an 8.7 mm neck, the shoulder between
    whole += cyl(rear_face + CARTRIDGE_OUT_MM - CARTRIDGE_LEN_MM - 0.5, rear_face - 0.3, CARTRIDGE_R_MM + 0.5)
    whole += cartridge(rear_face)
    whole += parts["nose"]
    # the wheel supports (support_engine.py): the saddle is the body's own underside, so it adds only where
    # the arm and the transition leave the body; inside the wheel's cylinder the CFD wheel is solid, so the
    # disc, the axle and the hubcap are inside it and the transition's end section runs on into the wheel
    if (out / "parts" / "supports.stl").exists():
        sup = cfd_supports(full(out / "parts" / "supports.stl"), asm)
        rep["supports_mm3"] = round(sup.volume(), 1)
        whole += sup
    for n in ("fwing", "rwing", "tethers"):
        if n not in parts:
            continue
        for piece in parts[n].decompose():
            if piece.volume() < 0.5:
                continue
            if (piece ^ whole).volume() < 0.5:                  # only touching, or apart: join it
                (x0, y0, z0, x1, y1, z1) = piece.bounding_box()
                above = (box(x0, x1, max(abs(y0), abs(y1)), z1 - 0.2, 70.0) ^ whole)
                below = (box(x0, x1, max(abs(y0), abs(y1)), 0.0, z0 + 0.2) ^ whole)
                if n == "tethers" and above.volume() > 0:       # a stalk up into the floor above the ring
                    zt = above.bounding_box()[2]
                    xc = 0.5 * (x0 + x1)
                    piece += box(xc - 1.0, xc + 1.0, 1.5, z1 - 0.5, zt + 1.0)
                    rep.setdefault("stalks_mm", []).append(round(zt - z1, 2))
                elif below.volume() > 0:                        # stands on the body: sink it in
                    gap = z0 - below.bounding_box()[5]
                    piece += piece.translate([0.0, 0.0, -(max(gap, 0.0) + 0.6)])
                    rep.setdefault("sunk_mm", {})[n] = round(float(gap), 2)
            whole += piece
    pcs = sorted(whole.decompose(), key=lambda q: -q.volume())
    rep["floating_mm3"] = round(sum(abs(q.volume()) for q in pcs[1:]), 2)
    # enclosed cavities are separate pieces with negative volume: dropped with the specks above
    rep["cavity_mm3"] = round(float(sum(-q.volume() for q in pcs[1:] if q.volume() < 0)), 2)
    M = pcs[0].simplify(0.02)                                   # no edge under 20 um (Boolean slivers)

    # The surface the CFD gets: the merged solid rebuilt from its own occupancy on a 0.3 mm grid
    # (layout_wrap.py). The union carries needle triangles, 0.03 mm ledges and edges where two surfaces
    # only touch; the CAD importer rejects those, and nothing thinner than a cell survives the rebuild,
    # so no slit or void is left either.
    import layout_wrap as LW
    w = LW.wrap(M, pitch_mm or WRAP_PITCH_MM)
    q = LW.quality(w)
    pts, _ = trimesh.sample.sample_surface(w, 20000, seed=0)
    _, dist, _ = trimesh.proximity.closest_point(CI.from_manifold(M), pts)              # mm, wrap -> merged solid
    rep["wrap_pitch_mm"] = pitch_mm or WRAP_PITCH_MM
    rep["cartridge"] = "official"
    rep.update(cfd_faces=q["faces"], car_volume_mm3=round(float(w.volume)), merged_volume_mm3=round(float(M.volume())),
               wrap_p99_mm=round(float(np.percentile(dist, 99)), 3), wrap_max_mm=round(float(dist.max()), 3),
               one_piece=rep["floating_mm3"] <= FLOATING_MAX_MM3,
               watertight=bool(q["closed"] and w.volume > 0 and abs(w.volume / M.volume() - 1) < 0.01))
    w.export(str(out / "body_wrap_mm.stl"))
    wm = w.copy(); wm.apply_scale(1e-3)
    wm.export(str(out / "body_full.stl"))
    wheels = {}
    for s in asm["extra_surfaces"]:
        if s["name"] in ("wheelF", "wheelR"):
            wheels[s["name"]] = tuple(s["rotating"]["origin"])
            src = Path(s["stl"]) if Path(s["stl"]).exists() else out / "parts" / Path(s["stl"]).name
            (out / f"{s['name']}.stl").write_bytes(src.read_bytes())
    (out / "wheels.json").write_text(json.dumps(wheels))
    rep["wheels"] = wheels
    return {"cfd": rep}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--x-json", default=None)
    a = ap.parse_args()
    x = np.array(json.loads(a.x_json)) if a.x_json else LG.x_default()
    r = build(x, Path(a.out))
    print(json.dumps({k: r.get(k) for k in ("legal", "layout", "failed", "mass", "build_s", "cfd", "why")}, indent=1, default=float))
