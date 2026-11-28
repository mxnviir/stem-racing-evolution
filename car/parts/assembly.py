"""
assembly.py -- every non-body part of the car, as CFD patches + a mass rollup.

    build(W, x_front, d_halo, body_half_stl, out_dir) -> dict

writes, in out_dir (right halves, metres, ASCII STL):
    wheelF.stl, wheelR.stl   closed wheels, sunk 0.3 mm into the track, ROTATING
    supports.stl             v2 CAD wheel supports (Part 1 hardware_cad)
    halo.stl                 halo + helmet (Part 1 hardware_cad)
    fwing.stl                front wing + its centreline mount
    rwing.stl                rear wing + pylon
    tethers.stl              T6 tether line guides
    assembly.json            everything below

assembly.json carries:
    extra_surfaces           exactly what Part 2's OpenFOAMRunConfig /
                             AdjointRunConfig.extra_surfaces take
    fixed_hardware_kwargs    exactly what Part 2's FixedHardwareSpec takes, so
                             Part 3's mass rollup sees the real parts
    wheel_moi_kg_m2          the designed wheel, for the race objective
    gates                    every regulation margin Part 4 can measure (mm)

Why each part is in the CFD: measured on GitHub Actions (rnd-cfd, 2026-09-25)
wheels are 64-75 % of the car's drag and the body only 13-22 %.
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in ("body", "race"):
    _d = HERE.parent / _p
    if _d.is_dir() and str(_d) not in sys.path:
        sys.path.append(str(_d))

import component_contract as cc   # noqa: E402
import wheel as wh                # noqa: E402
import wings as wg                # noqa: E402
import nose as ns                 # noqa: E402
import beam_support as bsm        # noqa: E402
import support_engine as se      # noqa: E402
import front_wing_engine as fwe   # noqa: E402

SPEED_MPS = 20.0
SINK_MM = 0.3
PART_MATERIAL = "PLA"           # wings, mounts, tether guides: printed, solid
# No separate axles (team spec 2026-09-27): 3 mm stubs printed as part of each
# support, bearings pressed on. Their mass is in the support's.


def _export_half(mesh, path: Path) -> None:
    m = mesh.copy()
    m.vertices[m.vertices[:, 1] < 0.0, 1] = 0.0
    if m.is_volume and len(m.split(only_watertight=False)) == 1:
        # one solid only: several in one file (the rear wing and its struts, overlapping) would be merged by
        # the cleanup into whatever manifold3d makes of them (the struts were lost, gallery 5 G00); those
        # are made safe one by one before they are put together
        m = _stl_safe(m)
    m.export(str(path), file_type="stl_ascii")


STL_EDGE_MM = 0.01       # mm: merged before a part is written (_stl_safe); above support_engine._weld's 2 um


def _stl_safe(m):
    """A closed part (metres) that is still closed once written to STL and read back: marching cubes
    leaves edges of 1e-6 mm, which the file's float32 collapses, and the part came back open (gallery G02's
    nose: the referee could not measure it). Edges under STL_EDGE_MM are merged by manifold3d's simplify,
    which keeps the solid closed. Unchanged where that would not give a closed solid."""
    import io
    import trimesh
    import support_engine as se
    q = m.copy(); q.apply_scale(1e3)
    try:
        t = se._to_trimesh(se._to_manifold(q).simplify(STL_EDGE_MM))
    except Exception:  # noqa: BLE001
        return m
    if not t.is_volume:
        return m
    t.apply_scale(1e-3)
    back = trimesh.load(io.BytesIO(t.export(file_type="stl_ascii").encode()), file_type="stl")
    return t if back.is_volume else m


def _half_of(mesh):
    """Right half (y >= 0) of a centreline part, capped on the symmetry plane."""
    import trimesh
    return trimesh.intersections.slice_mesh_plane(mesh, [0, 1, 0], [0, 0, 0], cap=True)


def tether_guide(x_mm: float, id_mm: float = cc.TETHER_GUIDE_ID, od_mm: float = cc.TETHER_GUIDE_OD,
                 length_mm: float = 2.0, z_bottom_mm: float = cc.TETHER_GUIDE_BOTTOM_Z):
    import trimesh
    ring = trimesh.creation.annulus(r_min=id_mm / 2e3, r_max=od_mm / 2e3,
                                    height=length_mm / 1e3, sections=48)
    ring.apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [0, 1, 0]))
    ring.apply_translation([x_mm / 1e3, 0.0, (z_bottom_mm + od_mm / 2) / 1e3])
    return _half_of(ring), dict(id_mm=id_mm, od_mm=od_mm, x_mm=x_mm,
                                bottom_mm=z_bottom_mm)


def _full_mass_com(half_meshes, density_g_cm3: float):
    """Mass (kg) and COM (m) of the FULL symmetric part from its right half."""
    vol = sum(abs(m.volume) for m in half_meshes)
    if vol <= 0:
        return 0.0, (0.0, 0.0, 0.0)
    cx = sum(abs(m.volume) * m.center_mass[0] for m in half_meshes) / vol
    cz = sum(abs(m.volume) * m.center_mass[2] for m in half_meshes) / vol
    return 2 * vol * density_g_cm3 * 1e3, (cx, 0.0, cz)


def fixed_hardware_kwargs(W_mm, x_front_mm, d_halo_mm, parts_mass: dict,
                          wheel_design, support_kg=(1.396e-3, 1.545e-3),
                          rear_face_mm: float | None = None) -> dict:
    """Kwargs for Part 2's FixedHardwareSpec from the REAL parts.

    The spec has a single "rear_wing" slot; it carries every aero part here
    (front wing + mount, rear wing + pylon, tether guides) at their combined
    COM, replacing the old 5 g placeholder.
    """
    import fixed_hardware as fh
    from geometry_contract import CO2_MASS_KG, R_WHEEL_M
    ref_a = cc.ref_plane_A(x_front_mm) / 1e3
    ref_b = cc.ref_plane_B(x_front_mm, W_mm) / 1e3
    rear_face = ref_b + 0.040 if rear_face_mm is None else rear_face_mm / 1e3
    hw = fh.compute_default_fixed_hardware_inputs(W_mm, x_front_mm, d_halo_mm,
                                                  ref_a, ref_b, rear_face_x_m=rear_face)
    halo = hw["halo_geometry"]
    hz = [z for _y, z in halo.cross_section_yz_m]
    f, r = wh.design(wheel_design)
    sup_f, sup_r = support_kg                       # each; default = v2 CAD (fixed_hardware)
    mf = 2 * (f.mass * 1e-3 + sup_f)
    mr = 2 * (r.mass * 1e-3 + sup_r)
    xf, xr = x_front_mm / 1e3, (x_front_mm + W_mm) / 1e3
    aero = [v for k, v in parts_mass.items() if k in ("fwing", "rwing", "tethers", "nose")]
    m_aero = sum(m for m, _c in aero)
    c_aero = tuple(sum(m * c[i] for m, c in aero) / m_aero for i in range(3)) if m_aero else (rear_face, 0, 0.05)
    return {
        "co2_cartridge_mass_kg": CO2_MASS_KG,
        "co2_cartridge_com": tuple(v / 1e3 for v in hw["canister_com_mm"]),
        "rear_wing_mass_kg": m_aero,
        "rear_wing_com": c_aero,
        "wheels_axles_mass_kg": mf + mr,
        "wheels_axles_com": ((mf * xf + mr * xr) / (mf + mr), 0.0, R_WHEEL_M),
        "wheels_front_mass_kg": mf, "wheels_front_com": (xf, 0.0, R_WHEEL_M),
        "wheels_rear_mass_kg": mr, "wheels_rear_com": (xr, 0.0, R_WHEEL_M),
        "halo_mass_kg": fh.HALO_MASS_KG,
        "halo_com": (0.5 * (halo.x_front_m + halo.x_rear_m), 0.0, 0.5 * (min(hz) + max(hz))),
    }


def _printed_minus_body(part_m, body_m):
    """A printed part (right half, metres) less what lies inside the foam body: it ends on the body's surface instead
    of running on through the foam to the centre plane, where it was cut open (a foam nose's wing mount)."""
    import manifold3d as mf
    import support_engine as se
    import wing_loft as wl
    p = part_m.copy(); p.apply_scale(1e3); b = body_m.copy(); b.apply_scale(1e3)
    pcs = sorted((se._to_manifold(p) - se._to_manifold(b)).decompose(), key=lambda q: -q.volume())
    t = wl.to_trimesh(mf.Manifold.batch_boolean([q for q in pcs if q.volume() > 1e-3 * pcs[0].volume()], mf.OpType.Add)) if pcs else p
    t.apply_scale(1e-3)
    return t


def _union_closed(parts_m):
    """Union of closed parts (metres) as ONE closed solid: in mm through manifold3d and support_engine's
    cleanup (Boolean seams leave zero-area flaps; a plain trimesh union of the nose and its pylons came out
    open, gallery G02, and the referee could not measure it). Not one closed solid: the car is refused."""
    import manifold3d as mf
    import support_engine as se
    import wing_loft as wl
    Ms = []
    for p in parts_m:
        q = p.copy(); q.apply_scale(1e3); Ms.append(se._to_manifold(q))
    # its own mesh if closed, else the cleaned one, else with its slivers merged (wing_loft.to_trimesh): the
    # nose and its swept struts came out open through the cleanup alone (gallery 7 G00)
    t = wl.to_trimesh(mf.Manifold.batch_boolean(Ms, mf.OpType.Add))
    if not t.is_volume:
        raise ValueError("illegal car: the nose and its wing pylons do not make one closed solid")
    t.apply_scale(1e-3)
    return t


def _legacy_front_wing(front, nm, x_front_mm, R, body, out, gates, info, surfaces, parts_mass):
    """Part 4's parametric front wing (wings.FrontWing) through the nose, on its centreline mount."""
    import trimesh
    fwm = wg.build_front_wing(front, x_front_mm, R, body)
    blockers = [body] + ([nm["nose"]] if nm is not None else [])
    y_cut = wg.front_span_cut_mm(fwm["wing"], blockers)
    need = y_cut + cc.FRONT_SPAN_SEGMENT_MIN + 0.05
    if y_cut > front.mount_t_mm / 2 and front.half_span_mm < need:
        front = replace(front, half_span_mm=need)
        fwm = wg.build_front_wing(front, x_front_mm, R, body)
        y_cut = wg.front_span_cut_mm(fwm["wing"], blockers)
    info["front_span_cut_mm"] = y_cut
    gates.update(wg.front_wing_gates(front, fwm, x_front_mm, R, nose_half_width_mm=y_cut))
    _export_half(trimesh.util.concatenate(list(fwm.values())), out / "fwing.stl")
    surfaces.append({"name": "fwing", "stl": str(out / "fwing.stl"), "rotating": None})
    parts_mass["fwing"] = _full_mass_com(list(fwm.values()), cc.DENSITY_G_CM3[PART_MATERIAL])


def build(W_mm: float, x_front_mm: float, d_halo_mm: float, body_half_stl: str,
          out_dir: str, wheel_design="team_stl",
          front: wg.FrontWing = wg.FrontWing(), rear: wg.RearWing = wg.RearWing(),
          rotate_wheels: bool = True, nose: "ns.NoseCone | None" = None,
          support: "bsm.BeamSupport | None" = None, support_modes: dict | None = None,
          supports: bool = True, tethers: bool = True, rear_face_mm: float | None = None,
          tether_offsets_mm: tuple = (5.0, 5.0), freeform=None) -> dict:
    """support: None = the team's v2 CAD supports (reference), or a BeamSupport.
    supports / tethers False: no wheel supports / tether guides at all (the layout search,
    search/layout/layout_car.py, adds them later); the wheels then float in the CFD.
    rear_face_mm: the body's rear face when it is not Ref B + 40 (sets the cartridge's mass)."""
    import trimesh
    import hardware_geometry as hg
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    body = trimesh.load(str(body_half_stl), force="mesh")
    f, r = wh.design(wheel_design)
    from geometry_contract import FRONT_WHEEL_INNER_Y_MM as FY, REAR_WHEEL_INNER_Y_MM as RY
    x_rear = x_front_mm + W_mm
    R = 0.5 * (f.R + r.R)

    surfaces, gates, parts_mass, info = [], {}, {}, {}
    wheel_meshes = {}

    # -- wheels (rotating closed cylinders) --------------------------------
    for name, x, y_in, w in (("wheelF", x_front_mm, FY, f.w), ("wheelR", x_rear, RY, r.w)):
        cyl = wh.cfd_surface(R, w, x, y_in, SINK_MM)
        wheel_meshes[name] = cyl
        _export_half(cyl, out / f"{name}.stl")
        origin = (x / 1e3, (y_in + w / 2) / 1e3, (R - SINK_MM) / 1e3)
        surfaces.append({"name": name, "stl": str(out / f"{name}.stl"),
                         "rotating": ({"origin": origin, "axis": (0, 1, 0),
                                       "omega": -SPEED_MPS / (R / 1e3)} if rotate_wheels else None)})
    gates["T7.5_diameter_min"] = 2 * R - cc.WHEEL_DIA_MIN
    gates["T7.5_diameter_max"] = cc.WHEEL_DIA_MAX - 2 * R
    gates["T7.4_front_contact"] = f.w - cc.FRONT_CONTACT_MIN
    gates["T7.4_rear_contact"] = r.w - cc.REAR_CONTACT_MIN
    gates["T7.2_front_gap"] = 2 * FY - cc.FRONT_GAP_MIN
    gates["T7.2_rear_gap"] = 2 * RY - cc.REAR_GAP_MIN

    # -- supports: v2 CAD, or the parametric beam -------------------------
    support_kg = (1.396e-3, 1.545e-3)
    if not supports:
        sup, support_kg = [], (0.0, 0.0)
    elif isinstance(support, se.SupportDesign):
        # the searches' own supports, generated from THIS body at each axle (support_engine.py):
        # a saddle flush with the underside, a lofted arm, a fillet into the disc, the axle stub
        hubs = wh.hub_span_mm()
        import manifold3d as mf
        from geometry_contract import CHAMBER_AXIS_Z_MM, MIN_MACHINED_WALL_MM, MODEL_BLOCK_BORE_DIA_MM, \
            MODEL_BLOCK_BORE_DEPTH_MM
        face = rear_face_mm if rear_face_mm is not None else cc.ref_plane_B(x_front_mm, W_mm) + 40.0
        rk = MODEL_BLOCK_BORE_DIA_MM / 2 + MIN_MACHINED_WALL_MM
        x0 = face - MODEL_BLOCK_BORE_DEPTH_MM - MIN_MACHINED_WALL_MM
        chamber = mf.Manifold.cylinder(face - x0, rk, rk, 96).rotate([0.0, 90.0, 0.0]).translate([x0, 0.0, CHAMBER_AXIS_Z_MM])
        res = se.build(support, body, [("front", x_front_mm, R - SINK_MM, R, FY, f.total_width, hubs[0]),
                                       ("rear", x_rear, R - SINK_MM, R, RY, r.total_width, hubs[1])], out_dir=out,
                       keep_out_mm=chamber)
        sup = [res["support"]]
        support_kg = tuple(res["net_kg_each"])               # per side: the part and its hubcap, less the foam it replaces
        gates.update(res["gates"])
        info["support_engine"] = dict(res["info"], design=se.params_dict(support))
        info["support_mode"] = {"front": "engine", "rear": "engine"}
        for tag, cap in res["hubcaps"].items():
            cap.export(str(out / f"hubcap_{tag}.stl"))
            info.setdefault("hubcap_stls", {})[tag] = str(out / f"hubcap_{tag}.stl")
    elif support is None:
        sup = [hg.build_wheel_assembly(a, x)[f"{a}_wheel_support_right"]
               for a, x in (("front", x_front_mm), ("rear", x_rear))]
        # Safety net: the rear CAD bottomed at 1.40 mm before it was fixed.
        sup = [trimesh.intersections.slice_mesh_plane(
            m, [0, 0, 1.0], [0, 0, (cc.TRACK_CLEARANCE_MIN + 0.01) / 1e3], cap=True) for m in sup]
    else:
        # The team's architecture, parametric: one PA12 beam per axle through
        # the body, stub axles, optional wheel-face discs.
        sup, support_kg = [], []
        for tag, x, y_in, w, disc, hub in zip(
                ("front", "rear"), (x_front_mm, x_rear), (FY, RY), (f.total_width, r.total_width),
                (support.disc_front, support.disc_rear), wh.hub_span_mm()):
            # the pod's channel comes from THIS body at this axle, so the
            # support follows the body through the search and the sizing
            ch = bsm.pod_channel(support, body, x, R - SINK_MM, R) if support.pod else {}
            bm = bsm.build(support, x, R - SINK_MM, y_in, w, disc, hub_mm=hub, R_mm=R, channel=ch,
                           mode=(support_modes or {}).get(tag))
            info.setdefault("support_mode", {})[tag] = bm["support_mode"]
            gates.update(bsm.gates(support, bm, x, R - SINK_MM, R, tag))
            if bm["_pod"]:
                info.setdefault("pod", {})[tag] = ch
                half = ch["length_mm"] / 2
                corner = max(math.hypot(half, R - SINK_MM - ch["z_floor_mm"]),
                             math.hypot(half, ch["z_arch_mm"] - (R - SINK_MM)))
                gates[f"T7.12.1_{tag}_pod_in_cylinder"] = R - corner
                gates[f"machining_{tag}_foam_over_pod"] = ch["z_body_top_mm"] - ch["z_arch_mm"] - 4.0
            sup.append(bm["support"])
            support_kg.append(bsm.mass_kg(support, bm) / 2)      # per side, hubcap included
            info.setdefault("support_structure", {})[tag] = bm["_struct"]
            if bm["hubcap"] is not None:
                # its own printed part (pressed on the stub's end): not a CFD
                # surface, it sits flush inside the wheel's disc
                bm["hubcap"].export(str(out / f"hubcap_{tag}.stl"))
                info.setdefault("hubcap_stls", {})[tag] = str(out / f"hubcap_{tag}.stl")
        support_kg = tuple(support_kg)
    if sup:
        _export_half(trimesh.util.concatenate(sup), out / "supports.stl")
        surfaces.append({"name": "supports", "stl": str(out / "supports.stl"), "rotating": None})
    halo = hg.build_halo(cc.ref_plane_A(x_front_mm) / 1e3, d_halo_mm)["halo_right"]
    _export_half(halo, out / "halo.stl")
    surfaces.append({"name": "halo", "stl": str(out / "halo.stl"), "rotating": None})

    # -- front wing --------------------------------------------------------
    # -- nose cone (printed, ahead of Ref A) --------------------------------
    nm = None
    foam_nose = freeform is not None and getattr(freeform, "nose_foam", False)
    if foam_nose:
        import numpy as np
        # the nose zone is foam, milled with the body (owner 2026-11-28): no printed nose; the front wing holds to the
        # body's own nose ahead of Ref A, and its pylons (if any) are the printed wing mount
        ra = cc.ref_plane_A(x_front_mm) * 1e-3
        b_ = body.slice_plane([ra + 2e-3, 0.0, 0.0], [-1.0, 0.0, 0.0], cap=True)
        nm = {"nose": b_, "_full": b_ + b_.copy().apply_transform(np.diag([1.0, -1.0, 1.0, 1.0]))} if len(b_.faces) else None
        info["foam_nose"] = True
    elif nose is not None:
        if freeform is not None:
            # the free-form car's nose: its own surfaces ahead of Ref A (car/body/freeform.py)
            import freeform as fr
            nm = fr.build_nose(freeform, freeform.x_tip_mm, freeform.x_rear_mm, cc.ref_plane_A(x_front_mm),
                               blend_after_mm=nose.blend_after_ref_a_mm)
            if nm is None:
                raise ValueError("illegal car: the free-form surfaces leave no nose ahead of Ref A")
        else:
            nm = ns.build_nose(nose, x_front_mm, body)
        gates.update(ns.nose_gates(nose, nm, x_front_mm))
        _export_half(nm["nose"], out / "nose.stl")
        surfaces.append({"name": "nose", "stl": str(out / "nose.stl"), "rotating": None})
        parts_mass["nose"] = ns.nose_mass_com(nose, nm["_full"])

    # -- front wing --------------------------------------------------------
    # A nose or body the wing passes through splits its span there (T8.6.1): each outer
    # segment must then be 25 mm, so the tip goes out to the crossing + 25 (Part 1 found
    # 24.23 mm segments on car 41518/1, 2026-11-12, when the gate assumed no crossing).
    if isinstance(front, fwe.FrontWingDesign):
        # the searches' wing (front_wing_engine.py): one piece or two, its own pylons to THIS nose, printed
        # with the nose (nose.stl), and every T8 rule measured
        if nm is None:
            raise ValueError("the generated front wing needs a nose to hold it")
        res = fwe.build(front, nm["nose"], x_front_mm, (x_front_mm, FY, R), body)
        gates.update({f"P4.{k}" if not k.startswith("P4.") else k: v for k, v in res["gates"].items()})
        info["front_wing_engine"] = {k: v for k, v in res["info"].items()}
        if res["pylons"] is not None and foam_nose:          # the printed wing mount on the foam nose
            mount = _printed_minus_body(res["pylons"], body)   # ends on the foam's surface: one closed solid a side
            _export_half(mount, out / "nose.stl")
            surfaces.append({"name": "nose", "stl": str(out / "nose.stl"), "rotating": None})
            parts_mass["nose"] = _full_mass_com([mount], cc.DENSITY_G_CM3[PART_MATERIAL])
        elif res["pylons"] is not None:
            nose_with = _union_closed([nm["nose"], res["pylons"]])
            _export_half(nose_with, out / "nose.stl")
            parts_mass["nose"] = _full_mass_com([nose_with], cc.DENSITY_G_CM3[PART_MATERIAL])
        _export_half(res["wing"], out / "fwing.stl")
        surfaces.append({"name": "fwing", "stl": str(out / "fwing.stl"), "rotating": None})
        parts_mass["fwing"] = _full_mass_com([res["wing"]], cc.DENSITY_G_CM3[PART_MATERIAL])
    else:
        _legacy_front_wing(front, nm, x_front_mm, R, body, out, gates, info, surfaces, parts_mass)

    # -- rear wing ---------------------------------------------------------
    forward = [body, halo] + sup
    import rear_wing_engine as rwe
    if isinstance(rear, rwe.RearWingDesign):
        # the searches' rear wing (rear_wing_engine.py): its struts cast to this body, wherever it is
        wheel_r = (x_rear, RY, R, r.w)
        res = rwe.build(rear, body, x_front_mm, W_mm, wheel_r, clear_of=[halo] + sup + [wheel_meshes["wheelR"]],
                        forward_parts=forward)
        gates.update({f"P4.{k}": v for k, v in res["gates"].items()})
        info["rear_wing_engine"] = res["info"]
        rparts = [res["wing"]] + res["struts"]
        info["rear_wing_slot_mm3"] = float(sum(abs(p.volume) for p in res["plugs"]) * 1e9)
    else:
        rwm = wg.build_rear_wing(rear, x_front_mm, W_mm, body, forward)
        rd = {k: rwm[k] for k in ("wing", "pylon")}
        gates.update(wg.rear_wing_gates(rear, rd, x_front_mm, W_mm, body, forward))
        rparts = list(rd.values())
        info["rear_wing_z_chord_mm"] = rwm["_z_chord_mm"]
        info["body_top_under_rear_wing_mm"] = rwm["_body_top_mm"]
    # the wing and its supports as separate solids in one file: the referee tells them apart (T9.5.1)
    _export_half(trimesh.util.concatenate([_stl_safe(p) if p.is_volume else p for p in rparts]), out / "rwing.stl")
    surfaces.append({"name": "rwing", "stl": str(out / "rwing.stl"), "rotating": None})
    parts_mass["rwing"] = _full_mass_com(rparts, cc.DENSITY_G_CM3[PART_MATERIAL])

    # -- tether guides (T6) -------------------------------------------------
    tg = []
    for x in ((x_front_mm - tether_offsets_mm[0], x_rear + tether_offsets_mm[1]) if tethers else ()):
        m, meta = tether_guide(x)
        # A guide must be held by something (T6.3): where the floor is above the ring,
        # a stalk takes it up into the body.
        ring_top = meta["bottom_mm"] + meta["od_mm"]
        up, _, _ = body.ray.intersects_location([[x / 1e3, 0.002, 0.0]], [[0.0, 0.0, 1.0]],
                                                multiple_hits=False)
        if len(up) and up[0, 2] * 1e3 > ring_top - 0.5:
            stalk = wg.box(x - 1.0, x + 1.0, 0.0, 1.5, ring_top - 0.5, up[0, 2] * 1e3 + 1.0)
            m = trimesh.boolean.union([m, stalk], engine="manifold")
        tg.append(m)
        info.setdefault("tethers", []).append(meta)
    if tg:
        _export_half(trimesh.util.concatenate(tg), out / "tethers.stl")
        surfaces.append({"name": "tethers", "stl": str(out / "tethers.stl"), "rotating": None})
        parts_mass["tethers"] = _full_mass_com(tg, cc.DENSITY_G_CM3[PART_MATERIAL])
        gates["T6.1_front_guide"] = cc.TETHER_FROM_AXLE_MAX - tether_offsets_mm[0]
        gates["T6.1_rear_guide"] = cc.TETHER_FROM_AXLE_MAX - tether_offsets_mm[1]
        gates["T6.2_id_min"] = cc.TETHER_GUIDE_ID - cc.TETHER_ID_MIN
        gates["T6.2_id_max"] = cc.TETHER_ID_MAX - cc.TETHER_GUIDE_ID
        gates["T3.7_tether_clearance"] = cc.TETHER_GUIDE_BOTTOM_Z - cc.TRACK_CLEARANCE_MIN

    fhk = fixed_hardware_kwargs(W_mm, x_front_mm, d_halo_mm, parts_mass, wheel_design,
                                support_kg, rear_face_mm=rear_face_mm)
    info["support_mass_g_each"] = [m * 1e3 for m in support_kg]
    result = {
        "W_mm": W_mm, "x_front_mm": x_front_mm, "d_halo_mm": d_halo_mm,
        "wheel_design": wh.summary(wheel_design),
        "wheel_moi_kg_m2": wh.mean_inertia_kg_m2(wheel_design),
        "extra_surfaces": surfaces,
        "parts_mass_g": {k: v[0] * 1e3 for k, v in parts_mass.items()},
        "fixed_hardware_kwargs": fhk,
        "gates": gates,
        "failed_gates": sorted(k for k, v in gates.items() if v < -1e-9),
        "front_wing": wg.params_dict(front), "rear_wing": wg.params_dict(rear),
        "nose": None if nose is None else wg.params_dict(nose),
        "support": None if support is None else (se.params_dict(support) if isinstance(support, se.SupportDesign)
                                                  else wg.params_dict(support)),
        "info": info,
    }
    (out / "assembly.json").write_text(json.dumps(result, indent=2, default=list))
    return result


CARTRIDGE_LENGTH_MM = 65.0       # Appendix vi, nominal (65.00 -0.5)


LEAD_G_CM3 = 11.34


def ballast_com_mm(asm: dict, lead_g_cm3: float = LEAD_G_CM3) -> list:
    """Centre of the lead in the T1.22 container under the halo, filled from the bottom
    (it was put at the halo's own centre before 2026-11-11, ~10 mm too high)."""
    import trimesh
    from halo_pocket import ballast_capsule_mesh_m
    cap = ballast_capsule_mesh_m((asm["x_front_mm"] - 16.0) / 1e3, asm["d_halo_mm"], above_floor_mm=0.0)
    cap.apply_scale(1e3)
    want = min(asm["ballast_g"] / lead_g_cm3 * 1e3, abs(cap.volume))
    lo, hi = cap.bounds[0, 2], cap.bounds[1, 2]
    fill = cap
    for _ in range(40):                              # the level that holds that volume
        z = 0.5 * (lo + hi)
        fill = trimesh.intersections.slice_mesh_plane(cap, [0, 0, -1], [0, 0, z], cap=True)
        lo, hi = (z, hi) if abs(fill.volume) < want else (lo, z)
    c = fill.center_mass
    return [float(c[0]), 0.0, float(c[2])]


def write_mass_report(asm: dict, path) -> dict:
    """mass_report.json: the car as it is made. Every manufactured file's own volume x its
    material at its own centre of mass, the wheels with their bearings on their axles, the
    halo and helmet, any ballast in its container; then the cartridge. Masses in g, centres
    in mm (x from the nose, z up from the track; y = 0, the car is symmetric)."""
    import json
    import trimesh
    rows = {}
    for k, v in asm["joints"]["manufactured"].items():
        if "file" not in v:
            continue
        m = trimesh.load(v["file"], force="mesh")
        c = m.center_mass if m.is_volume else m.bounds.mean(axis=0)
        rows[k] = (float(v["g"]), [float(c[0]), 0.0, float(c[2])])
    wd, fhk = asm["wheel_design"], asm["fixed_hardware_kwargs"]
    for tag, x in (("front", asm["x_front_mm"]), ("rear", asm["x_front_mm"] + asm["W_mm"])):
        R = wd[tag]["R_mm"]
        rows[f"wheels_{tag}_x2"] = (2.0 * wd[tag]["mass_g"], [x, 0.0, R - SINK_MM])
    hc = fhk["halo_com"]
    rows["halo_and_helmet"] = (fhk["halo_mass_kg"] * 1e3, [hc[0] * 1e3, 0.0, hc[2] * 1e3])
    if asm.get("ballast_g"):
        rows["ballast"] = (float(asm["ballast_g"]), ballast_com_mm(asm))
    car_g = sum(g for g, _c in rows.values())
    com = [sum(g * c[i] for g, c in rows.values()) / car_g for i in range(3)]
    # The cartridge's EMPTY shell (Appendix vi: 28.9-29.4 g full, less the 8 g charge): Part 2's
    # race model adds the gas itself.
    cart_g = fhk["co2_cartridge_mass_kg"] * 1e3
    cc_ = [v * 1e3 for v in fhk["co2_cartridge_com"]]
    tot = car_g + cart_g
    com_w = [(car_g * com[i] + cart_g * cc_[i]) / tot for i in range(3)]
    jr = asm["joints"]
    guides = asm.get("info", {}).get("tethers", [])
    report = {
        "car_weight_kg": tot * 1e-3,                       # car + empty cartridge shell, no CO2
        "com_mm": com_w,                                   # [x from the nose tip, y, z], same, no CO2
        "wheel_moi_kg_m2": asm["wheel_moi_kg_m2"],         # one wheel; the race model takes 4
        "mass_g_without_cartridge": car_g,                # T3.6, the scrutineer's scale
        "com_mm_without_cartridge": com,
        "cartridge_shell_g": cart_g,
        # A 65.0 mm cartridge (Appendix vi) home against the chamber's closed end.
        "cartridge_rear_centre_x_mm": (jr["bore_end_x_mm"] + CARTRIDGE_LENGTH_MM
                                       if "bore_end_x_mm" in jr else None),
        "cartridge_axis_z_mm": jr.get("bore_axis_z_mm"),
        # The tether line through the guides' centres (T6), as built.
        "tether_line_z_mm": (guides[0]["bottom_mm"] + guides[0]["od_mm"] / 2.0) if guides else None,
        "parts": {k: {"g": g, "com_mm": c, **{f: asm["joints"]["manufactured"][k][f]
                                              for f in ("cm3", "g_cm3", "density_source", "weighed")
                                              if f in asm["joints"]["manufactured"].get(k, {})}}
                  for k, (g, c) in rows.items()},
        "placeholders": "printed-part masses are volume x datasheet density (each row names its source); "
                        "none is weighed yet, nor the wheels: replace with the team's real print weights",
    }
    Path(path).write_text(json.dumps(report, indent=1))
    return report


def extra_surfaces_from(assembly_json: str) -> tuple:
    """Load the Part 2 extra_surfaces tuple back from assembly.json."""
    d = json.loads(Path(assembly_json).read_text())
    out = []
    for s in d["extra_surfaces"]:
        rot = s["rotating"]
        if rot:
            rot = {"origin": tuple(rot["origin"]), "axis": tuple(rot["axis"]),
                   "omega": rot["omega"]}
        out.append({"name": s["name"], "stl": s["stl"], "rotating": rot})
    return tuple(out)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--W", type=float, default=120.3)
    ap.add_argument("--x-front", type=float, default=46.0)
    ap.add_argument("--d-halo", type=float, default=43.72)
    ap.add_argument("--body", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--wheels", default="team_stl")
    a = ap.parse_args()
    res = build(a.W, a.x_front, a.d_halo, a.body, a.out, a.wheels)
    print(json.dumps({k: res[k] for k in ("parts_mass_g", "wheel_moi_kg_m2", "failed_gates")},
                     indent=2))
