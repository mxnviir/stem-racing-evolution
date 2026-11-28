"""The DAFoam drag (and lift) gradient -> the race-time gradient in the REAR WING's and its STRUTS' numbers
(carhybrid/METHOD.md section 3; owner 2026-11-28: the wing and strut numbers in our hybrid, with a DAFoam-based
gradient).

    python carhybrid/wing_grad.py <evolution checkout> <car folder> <dafoam_sens.npz>[,<frame 2 npz>,...] <out json>
         [--names rw_chord,rw_rise,...] [--delta 0.05] [--d20 N] [--l20 N] [--one-sided]

car folder: make_car.py's <out>: <out>/rec.json and <out>/evo/parts/rwing_ctx.pkl (evo_patches/rear_wing_shape.py
W8: the wing engine's inputs: the design, the body, the parts it keeps clear of).

The numbers s_j in [-1, 1] are each searched parameter's place in its range (search/layout/layout_params.json lo..hi):
theta_j = mid_j + s_j * half_j. Default: every continuous rear-wing and strut number of mode shape or evolve
(DEFAULT_NAMES; rwsup_root is not here: it picks the support's topology, which is evolution's and the layout
search's to change, not a gradient's).

How: for each number the wing and its struts are REBUILT by the engine itself (rear_wing_engine.build, with every rule
it applies: placement, T7.9, T9.4-T9.7, the strict span, the struts' stress sizing) at theta_j +- delta * half_j. Each
DAFoam surface point on the wing or a strut moves along the outward normal at V_ij = d sd_i / d s_j (mm per unit),
sd = the signed distance to the rebuilt parts (positive inside; union = max over parts). Then exactly as grad.py:
  dD/ds_j = 2 sum_i (G_i . n_i) V_ij (whole car, x 1e-3 for mm -> m), dL/ds_j the same from dL_dXs,
  dm/ds_j = d(parts volume x 2 x ABS density)/ds_j,
  dT/ds_j = dT/dD dD/ds_j + dT/dL dL/ds_j + dT/dm dm/ds_j (the race model's partials at this car, as grad.py).
A number whose +- builds change the wing's topology (a support moved to another root, the wing moved forward for its
supports, supported / unsupported) is a branch change: no derivative (null, listed in "branch_change"); the
trust-region loop's measured steps still move it, and it is reported.
"""
from __future__ import annotations

import dataclasses, json, pickle, sys, time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

DEFAULT_NAMES = ("rw_chord", "rw_t", "rw_half_span", "rw_x_le", "rw_lift", "rw_camber", "rw_camber_pos", "rw_rise",
                 "rw_rise_p", "rw_sweep", "rw_sweep_p", "rw_tip_chord_ratio", "rw_tip_t_ratio",
                 "rwsup_out", "rwsup_fwd", "rwsup_leave", "rwsup_arrive", "rwsup_chord", "rwsup_thick", "rwsup_flare",
                 "rwsup_round", "rwsup_mid_x", "rwsup_mid_y", "rwsup_mid_z", "rwsup_swell")
ON_WING_MM = 0.3          # a DAFoam point within this of the wing or a strut (and nearer it than the body) is on it


def ranges(evo):
    P = json.loads((Path(evo) / "search" / "layout" / "layout_params.json").read_text())["params"]
    return {q["name"]: (float(q["lo"]), float(q["hi"])) for q in P}


def field_of(name):
    """(where, field) on the RearWingDesign for a layout_params name: ("rw", field) or ("sup", field)."""
    from make_car import RW_KEYS, SUP_KEYS
    if name in RW_KEYS:
        return "rw", RW_KEYS[name]
    if name.startswith("rwsup_") and name[6:] in SUP_KEYS:
        return "sup", SUP_KEYS[name[6:]]
    raise SystemExit(f"wing_grad: unknown number {name!r}")


def get_theta(design, name):
    w, f = field_of(name)
    return float(getattr(design.support if w == "sup" else design, f))


def with_theta(design, name, value):
    w, f = field_of(name)
    if w == "sup":
        return dataclasses.replace(design, support=dataclasses.replace(design.support, **{f: float(value)}))
    return dataclasses.replace(design, **{f: float(value)})


def topology(info):
    st = info.get("struts") or []
    return (bool(info.get("supported")), tuple(round(float(p.get("root", 1.0)), 2) for p in st),
            round(float(info.get("moved_forward_for_supports_mm", 0.0)), 3))


class Built:
    """One engine build: its parts (mm, right halves) and the signed distance to their union."""

    def __init__(self, rwe, design, ctx):
        res = rwe.build(design, ctx["body"], ctx["x_front_mm"], ctx["W_mm"], ctx["wheel_rear"],
                        clear_of=ctx["clear_of"], forward_parts=ctx["forward_parts"])
        self.info, self.gates = res["info"], res["gates"]
        self.parts = []
        for m in [res["wing"]] + list(res["struts"]):
            q = m.copy(); q.apply_scale(1e3); self.parts.append(q)
        self.volume_mm3 = float(sum(abs(p.volume) for p in self.parts))
        self.topo = topology(self.info)

    def sd(self, P):
        import trimesh
        return np.max([trimesh.proximity.signed_distance(p, P) for p in self.parts], axis=0)


def main(evo, car, sens_p, out_p, names=None, delta=0.05, d20=None, l20=None, one_sided=False, quiet=False):
    t_start = time.time()
    import grad as G
    fr, SL, SS = G.evo_import(evo)
    import rear_wing_engine as rwe
    import trimesh
    import component_contract as cc
    import assembly as AS
    car = Path(car)
    rec = json.loads((car / "rec.json").read_text())
    ctx = pickle.loads((car / "evo" / "parts" / "rwing_ctx.pkl").read_bytes())
    design0 = ctx["design"]
    R = ranges(evo)
    names = list(names or DEFAULT_NAMES)
    rho = cc.DENSITY_G_CM3[AS.PART_MATERIAL] * 1e-3          # g per mm^3

    base = Built(rwe, design0, ctx)
    body = ctx["body"].copy(); body.apply_scale(1e3)
    paths = [p for p in str(sens_p).split(",") if p]
    frames = []
    for p in paths:
        S = np.load(p)
        X = S["Xs"] * 1e3
        Gs = [S["dD_dXs"]] + ([S["dL_dXs"]] if "dL_dXs" in S.files else [])
        dw = np.min([trimesh.proximity.closest_point(q, X)[1] for q in base.parts], axis=0)
        db = trimesh.proximity.closest_point(body, X)[1]
        on = (dw < ON_WING_MM) & (dw < db)
        Q = X[on]
        # outward normals of the parts at the points (the nearest part's face normal)
        best, nrm = np.full(len(Q), np.inf), np.zeros((len(Q), 3))
        for q in base.parts:
            _, dist, tri = trimesh.proximity.closest_point(q, Q)
            k = dist < best
            best[k], nrm[k] = dist[k], q.face_normals[tri[k]]
        Wn = np.stack([(g[on] * nrm).sum(1) for g in Gs], 1)       # N/m along the outward normal
        frames.append({"npz": str(p), "Q": Q, "Wn": Wn, "n_on_wing": int(on.sum()), "n_total": int(len(X)),
                       "D_N": float(S["D_N"]) if "D_N" in S.files else None,
                       "L_N": float(S["L_N"]) if "L_N" in S.files else None, "has_lift": len(Gs) > 1})
    have_lift = all(f["has_lift"] for f in frames)
    sd0 = [base.sd(f["Q"]) for f in frames] if one_sided else None

    rows, dD, dL, dm, branch = [], [], [], [], []
    for name in names:
        lo, hi = R[name]
        half = 0.5 * (hi - lo)
        th0 = get_theta(design0, name)
        s0 = (th0 - 0.5 * (lo + hi)) / half if half > 0 else 0.0
        hp = min(delta, 1.0 - s0) if not one_sided else min(delta, 1.0 - s0)
        hm = 0.0 if one_sided else min(delta, 1.0 + s0)
        if hp <= 1e-6 and hm <= 1e-6:
            rows.append({"name": name, "s": s0, "skipped": "range of zero width"}); dD.append(None); dL.append(None); dm.append(None)
            continue
        if hp <= 1e-6:                                      # at the top of its range: one-sided downwards
            hp, hm = 0.0, min(delta, 1.0 + s0)
        t1 = time.time()
        bp = Built(rwe, with_theta(design0, name, th0 + hp * half), ctx) if hp > 0 else base
        bm = Built(rwe, with_theta(design0, name, th0 - hm * half), ctx) if hm > 0 else base
        row = {"name": name, "theta": th0, "s": s0, "h_plus": hp, "h_minus": hm, "build_s": round(time.time() - t1, 1)}
        if bp.topo != base.topo or bm.topo != base.topo:
            row["branch_change"] = {"base": base.topo, "plus": bp.topo, "minus": bm.topo}
            branch.append(name); dD.append(None); dL.append(None); dm.append(None); rows.append(row)
            continue
        hs = hp + hm
        rD, rL = 0.0, 0.0
        for k, f in enumerate(frames):
            sp_ = bp.sd(f["Q"]) if hp > 0 else sd0[k]
            sm_ = bm.sd(f["Q"]) if hm > 0 else (sd0[k] if sd0 is not None else base.sd(f["Q"]))
            V = (sp_ - sm_) / hs                             # mm per unit s, outward
            rD += 2.0 * float(f["Wn"][:, 0] @ V) * 1e-3
            if have_lift:
                rL += 2.0 * float(f["Wn"][:, 1] @ V) * 1e-3
        dD.append(rD / len(frames)); dL.append(rL / len(frames) if have_lift else 0.0)
        dm.append(2.0 * rho * (bp.volume_mm3 - bm.volume_mm3) / hs)
        row.update(dD_ds_N=dD[-1], dL_ds_N=dL[-1], dm_ds_g=dm[-1])
        rows.append(row)
        if not quiet:
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)

    # the race model's partials, as grad.py
    m = rec["mass"]
    lead = float(m.get("lead_g", 0.0) or 0.0)
    D20 = float(d20) if d20 is not None else (float(rec["D20_N"]) if rec.get("D20_N") is not None else
          (2.0 * float(np.mean([f["D_N"] for f in frames])) if all(f["D_N"] is not None for f in frames) else 0.36))
    L20 = float(l20) if l20 is not None else (float(rec["L20_N"]) if rec.get("L20_N") is not None else
          (2.0 * float(np.mean([f["L_N"] for f in frames])) if all(f["L_N"] is not None for f in frames) else 0.0))
    RT = G.race_fn(SL, rec)
    T0 = RT(D20, L20)
    h = 1e-3
    dT_dD = (RT(D20 + h, L20) - RT(D20 - h, L20)) / (2 * h)
    dT_dL = (RT(D20, L20 + h) - RT(D20, L20 - h)) / (2 * h) if have_lift else 0.0
    dT_dm = (RT(D20, L20, 0.1) - T0) / 0.1 if lead <= 0 else 0.0
    dT = [None if a is None else dT_dD * a + dT_dL * b + dT_dm * c for a, b, c in zip(dD, dL, dm)]
    g = np.array([0.0 if v is None else v for v in dT])
    gN = float(np.linalg.norm(g))
    out = {"names": names, "s": [r.get("s") for r in rows], "theta": [r.get("theta") for r in rows],
           "ranges": {n: R[n] for n in names}, "delta": delta, "one_sided": one_sided,
           "n_frames": len(frames), "n_points_on_wing": int(np.mean([f["n_on_wing"] for f in frames])),
           "base_topology": base.topo, "branch_change": branch,
           "D20_N_used": D20, "L20_N_used": L20, "T0_s": T0, "dT_dD_s_per_N": dT_dD, "dT_dL_s_per_N": dT_dL,
           "dT_dm_s_per_g": dT_dm, "lead_g": lead,
           "dD_ds_N": dD, "dL_ds_N": dL if have_lift else None, "dm_ds_g": dm, "dT_ds_s": dT,
           "norm_dT_ds_s": gN, "direction": (-g / gN).tolist() if gN > 0 else None, "rows": rows,
           "run_s": round(time.time() - t_start, 1)}
    Path(out_p).write_text(json.dumps(out, indent=1, default=float))
    if not quiet:
        print(json.dumps({k: out[k] for k in ("n_points_on_wing", "base_topology", "branch_change", "T0_s",
                                              "norm_dT_ds_s", "run_s")}, default=float))
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) < 4 or a[0] in ("-h", "--help"):
        print(__doc__); sys.exit(0)
    opt = lambda k: a[a.index(k) + 1] if k in a else None
    main(*a[:4], names=(opt("--names") or "").split(",") if opt("--names") else None,
         delta=float(opt("--delta") or 0.05), d20=opt("--d20"), l20=opt("--l20"), one_sided="--one-sided" in a)
