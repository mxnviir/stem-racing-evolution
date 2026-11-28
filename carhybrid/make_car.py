"""Build one whole car with the pinned evolution builder and adapt it for OpenFOAM (carhybrid/DESIGN.md P0).

    python carhybrid/make_car.py <evolution checkout> <layout: "team" | genes json | Lxxxx (search/layout/results)>
                                 <stage-2 numbers u: "zero" | json list> <out folder> [--grid gu,gv,gw,gb]
                                 [--wing '{"rw_rise": 6, "rw_camber": -0.03, "rwsup_root": 0.9, ...}']

Writes <out>/evo/ (the evolution car folder: summary.json, layout_car.json, body_half.stl, parts/, body_full.stl,
wheels), <out>/rec.json (the build record: genes, decoded layout incl. the free-form seed, legal, mass, cfd report) and
<out>/cfd/ (the OpenFOAM half-car folder, carhybrid/adapt.py). Exit code 0 only if the car is legal and adapted.
The stage-2 numbers move the free-form body and nose exactly as stage 2 does (simscale_shape.ShapeCampaign.job:
freeform.from_unit(seed, u, AMP_MM, *grid)); the seed is the stage-1 car's own organic free-form surfaces.
grid: [--grid gu,gv,gw,gb] (carhybrid/grad.py's), else from the length of u: 45 = stage 2's COARSE (5,2,3,4), the
seed's own grid's count (team car: 720 = 24,9,12,9, carhybrid/FINE_ADJOINT.md) = that grid; any other length is
refused (exit 2: u is never truncated, so a fine-grid u can not silently build a coarse car).
--wing: the rear wing's shape and strut numbers by their layout_params.json names (rw_*: chord, t, half_span, x_le,
lift and evo_patches/rear_wing_shape.py's camber, camber_pos, rise, rise_p, sweep, sweep_p, tip_chord_ratio,
tip_t_ratio; rwsup_*: the strut's 13 numbers), over the layout's own (rec.json records them as "wing").
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


def grid_of(fr, SS, seed, n_u, grid=None):
    """The free-form grid u is on: explicit, else stage 2's COARSE or the seed's own grid, by the count."""
    if grid:
        g = tuple(int(v) for v in str(grid).split(","))
        return g if len(g) == 4 else g + (None,)
    for g in (tuple(SS.COARSE), (seed.nu, seed.nv, seed.nw, seed.nv)):
        if fr.n_genes(*g) == n_u:
            return g
    raise SystemExit(f"u has {n_u} numbers: neither stage 2's {fr.n_genes(*SS.COARSE)} nor the seed grid's; pass --grid")


RW_KEYS = {"rw_chord": "chord_mm", "rw_t": "t_mm", "rw_half_span": "half_span_mm", "rw_x_le": "x_le_mm",
           "rw_lift": "lift_mm", "rw_camber": "camber", "rw_camber_pos": "camber_pos", "rw_rise": "rise_mm",
           "rw_rise_p": "rise_p", "rw_sweep": "sweep_mm", "rw_sweep_p": "sweep_p",
           "rw_tip_chord_ratio": "tip_chord_ratio", "rw_tip_t_ratio": "tip_t_ratio"}
SUP_KEYS = {"out": "aim_out_deg", "fwd": "aim_fwd_deg", "leave": "leave", "arrive": "arrive", "chord": "chord_mm",
            "thick": "thick_mm", "flare": "flare", "round": "round", "mid_x": "mid_x", "mid_y": "mid_y",
            "mid_z": "mid_z", "swell": "swell", "root": "root"}


def rear_wing_of(d, over):
    """The layout's rear wing (decoded) with the --wing numbers over it."""
    rw = json.loads(json.dumps(d["rwing"]))
    sup = rw.setdefault("support", {})
    for k, v in (over or {}).items():
        if k in RW_KEYS:
            rw[RW_KEYS[k]] = float(v)
        elif k.startswith("rwsup_") and k[6:] in SUP_KEYS:
            sup[SUP_KEYS[k[6:]]] = float(v)
        else:
            raise SystemExit(f"--wing: unknown number {k!r} (rw_* or rwsup_*: {sorted(RW_KEYS) + ['rwsup_' + q for q in SUP_KEYS]})")
    return rw


def main(evo, layout, u_s, out, grid=None, wing=None):
    evo, out = Path(evo).resolve(), Path(out).resolve()
    sys.path.insert(0, str(evo)); import repo_paths; repo_paths.add()
    import layout_car as LC, layout_genes as LG, freeform as fr, simscale_shape as SS
    x = genes_of(evo, layout)
    d = LG.decode(x)
    freeform, g = None, None
    if u_s != "zero":
        u = np.asarray(json.loads(u_s), float)
        seed = fr.Freeform.from_dict(d["freeform"])
        g = grid_of(fr, SS, seed, len(u), grid)
        if fr.n_genes(*g) != len(u):
            print(f"u has {len(u)} numbers, grid {g} {fr.n_genes(*g)}", file=sys.stderr)
            return 2
        freeform = fr.from_unit(seed, u, SS.AMP_MM, *g).as_dict()
    car = out / "evo"
    over = json.loads(wing) if wing else None
    rec = LC.build(x, car, freeform=freeform, rwing=rear_wing_of(d, over) if over else None)
    rec["wing"] = over
    rec["u"] = None if u_s == "zero" else json.loads(u_s)
    rec["grid"] = None if g is None else [v for v in g]
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
    a = sys.argv[1:]
    sys.exit(main(*a[:4], grid=a[a.index("--grid") + 1] if "--grid" in a else None,
                  wing=a[a.index("--wing") + 1] if "--wing" in a else None))
