"""The whole-car layout search on SimScale.

The numbers that place the car's hard parts -- wheelbase, rear axle -> rear face, halo -> rear face,
cargo position and orientation, nose length / profile / tip height, front wing height, rear wing
position (layout_params.json, mode "evolve") -- are searched with a surrogate, not by blind evolution:

  1. the first batches are space-filling (Sobol) over the legal set;
  2. after each batch a Gaussian process (one length scale per number) is fitted to the race times;
  3. EVERY legal layout on the 1 mm grid (32 424 of them), each with several settings of the other
     numbers, is scored on the model (expected improvement); the 8 most promising are built and sent
     to CFD (the chosen ones are added with their predicted time before the next is picked, so a batch
     does not pile onto one spot).
  So all legal combinations are gone through every round, on the model; CFD time is only spent on 8.

Each car is built whole by layout_car.py (Parts 1-5, the referee, one closed solid) and timed by the
race model with ITS OWN wheelbase, thrust point, mass and COM. At equal drag the layout moves the race
time by 0.15 ms over the whole legal range (the car is always ballasted to 48.2 g), 1 % of drag by
0.72 ms: the layout matters through drag, which is why it needs CFD at all.

A process of its own, in its own folder (runs/layout), next to the sidepod campaign: the same
keys (key state seeded from the sidepod run, so a retired key stays retired), the same SimScale setup
(half car, 20 m/s, k-omega SST, 1000 iterations, manual mesh max 30 / min 0.8 mm, road 0.15 mm under
the wheels). The whole car is faceted (every part of it is generated); the reference (every number at
its default) goes through the same route at two facet tolerances: like for like, and their difference
is the noise a gain has to beat.

    python simscale_layout.py init runs/layout --budget 48
    python simscale_layout.py ref  runs/layout            # the reference, twice
    python simscale_layout.py run  runs/layout            # ticks every 60 s
    python simscale_layout.py status runs/layout
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
import layout_genes as LG  # noqa: E402
import layout_car as LC  # noqa: E402
import simscale_step_eval as SE  # noqa: E402

CAR_FACES = (24000, 32000, 40000, 28000, 20000, 48000)     # whole-car triangle counts tried in turn (16000 meshed too coarse, 2026-11-22)
SE.FACES = CAR_FACES                                         # the evaluator's retry ladder (this process only)
BOX_M = ((-0.01, 0.0, -0.001), (0.30, 0.06, 0.07))           # car <= 264 mm with the cartridge; rear wing <= 65 mm
WHEEL_FILLET_MM = 1.0                                        # rim rounding of the CFD wheels
# The mesher's small-feature suppression for a faceted whole car: every triangle edge is a CAD edge, and
# at SimScale's 0.05 mm default the 0.15-0.5 mm edges made it fail ("mesh too coarse", "size smaller than
# tolerance") on 5 of 9 cars and on the reference at every triangle count. At 0.8 mm (the smallest cell)
# the reference meshed first time (5.0 M cells, 2026-11-22).
SMALL_FEATURE_MM = 0.8


def cfd_setup(faces, pitch_mm, supports=None):
    """What a result depends on besides the car: compared only within one setup."""
    return {"small_feature_mm": SMALL_FEATURE_MM, "wheel_fillet_mm": WHEEL_FILLET_MM, "faces": int(faces),
            "wrap_pitch_mm": pitch_mm, "cartridge": CARTRIDGE, "supports": SUPPORTS if supports is None else supports}


CARTRIDGE = "official"          # the official cartridge's shape (layout_car.cartridge); "cylinder" before 2026-11-22 17:50
SUPPORTS = "engine1"            # the wheel supports in the CFD car (support_engine.py); "none": a car built before them


def setup_key(su):
    """The part of a CFD setup that shifts every car alike: results compare only within one key."""
    su = su or {}
    return (su.get("small_feature_mm"), su.get("wheel_fillet_mm"), su.get("cartridge", "cylinder"),
            su.get("supports", "none"))


CURRENT_SETUP = (SMALL_FEATURE_MM, WHEEL_FILLET_MM, CARTRIDGE, SUPPORTS)
# Setups that differ from the current one only in parts that every car carries alike (the cartridge, the
# supports) are brought over by the difference between their reference runs (data()); the mesh settings
# are not (they move cars differently).
OFFSET_OK = 2                   # setup_key[:OFFSET_OK] must match
BUILD_WORKERS = 3                                            # one Parts 1-5 build peaks at ~5.5 GB


DEV_P99_MAX_MM, DEV_MAX_MM = 0.25, 0.6     # the faceted (timed) surface's distance from the built car


def surface_deviation(a, b, n=20000) -> dict:
    """Two-way distance (mm) between two surfaces given in metres, from n points sampled on each."""
    import trimesh
    out = []
    for P, Q in ((a, b), (b, a)):
        pts, _ = trimesh.sample.sample_surface(P, n, seed=0)
        _, d, _ = trimesh.proximity.closest_point(Q, pts)
        out.append(d)
    d = np.concatenate(out) * 1e3
    return {"mean_mm": round(float(d.mean()), 4), "p99_mm": round(float(np.percentile(d, 99)), 4), "max_mm": round(float(d.max()), 4)}


def air_volume(car_mm, wheel_boxes, out, domain=None):
    """The half wind tunnel with the car cut out, as STEP (mm), faces named for the boundary conditions.
    The car is the faceted solid (it crosses the symmetry plane, so the cut there is clean); the wheels are
    exact cylinders the size and place of the car's own wheel files: wheels given as triangles made the
    mesher fail on their faces (both reference cars, 2026-11-22), and the sidepod route's wheels are exact too."""
    import simscale_step_domain as SD
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.OCP.collections import List_TopoDS_Shape
    from OCP.TopAbs import TopAbs_SOLID
    from OCP.TopExp import TopExp_Explorer
    from OCP.TopoDS import TopoDS
    from OCP.gp import gp_Ax2, gp_Dir, gp_Pnt
    from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
    from OCP.TopAbs import TopAbs_EDGE
    car_s = SD.stl_solid(car_mm)
    wheel_s, zmin, wvol = {}, [], 0.0
    for name, b in wheel_boxes.items():
        (x0, y0, z0), (x1, y1, _) = b
        R = 0.5 * (x1 - x0)
        cyl = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(0.5 * (x0 + x1), y0, z0 + R), gp_Dir(0, 1, 0)), R, y1 - y0).Shape()
        # rounded rims, as the real wheels: a sharp rim 0.15 mm over the road could not be meshed at the
        # 0.8 mm floor ("mesh too coarse" on the rim edge, both reference cars, 2026-11-22)
        fil = BRepFilletAPI_MakeFillet(cyl)
        ex = TopExp_Explorer(cyl, TopAbs_EDGE)
        while ex.More():
            fil.Add(WHEEL_FILLET_MM, TopoDS.Edge_s(ex.Current()) if hasattr(TopoDS, "Edge_s") else TopoDS.Edge(ex.Current()))
            ex.Next()
        fil.Build()
        if not fil.IsDone():
            raise SE.GeometryError(f"wheel fillet failed ({name})")
        wheel_s[name] = fil.Shape()
        zmin.append(z0); wvol += np.pi * R * R * (y1 - y0)
    (x0, y0, _), (x1, y1, z1) = domain or SD.DOMAIN_MM
    z0 = min(zmin) - SD.ROAD_GAP_MM
    dom = ((x0, y0, z0), (x1, y1, z1))
    box = BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape()
    cut = BRepAlgoAPI_Cut()
    args, tools = List_TopoDS_Shape(), List_TopoDS_Shape()
    args.Append(box)
    for sh in [car_s, *wheel_s.values()]:
        tools.Append(sh)
    cut.SetArguments(args); cut.SetTools(tools); cut.SetFuzzyValue(1e-5); cut.Build()
    if not cut.IsDone():
        raise SE.GeometryError("air-volume Boolean failed")
    fluid = cut.Shape()
    sols = []
    ex = TopExp_Explorer(fluid, TopAbs_SOLID)
    while ex.More():
        sols.append(TopoDS.Solid(ex.Current())); ex.Next()
    sealed = 0.0
    if len(sols) > 1:                                               # air sealed inside the car: no flow reaches it
        vols = [SD.volume(x) for x in sols]
        k = int(np.argmax(vols)); sealed = float(sum(vols) - vols[k]); fluid = sols[k]
    removed = (x1 - x0) * (y1 - y0) * (z1 - z0) - SD.volume(fluid) - sealed
    labels = SD.label_faces(fluid, cut, wheel_s, dom)
    SD.write_step(fluid, labels, out)
    counts = {}
    for n in labels.values():
        counts[n] = counts.get(n, 0) + 1
    return {"removed_share": removed / (0.5 * car_mm.volume + wvol), "sealed_mm3": sealed, "faces": counts}


class LayoutEvaluator(SE.SimScaleStepEvaluator):
    """The sidepod evaluator with the whole generated car as the geometry and per-car wheel centres."""

    def _build_step(self, cid, faces):
        import trimesh
        import simscale_step_domain as SD
        folder = Path(self.st["cases"][cid]["folder"])
        import layout_wrap as LW
        full = trimesh.load(folder / "body_full.stl", force="mesh")
        try:
            w = trimesh.load(folder / "body_wrap_mm.stl", force="mesh")    # the rebuilt surface (layout_car.export_cfd)
            dec = LW.reduce(w, faces)
        except Exception as e:                                       # noqa: BLE001
            raise SE.GeometryError(f"surface reduction: {e!r}") from e
        q = LW.quality(dec)
        q["self_intersecting_faces"] = LW.self_intersections(dec)
        if not q["closed"] or q["needles"] or q["self_intersecting_faces"] or dec.volume <= 0:
            raise SE.GeometryError(f"faceted car not clean: {q}")
        dec.apply_scale(1e-3)
        if abs(dec.volume / full.volume - 1) > 0.01:
            raise SE.GeometryError(f"faceting changed the volume by {dec.volume / full.volume - 1:+.2%}")
        # how far the timed surface is from the built one, mm (from the timed surface: the built one also
        # carries the hair-thin internal slits the rebuild is there to remove)
        pts, _ = trimesh.sample.sample_surface(dec, 20000, seed=0)
        _, d, _ = trimesh.proximity.closest_point(full, pts)
        dev = {"p99_mm": round(float(np.percentile(d, 99)) * 1e3, 3), "max_mm": round(float(d.max()) * 1e3, 3), **q}
        self.st["cases"][cid].setdefault("faceting", {})[str(faces)] = dev
        if dev["p99_mm"] > DEV_P99_MAX_MM or dev["max_mm"] > DEV_MAX_MM:
            raise SE.GeometryError(f"faceted car too far from the built one: {dev}")
        out = folder / f"domain_{faces}.step"
        dec.apply_scale(1e3)                                         # mm for the CAD kernel
        boxes = {}
        for wname in ("wheelF", "wheelR"):
            b = trimesh.load(folder / f"{wname}.stl", force="mesh").bounds * 1e3
            boxes[wname] = b
        try:
            rep = air_volume(dec, boxes, out)
        except SE.GeometryError:
            raise
        except Exception as e:                                       # noqa: BLE001
            raise SE.GeometryError(f"air volume build: {e!r}") from e
        self.st["cases"][cid]["faceting"][str(faces)].update(removed_share=round(rep["removed_share"], 4), sealed_air_mm3=round(rep["sealed_mm3"], 1))
        cfd = json.loads((folder / "layout_car.json").read_text()).get("cfd", {})
        # the car's own record says whether its solid carries the supports (a car built before them does not)
        self.st["cases"][cid]["setup"] = cfd_setup(faces, cfd.get("wrap_pitch_mm", 0.3),
                                                   SUPPORTS if cfd.get("supports_mm3") else "none")
        if not 0.97 < rep["removed_share"] < 1.03:
            raise SE.GeometryError(f"air volume: removed {rep['removed_share']:.3f} of the car's volume")
        return out

    def _assign(self, cid):
        ok = super()._assign(cid)
        c = self.st["cases"][cid]
        if ok and c.get("attempt0"):                                 # a case that starts further down the face ladder
            c["attempt"] = c["attempt0"]
        return ok

    def submit(self, cases):
        super().submit(cases)
        for c in cases:
            if c.get("attempt0"):
                self.st["cases"][c["id"]]["attempt0"] = c["attempt0"]
        self._save()

    def _advance(self, cid):
        c = self.st["cases"][cid]
        self._wheels = {k: tuple(v) for k, v in json.loads((Path(c["folder"]) / "wheels.json").read_text()).items()}
        return super()._advance(cid)


def prune(folder):
    """A finished car keeps a 60k-face copy of its CFD surface (body_full.stl, for pictures); the 840k-face
    rebuild and the STEP air volumes go (130 MB a car)."""
    import trimesh
    import layout_wrap as LW
    folder = Path(folder)
    wp = folder / "body_wrap_mm.stl"
    try:
        if wp.exists():
            d = LW.reduce(trimesh.load(wp, force="mesh"), 60000); d.apply_scale(1e-3)
            d.export(folder / "body_full.stl")
            wp.unlink()
        for f in list(folder.glob("domain_*.step")) + list(folder.glob("domain_*.zip")):
            f.unlink()
    except Exception as e:                                           # noqa: BLE001  (tidying must not stop the search)
        print("prune", folder.name, repr(e)[:200], flush=True)


SWITCHES = ("cargo_flip", "fw_pieces")   # genes that are a choice (below 0.5 / above), not a dial


def race_model(d):
    """objective.RaceModel for this car's own layout. RaceModel works in the team STEP's mirrored
    frame (x decreasing to the rear); x -> -x maps the car frame (x from the nose tip) onto it."""
    from objective import RaceModel
    rear_face = d["x_front"] + d["W"] + d["rear"]
    spec = SimpleNamespace(x_nose_mm=0.0, x_front_axle=-d["x_front"], x_rear_axle=-(d["x_front"] + d["W"]),
                           tail_x=-rear_face, cartridge_out_mm=LC.CARTRIDGE_OUT_MM, thrust_z=LC.CHAMBER_Z_MM)
    return RaceModel(spec)


def race_time(rec, D20, L20):
    m = rec["mass"]
    T, _ = race_model(rec["decoded"])(D20, L20, {"race_mass_g": m["race_mass_g"], "com_z_mm": m["com_z_mm"],
                                                  "com_x_car_mm": m["com_x_car_mm"]})
    return T


def timed_rows(rows, refs=False):
    """([(row, T on the current setup)], offsets ms, whether the current setup's reference is timed) of every
    car timed so far (the reference runs too if refs): a setup that differs from the current one only in the
    cartridge or the supports (parts every car carries) is brought over by the difference between its
    reference and the current one; any other setup, or one whose reference has no time yet, is left out.
    With the supports the shift is only approximately the same for every layout (rerun measures how far)."""
    ref = {}
    for row in rows:
        if row.get("kind") == "ref" and row.get("T") is not None:
            ref.setdefault(setup_key(row.get("setup")), []).append(row["T"])
    ref = {k: float(np.mean(v)) for k, v in ref.items()}
    off_ref = {k: ref[CURRENT_SETUP] - ref[k] for k in ref
               if CURRENT_SETUP in ref and k[:OFFSET_OK] == CURRENT_SETUP[:OFFSET_OK]}
    # Cars timed on an earlier setup and again on the current one (rerun): what the reference shift misses
    # for the supports. On 2026-11-22 the team layout's shift was 2.65 ms and three re-timed cars' 5.4-6.3:
    # the supports cost the team layout less than a typical car. The mean miss is added to every setup
    # without the supports.
    last = {}
    for row in rows:
        if row.get("T") is not None:
            last[row["id"]] = (row["T"], setup_key(row.get("setup")))
    miss = []
    for row in rows:
        old, new = row.get("rerun_of"), row.get("id")
        if old and old in last and new in last and last[new][1] == CURRENT_SETUP and last[old][1] in off_ref \
                and last[old][1][3:] != CURRENT_SETUP[3:]:
            miss.append(last[new][0] - last[old][0] - off_ref[last[old][1]])
    corr = float(np.mean(miss)) if miss else 0.0
    offsets, out = {}, []
    if miss:
        offsets["supports_pair_correction_ms"] = round(corr * 1e3, 3)
        offsets["supports_pair_spread_ms"] = round(float(np.std(miss)) * 1e3, 3)
        offsets["supports_pairs"] = len(miss)
    for row in rows:
        if row.get("T") is None or row.get("genes") is None or (row.get("kind") == "ref" and not refs):
            continue
        k = setup_key(row.get("setup"))
        if k == CURRENT_SETUP:
            off = 0.0
        elif k in off_ref:
            off = off_ref[k] + (corr if k[3:] != CURRENT_SETUP[3:] else 0.0)
            offsets[str(k)] = round(off * 1e3, 3)
        else:
            continue
        out.append((row, row["T"] + off))
    return out, offsets, CURRENT_SETUP in ref


def _build_one(args):
    x, folder, lattice, support, fwing, freeform, rwing = (tuple(args) + (None,) * 5)[:7]
    try:
        return LC.build(np.asarray(x), Path(folder), lattice=lattice, support=support, fwing=fwing, freeform=freeform,
                        rwing=rwing)
    except Exception as e:                                           # noqa: BLE001
        return {"legal": False, "why": f"builder crashed: {e!r}"[:500], "genes": list(map(float, x))}


class Campaign:
    REF = "REF"                         # the reference car's folder under cars/
    def __init__(self, root):
        self.root = Path(root)
        self.cfg = json.loads((self.root / "config.json").read_text())
        self.st = json.loads((self.root / "state.json").read_text())
        for k, v in (("sobol_used", 0), ("batch", []), ("inflight", {}), ("pending_build", {}), ("best", None)):
            self.st.setdefault(k, v)
        self.cfg.setdefault("batch", 8)
        self.cfg.setdefault("surrogate", {"n_initial": 16, "draws_per_layout": 3, "seed": 1})

    def save(self):
        p = self.root / "state.json"; t = p.with_suffix(".tmp")
        t.write_text(json.dumps(self.st, indent=1)); t.replace(p)

    def log(self, row):
        with open(self.root / "designs.jsonl", "a") as f:
            f.write(json.dumps(row, default=float) + "\n")

    # ------------------------------------------------------------------ the surrogate
    def timed(self, refs=False):
        """[(row, T on the current setup)] of every car timed so far (timed_rows)."""
        out, self.offsets, self.ref_ready = timed_rows(self.designs(), refs)
        return out

    def data(self):
        """(X, T) of every car timed so far on (or brought over to) the current setup (timed(), the
        reference runs included), and X of the cars that could not be built or broke a rule."""
        t = self.timed(refs=True)
        X = [r["genes"] for r, _ in t]
        y = [T for _, T in t]
        bad = [r["genes"] for r in self.designs() if r.get("legal") is False and r.get("genes") is not None]
        return np.array(X, float).reshape(-1, LG.DIM), np.array(y, float), np.array(bad, float).reshape(-1, LG.DIM)

    def best_adjusted(self):
        """[id, T] of the fastest car on the current setup's scale (raw times of different setups do not compare)."""
        t = self.timed()
        if not t:
            return None
        r, T = min(t, key=lambda q: q[1])
        return [r["id"], T]

    def designs(self):
        p = self.root / "designs.jsonl"
        return [json.loads(l) for l in p.read_text().splitlines() if l.strip()] if p.exists() else []

    def candidates(self, rng, best_x=None):
        """Gene vectors covering EVERY legal layout on the 1 mm grid, the other numbers drawn at random
        (and, round the best car so far, perturbed), all snapped to what decode() will build."""
        iW, iR, iH = (LG.NAMES.index(k) for k in ("W", "rear", "halo_rear"))
        L = np.array(LG.legal_layouts(1.0))
        per = self.cfg["surrogate"].get("draws_per_layout", 3)
        X = rng.random((len(L) * per, LG.DIM))
        Lr = np.repeat(L, per, axis=0)
        X[:, iW] = (Lr[:, 0] - LG.LO[iW]) / (LG.HI[iW] - LG.LO[iW])
        X[:, iR] = (Lr[:, 1] - LG.LO[iR]) / np.maximum(np.minimum(LG.HI[iR], LG.BLOCK_MM - 16.0 - Lr[:, 0]) - LG.LO[iR], 1e-9)
        lo = np.maximum(Lr[:, 1] + 50.0, 102.0); hi = Lr[:, 0] + Lr[:, 1]
        X[:, iH] = np.clip((Lr[:, 2] - lo - 0.5) / np.maximum(hi - lo - 0.5, 1e-9), 0, 1)
        if best_x is not None:                                   # the neighbourhood of the best car
            near = np.clip(np.asarray(best_x)[None] + 0.08 * rng.standard_normal((4000, LG.DIM)), 0, 1)
            X = np.vstack([X, near])
        for name in SWITCHES:                                    # switches, not dials
            if name in LG.NAMES:
                k = LG.NAMES.index(name); X[:, k] = (X[:, k] >= 0.5).astype(float)
        return X

    def ask(self, n):
        """The next n gene vectors: Sobol while there is little data, then the model's best over all
        legal layouts (expected improvement, kriging believer for the batch)."""
        from scipy.stats import norm, qmc
        sg = self.cfg["surrogate"]
        X, y, bad = self.data()
        rng = np.random.default_rng(sg.get("seed", 1) + len(y) + 1000 * self.st["g"])
        if len(y) < sg["n_initial"]:
            S = qmc.Sobol(LG.DIM, seed=sg.get("seed", 1)).random(64)[self.st["sobol_used"]:self.st["sobol_used"] + n]
            self.st["sobol_used"] += n
            for name in SWITCHES:
                if name in LG.NAMES:
                    k = LG.NAMES.index(name); S[:, k] = (S[:, k] >= 0.5).astype(float)
            return [x for x in S], {"mode": "sobol"}
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
        mu0, sd0 = y.mean(), y.std() + 1e-12
        kern = ConstantKernel(1.0, (0.05, 20.0)) * Matern(length_scale=np.full(LG.DIM, 0.5), length_scale_bounds=(0.05, 20.0), nu=2.5) \
            + WhiteKernel(0.05, (1e-4, 1.0))
        gp = GaussianProcessRegressor(kern, normalize_y=False, n_restarts_optimizer=6, random_state=0).fit(X, (y - mu0) / sd0)
        ls = gp.kernel_.k1.k2.length_scale
        C = self.candidates(rng, X[np.argmin(y)])
        if len(bad):                                             # keep away from cars that could not be built
            from scipy.spatial import cKDTree
            C = C[cKDTree(bad).query(C)[0] > 0.12]
        picks, Xa, ya = [], X.copy(), (y - mu0) / sd0
        fixed = GaussianProcessRegressor(gp.kernel_, optimizer=None, normalize_y=False)
        best = ya.min()
        info = {"mode": "model", "setup_offsets_ms": getattr(self, "offsets", {}), "n_fit": int(len(y)),
                "length_scales": {k: round(float(v), 2) for k, v in zip(LG.NAMES, ls)},
                "noise_ms": round(float(np.sqrt(gp.kernel_.k2.noise_level)) * sd0 * 1e3, 3), "candidates": int(len(C)), "predicted": []}
        for _ in range(n):
            fixed.fit(Xa, ya)
            m, sdev = fixed.predict(C, return_std=True)
            z = (best - m) / np.maximum(sdev, 1e-9)
            ei = (best - m) * norm.cdf(z) + sdev * norm.pdf(z)
            k = int(np.argmax(ei))
            picks.append(C[k].copy()); info["predicted"].append([round(float(m[k] * sd0 + mu0), 5), round(float(sdev[k] * sd0 * 1e3), 3)])
            Xa = np.vstack([Xa, C[k]]); ya = np.append(ya, m[k])       # believe the model for the rest of the batch
            C = np.delete(C, k, axis=0)
        mb, _ = gp.predict(self.candidates(np.random.default_rng(0)), return_std=True)
        info["model_best_T"] = round(float(mb.min() * sd0 + mu0), 5)
        return picks, info

    def evaluator(self):
        import simscale_template as ST
        ST.BOX_M = BOX_M
        ST.SMALL_FEATURE_MM = float(self.cfg["evaluator"].get("small_feature_mm", SMALL_FEATURE_MM))
        e = self.cfg["evaluator"]
        return LayoutEvaluator(None, self.root / "simscale_state.json", runs_per_key=None, speed_mps=20.0,
                               fineness=6, box_mm=e.get("box_mm", 1.5), iters=e.get("iters", 1000),
                               sizing=e.get("sizing"), max_cost=e.get("max_cost"))

    def job(self, x):
        """What to build for one asked vector (the shape search builds a lattice on a fixed layout)."""
        return {"x": np.asarray(x).tolist()}

    def batch_done(self, ids):
        """Called once a whole batch is back (the shape search tells CMA here)."""

    def new_id(self):
        cid = f"L{self.st['n']:04d}"; self.st["n"] += 1
        return cid

    def submit_reference(self, ev, ladder=(0, 1)):
        """The team's own layout (every number at its default) through this route, once per triangle
        count in `ladder` (indices into CAR_FACES): the like-for-like reference for the search, and the
        difference between the two runs is the noise a layout's gain has to beat."""
        import shutil
        src = self.root / "cars" / "REF"
        rec = json.loads((src / "layout_car.json").read_text()) if (src / "layout_car.json").exists() else LC.build(LG.x_default(), src)
        self.log({"id": "REF", "kind": "ref", "genes": None, "legal": rec.get("legal"), "layout": rec.get("layout"),
                  "failed": rec.get("failed"), "why": rec.get("why"), "mass": rec.get("mass"), "cfd": rec.get("cfd")})
        if not rec.get("legal"):
            raise SystemExit(f"reference car is not legal: {rec.get('failed')} {rec.get('why')}")
        cases, refs = [], {}
        for i in ladder:
            cid = f"REF{CAR_FACES[i]}"
            shutil.copytree(src, self.root / "cars" / cid, dirs_exist_ok=True)
            cases.append({"id": cid, "folder": str(self.root / "cars" / cid), "attempt0": i})
            refs[cid] = None
        ev.submit(cases)
        self.st["refs"] = refs
        self.save()

    def submit_noise_reference(self, ev, pitch_mm=0.25, src="REF", name=None):
        """The reference car again with its CFD surface rebuilt on another grid: a different faceting of
        the same car, so the two results differ by the route's noise alone (the two reference copies at
        two triangle counts can both fall back onto the same rung of the ladder, as on 2026-11-22)."""
        import shutil
        cid = name or f"{src}P{int(round(pitch_mm * 100)):02d}"
        dst = self.root / "cars" / cid
        shutil.copytree(self.root / "cars" / src, dst, dirs_exist_ok=True)
        rec = json.loads((dst / "layout_car.json").read_text())
        rec.update(LC.export_cfd(dst.resolve(), rec["decoded"], pitch_mm=pitch_mm))
        (dst / "layout_car.json").write_text(json.dumps(rec, indent=1, default=float))
        if not (rec["cfd"]["one_piece"] and rec["cfd"]["watertight"]):
            raise SystemExit(f"noise reference not one closed solid: {rec['cfd']}")
        ev.submit([{"id": cid, "folder": str(dst.resolve())}])
        self.st.setdefault("refs", {})[cid] = None
        self.save()
        return cid

    def submit_fresh_reference(self, ev, name):
        """The team's layout built afresh with the current code (a new CFD setup: the supports) and timed
        as a reference: what brings the earlier setup's results over (data())."""
        dst = self.root / "cars" / name
        rec = LC.build(LG.x_default(), dst.resolve())
        self.log({"id": name, "kind": "ref", "genes": None, "legal": rec.get("legal"), "layout": rec.get("layout"),
                  "failed": rec.get("failed"), "why": rec.get("why"), "mass": rec.get("mass"), "cfd": rec.get("cfd"),
                  "support": rec.get("support")})
        if not rec.get("legal"):
            raise SystemExit(f"reference car is not legal: {rec.get('failed')} {rec.get('why')}")
        ev.submit([{"id": name, "folder": str(dst.resolve())}])
        self.st.setdefault("refs", {})[name] = None
        self.save()
        return name

    def rerun(self, ids):
        """Cars timed on an earlier setup built again with the current code and timed (new ids, budget
        counted): how far the setup's shift differs from car to car."""
        rows = {r["id"]: r for r in self.designs() if r.get("genes") is not None and r.get("kind") != "ref"}
        pb = self.st.setdefault("pending_build", {}) or {}
        for i in ids:
            pb[self.new_id()] = {"x": rows[i]["genes"], "rerun_of": i}
        self.st["pending_build"] = pb
        self.save()
        return list(pb)

    def collect(self, ev):
        """Results back from the CFD: race time with each car's own layout and mass, logged, files pruned."""
        st = self.st
        for cid, c in ev.st["cases"].items():                    # setups survive the case leaving the evaluator
            if c.get("setup"):
                st.setdefault("setups", {})[cid] = c["setup"]
        res = ev.collect()
        for cid, r in res.items():
            is_ref = cid in st.get("refs", {}) and st["refs"][cid] is None
            v = st["inflight"].pop(cid, None)
            if v is None and not is_ref:
                continue
            rec = json.loads((self.root / "cars" / cid / "layout_car.json").read_text())
            T = None
            if r.get("ok") and r.get("D20_N") is not None:
                T = race_time(rec, r["D20_N"], r["L20_N"] or 0.0)
            setup = (ev.st["cases"].get(cid) or {}).get("setup") or self.st.get("setups", {}).get(cid)
            self.log({"id": cid, "kind": "ref" if is_ref else "car", "genes": rec["genes"], "D20_N": r.get("D20_N"),
                      "L20_N": r.get("L20_N"), "T": T, "layout": rec.get("layout"), "eval": r.get("info"),
                      "setup": setup, "lattice": rec.get("lattice")})
            prune(self.root / "cars" / cid)
            if is_ref:
                st["refs"][cid] = {"D20_N": r.get("D20_N"), "L20_N": r.get("L20_N"), "T": T, "ok": r.get("ok")}
        st["best"] = self.best_adjusted() or st.get("best")      # on the current setup's scale
        self.save()
        return res

    def tick(self, ev):
        st = self.st
        if "refs" not in st:
            self.submit_reference(ev, ladder=(0,))
            self.submit_noise_reference(ev, src=self.REF)
        # 1) the next batch, when the last one is all back
        # (not before the current CFD setup's reference is timed: until then nothing earlier is comparable
        # and the model would have no data)
        self.timed()
        if not st["batch"] and not st.get("pending_build") and st["used"] < self.cfg["budget"] and ev.capacity() > 0 \
                and getattr(self, "ref_ready", True):
            xs, info = self.ask(min(self.cfg["batch"], self.cfg["budget"] - st["used"]))
            st["pending_build"] = {}
            for x in xs:
                st["pending_build"][self.new_id()] = self.job(x)
            st["batch"] = list(st["pending_build"])
            st["ask"] = info
            print(time.strftime("%H:%M:%S"), "batch", st["g"], json.dumps(info), flush=True)
            self.save()
        # 2) build what is pending (in parallel); the legal cars go to CFD
        pb = st.get("pending_build") or {}
        if pb:
            jobs = [(v["x"], str((self.root / "cars" / cid).resolve()), v.get("lattice"), v.get("support"), v.get("fwing"),
                     v.get("freeform"), v.get("rwing"))
                    for cid, v in pb.items()]
            with ProcessPoolExecutor(BUILD_WORKERS) as pool:
                recs = list(pool.map(_build_one, jobs))
            cases = []
            for (cid, v), rec in zip(pb.items(), recs):
                self.log({"id": cid, "batch": st["g"], "genes": rec.get("genes") or v["x"], "lattice": v.get("lattice"),
                          "v": v.get("v"), "support": rec.get("support"), "rerun_of": v.get("rerun_of"),
                          "legal": bool(rec.get("legal")),
                          "layout": rec.get("layout"), "failed": rec.get("failed"), "why": rec.get("why"),
                          "mass": rec.get("mass"), "build_s": rec.get("build_s"), "cfd": rec.get("cfd")})
                if not rec.get("legal"):
                    prune(self.root / "cars" / cid)
                if rec.get("legal"):
                    cases.append({"id": cid, "folder": str((self.root / "cars" / cid).resolve())})
                    st["inflight"][cid] = v
                    st["used"] += 1
            if cases:
                ev.submit(cases)
            st["pending_build"] = {}
            self.save()
        # 3) results
        res = self.collect(ev)
        # 4) the batch is done when none of it is still out
        if st["batch"] and not st.get("pending_build") and not any(i in st["inflight"] for i in st["batch"]):
            self.batch_done(st["batch"])
            st["g"] += 1; st["batch"] = []
        self.save()
        return res


def reference(root):
    """Only the reference cars (see Campaign.submit_reference), then stop."""
    c = Campaign(root)
    ev = c.evaluator()
    if "refs" not in c.st:
        c.submit_reference(ev)
    while any(v is None for v in c.st["refs"].values()):
        c.collect(ev)
        print(time.strftime("%H:%M:%S"), "reference", {k: (v and v.get("T")) for k, v in c.st["refs"].items()},
              {cid: (cc.get("phase"), cc.get("attempt"), cc.get("tries"), cc.get("last_error")) for cid, cc in ev.st["cases"].items()}, flush=True)
        if any(v is None for v in c.st["refs"].values()):
            time.sleep(60)
    print("reference done:", json.dumps(c.st["refs"]))


def init(root, budget, seed_keys_from=str(repo_paths.ROOT / "runs" / "sidepod" / "simscale_state.json")):
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    (root / "config.json").write_text(json.dumps({
        "_comment": "Layout search (simscale_layout.py): a Gaussian-process model over the numbers marked 'evolve' in "
                    "layout_params.json; every legal layout is scored on it each round, 8 go to CFD.",
        "budget": budget, "batch": 8, "surrogate": {"n_initial": 16, "draws_per_layout": 3, "seed": 1},
        "evaluator": {"iters": 1000, "box_mm": 1.5, "max_cost": 80.0, "sizing": {"max_mm": 30, "min_mm": 0.8},
                      "small_feature_mm": SMALL_FEATURE_MM},
        "rules": "no wheel supports (added later, their mass counted); cargo and tether guides on; milled top and bottom, "
                 "3.175 mm ball, 3 mm walls (Part 1); sidepod top <= 27.175 mm; road 0.15 mm under the wheels; lead to 48.2 g"}, indent=1))
    (root / "state.json").write_text(json.dumps({"n": 0, "g": 0, "used": 0, "batch": [], "inflight": {}, "sobol_used": 0,
                                                 "best": None, "pending_build": {}}, indent=1))
    src = Path(seed_keys_from)
    if src.exists() and not (root / "simscale_state.json").exists():
        k = json.loads(src.read_text())["keys"]
        keys = {fp: {"runs": v.get("runs", 0), "project": v.get("project"), "box": None, "retired": v.get("retired")}
                for fp, v in k.items()}
        (root / "simscale_state.json").write_text(json.dumps({"keys": keys, "cases": {}, "done": {}}, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("init", "ref", "noise", "refbuild", "rerun", "run", "tick", "status"))
    ap.add_argument("root")
    ap.add_argument("--budget", type=int, default=48)
    ap.add_argument("--pitch", type=float, default=0.25, help="noise: the grid the copy's surface is rebuilt on, mm")
    ap.add_argument("--name", default=None, help="noise: the copy's id; refbuild: the new reference's id")
    ap.add_argument("--ids", default="", help="rerun: comma-separated car ids to build and time again")
    a = ap.parse_args()
    if a.cmd == "init":
        return init(a.root, a.budget)
    if a.cmd == "ref":
        return reference(a.root)
    c = Campaign(a.root)
    if a.cmd == "noise":                  # stop `run` first: both write the state files
        ev = c.evaluator()
        print("submitted", c.submit_noise_reference(ev, pitch_mm=a.pitch, src=c.REF, name=a.name))
        return
    if a.cmd == "refbuild":               # stop `run` first: both write the state files
        ev = c.evaluator()
        print("submitted", c.submit_fresh_reference(ev, a.name or "REFS"))
        return
    if a.cmd == "rerun":                  # stop `run` first; the next tick builds them
        print("queued", c.rerun([i for i in a.ids.split(",") if i]))
        return
    if a.cmd == "status":
        ev = c.evaluator()
        print(json.dumps({"state": {k: v for k, v in c.st.items() if k != "pending_build"}, "eval": ev.status()}, indent=1, default=str))
        return
    ev = c.evaluator()
    while True:
        c.tick(ev)
        print(time.strftime("%H:%M:%S"), "batch", c.st["g"], "used", c.st["used"], "best", c.st["best"],
              "inflight", len(c.st["inflight"]), flush=True)
        if a.cmd == "tick" or (c.st["used"] >= c.cfg["budget"] and not c.st["inflight"] and not c.st["batch"]):
            return
        time.sleep(60)


if __name__ == "__main__":
    main()
