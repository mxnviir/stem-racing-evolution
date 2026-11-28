"""Stage 2 of the whole-car search: free-form evolution of the body on one fixed layout.

Stage 1 (simscale_layout.py) finds the layout -- wheelbase, tail, halo, cargo, nose, wing height and
position -- with a surrogate model. This stage takes the best car of stage 1 (or any car's genes),
holds its layout, and lets CMA-ES design the body and nose freely (car/body/freeform.py: a top surface, a
bottom surface and a side outline from the nose tip to the rear face, user 2026-11-23: no fixed stations,
no constant-width sidepods), together with the front wing (one piece or two, its own supports) and the
wheel supports' plates. The surfaces start as the stage-1 car itself (fitted to its body and nose); CMA
moves them by a coarse grid of displacements (the number of free numbers the CFD runs can search, a
setting), before Part 1 applies the rule masks (cargo, halo pocket, chamber wall, wheel keep-clear zones,
the 27.175 mm sidepod top) and the 3.175 mm ball-end milling. So whatever the numbers, the car stays
legal and millable; a vector that breaks a rule anyway is told to CMA as the worst.

The search, the CFD route (one closed solid, shrink-wrapped surface, exact wheels, SimScale) and the race
model (each car's own mass and COM) are stage 1's (simscale_layout.Campaign); only what is asked differs.

    python simscale_shape.py init runs/shape --from runs/layout [--car L0012] --budget 40
    python simscale_shape.py run  runs/shape
    python simscale_shape.py status runs/shape
"""
from __future__ import annotations

import sys
import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
import pickle
import time
from pathlib import Path

import numpy as np

import simscale_layout as SL

COARSE = (5, 2, 3, 4)        # displacement control points: along the car, across the top, up the side, across the floor
AMP_MM = (6.0, 4.0, 6.0)     # mm per unit of CMA's numbers: top, bottom, side
FIT_GRID = (24, 9, 12)       # the surfaces' own control grid, fitted to the stage-1 car (IoU 0.85 on L0044)
NOSE_LEAD_MM = 2.0          # the seed nose starts this far ahead of the seed front wing's leading edge


class ShapeCampaign(SL.Campaign):
    REF = "BASE"
    def __init__(self, root):
        super().__init__(root)
        self.base = np.asarray(self.cfg["base_genes"], float)
        f = self.cfg["freeform"]
        import freeform as fr
        self.fr = fr
        self.seed = fr.Freeform.from_dict(f["seed"])
        self.coarse = tuple(f["coarse"])
        self.amp = tuple(f["amp_mm"])

    # ------------------------------------------------------------------ CMA
    def es(self):
        p = self.root / "cma.pkl"
        if p.exists():
            return pickle.loads(p.read_bytes())
        import cma
        c = self.cfg["cma"]
        es = cma.CMAEvolutionStrategy(np.zeros(self.dim()), c["sigma"],
                                      {"popsize": c["popsize"], "seed": c.get("seed", 1), "verbose": -9,
                                       "bounds": [-2.0, 2.0]})
        self._save_es(es)
        return es

    def _save_es(self, es):
        p = self.root / "cma.pkl"; t = p.with_suffix(".tmp")
        t.write_bytes(pickle.dumps(es)); t.replace(p)

    def ask(self, n):
        es = self.es()
        xs = es.ask(n)
        self.st["asked"] = [np.asarray(x).tolist() for x in xs]
        self._save_es(es)
        return xs, {"mode": "cma", "sigma": round(float(es.sigma), 4), "generation": int(es.countiter)}

    def dim(self):
        """CMA's vector: the surfaces' coarse displacements, then the shape genes (layout_genes.SHAPE_GENES:
        the wheel supports and the front wing; 0 = the catalogue's value, -1 / +1 its lo / hi)."""
        return self.fr.n_genes(*self.coarse) + len(SL.LG.SHAPE_GENES)

    def job(self, v):
        v = np.asarray(v, float)
        n = self.fr.n_genes(*self.coarse)
        ff = self.fr.from_unit(self.seed, v[:n], self.amp, *self.coarse)
        sh = SL.LG.shape_from_unit(v[n:], stage1=self.base)      # the layout search's wing numbers kept
        return {"x": self.base.tolist(), "v": v.tolist(), "freeform": ff.as_dict(),
                "support": sh["support"], "fwing": sh["fwing"], "rwing": sh["rwing"]}

    def batch_done(self, ids):
        """Every car of the generation is back: tell CMA its race times (a car that could not be built,
        broke a rule or failed in the CFD as the worst of the generation plus 5 ms)."""
        es = self.es()
        T, V = {}, {}
        for row in self.designs():
            if row.get("id") in ids:
                if "T" in row:
                    T[row["id"]] = row["T"]
                if row.get("v") is not None:                  # the whole CMA vector (the build's row)
                    V[row["id"]] = row["v"]
                elif row.get("lattice") and row["id"] not in V:
                    V[row["id"]] = row["lattice"]["values"]
        vs, ys = [], []
        for cid in ids:
            if cid not in V:                         # nothing to tell (should not happen: every build row carries it)
                continue
            vs.append(np.asarray(V[cid], float))
            ys.append(T.get(cid))
        good = [y for y in ys if y is not None]
        worst = (max(good) if good else 1.0) + 0.005
        es.tell(vs, [y if y is not None else worst for y in ys])
        self._save_es(es)
        print(time.strftime("%H:%M:%S"), "told CMA generation", int(es.countiter), "best this generation",
              min(good) if good else None, flush=True)

    def submit_reference(self, ev, ladder=(0,)):
        """The base car itself (every number of CMA's at 0: the fitted surfaces, the catalogue's wing and
        supports) through the CFD once: what the shape changes are measured against."""
        import shutil
        src = self.root / "cars" / "BASE"
        if not (src / "layout_car.json").exists():
            j = self.job(np.zeros(self.dim()))
            rec = SL.LC.build(self.base, src.resolve(), freeform=j["freeform"], support=j["support"], fwing=j["fwing"],
                              rwing=j.get("rwing"))
        else:
            rec = json.loads((src / "layout_car.json").read_text())
        self.log({"id": "BASE", "kind": "ref", "genes": None, "legal": rec.get("legal"), "layout": rec.get("layout"),
                  "failed": rec.get("failed"), "why": rec.get("why"), "mass": rec.get("mass")})
        if not rec.get("legal"):
            raise SystemExit(f"base car is not legal: {rec.get('failed')} {rec.get('why')}")
        cid = f"BASE{SL.CAR_FACES[0]}"
        shutil.copytree(src, self.root / "cars" / cid, dirs_exist_ok=True)
        ev.submit([{"id": cid, "folder": str((self.root / "cars" / cid).resolve()), "attempt0": 0}])
        self.st["refs"] = {cid: None}
        self.save()

    def timed(self):                     # CMA does not need the base car's time to ask a generation
        out = super().timed()
        self.ref_ready = True
        return out

    def data(self):                      # the GP of stage 1 is not used here
        raise NotImplementedError


def best_of(layout_root: Path, car: str | None = None):
    """(id, genes, T) of a stage-1 car: the named one, or the fastest that came back."""
    rows = {}
    for line in (layout_root / "designs.jsonl").read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            rows.setdefault(r["id"], {}).update(r)
    if car:
        r = rows[car]
        return car, r["genes"], r.get("T")
    # the fastest on the current CFD setup's scale (raw times of different setups do not compare)
    lines = [json.loads(l) for l in (layout_root / "designs.jsonl").read_text().splitlines() if l.strip()]
    t, _, _ = SL.timed_rows(lines)
    if not t:
        raise SystemExit(f"no timed car in {layout_root}")
    r, T = min(t, key=lambda q: q[1])
    return r["id"], rows[r["id"]]["genes"], T


def fit_seed(car_dir: Path, genes):
    """The free-form surfaces of a built car's BODY, on FIT_GRID, with a nose made from the body itself
    (freeform.seed_nose: its Ref A section carried forward and rounded off just ahead of the front wing
    stage 2 starts from), not the stage-1 car's nose cone (user 2026-11-23: no cone)."""
    import trimesh
    import freeform as fr
    d = SL.LG.decode(genes)
    body = trimesh.load(str(car_dir / "body_half.stl"), force="mesh"); body.apply_scale(1e3)
    x_tip = d["x_front"] - 16.0 - d["nose"]["length_mm"]
    x_rear = d["x_front"] + d["W"] + d["rear"]
    ref_a = d["x_front"] - 16.0
    fw = SL.LG.shape_from_unit(np.zeros(len(SL.LG.SHAPE_GENES)), stage1=genes)["fwing"]
    x_le = ref_a - fw["te_gap_mm"] - fw["chord_mm"]
    ff = fr.fit([body], x_tip, x_rear, *FIT_GRID)
    return fr.seed_nose(ff, body, ref_a, max(x_le - NOSE_LEAD_MM, x_tip)).as_dict()


def init(root, layout_root, car=None, budget=40, popsize=8, sigma=0.5, seed_ff=None):
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    cid, genes, T = best_of(Path(layout_root), car)
    if seed_ff is None:
        src = Path(layout_root) / "cars" / cid
        if not (src / "body_half.stl").exists():
            SL.LC.build(np.asarray(genes), (root / "cars" / "SEED").resolve())
            src = root / "cars" / "SEED"
        seed_ff = fit_seed(src, genes)
    (root / "config.json").write_text(json.dumps({
        "_comment": "Stage 2 (simscale_shape.py): CMA-ES on the free-form body and nose (fitted to the stage-1 car), "
                    "the front wing and the wheel supports, layout held at base_genes.",
        "base_from": {"run": str(layout_root), "car": cid, "T": T}, "base_genes": genes,
        "freeform": {"seed": seed_ff, "coarse": list(COARSE), "amp_mm": list(AMP_MM)},
        "budget": budget, "batch": popsize, "cma": {"popsize": popsize, "sigma": sigma, "seed": 1},
        "evaluator": {"iters": 1000, "box_mm": 1.5, "max_cost": 80.0, "sizing": {"max_mm": 30, "min_mm": 0.8},
                      "small_feature_mm": SL.SMALL_FEATURE_MM}}, indent=1))
    (root / "state.json").write_text(json.dumps({"n": 0, "g": 0, "used": 0, "batch": [], "inflight": {}, "best": None,
                                                 "pending_build": {}}, indent=1))
    src = Path(layout_root) / "simscale_state.json"
    if src.exists() and not (root / "simscale_state.json").exists():
        k = json.loads(src.read_text())["keys"]
        keys = {fp: {"runs": v.get("runs", 0), "project": v.get("project"), "box": None, "retired": v.get("retired")}
                for fp, v in k.items()}
        (root / "simscale_state.json").write_text(json.dumps({"keys": keys, "cases": {}, "done": {}}, indent=1))
    print("stage 2 on", cid, "T", T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("init", "run", "tick", "status"))
    ap.add_argument("root")
    ap.add_argument("--from", dest="src", default="runs/layout")
    ap.add_argument("--car", default=None)
    ap.add_argument("--budget", type=int, default=40)
    a = ap.parse_args()
    if a.cmd == "init":
        return init(a.root, a.src, a.car, a.budget)
    c = ShapeCampaign(a.root)
    ev = c.evaluator()
    if a.cmd == "status":
        print(json.dumps({"state": {k: v for k, v in c.st.items() if k not in ("pending_build", "asked")}, "eval": ev.status()},
                         indent=1, default=str))
        return
    while True:
        c.tick(ev)
        print(time.strftime("%H:%M:%S"), "generation", c.st["g"], "used", c.st["used"], "best", c.st["best"],
              "inflight", len(c.st["inflight"]), flush=True)
        if a.cmd == "tick" or (c.st["used"] >= c.cfg["budget"] and not c.st["inflight"] and not c.st["batch"]):
            return
        time.sleep(60)


if __name__ == "__main__":
    main()
