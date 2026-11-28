"""Stage 2 (CMA-ES on the free-form body, simscale_shape.ShapeCampaign) with its CFD on GitHub Actions and the team's
OpenFOAM pipeline (search/layout/openfoam_eval.py), and -- when the campaign's `adjoint.mode` is "guided" -- an adjoint
at CMA's mean between generations (guided.py; search/adjoint: DAFoam on AWS Batch or a coarse one on a runner).

    python actions_shape.py init   <run dir> (--layout-run <stage-1 run dir> [--car ID] | --genes '[...]')
                                   [--budget 200] [--popsize 20] [--sigma 0.5] [--adjoint guided|off] [--every 1]
    python actions_shape.py ingest <run dir> <results dir>       # every result.json under it; tells CMA
    python actions_shape.py next   <run dir>                     # adjoint_request.json, or request.json, or STOP
    python actions_shape.py gradient <run dir> <adjoint result.json>   # the adjoint's result for this generation
    python actions_shape.py status <run dir>

One generation:
  next   -> (guided, and an adjoint due) writes <run>/adjoint_request.json {"gen", "job": fine_adjoint job for CMA's
            mean}; the adjoint workflow runs it and hands its result.json to `gradient`;
  next   -> writes <run>/request.json: popsize cars (guided: 4 of them on the adjoint's line) with their full specs
            in <run>/specs/<id>.json; the forward workflow runs them and hands the results to `ingest`.
A generation whose adjoint failed or was not usable is a plain CMA generation (the bench's rule).
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
import time
from pathlib import Path

import numpy as np

import guided
import openfoam_eval as OE
import simscale_shape as SS


class ActionsShape(SS.ShapeCampaign):
    def ask_guided(self, n, grad=None):
        es = self.es()
        info = {"mode": "cma", "generation": int(es.countiter)}
        if grad is not None:
            info.update(guided.inject_line(es, grad))
        xs = es.ask(n)
        info["sigma"] = round(float(es.sigma), 4)
        self._save_es(es)
        return xs, info

    def timed(self, refs=False):
        out = []
        for r in self.designs():
            if r.get("route") == OE.ROUTE and r.get("T") is not None and (refs or r.get("kind") != "ref"):
                out.append((r, float(r["T"])))
        self.offsets, self.ref_ready = {}, True
        return out

    def adjoint_job(self, v, gen):
        """The fine_adjoint job for design vector v (CMA's mean)."""
        j = self.job(v)
        a = self.cfg["adjoint"]
        return {"id": f"{self.root.name}-g{gen}-mean",
                "car": {"genes": j["x"], "freeform": j["freeform"], "support": j["support"], "fwing": j["fwing"],
                        "rwing": j["rwing"]},
                "parametrisation": {"freeform_seed": self.cfg["freeform"]["seed"], "v": [float(t) for t in v],
                                    "coarse": list(self.coarse), "amp_mm": list(self.amp)},
                "mesh": a.get("mesh", {"res": "resolved", "wall_resolved": True, "frame_k": 0}),
                "steady": a.get("steady", {"extra": 2000, "span_mN": 0.10}),
                "dafoam": a.get("dafoam", {}), "link": a.get("link", {"tol_mm": 1.0})}


def init(root, layout_run=None, car=None, genes=None, budget=200, popsize=20, sigma=0.5, adjoint="guided", every=1,
         seed_ff=None):
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    if genes is not None:                                   # a given layout (no stage-1 run folder)
        lr = root / "_layout"; (lr / "cars").mkdir(parents=True, exist_ok=True)
        (lr / "designs.jsonl").write_text(json.dumps({"id": "GIVEN", "genes": list(genes), "T": None}) + "\n")
        layout_run, car = lr, "GIVEN"
    SS.init(root, layout_run, car=car, budget=budget, popsize=popsize, sigma=sigma, seed_ff=seed_ff)
    cfg = json.loads((root / "config.json").read_text())
    cfg.update(route=OE.ROUTE, batch=popsize,
               adjoint={"mode": adjoint, "every": every, "backend": "aws",
                        "mesh": {"res": "resolved", "wall_resolved": True, "frame_k": 0},
                        "steady": {"extra": 2000, "span_mN": 0.10},
                        "dafoam": {"primal_iters": 300, "env": {"DAF_ADJ": "Krylov", "DAF_PC_FILL": "1", "DAF_PC_LEVEL": "2",
                                                                 "DAF_GMRES_RESTART": "400", "DAF_GMRES_ITERS": "3000"}},
                        "link": {"tol_mm": 1.0}})
    cfg["_comment"] = ("Stage 2 on GitHub Actions + OpenFOAM (actions_shape.py); adjoint.mode guided = an adjoint at "
                       "CMA's mean every `every` generations (guided.py, search/adjoint)")
    (root / "config.json").write_text(json.dumps(cfg, indent=1))
    (root / "state.json").write_text(json.dumps({"n": 0, "g": 0, "batch": [], "inflight": {}, "best": None,
                                                 "gradients": {}}, indent=1))
    print(json.dumps({"init": str(root), "dim": ActionsShape(root).dim(), "adjoint": adjoint}))


def _adjoint_due(c) -> bool:
    a = c.cfg.get("adjoint") or {}
    if a.get("mode") != "guided":
        return False
    g = c.st["g"]
    return g % max(int(a.get("every", 1)), 1) == 0 and str(g) not in c.st.get("gradients", {})


def next_step(root):
    c = ActionsShape(root)
    cars = {r["id"] for r in c.designs() if r.get("kind") == "car"}
    if len(cars) >= c.cfg["budget"]:
        (c.root / "STOP").write_text(f"budget {c.cfg['budget']} reached\n")
        print(json.dumps({"stopped": True, "cars": len(cars)}))
        return None
    (c.root / "adjoint_request.json").unlink(missing_ok=True)
    if _adjoint_due(c):
        es = c.es()
        v = guided.mean_phenotype(es)
        req = {"gen": c.st["g"], "job": c.adjoint_job(v, c.st["g"]), "backend": c.cfg["adjoint"].get("backend", "aws"),
               "written": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        (c.root / "adjoint_request.json").write_text(json.dumps(req, indent=1))
        c.st["phase"] = "adjoint"; c.save()
        print(json.dumps({"phase": "adjoint", "gen": c.st["g"]}))
        return req
    gr = (c.st.get("gradients") or {}).get(str(c.st["g"]))
    grad = None
    if gr and gr.get("usable"):
        n = c.fr.n_genes(*c.coarse)
        grad = np.zeros(c.dim())
        grad[:n] = gr["dD_dv_car_mN"]
    n = min(c.cfg["batch"], c.cfg["budget"] - len(cars))
    xs, info = c.ask_guided(n, grad)
    cases = []
    (c.root / "specs").mkdir(exist_ok=True)
    for x in xs:
        cid = c.new_id()
        spec = c.job(np.asarray(x, float))
        (c.root / "specs" / f"{cid}.json").write_text(json.dumps(spec))
        cases.append({"id": cid, "k": 0, "spec": f"specs/{cid}.json", "kind": "car"})
    c.st["inflight"] = {f"{q['id']}~{q['k']}": q for q in cases}
    c.st["batch"] = [q["id"] for q in cases]
    c.st["phase"] = "forward"
    req = {"gen": c.st["g"], "cases": [{"id": q["id"], "k": q["k"], "spec": q["spec"]} for q in cases], "info": info}
    (c.root / "request.json").write_text(json.dumps(req, indent=1, default=float))
    c.save()
    print(json.dumps({"phase": "forward", "gen": c.st["g"], "cases": len(cases), "guided": info.get("guided", False)}))
    return req


def gradient(root, result_json):
    """The adjoint's result for this generation: stored; a failed or unusable one is stored as such (the generation
    then runs plain)."""
    c = ActionsShape(root)
    p = Path(result_json)
    r = json.loads(p.read_text()) if p.exists() else {"usable": False, "error": "no result"}
    g = r.get("gradient") or {}
    rec = {"usable": bool(r.get("usable")) and "coarse_car_mN_per_unit" in g, "error": r.get("error"),
           "production": r.get("production"), "steady": r.get("steady"), "dafoam": r.get("dafoam"),
           "n_active": g.get("n_active"), "active_share_of_dD": g.get("active_share_of_dD")}
    if rec["usable"]:
        rec["dD_dv_car_mN"] = g["coarse_car_mN_per_unit"]
    c.st.setdefault("gradients", {})[str(c.st["g"])] = rec
    c.save()
    print(json.dumps({"gen": c.st["g"], "usable": rec["usable"], "error": rec["error"]}))


def ingest(root, results):
    c = ActionsShape(root)
    n = 0
    for p in sorted(Path(results).rglob("result.json")):
        r = json.loads(p.read_text())
        key = f"{r['id']}~{r['k']}"
        q = c.st["inflight"].pop(key, None)
        if q is None:
            continue
        c.log({"id": r["id"], "k": r["k"], "kind": "car", "route": OE.ROUTE, "gen": c.st["g"], "v": r.get("v"),
               "genes": r.get("genes"), "legal": r.get("legal"), "failed": r.get("failed"), "why": r.get("why"),
               "D20_N": r.get("D20_N"), "L20_N": r.get("L20_N"), "D_half_se_mN": r.get("D_half_se_mN"), "T": r.get("T"),
               "cells": r.get("cells"), "error": r.get("error"), "mass": r.get("mass"), "cfd_s": r.get("cfd_s")})
        n += 1
    for key, q in list(c.st["inflight"].items()):            # no result (runner lost): told to CMA as the worst
        spec = json.loads((c.root / q["spec"]).read_text())
        c.log({"id": q["id"], "k": q["k"], "kind": "car", "route": OE.ROUTE, "gen": c.st["g"], "v": spec["v"],
               "legal": None, "T": None, "error": "no result"})
        c.st["inflight"].pop(key)
    c.batch_done(c.st["batch"])
    c.st["g"] += 1
    t = c.timed()
    c.st["best"] = min(([r["id"], T] for r, T in t), key=lambda q: q[1]) if t else None
    c.save()
    print(json.dumps({"ingested": n, "best": c.st["best"]}))


def status(root):
    c = ActionsShape(root)
    rows = [r for r in c.designs() if r.get("route") == OE.ROUTE]
    gr = c.st.get("gradients", {})
    print(json.dumps({"gen": c.st["g"], "phase": c.st.get("phase"), "rows": len(rows), "timed": len(c.timed()),
                      "best": c.st.get("best"), "gradients": {k: v.get("usable") for k, v in gr.items()},
                      "sigma": round(float(c.es().sigma), 4)}, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("root"); ap.add_argument("args", nargs="*")
    ap.add_argument("--layout-run"); ap.add_argument("--car"); ap.add_argument("--genes")
    ap.add_argument("--budget", type=int, default=200); ap.add_argument("--popsize", type=int, default=20)
    ap.add_argument("--sigma", type=float, default=0.5); ap.add_argument("--adjoint", default="guided")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--seed-car", help="a built car folder whose free-form body is the seed (no rebuild, no fit)")
    a = ap.parse_args()
    if a.cmd == "init":
        sf = json.loads((Path(a.seed_car) / "layout_car.json").read_text())["decoded"]["freeform"] if a.seed_car else None
        init(a.root, a.layout_run, a.car, json.loads(a.genes) if a.genes else None, a.budget, a.popsize, a.sigma,
             a.adjoint, a.every, seed_ff=sf)
    elif a.cmd == "ingest":
        ingest(a.root, a.args[0])
    elif a.cmd == "next":
        next_step(a.root)
    elif a.cmd == "gradient":
        gradient(a.root, a.args[0])
    elif a.cmd == "status":
        status(a.root)
