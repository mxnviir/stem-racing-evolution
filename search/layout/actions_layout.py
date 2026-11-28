"""The layout search (Bayesian, simscale_layout.Campaign's model) with its CFD on GitHub Actions and the team's
OpenFOAM pipeline (openfoam_eval.py) instead of SimScale. One batch per Actions run; the run's collect job ingests
the results, asks the next batch and starts it (.github/workflows/layout_openfoam.yml).

    python actions_layout.py init   <run dir> [--budget 100] [--batch 20]
    python actions_layout.py ingest <run dir> <results dir>      # every result.json under it
    python actions_layout.py next   <run dir>                    # writes <run dir>/request.json (or STOP)
    python actions_layout.py status <run dir>

Only results of this route are fitted (openfoam_eval.ROUTE): the SimScale times are another setup and do not mix.
The first batch: the reference (the team layout, every gene at its value) on meshing frames 0 and 1 -- its two
runs are the route's noise for this car -- and Sobol points. After it, every batch is the model's
(expected improvement, kriging believer), plus a re-run on a new meshing frame of the best car while it has
fewer than RERUN_TO runs, so the best is never one lucky mesh.
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import argparse
import json
from pathlib import Path

import numpy as np

import layout_genes as LG
import openfoam_eval as OE
import simscale_layout as SL

RERUN_TO = 2            # runs (meshing frames) the best car gets before it is trusted


class ActionsCampaign(SL.Campaign):
    def timed(self, refs=False):
        """[(row, T)] of every car timed on this route (the reference runs too if refs)."""
        out = []
        for r in self.designs():
            if r.get("route") != OE.ROUTE or r.get("T") is None:
                continue
            if r.get("kind") == "ref" and not refs:
                continue
            out.append((r, float(r["T"])))
        self.offsets, self.ref_ready = {}, any(r.get("kind") == "ref" for r, _ in out)
        return out

    def runs_of(self, cid):
        return [r for r in self.designs() if r.get("id") == cid and r.get("T") is not None]


def init(root, budget=100, batch=20, seed=1):
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    (root / "config.json").write_text(json.dumps({
        "_comment": "Layout search on GitHub Actions + the team's OpenFOAM pipeline (actions_layout.py, openfoam_eval.py).",
        "route": OE.ROUTE, "budget": budget, "batch": batch, "genes": LG.NAMES,
        "surrogate": {"n_initial": batch, "draws_per_layout": 3, "seed": seed}}, indent=1))
    (root / "state.json").write_text(json.dumps({"n": 0, "g": 0, "sobol_used": 0, "batch": [], "inflight": {}, "best": None}, indent=1))
    c = ActionsCampaign(root)
    ref = [float(v) for v in LG.x_default()]
    cases = [{"id": "REF", "k": 0, "x": ref, "kind": "ref"}, {"id": "REF", "k": 1, "x": ref, "kind": "ref"}]
    xs, info = c.ask(batch - len(cases))
    for x in xs:
        cases.append({"id": c.new_id(), "k": 0, "x": [float(v) for v in x], "kind": "car"})
    return _write(c, cases, info)


def _write(c, cases, info):
    c.st["inflight"] = {f"{q['id']}~{q['k']}": q for q in cases}
    req = {"gen": c.st["g"], "cases": [{"id": q["id"], "k": q["k"], "x": q["x"]} for q in cases], "info": info}
    (c.root / "request.json").write_text(json.dumps(req, indent=1, default=float))
    c.save()
    print(json.dumps({"gen": c.st["g"], "cases": len(cases), "mode": info.get("mode")}))
    return req


def ingest(root, results):
    c = ActionsCampaign(root)
    n = 0
    for p in sorted(Path(results).rglob("result.json")):
        r = json.loads(p.read_text())
        key = f"{r['id']}~{r['k']}"
        q = c.st["inflight"].pop(key, None)
        if q is None:
            continue
        c.log({"id": r["id"], "k": r["k"], "kind": q.get("kind", "car"), "route": OE.ROUTE, "gen": c.st["g"],
               "genes": r["genes"], "legal": r.get("legal"), "failed": r.get("failed"), "why": r.get("why"),
               "D20_N": r.get("D20_N"), "L20_N": r.get("L20_N"), "D_half_se_mN": r.get("D_half_se_mN"), "T": r.get("T"),
               "cells": r.get("cells"), "y_plus_max": r.get("y_plus_max"), "error": r.get("error"),
               "layout": r.get("layout"), "mass": r.get("mass"), "build_s": r.get("build_s"), "cfd_s": r.get("cfd_s")})
        n += 1
    for key, q in list(c.st["inflight"].items()):          # a case with no result (runner lost): logged as failed
        c.log({"id": q["id"], "k": q["k"], "kind": q.get("kind", "car"), "route": OE.ROUTE, "gen": c.st["g"], "genes": q["x"],
               "legal": None, "T": None, "error": "no result"})
        c.st["inflight"].pop(key)
    c.st["g"] += 1
    c.st["best"] = c.best_adjusted()
    c.save()
    print(json.dumps({"ingested": n, "best": c.st["best"]}))


def next_batch(root):
    c = ActionsCampaign(root)
    cfg = c.cfg
    cars = {r["id"] for r in c.designs() if r.get("kind") == "car"}
    if len(cars) >= cfg["budget"]:
        (c.root / "STOP").write_text(f"budget {cfg['budget']} reached\n")
        print(json.dumps({"stopped": True, "cars": len(cars)}))
        return None
    cases = []
    best = c.best_adjusted()
    if best:
        runs = c.runs_of(best[0])
        if len(runs) < RERUN_TO:
            cases.append({"id": best[0], "k": max(r.get("k", 0) for r in runs) + 1, "x": runs[0]["genes"], "kind": "car"})
    n = min(cfg["batch"] - len(cases), cfg["budget"] - len(cars))
    xs, info = c.ask(n)
    for x in xs:
        cases.append({"id": c.new_id(), "k": 0, "x": [float(v) for v in x], "kind": "car"})
    return _write(c, cases, info)


def status(root):
    c = ActionsCampaign(root)
    rows = [r for r in c.designs() if r.get("route") == OE.ROUTE]
    t = c.timed(refs=True)
    print(json.dumps({"gen": c.st["g"], "rows": len(rows), "timed": len(t), "illegal": sum(1 for r in rows if r.get("legal") is False),
                      "best": c.best_adjusted(), "inflight": len(c.st["inflight"])}, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd"); ap.add_argument("root"); ap.add_argument("args", nargs="*")
    ap.add_argument("--budget", type=int, default=100); ap.add_argument("--batch", type=int, default=20)
    a = ap.parse_args()
    if a.cmd == "init":
        init(a.root, a.budget, a.batch)
    elif a.cmd == "ingest":
        ingest(a.root, a.args[0])
    elif a.cmd == "next":
        next_batch(a.root)
    elif a.cmd == "status":
        status(a.root)
