"""P3: the stage-1 layout search (Bayesian: Gaussian process + expected improvement over the 12 layout numbers) on
OpenFOAM (carhybrid/DESIGN.md P3; carhybrid/BO.md).

The search itself is the evolution builder's own (search/layout/simscale_layout.Campaign.ask: Sobol while there is
little data, then a GP with one length scale per number fitted to the race times, every legal layout on the 1 mm grid
x a few draws of the other numbers scored on the model by expected improvement, the batch picked one by one with the
kriging believer, kept away from layouts that could not be built). Only the data come from a different place: each
car is built by carhybrid/make_car.py (the patched builder, the stage-1 car's organic free-form body, u = 0) and
timed by the audit protocol (fresh wall-resolved production mesh, 1000 iterations, frame k) on the GitHub workflow;
its race time is the evolution race model at the car's own mass and layout, with that frame's drag and lift
(carhybrid/loop.py frame_score), so search and loop use one race time.

    python carhybrid/bo.py init [--budget 48] [--batch 8] [--initial 16] [--seed 1]
    python carhybrid/bo.py advance <results folder> <run id> [--evo <builder checkout>]   (0 next, 3 done, 2 refused)
    python carhybrid/bo.py status

State: carhybrid/bo/state.json (config, the batches, what is pending), carhybrid/bo/designs.jsonl (one row per car:
id, genes, batch, frame, legal, why, drag/lift, race time). Batches: bo_<n> ("search", one frame per car: frame =
car number mod 4, so the GP's noise term sees the frame-to-frame scatter) and, once the budget is spent, bo_confirm
(the three fastest layouts on frames 0-3, P3's mark). Each advance writes the next batch into carhybrid/request.json
and sets "launch"; the workflow's collect hook dispatches it (as for loop batches).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE / "bo"
REQ = HERE / "request.json"
EXIT_OK, EXIT_REFUSED, EXIT_DONE = 0, 2, 3
CONFIRM_TOP, CONFIRM_FRAMES = 3, (0, 1, 2, 3)
# re-run alongside the best layouts in the confirm stage, paired by frame (P3's mark: does the search's best beat the
# layout P2 starts from by >= 3 SE): L0039 (the fastest SimScale-timed layout, P2's start) and the team car
CONFIRM_REFS = (("REF_L0039", "L0039"), ("REF_team", "team"))


def load(p, default=None):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else default


def save(p, d):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    t = p.with_suffix(".tmp"); t.write_text(json.dumps(d, indent=1, default=float) + "\n"); t.replace(p)


def rows(root=ROOT):
    p = Path(root) / "designs.jsonl"
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()] if p.exists() else []


def append_rows(new, root=ROOT):
    p = Path(root) / "designs.jsonl"; p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a") as f:
        for r in new:
            f.write(json.dumps(r, default=float) + "\n")


# ------------------------------------------------------------------ the evolution search, fed with OpenFOAM data
def evo_modules(evo):
    """(layout_genes, simscale_layout) from the patched builder checkout."""
    sys.path.insert(0, str(HERE))
    from grad import evo_import
    evo_import(evo)
    import layout_genes as LG
    import simscale_layout as SL
    return LG, SL


def campaign(evo, st, data_rows):
    """simscale_layout.Campaign with its data from data_rows (this module's designs) instead of SimScale's."""
    LG, SL = evo_modules(evo)

    class OFCampaign(SL.Campaign):
        def __init__(self):                               # no folder: config and state are this module's
            self.root = ROOT
            self.cfg = {"batch": st["batch"], "surrogate": {"n_initial": st["initial"], "draws_per_layout": 3,
                                                            "seed": st["seed"]}}
            self.st = {"sobol_used": st["sobol_used"], "g": st["n_batches"], "n": st["n_cars"]}
            self.offsets = {}

        def data(self):
            t = [(r["genes"], r["T"]) for r in data_rows if r.get("T") is not None]
            bad = [r["genes"] for r in data_rows if r.get("legal") is False and r.get("genes") is not None]
            X = np.array([g for g, _ in t], float).reshape(-1, LG.DIM)
            y = np.array([T for _, T in t], float)
            return X, y, np.array(bad, float).reshape(-1, LG.DIM)

    return OFCampaign(), LG, SL


# ------------------------------------------------------------------ batches
def search_cases(st, genes_list):
    cases, ids = [], []
    for x in genes_list:
        cid = f"B{st['n_cars']:03d}"
        k = st["n_cars"] % 4
        st["n_cars"] += 1
        ids.append({"id": cid, "genes": [round(float(v), 6) for v in x], "k": k})
        cases.append({"name": f"bo_{cid}_k{k}", "type": "score", "layout": json.dumps([round(float(v), 6) for v in x]),
                      "k": k})
    return cases, ids


def write_batch(st, name, cases, meta):
    R = load(REQ, {}) or {}
    R.setdefault("batches", {})[name] = cases
    R["launch"] = name
    R.setdefault("_bo_note", "P3 (carhybrid/bo.py, carhybrid/BO.md): batches bo_<n> and bo_confirm are written by "
                             "bo.py advance; its state is carhybrid/bo/state.json and designs.jsonl")
    save(REQ, R)
    st["pending"] = {"batch": name, **meta}
    st.setdefault("batches", []).append({"batch": name, **{k: v for k, v in meta.items() if k != "cars"},
                                         "n": len(cases)})


def ask_next(st, evo, n):
    C, LG, SL = campaign(evo, st, rows())
    picks, info = C.ask(n)
    st["sobol_used"] = C.st["sobol_used"]
    return picks, info


# ------------------------------------------------------------------ commands
def init(budget=48, batch=8, initial=16, seed=1, evo="evo", force=False):
    if load(ROOT / "state.json") and not force:
        print("bo: already initialised (carhybrid/bo/state.json); --force starts again", file=sys.stderr)
        return EXIT_REFUSED
    if force and (ROOT / "designs.jsonl").exists():
        (ROOT / "designs.jsonl").unlink()
    st = {"budget": budget, "batch": batch, "initial": initial, "seed": seed, "sobol_used": 0, "n_cars": 0,
          "n_batches": 0, "status": "search", "runs": []}
    picks, info = ask_next(st, evo, initial)              # the first batch: the whole Sobol start at once
    st["n_batches"] = 1
    cases, ids = search_cases(st, picks)
    write_batch(st, "bo_1", cases, {"kind": "search", "mode": info.get("mode"), "cars": ids})
    save(ROOT / "state.json", st)
    print(f"bo: init, batch bo_1 ({len(cases)} cars, {info.get('mode')})")
    return EXIT_OK


def frame_T(results, name, evo):
    """(T s, block SE s, info) for one scored car: loop.frame_score (the race model at the car's own mass)."""
    import loop
    return loop.frame_score(results, name, evo)


def ingest_search(st, results, run_id, evo):
    P = st["pending"]
    new = []
    for c in P["cars"]:
        name = f"bo_{c['id']}_k{c['k']}"
        rec = load(Path(results) / name / "rec.json")
        res = load(Path(results) / name / "result.json")
        row = {"id": c["id"], "genes": c["genes"], "batch": P["batch"], "k": c["k"], "run": str(run_id)}
        if rec is None:
            row.update(legal=None, why="no build record (job failed before the build finished)")
        elif not rec.get("legal") or rec.get("failed"):
            row.update(legal=False, why="; ".join(map(str, rec.get("failed") or [rec.get("why")]))[:300],
                       decoded={k: rec.get("decoded", {}).get(k) for k in ("W", "rear", "halo_rear", "x_front")})
        else:
            T, se, info = frame_T(results, name, evo)
            row.update(legal=True, decoded={k: rec.get("decoded", {}).get(k) for k in ("W", "rear", "halo_rear", "x_front")},
                       mass_g=(rec.get("mass") or {}).get("race_mass_g"))
            if T is None:
                row.update(T=None, why=info.get("why"))
            else:
                row.update(T=float(T), T_se=float(se), drag_mN=info.get("drag_mN"), lift_mN=info.get("lift_mN"),
                           drag_se_mN=info.get("se_mN"))
        new.append(row)
    append_rows(new)
    return new


def ingest_confirm(st, results, run_id, evo):
    out = {}
    for lay in st["pending"]["layouts"]:
        Ts = {}
        for k in CONFIRM_FRAMES:
            T, se, info = frame_T(results, f"bo_confirm_{lay['id']}_k{k}", evo)
            if T is not None:
                Ts[k] = {"T": float(T), "se": float(se), "drag_mN": info.get("drag_mN"), "lift_mN": info.get("lift_mN")}
        out[lay["id"]] = {"genes": lay["genes"], "frames": Ts,
                          "T_mean": float(np.mean([v["T"] for v in Ts.values()])) if Ts else None}
    # P3's mark: the search's best (lowest mean over the frames it has) against each reference, paired by frame
    best = min((k for k in out if not k.startswith("REF_") and out[k]["T_mean"] is not None),
               key=lambda k: out[k]["T_mean"], default=None)
    marks = {}
    for ref in (k for k in out if k.startswith("REF_")):
        if best is None:
            break
        ks = sorted(set(out[best]["frames"]) & set(out[ref]["frames"]))
        if len(ks) >= 2:
            d = np.array([out[best]["frames"][k]["T"] - out[ref]["frames"][k]["T"] for k in ks])
            se = float(d.std(ddof=1) / np.sqrt(len(d)))
            marks[f"{best}-{ref}"] = {"mean_ms": float(d.mean() * 1e3), "se_ms": se * 1e3, "frames": ks,
                                      "beats_by_3SE": bool(d.mean() <= -3 * se)}
    st["confirm"] = {"run": str(run_id), "layouts": out, "best": best, "marks": marks}
    return out


def confirm_cases(st):
    best = {}
    for r in rows():
        if r.get("T") is not None and (r["id"] not in best or r["T"] < best[r["id"]]["T"]):
            best[r["id"]] = r
    top = sorted(best.values(), key=lambda r: r["T"])[:CONFIRM_TOP]
    lays = [{"id": r["id"], "layout": json.dumps(r["genes"]), "genes": r["genes"], "T_search": r["T"]} for r in top]
    lays += [{"id": rid, "layout": lay, "genes": None, "T_search": None} for rid, lay in CONFIRM_REFS]
    cases = [{"name": f"bo_confirm_{L['id']}_k{k}", "type": "score", "layout": L["layout"], "k": k}
             for L in lays for k in CONFIRM_FRAMES]
    return cases, lays


def advance(results, run_id, evo="evo"):
    st = load(ROOT / "state.json")
    if not st or not st.get("pending"):
        print("bo: nothing pending", file=sys.stderr)
        return EXIT_REFUSED
    P = st["pending"]
    if Path(results).name != P["batch"]:
        print(f"bo: results {Path(results).name} are not the pending batch {P['batch']}", file=sys.stderr)
        return EXIT_REFUSED
    st["runs"].append({"batch": P["batch"], "run": str(run_id)})
    if P["kind"] == "confirm":
        out = ingest_confirm(st, results, run_id, evo)
        st["pending"], st["status"] = None, "done"
        save(ROOT / "state.json", st)
        print("bo: confirm ingested: " + ", ".join(f"{k} {v['T_mean']}" for k, v in out.items()))
        return EXIT_DONE
    new = ingest_search(st, results, run_id, evo)
    timed = [r for r in rows() if r.get("T") is not None]
    print(f"bo: {P['batch']}: {sum(r.get('T') is not None for r in new)} timed, "
          f"{sum(r.get('legal') is False for r in new)} illegal, {len(new)} cars; {len(timed)} timed in all")
    if st.get("paused"):                # owner, 2026-11-28: the search waits for the newer car; keep the data, ask nothing
        st["pending"], st["status"] = None, "paused"
        save(ROOT / "state.json", st)
        print(f"bo: paused ({st['paused']}); {P['batch']} ingested, no next batch")
        return EXIT_DONE
    spent = st["n_cars"]
    if spent >= st["budget"]:
        cases, lays = confirm_cases(st)
        write_batch(st, "bo_confirm", cases, {"kind": "confirm", "layouts": lays})
        st["status"] = "confirm"
        save(ROOT / "state.json", st)
        print(f"bo: budget spent ({spent} cars): bo_confirm, {len(lays)} layouts x {len(CONFIRM_FRAMES)} frames")
        return EXIT_OK
    n = min(st["batch"], st["budget"] - spent)
    picks, info = ask_next(st, evo, n)
    st["n_batches"] += 1
    name = f"bo_{st['n_batches']}"
    cases, ids = search_cases(st, picks)
    meta = {"kind": "search", "mode": info.get("mode"), "cars": ids}
    for k in ("length_scales", "noise_ms", "model_best_T", "n_fit", "predicted", "candidates"):
        if k in info:
            meta[k] = info[k]
    write_batch(st, name, cases, meta)
    save(ROOT / "state.json", st)
    print(f"bo: next {name} ({len(cases)} cars, {info.get('mode')}"
          + (f", model best {info['model_best_T']:.4f} s, noise {info.get('noise_ms')} ms" if "model_best_T" in info else "") + ")")
    return EXIT_OK


def status():
    st = load(ROOT / "state.json")
    if not st:
        print("bo: not initialised"); return EXIT_OK
    R = rows()
    timed = sorted((r for r in R if r.get("T") is not None), key=lambda r: r["T"])
    print(f"bo: {st['status']}; {st['n_cars']}/{st['budget']} cars asked, {len(R)} back, {len(timed)} timed, "
          f"{sum(r.get('legal') is False for r in R)} illegal; pending {(st.get('pending') or {}).get('batch')}")
    for r in timed[:5]:
        print(f"  {r['id']} T {r['T']:.5f} s (k{r['k']}), drag {r.get('drag_mN')} mN, layout {r.get('decoded')}")
    return EXIT_OK


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__); return EXIT_OK
    opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d
    evo = opt("--evo", os.environ.get("CARHYBRID_EVO", "evo"))
    if argv[0] == "init":
        return init(int(opt("--budget", 48)), int(opt("--batch", 8)), int(opt("--initial", 16)), int(opt("--seed", 1)),
                    evo, force="--force" in argv)
    if argv[0] == "advance":
        return advance(argv[1], argv[2], evo)
    if argv[0] == "status":
        return status()
    print(__doc__); return EXIT_REFUSED


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
