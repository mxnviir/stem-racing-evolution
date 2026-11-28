"""Seeded hybrid lineages: the driver (seedhybrid/README.md). Runs in the workflow's collect job (or by hand).

    python seedhybrid/driver.py init                       -> state.json + request.json (stage "setup")
    python seedhybrid/driver.py advance <results dir> <run id> <adjoint checkout>
        ingest the finished stage's results, then write the next stage's request.json (or stop)

The generation rule is the original loop's (stem-racing-adjoint research/sidepod_opt/hybrid/loop.py), per lineage:
- every child of generation k is the parent plus ONE move, chained morph of the lineage's own base mesh;
- the parent is re-run as <L>_g<k>_zero; the best child becomes the parent if it is below the re-run by >= 3 SE
  (whole-signal SE of both) AND by >= FLOOR mN; else the parent stays;
- directions: the newest DAFoam gradient at (or near) the parent, the newest ESI map, steepest descent in the span of
  the previous generation's measured moves (all: 'comb'; random only: 'es'), N_RND new random smooth fields.
Seeded version's additions: per-lineage base surface/mesh/state; moves confined to the legal sidepod region
(STEP_REGION=1); AMP halves when a generation accepts nothing; a lineage stops at max_gen, at AMP < AMP_MIN, or when
every child folds the morph (the morph's range is spent: re-meshing, "phase B", is not built yet).
"""
from __future__ import annotations

import json, re, shutil, subprocess, sys, tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CFG, ST, REQ = HERE / "config.json", HERE / "state.json", HERE / "request.json"
LIN = HERE / "lineages"
REL = "seedhybrid_lineages"          # the lineages folder as seen from the adjoint checkout's research/sidepod_opt


def load(p, default=None):
    return json.loads(p.read_text()) if p.exists() else default


def save(p, d):
    p.write_text(json.dumps(d, indent=1))


def init():
    cfg = load(CFG)
    S = {"stage": "setup", "gen": 0, "setup_run": None, "lineages": {}}
    for L in cfg["lineages"]:
        S["lineages"][L["name"]] = {"seed": L["seed"], "status": "setup", "parent": "", "amp": cfg["amp0"],
                                    "dafoam_latest": None, "dafoam_at": "", "esi_latest": None, "next_seed": 1000,
                                    "generations": {}, "base": None}
        (LIN / L["name"]).mkdir(parents=True, exist_ok=True)
    save(ST, S)
    cases = []
    for L in cfg["lineages"]:
        cases.append({"type": "base", "lineage": L["name"], "name": f"{L['name']}_base"})
        cases.append({"type": "dafoam", "lineage": L["name"], "name": f"{L['name']}_daf_g0", "parent": "", "tag": "g0"})
    save(REQ, {"stage": "setup", "gen": 0, "setup_run": None, "lineages": cfg["lineages"], "cases": cases})
    print(f"setup request: {len(cases)} cases")


def child_result(R, name):
    """(mean of the last 404 iterations, whole-signal SE) of a child, as hybrid/loop.py reads it, or None."""
    j, f = R / f"bump_{name}.json", R / "forces" / f"bump_{name}" / "forces.dat"
    if not (j.exists() and f.exists()):
        return None
    J = json.loads(j.read_text())
    if J.get("morph_fallback"):
        return None
    rows = [[float(x) for x in re.sub(r"[()]", " ", l).split()] for l in open(f) if l.strip() and not l.startswith("#")]
    if len(rows) < 404:
        return None
    return float(np.mean([r[1] for r in rows][-404:]) * 1e3), float(J["force_mean_stderr"]) * 1e3


def with_geom(adj, L, setup_run_dir):
    """Put lineage L's base surfaces at gate/original(_rm) of the adjoint checkout (combine.py reads them)."""
    sp = Path(adj) / "research" / "sidepod_opt"
    for d in ("original", "original_rm"):
        shutil.rmtree(sp / "gate" / d, ignore_errors=True)
    with tarfile.open(Path(setup_run_dir) / f"geom-{L}" / "geom.tgz") as tf:
        tf.extractall(sp, filter="fully_trusted") if sys.version_info >= (3, 12) else tf.extractall(sp)
    return sp


def advance(results, run_id, adj):
    cfg, S, req = load(CFG), load(ST), load(REQ)
    results = Path(results)
    stage, k = req["stage"], int(req["gen"])
    log = {"stage": stage, "gen": k, "run": run_id, "lineages": {}}
    if stage == "setup":
        S["setup_run"] = int(run_id)
        for name, Ls in S["lineages"].items():
            b = load(results / f"{name}_base" / "base.json")
            g = load(results / f"{name}_geom" / "geom.json")
            daf = results / f"{name}_daf_g0" / "dafoam_sens.npz"
            Ls["base"], Ls["geom"] = b, g
            if daf.exists():
                dst = LIN / name / "dafoam_g0.npz"; shutil.copy(daf, dst)
                Ls["dafoam_latest"] = f"{REL}/{name}/dafoam_g0.npz"
            ok = bool(b and b.get("run_dir_kept"))
            Ls["status"] = "active" if ok else "failed_setup"
            log["lineages"][name] = {"base_ok": ok, "steady_pass": (b or {}).get("steady_pass"),
                                     "base_D20_mN": (b or {}).get("D20_half_N", 0) * 1e3, "dafoam": daf.exists()}
    else:
        for name, Ls in S["lineages"].items():
            gen = Ls["generations"].get(str(k))
            if Ls["status"] != "active" or gen is None:
                continue
            R = results / "bumps"
            z = child_result(R, f"{name}_g{k}_zero")
            res = {n: child_result(R, n) for n in gen["children"]}
            ok = {n: v for n, v in res.items() if v is not None}
            rec = {"zero": z, "children": {n: (v[0] - z[0], v[1]) for n, v in ok.items()} if z else {},
                   "failed": [n for n in res if res[n] is None], "accepted": None}
            keep = LIN / name / "results"                  # keep each child's numbers and force history with the lineage
            for n in list(gen["children"]) + [f"{name}_g{k}_zero"]:
                for src, dst in ((R / f"bump_{n}.json", keep / f"bump_{n}.json"),
                                 (R / "forces" / f"bump_{n}" / "forces.dat", keep / "forces" / f"bump_{n}.dat")):
                    if src.exists():
                        dst.parent.mkdir(parents=True, exist_ok=True); shutil.copy(src, dst)
            if z is None or not ok:
                rec["note"] = "parent re-run missing" if z is None else "no child survived (morph folded)"
                Ls["amp"] /= 2
            else:
                best = min(ok, key=lambda n: ok[n][0])
                d = ok[best][0] - z[0]
                if d <= -3 * np.hypot(ok[best][1], z[1]) and d <= -cfg["floor_mN"]:
                    f, a = gen["children"][best]
                    Ls["parent"] = (Ls["parent"] + "|" if Ls["parent"] else "") + f"{f}@{a}"
                    rec["accepted"] = best
                    Ls.setdefault("drag_steady_mN", []).append(ok[best][0])
                else:
                    Ls["amp"] /= 2
                rec["best"], rec["best_change_mN"] = best, d
            # fresh gradients measured at this generation's parent: guides for the next one
            daf = results / f"{name}_daf_g{k}" / "dafoam_sens.npz"
            if daf.exists():
                shutil.copy(daf, LIN / name / f"dafoam_g{k}.npz"); Ls["dafoam_latest"] = f"{REL}/{name}/dafoam_g{k}.npz"
            esi = results / f"{name}_g{k}_zero" / "esi_sens.npz"
            if esi.exists():
                shutil.copy(esi, LIN / name / f"esi_g{k}.npz"); Ls["esi_latest"] = f"{REL}/{name}/esi_g{k}.npz"
            gen["result"] = rec
            # measured directions of generation k (moves of single fields, not chains) -> combined directions for k+1
            gdir = LIN / name / f"g{k + 1}"; gdir.mkdir(parents=True, exist_ok=True)
            meas = {n: {"field": gen["children"][n][0], "amp_mm": gen["children"][n][1], "dF_mN": dv}
                    for n, (dv, s) in rec["children"].items() if not n.endswith("_zero")}
            save(gdir / "measured.json", meas)
            if len(meas) >= 2:
                # combine.py writes under the adjoint checkout (it resolves "vec:" paths there): run it on a copy of
                # this repo's lineages folder placed at <checkout>/research/sidepod_opt/seedhybrid_lineages, copy back
                sp = with_geom(adj, name, results / "_setup")
                shutil.copytree(LIN, sp / REL, dirs_exist_ok=True)
                rel = f"{REL}/{name}/g{k + 1}"
                r = subprocess.run([sys.executable, "-W", "ignore", "hybrid/combine.py", f"{rel}/measured.json", rel],
                                   cwd=sp, capture_output=True, text=True, env={**__import__("os").environ, "STEP_REGION": "1"})
                shutil.copytree(sp / rel, gdir, dirs_exist_ok=True)
                (gdir / "combine.log").write_text(r.stdout[-4000:] + r.stderr[-4000:])
            if k >= cfg["max_gen"]:
                Ls["status"] = "done_max_gen"
            elif Ls["amp"] < cfg["amp_min"]:
                Ls["status"] = "done_amp"
            log["lineages"][name] = rec
            print(f"{name} g{k}: zero {z}; best {rec.get('best')} {rec.get('best_change_mN')}; accepted {rec['accepted']}; "
                  f"parent {Ls['parent']}; amp {Ls['amp']}; status {Ls['status']}", flush=True)
    S.setdefault("log", []).append(log)
    # next generation
    k1 = k + 1
    cases = []
    for name, Ls in S["lineages"].items():
        if Ls["status"] != "active":
            continue
        amp = Ls["amp"]
        kids = {"zero": ("zero", 0.0)}
        if Ls["dafoam_latest"]:
            kids["daf"] = (f"dafnpz:{Ls['dafoam_latest']}", amp)
            if k1 == 1:
                kids["dafh"] = (f"dafnpz:{Ls['dafoam_latest']}", amp / 2)
        if Ls["esi_latest"]:
            kids["esi"] = (f"esinpz:{Ls['esi_latest']}", amp)
        gdir = LIN / name / f"g{k1}"
        for arm, tag in (("guided", "comb"), ("es", "es")):
            if (gdir / f"{arm}.npy").exists():
                kids[tag] = (f"vec:{REL}/{name}/g{k1}/{arm}.npy", amp)
        n_rnd = cfg["n_rnd"] + (2 if k1 == 1 and not Ls["dafoam_latest"] else 0)
        s0 = Ls["next_seed"]
        kids.update({f"rnd{s}": (f"rnd{s}", amp) for s in range(s0, s0 + n_rnd)})
        Ls["next_seed"] = s0 + n_rnd
        Ls["generations"][str(k1)] = {"parent": Ls["parent"], "amp": amp,
                                      "children": {f"{name}_g{k1}_{n}": [f, a] for n, (f, a) in kids.items() if n != "zero"}}
        for n, (f, a) in kids.items():
            c = {"type": "child", "lineage": name, "name": f"{name}_g{k1}_{n}", "chainparent": Ls["parent"],
                 "field": f, "amp": f"{a:g}"}
            if n == "zero" and cfg.get("esi_after"):
                c["esi_after"] = str(cfg["esi_after"])
            cases.append(c)
        if Ls.get("dafoam_at") != Ls["parent"]:          # a fresh gradient only where the parent has moved since the last try
            cases.append({"type": "dafoam", "lineage": name, "name": f"{name}_daf_g{k1}", "parent": Ls["parent"], "tag": f"g{k1}"})
            Ls["dafoam_at"] = Ls["parent"]
    S["stage"], S["gen"] = ("gen", k1) if cases else ("stopped", k)
    save(ST, S)
    if cases:
        save(REQ, {"stage": "gen", "gen": k1, "setup_run": S["setup_run"], "lineages": cfg["lineages"], "cases": cases})
        print(f"generation {k1}: {len(cases)} cases")
    else:
        save(REQ, {"stage": "stopped", "gen": k, "setup_run": S["setup_run"], "lineages": cfg["lineages"], "cases": []})
        print("all lineages stopped")
    return bool(cases)


if __name__ == "__main__":
    if sys.argv[1] == "init":
        init()
    elif sys.argv[1] == "advance":
        sys.exit(0 if advance(*sys.argv[2:5]) else 3)
