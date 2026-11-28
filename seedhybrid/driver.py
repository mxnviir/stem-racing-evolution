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
        S["lineages"][L["name"]] = {**new_lineage(cfg, L), "status": "setup"}
        (LIN / L["name"]).mkdir(parents=True, exist_ok=True)
    save(ST, S)
    cases = []
    for L in cfg["lineages"]:
        if L.get("base_of"):
            continue
        cases.append({"type": "base", "lineage": L["name"], "name": f"{L['name']}_base"})
        cases.append({"type": "dafoam", "lineage": L["name"], "name": f"{L['name']}_daf_g0", "parent": "", "tag": "g0"})
    save(REQ, {"stage": "setup", "gen": 0, "setup_run": None, "lineages": cfg["lineages"], "cases": cases})
    print(f"setup request: {len(cases)} cases")


def enable_ablation(S, log):
    """Ablation arms: enabled by the rule in PROTOCOL.md E2, from which lineages have an adjoint gradient."""
    for name, Ls in S["lineages"].items():
        if not Ls.get("base_of"):
            continue
        B = S["lineages"][Ls["base_of"]]
        Ls["base"], Ls["geom"] = B.get("base"), B.get("geom")
        cond = Ls.get("enable_if") or {}
        has = lambda n: bool(S["lineages"].get(n, {}).get("dafoam_latest"))
        ok = (B["status"] == "active" and (not cond.get("dafoam") or has(cond["dafoam"]))
              and (not cond.get("no_dafoam") or not has(cond["no_dafoam"])))
        Ls["status"] = "active" if ok else "skipped_by_rule"
        log["lineages"][name] = {"base_of": Ls["base_of"], "enabled": ok, "rule": cond}


def redaf_request():
    """Re-run only the setup's DAFoam gradients (the bases are kept), before generation 1 has run: setup run
    46113354235 produced no gradient on any seed, which would have left every lineage random-only and the
    adjoint-vs-no-adjoint ablation empty. Generation 1's unrun children are dropped and drawn again after."""
    cfg, S = load(CFG), load(ST)
    assert S["stage"] in ("gen", "redaf") and int(S["gen"]) <= 1, "only before generation 1 has run"
    cases = []
    for name, Ls in S["lineages"].items():
        Ls["generations"].pop("1", None)
        Ls["next_seed"] = 1000
        if Ls.get("base_of"):
            Ls["status"] = "pending"                     # re-decided by the rule once the gradients are in
            continue
        if Ls["status"] == "active" and Ls["guides"] and not Ls.get("dafoam_latest"):
            cases.append({"type": "dafoam", "lineage": name, "name": f"{name}_daf_g0", "parent": "", "tag": "g0"})
    S["stage"], S["gen"] = "redaf", 0
    save(ST, S)
    save(REQ, {"stage": "redaf", "gen": 0, "setup_run": S["setup_run"], "lineages": cfg["lineages"], "cases": cases})
    print(f"redaf request: {len(cases)} cases")


def new_lineage(cfg, L):
    """A lineage's state. base_of: shares that lineage's seed, geometry, base mesh and base state (an ablation arm);
    guides: false = no adjoint (random and recombined moves only)."""
    return {"seed": L.get("seed"), "base_of": L.get("base_of"), "art": L.get("base_of") or L["name"],
            "guides": L.get("guides", True), "enable_if": L.get("enable_if"), "status": "pending", "parent": "",
            "amp": cfg["amp0"], "dafoam_latest": None, "dafoam_at": "", "esi_latest": None, "next_seed": 1000,
            "generations": {}, "base": None}


def score_case(store, lineage, art, chain, legal, k, role, bg_div=12, turb="kOmegaSST"):
    """One fresh-mesh production score (the audit protocol); identical designs are scored once (key: art, chain,
    legal, frame, mesh, model) and credited to every lineage that reaches them."""
    key = (art, chain, legal, k, bg_div, turb)
    if key in store:
        if lineage not in store[key]["lineages"]:
            store[key]["lineages"].append(lineage)
        store[key]["role"] = sorted(set(store[key]["role"].split("+")) | {role})
        store[key]["role"] = "+".join(store[key]["role"])
        return
    import hashlib
    tag = f"{art}_{'L' if legal else 'U'}{hashlib.md5(chain.encode()).hexdigest()[:8]}_k{k}" + ("" if bg_div == 12 else f"_bg{bg_div:g}") + \
          ("" if turb == "kOmegaSST" else "_lm")
    store[key] = {"type": "score", "name": f"c_{tag}", "lineage": lineage, "lineages": [lineage], "art": art,
                  "chain": chain, "legal": legal, "k": k, "bg_div": bg_div, "turb": turb, "role": role}


def confirm_cases(S):
    """PROTOCOL.md E3-E6: every accepted parent (and the seed) unlegalised on frames 0-1 (E3); each final unlegalised
    and legalised on frames 0-3 (E4); each legal seed on frames 0-3 (E5, F1)."""
    store = {}
    for name, Ls in S["lineages"].items():
        if not Ls.get("generations") or Ls["status"] in ("skipped_by_rule", "failed_setup", "pending"):
            continue
        levels = [l for l in Ls["parent"].split("|") if l]
        prefixes = ["|".join(levels[:i]) for i in range(len(levels) + 1)]
        for ch in prefixes:
            for k in (0, 1):
                score_case(store, name, Ls["art"], ch, False, k, "E3")
        for k in (2, 3):                                  # the unlegalised seed on all four frames too (E4's gain kept)
            score_case(store, name, Ls["art"], "", False, k, "E4")
        for k in (0, 1, 2, 3):
            score_case(store, name, Ls["art"], Ls["parent"], False, k, "E4")
            score_case(store, name, Ls["art"], Ls["parent"], True, k, "E4+E5")
            score_case(store, name, Ls["art"], "", True, k, "seed")
    return list(store.values())


def robust_cases(S):
    """PROTOCOL.md E5 robustness: the best legal final and its legal seed on the finer mesh (frames 0-1) and with the
    transition model (frame 0)."""
    C = S.get("confirm", {})
    finals = {}
    for n, c in C.items():
        if c["legal"] and c["chain"] and c.get("drag_mN") and c["bg_div"] == 12 and c["turb"] == "kOmegaSST":
            finals.setdefault((c["art"], c["chain"]), []).append(c["drag_mN"])
    full = {key: np.mean(v) for key, v in finals.items() if len(v) == 4}
    if not full:
        return []
    art, chain = min(full, key=full.get)
    win = next(l for l, Ls in S["lineages"].items() if Ls.get("art") == art and Ls["parent"] == chain)
    S["winner"] = {"lineage": win, "art": art, "chain": chain, "legal_mean_mN": float(full[(art, chain)])}
    store = {}
    for ch, role in ((chain, "winner"), ("", "winner_seed")):
        for k in (0, 1):
            score_case(store, win, art, ch, True, k, role, bg_div=15.12)
        score_case(store, win, art, ch, True, 0, role, turb="kOmegaSSTLM")
    return list(store.values())


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
    for L in cfg["lineages"]:                                 # lineages added to the config after init (ablation arms)
        if L["name"] not in S["lineages"]:
            S["lineages"][L["name"]] = new_lineage(cfg, L)
        S["lineages"][L["name"]].setdefault("art", L.get("base_of") or L["name"])
        S["lineages"][L["name"]].setdefault("guides", L.get("guides", True))
    for name in S["lineages"]:
        (LIN / name).mkdir(parents=True, exist_ok=True)     # empty folders are not in git
    if stage in ("confirm", "robust"):
        C = S.setdefault("confirm", {})
        (HERE / "confirm").mkdir(exist_ok=True)
        for c in req["cases"]:
            r = load(results / c["name"] / "result.json")
            if r is not None:
                save(HERE / "confirm" / f"{c['name']}.json", r)
            C[c["name"]] = {**{k: c[k] for k in ("lineages", "art", "chain", "legal", "k", "bg_div", "turb", "role")},
                            "drag_mN": (r or {}).get("drag_mN"), "body_mN": (r or {}).get("body_drag_mN"),
                            "wheels_mN": (r or {}).get("wheels_drag_mN"), "lift_mN": (r or {}).get("lift_mN"),
                            "error": (r or {}).get("error", None if r else "no result")}
        S["confirm_runs"] = S.get("confirm_runs", []) + [int(run_id)]
        nxt = robust_cases(S) if stage == "confirm" else []
        S["stage"] = "robust" if nxt else "stopped"
        save(ST, S)
        save(REQ, {"stage": S["stage"], "gen": k, "setup_run": S["setup_run"], "lineages": cfg["lineages"], "cases": nxt})
        print(f"{stage} ingested; next: {S['stage']} ({len(nxt)} cases)")
        return bool(nxt)
    if stage == "setup":
        S["setup_run"] = int(run_id)
        for name, Ls in S["lineages"].items():
            if Ls.get("base_of"):
                continue
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
        enable_ablation(S, log)
    elif stage == "redaf":                              # the setup's adjoint gradients again (redaf_request); bases kept
        for name, Ls in S["lineages"].items():
            if Ls.get("base_of") or Ls["status"] != "active":
                continue
            daf = results / f"{name}_daf_g0" / "dafoam_sens.npz"
            if daf.exists():
                dst = LIN / name / "dafoam_g0.npz"; shutil.copy(daf, dst)
                Ls["dafoam_latest"] = f"{REL}/{name}/dafoam_g0.npz"
            log["lineages"][name] = {"dafoam": daf.exists()}
        # a gradient still missing (redaf run 46137619748: team's GMRES stopped at 2000 iterations, still falling)
        # is re-run once more for that lineage only, with a longer Krylov budget, before generation 1 is drawn:
        # without it the lineage's adjoint arm would run random-only and its ablation would compare nothing
        miss = [n for n, Ls in S["lineages"].items()
                if not Ls.get("base_of") and Ls["status"] == "active" and Ls["guides"] and not Ls.get("dafoam_latest")]
        S["redaf_runs"] = S.get("redaf_runs", []) + [int(run_id)]
        if miss and len(S["redaf_runs"]) < 2:
            cases = [{"type": "dafoam", "lineage": n, "name": f"{n}_daf_g0", "parent": "", "tag": "g0",
                      "gmres_iters": "6000", "gmres_restart": "300"} for n in miss]
            log["redaf_again"] = miss
            S.setdefault("log", []).append(log)
            save(ST, S)
            save(REQ, {"stage": "redaf", "gen": 0, "setup_run": S["setup_run"], "lineages": cfg["lineages"], "cases": cases})
            print(f"redaf again for {miss} (GMRES 6000 iterations, restart 300)")
            return True
        enable_ablation(S, log)
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
                sp = with_geom(adj, Ls["art"], results / "_setup")
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
        if Ls["guides"] and Ls["dafoam_latest"]:
            kids["daf"] = (f"dafnpz:{Ls['dafoam_latest']}", amp)
            if k1 == 1:
                kids["dafh"] = (f"dafnpz:{Ls['dafoam_latest']}", amp / 2)
        if cfg.get("use_esi") and Ls["esi_latest"]:
            kids["esi"] = (f"esinpz:{Ls['esi_latest']}", amp)
        gdir = LIN / name / f"g{k1}"
        for arm, tag in (("guided", "comb"), ("es", "es")):
            if (gdir / f"{arm}.npy").exists():
                kids[tag] = (f"vec:{REL}/{name}/g{k1}/{arm}.npy", amp)
        # equal children per generation in every arm (PROTOCOL.md E2): random moves fill up to the target
        target = cfg["children_gen1"] if k1 == 1 else cfg["children_after"]
        n_rnd = max(cfg["n_rnd"], target - len(kids))
        s0 = Ls["next_seed"]
        kids.update({f"rnd{s}": (f"rnd{s}", amp) for s in range(s0, s0 + n_rnd)})
        Ls["next_seed"] = s0 + n_rnd
        Ls["generations"][str(k1)] = {"parent": Ls["parent"], "amp": amp,
                                      "children": {f"{name}_g{k1}_{n}": [f, a] for n, (f, a) in kids.items() if n != "zero"}}
        for n, (f, a) in kids.items():
            c = {"type": "child", "lineage": name, "art": Ls["art"], "name": f"{name}_g{k1}_{n}", "chainparent": Ls["parent"],
                 "field": f, "amp": f"{a:g}"}
            if n == "zero" and cfg.get("use_esi") and cfg.get("esi_after"):
                c["esi_after"] = str(cfg["esi_after"])
            cases.append(c)
        if Ls["guides"] and Ls.get("dafoam_at") != Ls["parent"]:   # a fresh gradient only where the parent has moved
            cases.append({"type": "dafoam", "lineage": name, "art": Ls["art"], "name": f"{name}_daf_g{k1}",
                          "parent": Ls["parent"], "tag": f"g{k1}"})
            Ls["dafoam_at"] = Ls["parent"]
    if not cases and not S.get("confirm"):
        cc = confirm_cases(S)
        S["stage"], S["gen"] = "confirm", k
        save(ST, S)
        save(REQ, {"stage": "confirm", "gen": k, "setup_run": S["setup_run"], "lineages": cfg["lineages"], "cases": cc})
        print(f"all lineages stopped; confirm stage: {len(cc)} cases")
        return True
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
    elif sys.argv[1] == "redaf":
        redaf_request()
    elif sys.argv[1] == "advance":
        sys.exit(0 if advance(*sys.argv[2:5]) else 3)
