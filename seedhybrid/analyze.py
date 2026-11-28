"""Scores the campaign against PROTOCOL.md after every stage and writes seedhybrid/RESULTS.md.

    python seedhybrid/analyze.py <adjoint checkout> <setup geometry dir: contains geom-<lineage>/geom.tgz>

Inside-loop numbers (E1, E2 trajectories, E7) are computed here. Claims (E2 outcome, E3-E6) need the confirm stage's
fresh-mesh scores, which this page reports once they exist (seedhybrid/confirm/).
Per-child metrics (move RMS, DAFoam prediction) are cached in seedhybrid/lineages/<name>/child_metrics.json.
"""
from __future__ import annotations

import json, os, shutil, subprocess, sys, tarfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LIN, REL = HERE / "lineages", "seedhybrid_lineages"
rng = np.random.default_rng(20261128)


def kind(field):
    if field.startswith("dafnpz:"):
        return "daf"
    if field.startswith("vec:") and field.endswith("/guided.npy"):
        return "comb"
    if field.startswith("vec:") and field.endswith("/es.npy"):
        return "es"
    if field.startswith("rnd"):
        return "rnd"
    if field.startswith("esinpz:"):
        return "esi"
    return "other"


def metrics(adj, S, geomdir):
    """Move RMS and DAFoam prediction of every measured child (cached)."""
    sp = Path(adj) / "research" / "sidepod_opt"
    for name, Ls in S["lineages"].items():
        cache_p = LIN / name / "child_metrics.json"
        cache = json.loads(cache_p.read_text()) if cache_p.exists() else {}
        todo = {}
        for k, g in Ls.get("generations", {}).items():
            if "result" not in g:
                continue
            grad = next((f[7:] for f, a in g["children"].values() if f.startswith("dafnpz:")), None)
            if grad is None and not Ls.get("guides", True):             # no-adjoint arm: its twin's gradient, for E1 only
                twin = S["lineages"].get(Ls.get("base_of") or "", {}).get("generations", {}).get(k, {})
                grad = next((f[7:] for f, a in twin.get("children", {}).values() if f.startswith("dafnpz:")), None)
            for c, (f, a) in g["children"].items():
                if c not in cache:
                    todo[c] = {"field": f, "amp": a, "grad": grad}
        if not todo:
            continue
        for d in ("original", "original_rm"):
            shutil.rmtree(sp / "gate" / d, ignore_errors=True)
        with tarfile.open(Path(geomdir) / f"geom-{Ls.get('art', name)}" / "geom.tgz") as tf:
            tf.extractall(sp, filter="fully_trusted") if sys.version_info >= (3, 12) else tf.extractall(sp)
        shutil.copytree(LIN, sp / REL, dirs_exist_ok=True)
        (sp / "_cm_in.json").write_text(json.dumps(todo))
        r = subprocess.run([sys.executable, "-W", "ignore", str(HERE / "child_metrics.py"), "_cm_in.json", "_cm_out.json"],
                           cwd=sp, capture_output=True, text=True, env={**os.environ, "STEP_REGION": "1"})
        if r.returncode != 0:
            print(f"{name}: child metrics failed: {r.stderr[-800:]}", flush=True)
            continue
        cache.update(json.loads((sp / "_cm_out.json").read_text()))
        cache_p.write_text(json.dumps(cache, indent=1))


def paired(C, a, b, frames):
    """a, b: (art, chain, legal, bg_div, turb). Paired by frame: mean and SE of a - b over the frames both have."""
    def get(key):
        return {c["k"]: c["drag_mN"] for c in C.values() if (c["art"], c["chain"], c["legal"], c["bg_div"], c["turb"]) == key
                and c.get("drag_mN") is not None and c["k"] in frames}
    A, B = get(a), get(b)
    ks = sorted(set(A) & set(B))
    if not ks:
        return None
    d = np.array([A[k] - B[k] for k in ks])
    se = float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else float("nan")
    return float(d.mean()), se, len(d), float(np.mean([A[k] for k in ks])), float(np.mean([B[k] for k in ks]))


def claims(S, L):
    C = S.get("confirm") or {}
    if not C:
        L.append("\n## Claims (E2-E6)\n\nNot yet: the confirm stage runs when every lineage has stopped.")
        return
    U, LG = 12, "kOmegaSST"
    K4, K2 = (0, 1, 2, 3), (0, 1)
    L.append("\n## Claims (fresh meshes, production scheme; the audit protocol)\n")
    rows, finals = [], {}
    # E3: yardstick fidelity
    xs, ys, steps = [], [], []
    for name, Ls in S["lineages"].items():
        if not Ls.get("generations") or Ls["status"] in ("skipped_by_rule", "failed_setup", "pending"):
            continue
        art, levels = Ls["art"], [l for l in Ls["parent"].split("|") if l]
        cum, inside = 0.0, [0.0]
        for kk in sorted(Ls["generations"], key=int):
            r = Ls["generations"][kk].get("result") or {}
            if r.get("accepted"):
                cum += r["best_change_mN"]; inside.append(cum)
        prev = None
        for i in range(1, len(levels) + 1):
            ch = "|".join(levels[:i])
            fr = paired(C, (art, ch, False, U, LG), (art, "", False, U, LG), K2)
            if fr and i < len(inside):
                xs.append(inside[i]); ys.append(fr[0])
                if prev is not None or i == 1:
                    step_in = inside[i] - inside[i - 1]
                    step_fr = fr[0] - (prev if prev is not None else 0.0)
                    steps.append((step_in, step_fr))
                prev = fr[0]
        fin_u = paired(C, (art, Ls["parent"], False, U, LG), (art, "", False, U, LG), K4)
        fin_l = paired(C, (art, Ls["parent"], True, U, LG), (art, "", True, U, LG), K4)
        leg = paired(C, (art, Ls["parent"], True, U, LG), (art, Ls["parent"], False, U, LG), K4)
        finals[name] = (art, Ls["parent"], fin_l)
        f1 = fin_l and fin_l[0] <= -1.4 and abs(fin_l[0]) >= 3 * fin_l[1]
        kept = (fin_l[0] / fin_u[0]) if (fin_l and fin_u and fin_u[0] < 0) else None
        rows.append(f"| {name} | {fmt(fin_u[0] if fin_u else None)} | {fmt(fin_l[0] if fin_l else None)} +- "
                    f"{fmt(fin_l[1] if fin_l else None)} | {fmt(leg[0] if leg else None)} | "
                    f"{'-' if kept is None else f'{100 * kept:.0f} %'} | {'pass' if f1 else 'not met'} |")
    L.append("| lineage | final vs seed, unlegalised (frames 0-3) | final vs seed, legal (F1) | legalisation's own change | gain kept (E4, mark >= 70 %) | F1 |")
    L.append("|---|---|---|---|---|---|")
    L += rows
    if len(xs) >= 3:
        from scipy.stats import spearmanr
        rho = float(spearmanr(xs, ys)[0])
        big = [(a, b) for a, b in steps if a <= -1.0]
        agree = sum(np.sign(a) == np.sign(b) for a, b in big)
        L.append(f"\n- **E3a** Spearman (inside-loop cumulative change vs fresh-mesh change, {len(xs)} parents): {rho:+.2f} "
                 f"(mark >= 0.8): **{'pass' if rho >= 0.8 else 'fail'}**.")
        L.append(f"- **E3b** accepted steps <= -1.0 mN inside the loop that keep their sign on fresh meshes: {agree}/{len(big)} "
                 f"(mark >= 80 %): **{'pass' if big and agree / len(big) >= 0.8 else 'fail' if big else 'n/a'}**.")
    # E2
    for name, Ls in S["lineages"].items():
        if Ls.get("guides", True) or name not in finals:
            continue
        twin = Ls["base_of"]
        if twin not in finals:
            continue
        a, b = finals[twin], finals[name]
        d = paired(C, (a[0], a[1], True, U, LG), (b[0], b[1], True, U, LG), K4)
        ok = bool(d) and d[0] < 0 and d[1] > 0 and abs(d[0]) >= 3 * d[1]      # a real, clear difference
        L.append(f"- **E2** {twin} (guided) minus {name} (no adjoint), legal finals, frames 0-3: "
                 f"{fmt(d[0] if d else None)} +- {fmt(d[1] if d else None)} mN: **{'pass' if ok else 'not met'}**.")
    # F2 and E6
    best = min((v for v in finals.items() if v[1][2]), key=lambda kv: kv[1][2][3], default=None)
    if best:
        name, (art, ch, fl) = best
        d7 = paired(C, (art, ch, True, U, LG), ("d7", "", True, U, LG), K4)
        L.append(f"- **F2** best legal design: {name} ({fl[3]:.2f} mN mean). Against d7's legal seed re-scored here: "
                 f"{fmt(d7[0] if d7 else None)} +- {fmt(d7[1] if d7 else None)} mN; against 171.44 (before): {fl[3] - 171.44:+.2f}. "
                 f"**{'pass' if d7 and d7[0] < 0 and abs(d7[0]) >= 3 * d7[1] else 'not met'}**.")
    gl = [(n, v) for n, v in finals.items() if S["lineages"][n].get("guides", True) and v[2]]
    if len(gl) >= 2:
        L.append("- **E6** legal finals (mean of frames 0-3): " + ", ".join(f"{n} {v[2][3]:.2f}" for n, v in gl))
    # E5 robustness
    W = S.get("winner")
    if W:
        art, ch = W["art"], W["chain"]
        base = paired(C, (art, ch, True, U, LG), (art, "", True, U, LG), K2)
        fine = paired(C, (art, ch, True, 15.12, LG), (art, "", True, 15.12, LG), K2)
        lm = paired(C, (art, ch, True, U, "kOmegaSSTLM"), (art, "", True, U, "kOmegaSSTLM"), (0,))
        b0 = paired(C, (art, ch, True, U, LG), (art, "", True, U, LG), (0,))
        if base and fine:
            tol = 0.25 * abs(base[0]) if abs(base[0]) >= 10 else 2.65
            ok = np.sign(base[0]) == np.sign(fine[0]) and abs(fine[0] - base[0]) <= tol
            L.append(f"- **E5 finer mesh** (winner {W['lineage']} vs its seed, frames 0-1): scoring mesh {base[0]:+.2f}, finer "
                     f"{fine[0]:+.2f} (tolerance {tol:.2f}): **{'pass' if ok else 'fail'}**.")
        if lm and b0:
            ok = np.sign(lm[0]) == np.sign(b0[0]) and 0.5 <= lm[0] / b0[0] <= 2
            L.append(f"- **E5 transition model** (frame 0): SST {b0[0]:+.2f}, SSTLM {lm[0]:+.2f}: **{'same answer' if ok else 'different'}**.")
    # flow-state flags: a frame whose body drag sits > 1.0 mN from the design's median body drag
    by = {}
    for c in C.values():
        if c.get("body_mN") is not None and c["bg_div"] == U and c["turb"] == LG:
            by.setdefault((c["art"], c["chain"], c["legal"]), []).append((c["k"], c["body_mN"]))
    flags = [f"{a}/{'legal' if lg else 'morphed'}/{(ch[-30:] or 'seed')} frame {k}" for (a, ch, lg), v in by.items() if len(v) >= 3
             for k, b in v if abs(b - np.median([x for _, x in v])) > 1.0]
    L.append(f"- **Flow-state flags** (body drag > 1.0 mN from the design's median): {', '.join(flags) or 'none'}.")


def boot_mean(groups, n=4000):
    """Mean over all values, with a 95 % interval from resampling the groups (generations), not the values."""
    groups = [np.asarray(g, float) for g in groups if len(g)]
    if not groups:
        return None
    allv = np.concatenate(groups)
    bs = [np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))]).mean() for _ in range(n)]
    return float(allv.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), int(len(allv))


def fmt(x, nd=2):
    return "-" if x is None else f"{x:+.{nd}f}" if isinstance(x, float) else str(x)


def main(adj, geomdir):
    S = json.loads((HERE / "state.json").read_text())
    cfg = json.loads((HERE / "config.json").read_text())
    metrics(adj, S, geomdir)
    L = []
    L.append("# Seeded hybrid: results\n")
    L.append("Written by `analyze.py` after every stage. Marks: [PROTOCOL.md](PROTOCOL.md). **Inside-loop numbers are the "
             "steady model on each lineage's morphed base mesh: the loop's yardstick, not a claim** (claims: E3-E6, "
             "after the confirm stage).\n")
    L.append(f"State: stage `{S.get('stage')}`, generation {S.get('gen')}, setup run {S.get('setup_run')}.\n")
    # E0
    L.append("## E0. Setup\n")
    L.append("| lineage | arm | status | base drag (steady, mN) | steady_pass | DAFoam at seed |")
    L.append("|---|---|---|---|---|---|")
    for name, Ls in S["lineages"].items():
        b = Ls.get("base") or {}
        arm = "no adjoint" if not Ls.get("guides", True) else "guided"
        daf = "-" if not Ls.get("guides", True) else ("converged" if (LIN / name / "dafoam_g0.npz").exists() else "no gradient")
        L.append(f"| {name} | {arm} (base: {Ls.get('art', name)}) | {Ls['status']} | "
                 f"{b.get('D20_half_N', 0) * 1e3:.2f} | {b.get('steady_pass')} | {daf} |")
    # per-generation table
    L.append("\n## Generations (inside-loop yardstick)\n")
    L.append("| lineage | gen | amp mm | parent re-run mN | children: change vs re-run (mN) | accepted |")
    L.append("|---|---|---|---|---|---|")
    E1 = {t: [] for t in ("daf", "comb", "es", "rnd")}       # per generation lists of change per mm RMS
    winners, acc = {}, {}
    e1a_hits, e1a_n = 0, 0
    sign_ok, sign_n, rs = 0, 0, []
    traj = {}
    for name, Ls in S["lineages"].items():
        cm = json.loads((LIN / name / "child_metrics.json").read_text()) if (LIN / name / "child_metrics.json").exists() else {}
        cum = 0.0
        traj[name] = []
        for k in sorted(Ls.get("generations", {}), key=int):
            g = Ls["generations"][k]
            res = g.get("result")
            if not res:
                continue
            ch = res.get("children", {})
            z = res.get("zero")
            parts = []
            for c, (dv, se) in sorted(ch.items(), key=lambda t: t[1][0]):
                t = kind(g["children"][c][0])
                parts.append(f"{c.split('_g' + k + '_')[-1]} {dv:+.2f}")
            L.append(f"| {name} | {k} | {g['amp']:g} | {fmt(z[0] if z else None)} | {', '.join(parts) or res.get('note', '-')} "
                     f"| {(res.get('accepted') or '-').split('_g' + k + '_')[-1]} |")
            if res.get("accepted"):
                cum += res["best_change_mN"]
                t = kind(g["children"][res["accepted"]][0])
                acc.setdefault(name, []).append((t, res["best_change_mN"]))
            traj[name].append((int(k), cum))
            if not Ls.get("guides", True) or not ch:
                continue
            best = min(ch, key=lambda c: ch[c][0])
            winners[kind(g["children"][best][0])] = winners.get(kind(g["children"][best][0]), 0) + 1
            per = {t: [] for t in E1}
            for c, (dv, se) in ch.items():
                t = kind(g["children"][c][0])
                m = cm.get(c, {})
                if t in per and m.get("rms_mm"):
                    per[t].append(dv / m["rms_mm"])
            for t in E1:
                E1[t].append(per[t])
            dk = [ch[c][0] for c in ch if kind(g["children"][c][0]) == "daf" and float(g["children"][c][1]) == float(g["amp"])]
            rk = [ch[c][0] for c in ch if kind(g["children"][c][0]) == "rnd"]
            if dk and rk:
                e1a_n += 1
                e1a_hits += int(dk[0] < np.median(rk))
            pm = [(ch[c][0], cm[c]["pred_mN"], ch[c][1]) for c in ch if c in cm and cm[c].get("pred_mN") is not None]
            for meas, pred, se in pm:
                if abs(meas) > 3 * np.hypot(se, z[1] if z else 0.0):       # clear: beyond 3 SE of child and re-run
                    sign_n += 1
                    sign_ok += int(np.sign(meas) == np.sign(pred))
            if len(pm) >= 3:
                a = np.array(pm)
                rs.append(float(np.corrcoef(a[:, 0], a[:, 1])[0, 1]))
    # E1
    L.append("\n## E1. Direction credit (guided lineages; paired within each generation)\n")
    L.append("| direction | generations it won | mean change per mm RMS (95 % interval over generations) | children |")
    L.append("|---|---|---|---|")
    stats = {t: boot_mean(E1[t]) for t in E1}
    for t in E1:
        s = stats[t]
        L.append(f"| {t} | {winners.get(t, 0)} | " + (f"{s[0]:+.2f} ({s[1]:+.2f} to {s[2]:+.2f}) mN/mm | {s[3]}" if s else "- | 0") + " |")
    if e1a_n:
        # difference of pooled means (daf minus random, per mm RMS), 95 % interval from resampling generations
        pairs = [(d, r) for d, r in zip(E1["daf"], E1["rnd"]) if d and r]
        excl = None
        if pairs:
            diff = lambda P: np.mean(np.concatenate([p[0] for p in P])) - np.mean(np.concatenate([p[1] for p in P]))
            bs = [diff([pairs[i] for i in rng.integers(0, len(pairs), len(pairs))]) for _ in range(4000)]
            d0, lo, hi = diff(pairs), np.percentile(bs, 2.5), np.percentile(bs, 97.5)
            excl = bool(hi < 0)
            L.append(f"\n- daf minus random, mean change per mm RMS: {d0:+.2f} mN/mm (95 % interval {lo:+.2f} to {hi:+.2f}, "
                     f"{len(pairs)} generations)")
        L.append(f"\n- **E1a**: `daf` below the median random child in {e1a_hits}/{e1a_n} generations "
                 f"({100 * e1a_hits / e1a_n:.0f} %; mark >= 70 %); interval of the difference below 0: {excl}. Verdict: **{'pass' if e1a_hits / e1a_n >= 0.7 and excl else 'not met (yet)'}**.")
    if sign_n:
        L.append(f"- **E1b**: DAFoam's predicted sign matches the measured sign on {sign_ok}/{sign_n} clear children "
                 f"({100 * sign_ok / sign_n:.0f} %; mark >= 80 %); Pearson r per generation: "
                 f"{', '.join(f'{r:+.2f}' for r in rs) or '-'}.")
    if acc:
        L.append("\n**Accepted steps by direction type** (inside-loop change, mN):\n")
        for name, steps in acc.items():
            L.append(f"- {name}: " + ", ".join(f"{t} {d:+.2f}" for t, d in steps) + f"; total {sum(d for _, d in steps):+.2f}")
    # E2 trajectories
    L.append("\n## E2. Ablation pairs (inside-loop trajectory; the verdict needs the confirm stage)\n")
    for name, Ls in S["lineages"].items():
        if Ls.get("guides", True) or Ls["status"] == "skipped_by_rule":
            continue
        twin = Ls.get("base_of")
        a, b = dict(traj.get(twin, [])), dict(traj.get(name, []))
        ks = sorted(set(a) | set(b))
        L.append(f"- **{twin}** (guided) vs **{name}** (no adjoint), cumulative accepted change after each generation: " +
                 "; ".join(f"g{k}: {fmt(a.get(k))} vs {fmt(b.get(k))}" for k in ks))
    # E7
    L.append("\n## E7. Cost so far (CFD jobs)\n")
    L.append("| lineage | children solved | DAFoam gradients |")
    L.append("|---|---|---|")
    for name, Ls in S["lineages"].items():
        n_ch = sum(len(g["children"]) + 1 for g in Ls.get("generations", {}).values() if g.get("result"))
        n_d = len(list((LIN / name).glob("dafoam_g*.npz")))
        L.append(f"| {name} | {n_ch} | {n_d} |")
    L.append("\nSetup: one base state (fresh mesh, ~4.5 h) per seed, shared by an ablation arm and its guided twin.")
    claims(S, L)
    (HERE / "RESULTS.md").write_text("\n".join(L) + "\n")
    print("RESULTS.md written", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
