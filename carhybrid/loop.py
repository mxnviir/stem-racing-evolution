"""P2F, the fine-adjoint loop (carhybrid/FINE_ADJOINT.md "The loop"): a stage machine, as seedhybrid/driver.py.

    python carhybrid/loop.py init <layout> <grid> [--arm adjoint|random] [--radius 0.25] [--seed 20261128]
                                  [--max-iter 20] [--bg-div 12.889] [--runner bigmem] [--np N] [--force]
                                  [--wing default|rw_chord,rwsup_out,... --evo <evolution checkout>]
        --wing: the rear wing's and its struts' numbers join the body's (x = [u; s], s in [-1, 1] over each number's
        layout_params.json range, started at the layout's own values); the gradient stage adds wing_grad.py's dT/ds
        a new arm at u = 0 (state carhybrid/loop/state.json) and its first stage in carhybrid/request.json
    python carhybrid/loop.py advance <results dir> <run id> [--evo <evolution checkout>] [--arm <arm>]
        ingest the finished stage (results dir = carhybrid/results/<batch>, the batch this arm is waiting on), write the
        next stage's cases into request.json; exit 0 if there is a next stage, 3 if the arm stopped (2: refused)
    python carhybrid/loop.py reissue <arm>      the arm's pending stage again (a new batch name: <batch>_r<n>)
    python carhybrid/loop.py status

<layout>: as make_car.py ("team" | genes json | Lxxxx). <grid>: the free-form grid "gu,gv,gw,gb" (5,2,3,4 = stage 2's
45 numbers; 24,9,12,9 = the seed's own 720). Every case carries "grid", so the build and grad.py use the same one.

Per iteration i of an arm, at its current numbers u_i:
  grad  (adjoint arm only) one case {"type": "fine_adjoint", "u": u_i, "k": i mod 4, "grid", "bg_div", "runner"} ->
        <case>/grad.json (carhybrid/grad.py). g_i = its race-time gradient dT_du_s (with lead on board grad.py has
        already zeroed dT/dm: drag + lift only). Gates: no grad.json, a non-finite or wrong-length gradient, or a
        gate known to be missed (fine_adjoint.json F1_pass / F2_pass / F3_pass false, steady.json "pass" false,
        dafoam_sens.json "gmres": {f: {"converged": false}}) stops the arm (a gate missed stops the line it gates).
  line  candidates u_i + a * r * d_i (a = 0.5, 1, 2; r the trust radius) clipped to the unit box [-1, 1], plus u_i
        itself ("base"), each {"type": "score"} on frames 0 and 1. d_i = -g / |g| projected on the box (components
        pushing out of a face u_j = +-1 dropped, then unit length); the random arm draws d_i as a unit normal vector
        (numpy default_rng([seed, i]), recorded) and skips the grad stage.
Race time of a candidate on frame k: simscale_layout.race_time(its own rec.json, 2 drag_mN / 1e3, 2 lift_mN / 1e3)
(whole car, lift + up: carhybrid/grad.py's docstring). A candidate is failed if any of its frames has no result,
an "error", or an illegal / failed build. Paired by frame: Delta_k = T_c,k - T_base,k, Delta = mean over k;
SE = max(SE from the block SEs: sqrt(sum_k se_c,k^2 + se_base,k^2) / K with se = dT/dD20 x 2 x drag_block_se_mN,
half the frame difference |Delta_0 - Delta_1| / 2 (for K frames: std(Delta_k) / sqrt(K))).
Accept the best (lowest Delta) non-failed candidate if Delta <= -3 SE and Delta <= -0.1 ms: u_{i+1} = its u. The
radius grows x 1.5 when the largest step (a = 2) is accepted, halves when nothing is accepted. Stop at r < 0.02 or
after 20 iterations (line stages).
Decisions this file makes that FINE_ADJOINT.md does not fix: carhybrid/FINE_ADJOINT.md "Decisions".
"""
from __future__ import annotations

import json, os, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = Path(os.environ.get("CARHYBRID_LOOP_DIR", HERE / "loop"))
REQ = Path(os.environ.get("CARHYBRID_REQUEST", HERE / "request.json"))

STEPS = (("a05", 0.5), ("a1", 1.0), ("a2", 2.0))
FRAMES = (0, 1)
N_SE, FLOOR_S = 3.0, 0.1e-3
GROW, SHRINK, R_MIN, R_MAX = 1.5, 0.5, 0.02, 1.0
RADIUS0, MAX_ITER, SEED = 0.25, 20, 20261128
BG_DIV = 12.889                    # production background division for the whole car (FINE_ADJOINT.md interfaces)
EXIT_NEXT, EXIT_STOP, EXIT_REFUSED = 0, 3, 2
DEFAULT_WING = ("rw_chord", "rw_t", "rw_half_span", "rw_x_le", "rw_lift", "rw_camber", "rw_camber_pos", "rw_rise",
                "rw_rise_p", "rw_sweep", "rw_sweep_p", "rw_tip_chord_ratio", "rw_tip_t_ratio",
                "rwsup_out", "rwsup_fwd", "rwsup_leave", "rwsup_arrive", "rwsup_chord", "rwsup_thick", "rwsup_flare",
                "rwsup_round", "rwsup_mid_x", "rwsup_mid_y", "rwsup_mid_z", "rwsup_swell")   # = wing_grad.DEFAULT_NAMES


# ------------------------------------------------------------------ files
def load(p, default=None):
    p = Path(p)
    return json.loads(p.read_text()) if p.exists() else default


def save(p, d):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(d, indent=1) + "\n")


def state_path():
    return LOOP / "state.json"


def grid_n(grid):
    g = [int(v) for v in str(grid).split(",")]
    gu, gv, gw = g[:3]
    gb = g[3] if len(g) > 3 else gv
    return gu * gv + gu * gb + gu * gw          # freeform.n_genes: top gu x gv, bottom gu x gb, side gu x gw


def u_str(u):
    u = np.asarray(u, float)
    return "zero" if not np.any(u) else json.dumps([round(float(v), 6) for v in u])


# ------------------------------------------------------------------ the rear wing's and struts' numbers (optional)
# An arm started with --wing carries s (one number per name in [-1, 1]: its place in layout_params.json's lo..hi)
# beside the body's u. The line works on x = [u; s]; the gradient stage adds carhybrid/wing_grad.py's dT/ds.
def nw(A):
    return len(A.get("wing_names") or [])


def xvec(A):
    return [float(v) for v in A["u"]] + [float(v) for v in (A.get("s") or [])]


def split_x(A, x):
    x = [float(v) for v in x]
    return x[:A["n"]], x[A["n"]:]


def wing_json(A, s):
    out = {}
    for nm, v in zip(A["wing_names"], s):
        lo, hi = A["wing_ranges"][nm]
        out[nm] = round(0.5 * (lo + hi) + float(v) * 0.5 * (hi - lo), 6)
    return json.dumps(out)


def wing_start(evo, layout, names):
    """(ranges, s0) of the named numbers on this layout: layout_params.json's lo..hi, the decoded layout's own values."""
    sys.path.insert(0, str(HERE))
    import make_car as MC
    evo = Path(evo).resolve()
    if str(evo) not in sys.path:
        sys.path.insert(0, str(evo))
    import repo_paths; repo_paths.add()
    import layout_genes as LG
    P = {q["name"]: q for q in json.loads((evo / "search" / "layout" / "layout_params.json").read_text())["params"]}
    rw = LG.decode(MC.genes_of(evo, layout))["rwing"]
    ranges, s0 = {}, []
    for nm in names:
        lo, hi = float(P[nm]["lo"]), float(P[nm]["hi"])
        if hi <= lo:
            raise SystemExit(f"loop: {nm} has a range of zero width (mode {P[nm]['mode']})")
        th = rw[MC.RW_KEYS[nm]] if nm in MC.RW_KEYS else rw.get("support", {}).get(MC.SUP_KEYS[nm[6:]], P[nm]["value"])
        ranges[nm] = (lo, hi)
        s0.append(round(float(np.clip((float(th) - 0.5 * (lo + hi)) / (0.5 * (hi - lo)), -1.0, 1.0)), 6))
    return ranges, s0


def batch_name(A, stage, retry=0):
    return f"loop_{A['arm']}_{A['i']}_{stage}" + (f"_r{retry}" if retry else "")


# ------------------------------------------------------------------ race time (the evolution race model)
_SL = None


def race_eval(rec, D20, L20, evo):
    """(race time s, dT/dD20 s per N) of a built car (its own mass and layout) at whole-car D20, L20 (N, lift + up).
    Replaced in the tests by a stand-in."""
    global _SL
    if _SL is None:
        sys.path.insert(0, str(HERE))
        from grad import evo_import
        _SL = evo_import(evo)[1]
    from grad import race_fn
    RT = race_fn(_SL, rec)                                   # simscale_layout.race_time on one jitted model
    h = 1e-3
    return RT(D20, L20), (RT(D20 + h, L20) - RT(D20 - h, L20)) / (2 * h)


_race_eval_evo = race_eval


# ------------------------------------------------------------------ stages
def project(u, d):
    """d with the components that push out of an active face of the unit box dropped, at unit length (None if 0)."""
    d = np.asarray(d, float).copy()
    d[(u >= 1 - 1e-9) & (d > 0)] = 0.0
    d[(u <= -1 + 1e-9) & (d < 0)] = 0.0
    nd = np.linalg.norm(d)
    return d / nd if nd > 0 and np.isfinite(nd) else None


def grad_cases(A, batch):
    c = {"name": f"{batch}_g", "type": "fine_adjoint", "layout": A["layout"], "u": u_str(A["u"]), "k": A["i"] % 4,
         "bg_div": A["bg_div"], "grid": A["grid"], "runner": A["runner"]}
    if A.get("np"):
        c["np"] = int(A["np"])
    if nw(A):
        c["wing"], c["wing_grad"] = wing_json(A, A["s"]), ",".join(A["wing_names"])
    return [c]


def line_cases(A, batch):
    x, d, r = np.asarray(xvec(A), float), np.asarray(A["direction"], float), A["radius"]
    cands = {"base": x}
    for tag, a in STEPS:
        cands[tag] = np.clip(x + a * r * d, -1.0, 1.0)
    A["candidates"] = {tag: [round(float(v), 6) for v in xc] for tag, xc in cands.items()}   # x = [u; s]
    cases = []
    for tag, xc in A["candidates"].items():
        uc, sc = split_x(A, xc)
        for k in FRAMES:
            c = {"name": f"{batch}_{tag}_k{k}", "type": "score", "layout": A["layout"], "k": k, "u": u_str(uc),
                 "grid": A["grid"]}
            if nw(A):
                c["wing"] = wing_json(A, sc)
            cases.append(c)
    return cases


def random_direction(A):
    rng = np.random.default_rng([int(A["seed"]), int(A["i"])])
    return rng.standard_normal(A["n"] + nw(A))


def write_request(A, stage, retry=0):
    batch = batch_name(A, stage, retry)
    cases = grad_cases(A, batch) if stage == "grad" else line_cases(A, batch)
    R = load(REQ, {}) or {}
    R.setdefault("batches", {})[batch] = cases
    R["launch"] = batch
    R.setdefault("_loop_note", "P2F (carhybrid/FINE_ADJOINT.md, carhybrid/loop.py): batches loop_<arm>_<i>_<stage> "
                               "are written by loop.py advance; its state is carhybrid/loop/state.json")
    save(REQ, R)
    A["stage"], A["pending"], A["pending_cases"], A["retry"] = stage, batch, [c["name"] for c in cases], retry
    return batch, cases


def next_stage(A):
    """Set up the stage after a finished line (or at init): grad (adjoint) or line with a new random direction."""
    if A["arm"] == "adjoint":
        return write_request(A, "grad")
    d0 = random_direction(A)
    d = project(np.asarray(xvec(A), float), d0)
    if d is None:
        return stop(A, "no feasible direction (every component at a face of the box)")
    A["direction"] = d.tolist()
    A["direction_info"] = {"kind": "random", "rng": [int(A["seed"]), int(A["i"])]}
    return write_request(A, "line")


def stop(A, why):
    A["stage"], A["status"], A["stop_reason"], A["pending"], A["pending_cases"] = "stopped", "stopped", why, None, []
    return None, []


# ------------------------------------------------------------------ ingest
def ingest_grad(A, R):
    name = A["pending_cases"][0]
    D = Path(R) / name
    g = load(D / "grad.json")
    st, ds, res = load(D / "steady.json"), load(D / "dafoam_sens.json"), load(D / "result.json")
    fa = load(D / "fine_adjoint.json") or {}                 # carhybrid/fine_adjoint.sh's summary: F1_pass, F2_pass, F3_pass
    rec = {"i": A["i"], "batch": A["pending"], "case": name, "k": A["i"] % 4, "u": A["u"],
           "score": {k: (res or {}).get(k) for k in ("drag_mN", "drag_block_se_mN", "lift_mN")} if res else None,
           "gates": {k: fa.get(k) for k in ("F1_dafoam_vs_pipeline_drag", "F1_pass", "F2_span_mN", "F2_pass",
                                            "F3_gmres_converged", "F3_pass")}}
    why = None
    # a gate stops the arm only when it is known to be missed (False); unknown (None, absent) does not
    if fa.get("F1_pass") is False:
        why = f"F1 (DAFoam's primal vs the pipeline's drag) missed: {fa.get('F1_dafoam_vs_pipeline_drag')}"
    elif fa.get("F2_pass") is False or (st is not None and st.get("pass") is False):
        why = f"F2 (steadiness) missed: span {fa.get('F2_span_mN', (st or {}).get('span_mN'))} mN"
    if why is None:
        conv = (ds or {}).get("converged")
        gm = {f: v.get("converged") for f, v in ((ds or {}).get("gmres") or {}).items() if isinstance(v, dict)}
        if fa.get("F3_pass") is False or conv is False or any(v is False for v in gm.values()) or \
                (isinstance(conv, dict) and not all(bool(v) for v in conv.values())):
            why = f"F3 (adjoint convergence) missed: {fa.get('F3_gmres_converged') or gm or conv}"
    if why is None and g is None:
        why = "no grad.json (the fine_adjoint job or its map failed)"
    if why is None:
        dT = np.asarray(g.get("dT_du_s") or [], float)
        if len(dT) != A["n"] or not np.all(np.isfinite(dT)) or not np.any(dT):
            why = f"unusable gradient (length {len(dT)}, grid {A['n']})"
    if why is None and nw(A):
        wg = load(D / "wing_grad.json")
        if wg is None or list(wg.get("names") or []) != list(A["wing_names"]):
            why = "no wing_grad.json for this arm's wing numbers (carhybrid/wing_grad.py failed or was not run)"
        else:
            ds_ = np.asarray([0.0 if v is None else float(v) for v in wg["dT_ds_s"]], float)
            if not np.all(np.isfinite(ds_)):
                why = "unusable wing gradient"
            else:
                rec["wing"] = {"branch_change": wg.get("branch_change"), "norm_dT_ds_s": wg.get("norm_dT_ds_s"),
                               "base_topology": wg.get("base_topology")}
                dT = np.r_[dT, ds_]
    if why is not None:
        rec["gate"] = why
        A["history"].append({"stage": "grad", **rec})
        return stop(A, f"iteration {A['i']}: {why}")
    # gradients taken at this same u on earlier frames (after a rejected line) are measurements of the same
    # derivative: the line uses their mean (FINE_ADJOINT.md Decisions)
    same = [np.asarray(x, float) for x in A.get("grads_at_u", [])]
    if same:
        rec["cos_with_previous_at_u"] = float(dT @ same[-1] / np.linalg.norm(dT) / np.linalg.norm(same[-1]))
    same.append(dT)
    A["grads_at_u"] = [x.tolist() for x in same]
    gm = np.mean(same, 0)
    rec.update(norm_dT_du_s=float(np.linalg.norm(dT)), n_averaged=len(same),
               **{k: g.get(k) for k in ("D20_N_used", "L20_N_used", "lead_g", "norm_terms_s", "largest_share")})
    d = project(np.asarray(xvec(A), float), -gm)
    A["history"].append({"stage": "grad", **rec})
    if d is None:
        return stop(A, f"iteration {A['i']}: no feasible descent direction in the box")
    A["direction"] = d.tolist()
    A["direction_info"] = {"kind": "adjoint", "n_averaged": len(same), "norm_mean_dT_du_s": float(np.linalg.norm(gm)),
                           "g": gm.tolist()}
    return write_request(A, "line")


def frame_score(Rd, name, evo):
    """(T s, block SE s, record) of one scored case, or (None, None, why)."""
    D = Path(Rd) / name
    res, rec = load(D / "result.json"), load(D / "rec.json")
    if res is None:
        return None, None, {"why": "no result"}
    if res.get("error") or res.get("drag_mN") is None:
        return None, None, {"why": f"score failed: {str(res.get('error'))[:200]}"}
    if rec is None or not rec.get("legal") or rec.get("failed"):
        return None, None, {"why": "illegal or failed build"}
    D20 = 2.0 * float(res["drag_mN"]) / 1e3
    L20 = 2.0 * float(res.get("lift_mN") or 0.0) / 1e3
    if (Path(evo) / "repo_paths.py").exists() or race_eval is not _race_eval_evo:
        T, dTdD = race_eval(rec, D20, L20, evo)
    elif res.get("race_time_s") is not None and res.get("dT_dD20_s_per_N") is not None:
        T, dTdD = float(res["race_time_s"]), float(res["dT_dD20_s_per_N"])   # computed by the case job
    else:
        raise RuntimeError(f"no evolution checkout at {evo} (--evo / CARHYBRID_EVO) for the race model, and "
                           f"{name}/result.json carries no race_time_s / dT_dD20_s_per_N")
    se = abs(dTdD) * 2.0 * float(res.get("drag_block_se_mN") or 0.0) / 1e3
    return T, se, {"drag_mN": res["drag_mN"], "lift_mN": res.get("lift_mN"), "se_mN": res.get("drag_block_se_mN"),
                   "mass_g": rec.get("mass", {}).get("race_mass_g"), "T_s": T, "se_s": se}


def ingest_line(A, R, evo):
    batch = A["pending"]
    rows = {}
    for tag in A["candidates"]:
        fr = {k: frame_score(R, f"{batch}_{tag}_k{k}", evo) for k in FRAMES}
        ok = all(v[0] is not None for v in fr.values())
        rows[tag] = {"ok": ok, "frames": {str(k): v[2] for k, v in fr.items()},
                     "T": [fr[k][0] for k in FRAMES] if ok else None, "se": [fr[k][1] for k in FRAMES] if ok else None}
    log = {"stage": "line", "i": A["i"], "batch": batch, "radius": A["radius"], "u": A["u"], "s": A.get("s"),
           "direction_info": {k: v for k, v in (A.get("direction_info") or {}).items() if k != "g"},
           "candidates": {}, "accepted": None}
    g = (A.get("direction_info") or {}).get("g")
    if not rows["base"]["ok"]:
        if A.get("retry", 0) < 1:                     # the base is the reference: without it nothing can be judged
            log["note"] = "base re-run missing: the line stage is issued again"
            A["history"].append(log)
            return write_request(A, "line", retry=A.get("retry", 0) + 1)
        log["note"] = "base re-run missing twice: counted as nothing accepted"
    else:
        Tb, sb = np.asarray(rows["base"]["T"]), np.asarray(rows["base"]["se"])
        best = None
        for tag, a in STEPS:
            r = rows[tag]
            c = {"a": a, "ok": r["ok"], "frames": r["frames"], "u": A["candidates"][tag]}
            if g is not None:
                c["predicted_s"] = float(np.dot(g, np.asarray(A["candidates"][tag]) - np.asarray(xvec(A))))
            if r["ok"]:
                Dk = np.asarray(r["T"]) - Tb
                se_b = float(np.sqrt(np.sum(np.asarray(r["se"]) ** 2 + sb ** 2)) / len(Dk))
                se_f = float(np.std(Dk, ddof=1) / np.sqrt(len(Dk))) if len(Dk) > 1 else 0.0
                c.update(delta_k_s=Dk.tolist(), delta_s=float(Dk.mean()), se_block_s=se_b, se_frames_s=se_f,
                         se_s=max(se_b, se_f))
                if best is None or c["delta_s"] < log["candidates"][best]["delta_s"]:
                    best = tag
            log["candidates"][tag] = c
        log["base"] = rows["base"]["frames"]
        if best is not None:
            c = log["candidates"][best]
            log["best"], log["best_delta_s"], log["best_se_s"] = best, c["delta_s"], c["se_s"]
            if c["delta_s"] <= -N_SE * c["se_s"] and c["delta_s"] <= -FLOOR_S:
                log["accepted"] = best
    r0 = A["radius"]
    if log["accepted"]:
        A["u"], s_new = split_x(A, A["candidates"][log["accepted"]])
        if nw(A):
            A["s"] = s_new
        A["grads_at_u"] = []
        if log["accepted"] == STEPS[-1][0]:
            A["radius"] = min(r0 * GROW, R_MAX)
        A.setdefault("accepted", []).append({"i": A["i"], "tag": log["accepted"], "delta_s": log["best_delta_s"],
                                             "se_s": log["best_se_s"]})
    else:
        A["radius"] = r0 * SHRINK
    log["radius_after"] = A["radius"]
    A["history"].append(log)
    A["i"] += 1
    if A["radius"] < R_MIN:
        return stop(A, f"trust radius {A['radius']:.4g} < {R_MIN}")
    if A["i"] >= A["max_iter"]:
        return stop(A, f"{A['max_iter']} iterations")
    return next_stage(A)


# ------------------------------------------------------------------ commands
def init(layout, grid, arm="adjoint", radius=RADIUS0, seed=SEED, max_iter=MAX_ITER, bg_div=BG_DIV, runner="bigmem",
         np_=None, force=False, wing=None, evo="evo"):
    assert arm in ("adjoint", "random"), arm
    S = load(state_path(), {"arms": {}})
    if arm in S["arms"] and not force:
        print(f"arm {arm} exists (status {S['arms'][arm]['status']}); --force to start it again", file=sys.stderr)
        return EXIT_REFUSED
    n = grid_n(grid)
    A = {"arm": arm, "layout": layout, "grid": str(grid), "n": n, "u": [0.0] * n, "radius": float(radius),
         "radius0": float(radius), "i": 0, "seed": int(seed), "max_iter": int(max_iter), "bg_div": float(bg_div),
         "runner": runner, "np": np_, "status": "active", "stage": None, "pending": None, "pending_cases": [],
         "retry": 0, "grads_at_u": [], "history": [], "runs": [], "accepted": []}
    if wing:
        names = list(DEFAULT_WING) if wing == "default" else [w for w in wing.split(",") if w]
        A["wing_names"] = names
        A["wing_ranges"], A["s"] = wing_start(evo, layout, names)
    batch, cases = next_stage(A)
    S["arms"][arm] = A
    save(state_path(), S)
    print(f"{arm}: {A['stage']} -> {batch} ({len(cases)} cases)")
    return EXIT_NEXT if batch else EXIT_STOP


def find_arm(S, results, arm=None):
    if arm:
        return arm if arm in S["arms"] else None
    b = Path(results).name
    hit = [a for a, A in S["arms"].items() if A.get("pending") == b]
    if hit:
        return hit[0]
    hit = [a for a, A in S["arms"].items() if A["status"] != "active" and A["runs"] and A["runs"][-1]["batch"] == b]
    if hit:                                                  # the last batch of an arm that has stopped since
        return hit[0]
    hit = [a for a, A in S["arms"].items() if A.get("pending_cases") and
           any((Path(results) / c).is_dir() for c in A["pending_cases"])]
    return hit[0] if len(hit) == 1 else None


def advance(results, run_id, evo="evo", arm=None):
    S = load(state_path())
    if not S:
        print("no loop state (run init first)", file=sys.stderr)
        return EXIT_REFUSED
    a = find_arm(S, results, arm)
    if a is None:
        print(f"no arm is waiting on {results}", file=sys.stderr)
        return EXIT_REFUSED
    A = S["arms"][a]
    if A["status"] != "active" or not A.get("pending"):
        print(f"arm {a} is {A['status']}", file=sys.stderr)
        return EXIT_STOP
    A["runs"].append({"batch": A["pending"], "run": str(run_id)})
    try:
        if A["stage"] == "grad":
            batch, cases = ingest_grad(A, results)
        else:
            batch, cases = ingest_line(A, results, evo)
    except Exception as e:                                   # noqa: BLE001  nothing is saved: the stage can be advanced again
        print(f"{a}: advance failed, state unchanged: {e!r}"[:2000], file=sys.stderr)
        return EXIT_REFUSED
    save(state_path(), S)
    if batch:
        print(f"{a}: iteration {A['i']} {A['stage']} -> {batch} ({len(cases)} cases); radius {A['radius']:.4g}")
        return EXIT_NEXT
    print(f"{a}: stopped ({A['stop_reason']}); accepted {len(A['accepted'])}")
    return EXIT_STOP


def reissue(arm):
    S = load(state_path())
    A = (S or {}).get("arms", {}).get(arm)
    if not A or A["status"] != "active" or not A.get("pending"):
        print(f"arm {arm} has no pending stage", file=sys.stderr)
        return EXIT_REFUSED
    batch, cases = write_request(A, A["stage"], retry=A.get("retry", 0) + 1)
    save(state_path(), S)
    print(f"{arm}: {A['stage']} reissued -> {batch} ({len(cases)} cases)")
    return EXIT_NEXT


def status():
    S = load(state_path(), {"arms": {}})
    for a, A in S["arms"].items():
        print(json.dumps({"arm": a, "status": A["status"], "stage": A["stage"], "i": A["i"], "radius": A["radius"],
                          "pending": A["pending"], "accepted": len(A["accepted"]), "stop_reason": A.get("stop_reason")}))
    return 0


def main(argv):
    opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d
    cmd = argv[0] if argv else ""
    if cmd == "init":
        return init(argv[1], argv[2], arm=opt("--arm", "adjoint"), radius=float(opt("--radius", RADIUS0)),
                    seed=int(opt("--seed", SEED)), max_iter=int(opt("--max-iter", MAX_ITER)),
                    bg_div=float(opt("--bg-div", BG_DIV)), runner=opt("--runner", "bigmem"),
                    np_=int(opt("--np")) if opt("--np") else None, force="--force" in argv,
                    wing=opt("--wing"), evo=opt("--evo", os.environ.get("CARHYBRID_EVO", "evo")))
    if cmd == "advance":
        return advance(argv[1], argv[2], evo=opt("--evo", os.environ.get("CARHYBRID_EVO", "evo")), arm=opt("--arm"))
    if cmd == "reissue":
        return reissue(argv[1])
    if cmd == "status":
        return status()
    print(__doc__, file=sys.stderr)
    return EXIT_REFUSED


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
