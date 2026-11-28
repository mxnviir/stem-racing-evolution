"""carhybrid/loop.py: dry runs of the P2F stage machine on synthetic results directories (fake result.json, rec.json,
grad.json), checking state.json and request.json at every step. The race model is replaced by a linear stand-in
(dT/dD20 = 0.185 s/N, as the team car) except in test_real_race_time (needs CARHYBRID_EVO, else skipped).
"""
import json, os, sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
CH = HERE.parent
sys.path.insert(0, str(CH))
import loop as LP  # noqa: E402

DTDD = 0.185
GRID, N = "5,2,3,4", 45


def fake_race(rec, D20, L20, evo):
    return 0.80 + DTDD * D20 - 0.002 * L20 + 0.0075 * (rec["mass"]["race_mass_g"] - 75.0), DTDD


@pytest.fixture
def env(tmp_path, monkeypatch):
    monkeypatch.setattr(LP, "LOOP", tmp_path / "loop")
    monkeypatch.setattr(LP, "REQ", tmp_path / "request.json")
    monkeypatch.setattr(LP, "race_eval", fake_race)
    (tmp_path / "request.json").write_text(json.dumps({"note": "keep me", "launch": "p1",
                                                       "batches": {"p1": [{"name": "x", "type": "score"}]}}))
    return tmp_path


def req(env):
    return json.loads((env / "request.json").read_text())


def state(env):
    return json.loads((env / "loop" / "state.json").read_text())


def arm(env, a="adjoint"):
    return state(env)["arms"][a]


def results_dir(env, batch):
    d = env / "results" / batch
    d.mkdir(parents=True, exist_ok=True)
    return d


def write_grad(env, A, g, dafoam=None, steady=None, summary=None):
    batch = A["pending"]
    R = results_dir(env, batch)
    c = R / A["pending_cases"][0]
    c.mkdir()
    if g is not None:
        (c / "grad.json").write_text(json.dumps({"dT_du_s": list(map(float, g)), "D20_N_used": 0.466, "lead_g": 0.0}))
    (c / "result.json").write_text(json.dumps({"drag_mN": 233.0, "drag_block_se_mN": 0.2, "lift_mN": 95.0}))
    if dafoam is not None:
        (c / "dafoam_sens.json").write_text(json.dumps(dafoam))
    if steady is not None:
        (c / "steady.json").write_text(json.dumps(steady))
    if summary is not None:
        (c / "fine_adjoint.json").write_text(json.dumps(summary))
    return R


def write_line(env, A, drag, legal=None, missing=()):
    """drag(tag, k) -> half-car drag_mN; legal(tag) -> bool; missing: case names with no folder (failed jobs)."""
    batch = A["pending"]
    R = results_dir(env, batch)
    for tag in A["candidates"]:
        for k in (0, 1):
            name = f"{batch}_{tag}_k{k}"
            if name in missing:
                continue
            c = R / name
            c.mkdir()
            ok = legal(tag) if legal else True
            (c / "rec.json").write_text(json.dumps({"legal": ok, "failed": [] if ok else ["T5.6"],
                                                    "mass": {"race_mass_g": 75.0}}))
            if ok:
                (c / "result.json").write_text(json.dumps({"drag_mN": drag(tag, k), "drag_block_se_mN": 0.2,
                                                           "lift_mN": 95.0}))
    return R


def advance(env, R, run="1"):
    return LP.main(["advance", str(R), run])


def check_request(env, A, stage, n_cases):
    r = req(env)
    assert r["note"] == "keep me" and "p1" in r["batches"]          # the rest of request.json is kept
    assert r["launch"] == A["pending"] == f"loop_{A['arm']}_{A['i']}_{stage}" + (f"_r{A['retry']}" if A["retry"] else "")
    cases = r["batches"][r["launch"]]
    assert len(cases) == n_cases and [c["name"] for c in cases] == A["pending_cases"]
    return cases


def line_u(cases, tag):
    c = next(c for c in cases if c["name"].endswith(f"_{tag}_k0"))
    return np.zeros(N) if c["u"] == "zero" else np.asarray(json.loads(c["u"]))


def test_adjoint_dry_run(env):
    rng = np.random.default_rng(3)
    g0 = rng.standard_normal(N) * 1e-3
    # ---- init -> grad 0
    assert LP.main(["init", "team", GRID]) == 0
    A = arm(env)
    c = check_request(env, A, "grad", 1)[0]
    assert c == {"name": "loop_adjoint_0_grad_g", "type": "fine_adjoint", "layout": "team", "u": "zero", "k": 0,
                 "bg_div": 12.889, "grid": GRID, "runner": "bigmem"}
    assert A["radius"] == 0.25 and A["n"] == N and A["stage"] == "grad"
    # ---- grad 0 -> line 0
    assert advance(env, write_grad(env, A, g0)) == 0
    A = arm(env)
    cases = check_request(env, A, "line", 8)
    assert {c["type"] for c in cases} == {"score"} and {c["k"] for c in cases} == {0, 1}
    assert all(c["grid"] == GRID and c["layout"] == "team" for c in cases)
    d = -g0 / np.linalg.norm(g0)
    assert np.allclose(A["direction"], d)
    assert line_u(cases, "base").tolist() == [0.0] * N
    for tag, a in (("a05", 0.5), ("a1", 1.0), ("a2", 2.0)):
        assert np.allclose(line_u(cases, tag), np.clip(a * 0.25 * d, -1, 1), atol=1e-6)
    # ---- line 0: a1 wins clearly (-2 mN half car on both frames) -> accepted, radius kept
    drag = {"base": 233.0, "a05": 232.5, "a1": 231.0, "a2": 234.0}
    assert advance(env, write_line(env, A, lambda t, k: drag[t] - 0.1 * k * (t == "a1"))) == 0
    A = arm(env)
    h = A["history"][-1]
    assert h["accepted"] == "a1" and A["radius"] == 0.25 and A["i"] == 1
    assert h["candidates"]["a1"]["delta_s"] == pytest.approx(DTDD * 2 * (-2.05) / 1e3)
    assert h["candidates"]["a1"]["se_s"] == pytest.approx(max(DTDD * 2 * 0.2e-3, DTDD * 2 * 0.1e-3 / 2))
    u1 = np.asarray(A["u"])
    assert np.allclose(u1, np.clip(0.25 * d, -1, 1), atol=1e-6)
    c = check_request(env, A, "grad", 1)[0]
    assert c["k"] == 1 and np.allclose(json.loads(c["u"]), u1)
    # ---- grad 1 -> line 1: nothing better -> rejected, radius halves, u kept
    g1 = rng.standard_normal(N) * 1e-3
    assert advance(env, write_grad(env, A, g1)) == 0
    A = arm(env)
    check_request(env, A, "line", 8)
    assert np.allclose(line_u(req(env)["batches"][A["pending"]], "base"), u1)
    assert advance(env, write_line(env, A, lambda t, k: 233.0 + 0.5 * (t != "base"))) == 0
    A = arm(env)
    assert A["history"][-1]["accepted"] is None and A["radius"] == 0.125 and np.allclose(A["u"], u1)
    c = check_request(env, A, "grad", 1)[0]
    assert c["k"] == 2
    # ---- grad 2 at the same u: the line uses the mean of the two gradients taken there
    g2 = rng.standard_normal(N) * 1e-3
    assert advance(env, write_grad(env, A, g2)) == 0
    A = arm(env)
    gm = (g1 + g2) / 2
    assert A["direction_info"]["n_averaged"] == 2 and np.allclose(A["direction"], -gm / np.linalg.norm(gm))
    # ---- line 2: the largest step wins -> accepted, radius x 1.5
    drag = {"base": 233.0, "a05": 232.6, "a1": 232.0, "a2": 230.0}
    assert advance(env, write_line(env, A, lambda t, k: drag[t])) == 0
    A = arm(env)
    assert A["history"][-1]["accepted"] == "a2" and A["radius"] == pytest.approx(0.1875) and A["grads_at_u"] == []
    # ---- grad 3 -> line 3: a1, a2 illegal; a05 better but within 3 SE -> rejected
    assert advance(env, write_grad(env, A, rng.standard_normal(N) * 1e-3)) == 0
    A = arm(env)
    assert advance(env, write_line(env, A, lambda t, k: 233.0 - 0.3 * (t == "a05"), legal=lambda t: t in ("base", "a05"))) == 0
    A = arm(env)
    h = A["history"][-1]
    assert h["candidates"]["a1"]["ok"] is False and h["candidates"]["a2"]["ok"] is False
    assert h["best"] == "a05" and h["accepted"] is None and A["radius"] == pytest.approx(0.09375)
    # ---- grad 4 -> line 4: the base's frame 1 job failed -> the same line stage again (r1)
    assert advance(env, write_grad(env, A, rng.standard_normal(N) * 1e-3)) == 0
    A = arm(env)
    first = req(env)["batches"][A["pending"]]
    assert advance(env, write_line(env, A, lambda t, k: 233.0, missing=(f"{A['pending']}_base_k1",))) == 0
    A = arm(env)
    again = check_request(env, A, "line", 8)
    assert A["pending"] == "loop_adjoint_4_line_r1" and A["i"] == 4 and A["radius"] == pytest.approx(0.09375)
    assert [c["u"] for c in again] == [c["u"] for c in first]
    assert advance(env, write_line(env, A, lambda t, k: 233.4 if t != "base" else 233.0)) == 0
    A = arm(env)
    assert A["radius"] == pytest.approx(0.046875) and A["i"] == 5
    # ---- rejections until the radius falls under 0.02 -> stop (exit 3)
    codes = []
    while arm(env)["status"] == "active":
        A = arm(env)
        R = write_grad(env, A, rng.standard_normal(N) * 1e-3) if A["stage"] == "grad" else \
            write_line(env, A, lambda t, k: 233.0 + 0.2 * (t != "base"))
        codes.append(advance(env, R))
    A = arm(env)
    assert codes == [0, 0, 0, 3] and A["i"] == 7 and A["radius"] < 0.02 and "trust radius" in A["stop_reason"]
    assert len(A["accepted"]) == 2 and len(A["runs"]) == 15
    assert req(env)["launch"] == "loop_adjoint_6_line"                # a stop leaves the last batch in place
    assert LP.main(["advance", str(env / "results" / "loop_adjoint_6_line"), "9"]) == 3


def test_random_arm(env):
    assert LP.main(["init", "team", GRID, "--arm", "random", "--max-iter", "2"]) == 0
    A = arm(env, "random")
    cases = check_request(env, A, "line", 8)                          # no grad stage
    d0 = np.asarray(A["direction"])
    assert A["direction_info"] == {"kind": "random", "rng": [20261128, 0]}
    rr = np.random.default_rng([20261128, 0]).standard_normal(N)
    assert np.allclose(d0, rr / np.linalg.norm(rr)) and np.allclose(line_u(cases, "a1"), 0.25 * d0, atol=1e-6)
    # an adjoint arm alongside: advance finds the arm by the results folder's batch name
    assert LP.main(["init", "team", GRID]) == 0
    assert advance(env, write_line(env, A, lambda t, k: 233.0 - 2.0 * (t == "a05"))) == 0
    S = state(env)
    assert S["arms"]["adjoint"]["i"] == 0 and S["arms"]["adjoint"]["stage"] == "grad"
    A = S["arms"]["random"]
    assert A["history"][-1]["accepted"] == "a05" and A["stage"] == "line" and A["i"] == 1
    assert not np.allclose(A["direction"], d0) and A["direction_info"]["rng"] == [20261128, 1]
    assert req(env)["launch"] == "loop_random_1_line"
    assert advance(env, write_line(env, A, lambda t, k: 233.0)) == 3   # max-iter 2
    assert "2 iterations" in arm(env, "random")["stop_reason"]
    # a stale or unknown results folder is refused
    assert LP.main(["advance", str(env / "results" / "nothing"), "1"]) == 2


@pytest.mark.parametrize("dafoam,steady,summary,g,why", [
    ({"gmres": {"D": {"converged": True}, "L": {"converged": False}}}, None, None, "ok", "F3"),
    ({"converged": {"D": True, "L": False}}, None, None, "ok", "F3"),
    (None, {"pass": False}, None, "ok", "F2"),
    (None, None, {"F1_pass": False, "F1_dafoam_vs_pipeline_drag": 0.05, "F2_pass": True, "F3_pass": True}, "ok", "F1"),
    (None, None, {"F1_pass": None, "F2_pass": True, "F3_pass": False}, "ok", "F3"),
    (None, None, None, None, "no grad.json"),
    (None, None, None, "short", "unusable"),
])
def test_grad_gates_stop(env, dafoam, steady, summary, g, why):
    assert LP.main(["init", "team", GRID]) == 0
    A = arm(env)
    vec = {"ok": np.ones(N) * 1e-3, "short": np.ones(3), None: None}[g]
    assert advance(env, write_grad(env, A, vec, dafoam=dafoam, steady=steady, summary=summary)) == 3
    A = arm(env)
    assert A["status"] == "stopped" and why in A["stop_reason"] and req(env)["launch"] == "loop_adjoint_0_grad"


def test_grad_gates_pass(env):
    """Gates passed or unknown (None): the line stage follows."""
    assert LP.main(["init", "team", GRID]) == 0
    A = arm(env)
    summ = {"F1_pass": True, "F2_pass": True, "F3_pass": None, "F3_gmres_converged": {}}
    dafoam = {"gmres": {"D": {"converged": True}, "L": {"converged": True}}}
    assert advance(env, write_grad(env, A, np.ones(N) * 1e-3, dafoam=dafoam, steady={"pass": True}, summary=summ)) == 0
    assert arm(env)["stage"] == "line" and arm(env)["history"][-1]["gates"]["F1_pass"] is True


def test_box_projection(env):
    u = np.array([1.0, -1.0, 0.3, 1.0])
    d = LP.project(u, np.array([1.0, -1.0, 1.0, -1.0]))
    assert np.allclose(d, np.array([0, 0, 1, -1]) / np.sqrt(2))
    assert LP.project(u[:2], np.array([1.0, -1.0])) is None
    assert LP.grid_n("24,9,12,9") == 720 and LP.grid_n("5,2,3,4") == 45 and LP.grid_n("5,2,3") == 35


def test_real_race_time(tmp_path):
    """frame_score with the evolution race model on P0's frame-0 score: T = race_time(rec, 2 D, 2 L) (whole car)."""
    evo = Path(os.environ.get("CARHYBRID_EVO", CH.parent / "evo"))
    if not (evo / "repo_paths.py").exists():
        pytest.skip("no evolution checkout (CARHYBRID_EVO)")
    src = CH / "results" / "p0" / "p0_team_score_k0"
    T, se, info = LP.frame_score(src.parent, src.name, evo)
    import grad as G
    SL = G.evo_import(evo)[1]
    res, rec = (json.loads((src / f).read_text()) for f in ("result.json", "rec.json"))
    assert T == pytest.approx(SL.race_time(rec, 2 * res["drag_mN"] / 1e3, 2 * res["lift_mN"] / 1e3), abs=1e-12)
    assert se == pytest.approx(0.185 * 2 * res["drag_block_se_mN"] / 1e3, rel=0.02)
    print(f"\nP0 frame 0: T {T:.5f} s, block SE {se * 1e3:.4f} ms (drag {res['drag_mN']:.1f}, lift {res['lift_mN']:.1f} mN half car)")


def test_race_time_without_checkout(env, monkeypatch):
    """No evolution checkout: the case job's race_time_s / dT_dD20_s_per_N are used; with neither, advance refuses
    (exit 2) and leaves the state as it was, so the stage can be advanced again once the checkout is there."""
    monkeypatch.setattr(LP, "race_eval", LP._race_eval_evo)
    monkeypatch.setenv("CARHYBRID_EVO", str(env / "no_checkout"))
    assert LP.main(["init", "team", GRID, "--arm", "random"]) == 0
    A = arm(env, "random")
    R = write_line(env, A, lambda t, k: 233.0)
    before = (env / "loop" / "state.json").read_text()
    assert advance(env, R) == 2 and (env / "loop" / "state.json").read_text() == before
    for c in R.iterdir():
        r = json.loads((c / "result.json").read_text())
        r.update(race_time_s=0.81 - 0.001 * c.name.endswith(("a1_k0", "a1_k1")), dT_dD20_s_per_N=DTDD)
        (c / "result.json").write_text(json.dumps(r))
    assert advance(env, R) == 0
    h = arm(env, "random")["history"][-1]
    assert h["accepted"] == "a1" and h["candidates"]["a1"]["delta_s"] == pytest.approx(-0.001)
