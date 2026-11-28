"""carhybrid/grad.py on P0's DAFoam gradient (carhybrid/results/p0/p0_team_dafoam).

Needs the pinned evolution checkout with this branch's evo_patches applied: CARHYBRID_EVO=<checkout> (else
<repo>/evo); skipped without it.
    CARHYBRID_EVO=/path/to/evo python -m pytest carhybrid/tests/test_grad.py -q -s
"""
import json, os, sys, time
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
CH = HERE.parent
sys.path.insert(0, str(CH))
import grad as G  # noqa: E402

P0 = CH / "results" / "p0" / "p0_team_dafoam"
EVO = Path(os.environ.get("CARHYBRID_EVO", CH.parent / "evo"))
pytestmark = pytest.mark.skipif(not (EVO / "repo_paths.py").exists(), reason="no evolution checkout (CARHYBRID_EVO)")
FINE = "24,9,12,9"


def run(tmp, npz, name, **kw):
    out = tmp / f"{name}.json"
    t = time.time()
    r = G.main(str(EVO), str(P0 / "rec.json"), npz, "zero", str(out), quiet=True, **kw)
    r["_wall_s"] = time.time() - t
    return r


@pytest.fixture(scope="module")
def tmp(tmp_path_factory):
    return tmp_path_factory.mktemp("grad")


@pytest.fixture(scope="module")
def g45(tmp):
    return run(tmp, str(P0 / "dafoam_sens.npz"), "g45", d20=0.36)


@pytest.fixture(scope="module")
def g720(tmp):
    r = run(tmp, str(P0 / "dafoam_sens.npz"), "g720", d20=0.36, grid=FINE)
    print(f"\n720-number map: {r['run_s']:.1f} s ({r['n_points_on_body']} points on the body)")
    return r


def test_reproduces_p0_drag(g45):
    """The default grid reproduces P0's stored dD/du (the field moved slightly since: evo_patches' soft cap v2)."""
    ref = np.asarray(json.loads((P0 / "grad.json").read_text())["dD_du_N"])
    c = G.cos(g45["dD_du_N"], ref)
    print(f"\ncos(dD/du, P0's) = {c:.5f}; points on body {g45['n_points_on_body']} (P0: 2203)")
    assert c > 0.99
    assert g45["n"] == 45 and g45["has_lift"] is False and g45["dL_du_N"] is None
    assert abs(g45["dT_dD_s_per_N"] - 0.17832) < 1e-3          # P0's partial at D20 = 0.36 N


def test_lin_matches_field_finite_difference(tmp, g45):
    """method lin (chain rule through the surfaces) == the whole-field finite difference (method fd, one-sided)."""
    fd = run(tmp, str(P0 / "dafoam_sens.npz"), "g45fd", d20=0.36, method="fd")
    a, b = np.asarray(g45["dD_du_N"]), np.asarray(fd["dD_du_N"])
    print(f"\nlin vs fd: cos {G.cos(a, b):.6f}, |a-b|/|b| {np.linalg.norm(a - b) / np.linalg.norm(b):.4f}")
    assert G.cos(a, b) > 0.999


def test_lift_equal_to_drag(tmp, g45):
    """A synthetic dL_dXs equal to dD_dXs maps to dL/du == dD/du; the race-time gradient adds dT/dL * dL/du."""
    S = np.load(P0 / "dafoam_sens.npz")
    p = tmp / "lift_eq_drag.npz"
    np.savez(p, Xs=S["Xs"], dD_dXs=S["dD_dXs"], dL_dXs=S["dD_dXs"], D_N=S["D_N"], L_N=0.0974)
    r = run(tmp, str(p), "glift", d20=0.465, l20=0.195)
    dD, dL = np.asarray(r["dD_du_N"]), np.asarray(r["dL_du_N"])
    assert r["has_lift"] and np.allclose(dL, dD, rtol=0, atol=1e-15)
    assert np.allclose(dD, g45["dD_du_N"], rtol=1e-12, atol=1e-15)
    # lift + up: more upforce unloads the wheels -> faster (race_objective's lift_20_n, + up)
    assert r["dT_dL_s_per_N"] < 0 and r["dT_dD_s_per_N"] > 0
    assert r["L20_N_used"] == 0.195 and r["L20_source"] == "cli"
    exp = (r["dT_dD_s_per_N"] + r["dT_dL_s_per_N"]) * dD     # no --body: dm/du = 0
    assert np.allclose(r["dT_du_s"], exp, rtol=1e-12, atol=1e-18)
    print(f"\ndT/dD {r['dT_dD_s_per_N'] * 1e3:.1f} ms/N, dT/dL {r['dT_dL_s_per_N'] * 1e3:.2f} ms/N at D20 0.465, L20 0.195")


def test_frames_averaged(tmp, g45):
    """Several npz files: the mapped gradients are averaged and each frame's cosine is reported."""
    S = np.load(P0 / "dafoam_sens.npz")
    rng = np.random.default_rng(1)
    p = tmp / "frame_b.npz"
    np.savez(p, Xs=S["Xs"], dD_dXs=S["dD_dXs"] * (1 + 0.5 * rng.standard_normal((len(S["Xs"]), 1))), D_N=S["D_N"])
    r = run(tmp, f"{P0 / 'dafoam_sens.npz'},{p}", "g2frames", d20=0.36)
    assert r["n_frames"] == 2 and len(r["frames"]) == 2 and len(r["cos_dD_frames"]) == 2
    one = np.asarray(g45["dD_du_N"])
    fb = 2 * np.asarray(r["dD_du_N"]) - one                    # the second frame's own map
    assert G.cos(fb, one) < 1 - 1e-6
    assert abs(r["cos_dD_frames"][0][1] - G.cos(fb, one)) < 1e-9
    print("\nper-frame cos with the mean:", [round(f["cos_dD_vs_mean"], 4) for f in r["frames"]],
          "pairwise", round(r["cos_dD_frames"][0][1], 4))


def test_smoothing(tmp):
    fr, _, SS = G.evo_import(EVO)
    seed = fr.Freeform.from_dict(json.loads((P0 / "rec.json").read_text())["decoded"]["freeform"])
    grid = G.parse_grid(FINE, SS)
    sm = G.smoother(fr, seed, grid, 15.0)
    n = fr.n_genes(*grid)
    assert np.allclose(sm(np.ones(n)), 1.0)                    # a constant is kept
    x = np.random.default_rng(0).standard_normal(n)
    assert np.linalg.norm(sm(x)) < 0.6 * np.linalg.norm(x)     # point-to-point noise is damped


def test_fine_grid_vs_coarse(g45, g720):
    """The 720-number gradient restricted to the coarse space (P^T g, P = the coarse numbers' shapes in fine numbers)
    is the 45-number gradient (the chain rule through freeform.from_unit); what share of it the coarse space holds is
    reported, not asserted."""
    fr, _, SS = G.evo_import(EVO)
    seed = fr.Freeform.from_dict(json.loads((P0 / "rec.json").read_text())["decoded"]["freeform"])
    P = G.prolongation(fr, seed, tuple(SS.COARSE), G.parse_grid(FINE, SS), SS.AMP_MM)
    gf, gc = np.asarray(g720["dD_du_N"]), np.asarray(g45["dD_du_N"])
    assert g720["n"] == 720 and np.all(np.isfinite(gf))
    c = G.cos(P.T @ gf, gc)
    Q, _ = np.linalg.qr(P)
    share = np.linalg.norm(Q.T @ gf) / np.linalg.norm(gf)
    print(f"\n720 vs 45: cos(P^T g720, g45) {c:.6f}, |P^T g720| / |g45| {np.linalg.norm(P.T @ gf) / np.linalg.norm(gc):.4f};"
          f" share of |g720| in the coarse space {share:.3f}; cos(P g45, g720) {G.cos(P @ gc, gf):.3f};"
          f" 720 run {g720['run_s']:.1f} s")
    assert c > 0.99
