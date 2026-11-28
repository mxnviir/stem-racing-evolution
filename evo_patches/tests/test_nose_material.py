"""Tests of evo_patches/nose_material.py (the nose's foam / printed split and the flush junction).

    EVO_CHECKOUT=<evolution checkout> [CAR_DIR=<a car built with the patch: make_car.py's <out>/evo>] \
        python -m pytest evo_patches/tests/test_nose_material.py -q

EVO_CHECKOUT is copied (only the files the patch touches are written, in a temporary copy), so it is never changed;
evo_patches/soft_cap.py is applied to the copy first, as in every run. The plan tests build the team's front wing
(about 10-20 s each). CAR_DIR: the checks on a built car (legality, masses, the junction, the files).
"""
from __future__ import annotations

import importlib
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
PATCHES = HERE.parent
sys.path.insert(0, str(PATCHES))
import nose_material as NM  # noqa: E402

EVO = os.environ.get("EVO_CHECKOUT")
CAR = os.environ.get("CAR_DIR")
need_evo = pytest.mark.skipif(not EVO, reason="set EVO_CHECKOUT to an evolution checkout")
need_car = pytest.mark.skipif(not CAR, reason="set CAR_DIR to a car built with the patch (<out>/evo)")

SHORT = [0.0, 0.0, 0.5, 0.1304, 0.5, 1.0, 0.6562, 0.5172, 0.0, 0.0, 0.175, 0.0]   # W 120, rear 36, wing 20 mm ahead


@pytest.fixture(scope="session")
def checkout(tmp_path_factory):
    """A patched copy of EVO_CHECKOUT (soft_cap, then nose_material), importable."""
    src = Path(EVO)
    dst = tmp_path_factory.mktemp("evo") / "evo"
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__", "results", ".git"))
    import soft_cap
    soft_cap.apply(dst)
    NM.apply(dst)
    sys.path.insert(0, str(dst))
    import repo_paths
    importlib.reload(repo_paths)
    repo_paths.add()
    sys.path.insert(0, str(dst / "car" / "parts"))
    for m in [m for m in list(sys.modules) if m in ("freeform", "nose_split", "layout_genes", "front_wing_engine")]:
        del sys.modules[m]
    return dst


@pytest.fixture(scope="session")
def plans(checkout):
    import layout_genes as LG
    import freeform as fr
    import nose_split as NS
    out = {}
    for name, x in (("team", LG.x_default()), ("short", np.array(SHORT))):
        d = LG.decode(x)
        ff = fr.Freeform.from_dict(d["freeform"])
        out[name] = (d, ff, NS.plan(ff, d["x_front"], d["W"], d["rear"], d["fwing"], d["nose"], d["skin_mm"]))
    return out


# ------------------------------------------------------------------------------------------- the patch itself
@need_evo
def test_apply_is_idempotent_and_exact(tmp_path):
    """Every anchor is found once; a second apply changes nothing; the module is written once."""
    dst = tmp_path / "evo"
    files = sorted({f for f, _o, _n in NM.PATCHES} | {"car/body/freeform.py"})
    for f in files:
        (dst / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(Path(EVO) / f, dst / f)
    import soft_cap
    soft_cap.apply(dst)
    first = NM.apply(dst)
    assert {s for _f, s in first} <= {"applied", "already applied", "written", "already written"}
    snap = {f: (dst / f).read_text() for f in files + [NM.MODULE]}
    second = NM.apply(dst)
    assert all(s.startswith("already") for _f, s in second), second
    assert snap == {f: (dst / f).read_text() for f in files + [NM.MODULE]}
    for f in files + [NM.MODULE]:
        compile((dst / f).read_text(), f, "exec")


@need_evo
def test_apply_refuses_a_changed_builder(tmp_path):
    dst = tmp_path / "evo"
    for f in {f for f, _o, _n in NM.PATCHES}:
        (dst / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(Path(EVO) / f, dst / f)
    p = dst / "car/parts/assembly.py"
    s = p.read_text()
    if "evo_patches/nose_material.py" in s:
        pytest.skip("EVO_CHECKOUT is already patched")
    p.write_text(s.replace('        res = fwe.build(front, nm["nose"]', '        res = fwe.build(front, nm["NOSE"]'))
    with pytest.raises(AssertionError):
        NM.apply(dst)


# ------------------------------------------------------------------------------------------- the decision
@need_evo
def test_team_plan(plans):
    """Team layout: the block's front is 5.2 mm ahead of Ref A and the wing's struts land up to 35.7 mm: the
    print stops just behind them (plus the collar), the foam carries on to Ref A, the glued face holds the loads."""
    import nose_split as NS
    d, ff, p = plans["team"]
    ref_a = d["x_front"] - 16.0
    assert p["ref_a_mm"] == pytest.approx(ref_a)
    assert p["block_front_mm"] == pytest.approx(d["x_front"] + d["W"] + d["rear"] - 223.0)
    assert p["x_cut_mm"] >= max(p["limits_mm"].values()) - 1e-9 or p["x_cut_mm"] == pytest.approx(ref_a)
    assert p["x_cut_mm"] >= p["block_front_mm"] + NS.BLOCK_MARGIN_MM - 1e-9
    assert p["x_cut_mm"] >= p["attachment"]["x_mm"] + NS.COLLAR_MM - 1e-9
    assert p["x_cut_mm"] <= ref_a + 1e-9
    assert p["decided_by"] in p["limits_mm"]
    assert 0.0 < p["foam_len_mm"] < 5.2                 # inside the block's 5.2 mm ahead of Ref A
    assert p["bond"]["util_mpa"] <= NS.FOAM_BOND_MPA
    assert p["tenon"] and min(p["tenon_walls_mm"]["top_mm"], p["tenon_walls_mm"]["side_mm"]) >= NS.FOAM_WALL_MM - 1e-9
    assert all(v >= -1e-9 for v in NS.gates(p).values()), NS.gates(p)


@need_evo
def test_short_plan_is_mostly_foam(plans):
    """A short layout (the block reaches 51 mm ahead of Ref A) with its wing 20 mm ahead: foam from just behind
    the struts' feet to Ref A, more than 10 mm of nose."""
    import nose_split as NS
    d, ff, p = plans["short"]
    assert p["block_front_mm"] < p["ref_a_mm"] - NS.NOSE_OVERHANG_MM      # the block is not what limits it
    assert p["decided_by"] in ("attachment", "bond", "thin")
    assert p["foam_len_mm"] > 10.0
    assert p["foam_nose_drawn_cm3"] > 2.0
    assert all(v >= -1e-9 for v in NS.gates(p).values())


@need_evo
def test_plan_is_deterministic(plans):
    import nose_split as NS
    d, ff, p = plans["team"]
    q = NS.plan(ff, d["x_front"], d["W"], d["rear"], d["fwing"], d["nose"], d["skin_mm"])
    assert json.dumps(p, sort_keys=True, default=float) == json.dumps(q, sort_keys=True, default=float)


@need_evo
def test_block_limits_a_long_car(plans):
    """The same nose and wing on a layout whose block ends at the wing's struts: no foam ahead of the block."""
    import nose_split as NS
    d, ff, _p = plans["team"]
    rear = 223.0 - 16.0 - d["W"] - 0.5           # the block's front 0.5 mm ahead of Ref A
    p = NS.plan(ff, d["x_front"], d["W"], rear, d["fwing"], d["nose"], d["skin_mm"])
    assert p["x_cut_mm"] == pytest.approx(p["ref_a_mm"])        # 0.5 mm of block, 1.0 mm margin: all printed
    assert p["foam_len_mm"] == pytest.approx(0.0)
    assert p["decided_by"] in ("block", "attachment")
    assert "nose_foam_inside_block" not in NS.gates(p)


@need_evo
def test_joint_stress_scales_with_lever(plans):
    import nose_split as NS
    sec = dict(A=300.0, zc=14.0, Iy=5000.0, Iz=9000.0, cz=8.0, cy=10.0)
    a = NS.joint_stress(sec, 30.0, (25.0, 8.0))
    b = NS.joint_stress(sec, 30.0, (15.0, 8.0))
    assert b["sigma_mpa"] > a["sigma_mpa"] > 0
    assert a["load_n"] == [2.0, 1.0, -4.0]


# ------------------------------------------------------------------------------------------- geometry and mass
class _Region:
    def __init__(self, origin_m, shape):
        self.origin_m, self.shape = origin_m, shape


class _Geom:
    def __init__(self, origin_m, shape):
        self.region = _Region(origin_m, shape)


@need_evo
def test_foam_weight(checkout):
    import nose_split as NS
    p = {"x_cut_mm": 10.0, "tenon": {"x0": 10.0, "x1": 18.0, "hw": 4.0, "z0": 2.0, "z1": 12.0}}
    g = _Geom((0.0, -0.02, 0.0), (41, 41, 41))          # 1 mm cells, x 0..40, y -20..20, z 0..40
    w = NS.foam_weight(g, p, 1.0)
    assert np.all(w[:10] == 0.0) and np.all(w[11:, :, 30] == 1.0)
    assert w[10, 20, 30] == pytest.approx(0.5)          # the cell on the face: half each side
    assert np.all(w[12:17, 17:24, 4:11] == 0.0)         # the tenon's pocket
    assert np.all(w[25:, :, :] == 1.0)


@need_evo
def test_printed_piece_meets_the_cut_body_flush(plans):
    """The print and a body made from the same field, cut at x_cut: the same section either side of the face
    (within the marching-cubes grid), one closed solid each, the tenon in the pocket with nothing overlapping."""
    import manifold3d as mf
    import trimesh
    from skimage import measure
    import nose_split as NS
    d, ff, p = plans["short"]
    xs, ys, zs = NS._grid(p["x_field_front_mm"], p["ref_a_mm"] + 12.0)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    phi = np.maximum(NS.drawn_phi(ff, p["ref_a_mm"], X, Y, Z), np.maximum(p["x_field_front_mm"] - X, X - (p["ref_a_mm"] + 12.0)))
    v, f, _, _ = measure.marching_cubes(np.pad(phi, 1, constant_values=1.0), 0.0, spacing=(NS.H_MM,) * 3)
    v = v - NS.H_MM + np.array([xs[0], ys[0], zs[0]])
    full = trimesh.Trimesh(v, f[:, ::-1], process=True)
    if full.volume < 0:
        full.invert()
    M = mf.Manifold(mf.Mesh(vert_properties=np.asarray(full.vertices, np.float32), tri_verts=np.asarray(full.faces, np.uint32)))
    half = NS._tm(M ^ NS._cube(-500, 500, 0.0, 500, -500, 500))
    body, ten = NS.cut_body(half, p)
    assert body.is_volume and ten is not None and ten.is_volume
    assert body.bounds[0, 0] * 1e3 == pytest.approx(p["x_cut_mm"], abs=1e-3)
    pn = NS.printed_nose(ff, p, ten)
    assert pn["nose"].is_volume and pn["_full"].is_volume
    b, n = body.copy(), pn["nose"].copy()
    b.apply_scale(1e3); n.apply_scale(1e3)
    xc = p["x_cut_mm"]

    def top(m, x, y):
        loc, _, _ = m.ray.intersects_location([[x, y, 100.0]], [[0, 0, -1.0]], multiple_hits=True)
        return float(loc[:, 2].max())
    for y in (0.3, 4.0, 7.0):
        assert abs(top(n, xc - 0.1, y) - top(b, xc + 0.1, y)) < 0.35
    over = NS._man(pn["nose"]) ^ NS._man(body)            # only the 0.3 mm the print runs into the foam
    assert over.volume() <= NS.OVERLAP_MM * 0.5 * p["bond"]["area_mm2"] * 1.3 + 1.0
    kg, com = NS.printed_mass_com(pn["nose"], p, 1.04)
    ahead = (NS._man(pn["nose"]) ^ NS._cube(-1e3, xc, -1e3, 1e3, -1e3, 1e3)).volume()
    assert kg * 1e3 == pytest.approx(2 * (ahead + NS._man(ten).volume()) / 1e3 * 1.04, rel=0.02)
    assert com[0] * 1e3 < xc + NS.TENON_LEN_MM


# ------------------------------------------------------------------------------------------- a built car
def _car():
    car = Path(CAR)
    S = json.loads((car / "summary.json").read_text())
    asm = json.loads((car / "parts" / "assembly.json").read_text())
    rec = json.loads((car / "layout_car.json").read_text())
    return car, S, asm, rec


@need_car
def test_car_is_legal_and_one_closed_solid():
    car, S, asm, rec = _car()
    assert rec["legal"], rec.get("failed")
    assert rec["cfd"]["one_piece"] and rec["cfd"]["watertight"]
    chk = S["legality"]["checks"]
    for k in ("T8.2_nose_overhang", "T8.5.1_nose_z", "T8.5.1_nose_half_width", "nose_joint_bond_stress",
              "T8.7_front_wing_clear_air_all_parts", "T8.8_front_wing_seen_from_front", "T4.1_single_body",
              "attached_nose"):
        assert chk[k]["pass"], (k, chk[k])


@need_car
def test_car_nose_split_and_masses():
    import trimesh
    car, S, asm, rec = _car()
    ns = asm["info"]["nose_split"]
    body = trimesh.load(str(car / "body_half.stl"), force="mesh")
    nose = trimesh.load(str(car / "parts" / "nose.stl"), force="mesh")
    assert body.is_volume and nose.is_volume
    assert body.bounds[0, 0] * 1e3 == pytest.approx(ns["x_cut_mm"], abs=0.05)              # the foam's face
    assert body.bounds[0, 0] * 1e3 >= ns["block_front_mm"] + 0.99                          # inside the block
    assert nose.bounds[1, 0] * 1e3 <= (ns["tenon"]["x1"] if ns.get("tenon") else ns["x_cut_mm"] + 0.31) + 0.05
    # masses: the print at its density in parts_mass_g and the race model's "printed"; the foam nose in the body
    assert S["parts_mass_g"]["nose"] == pytest.approx(ns["printed_g"], rel=1e-6)
    m = rec["mass"]
    assert m["nose_split"]["foam_g"] == pytest.approx(ns["foam_g"], rel=1e-6)
    assert m["items"]["printed"][0] == pytest.approx(
        asm["fixed_hardware_kwargs"]["rear_wing_mass_kg"] * 1e3, abs=2e-3)
    foam_ahead = 2 * abs(trimesh.intersections.slice_mesh_plane(
        body, [-1.0, 0, 0], [ns["ref_a_mm"] / 1e3, 0, 0], cap=True).volume) * 1e6 * 0.163
    assert ns["foam_g"] == pytest.approx(foam_ahead, rel=0.02, abs=0.01)
    # the referee's body mass (the field, unified_phi.compute_mass_com, in T3.6's mass as manufactured) and the race
    # model's (layout_car.mass: the mesh) agree: the field's cells ahead of the face and in the pocket are not foam
    wd = S["wheels"]
    field_g = (S["mass"]["manufactured_mass_g"] - 2 * wd["front"]["mass_g"] - 2 * wd["rear"]["mass_g"]
               - float(S["mass"].get("ballast_g") or 0.0) - asm["fixed_hardware_kwargs"]["halo_mass_kg"] * 1e3
               - sum(S["parts_mass_g"].values()))
    sup = 2 * sum(asm["info"].get("support_mass_g_each", []))        # (evo_patches/apply.py P1 adds the supports)
    assert min(abs(field_g - m["foam_g"]), abs(field_g - sup - m["foam_g"])) < 0.15


@need_car
def test_car_junction_is_flush():
    """The outer surface across the print-to-foam face: no jump bigger than the 1 mm grid's own error."""
    import trimesh
    car, S, asm, rec = _car()
    ns = asm["info"]["nose_split"]
    b = trimesh.load(str(car / "body_half.stl"), force="mesh"); b.apply_scale(1e3)
    n = trimesh.load(str(car / "parts" / "nose.stl"), force="mesh"); n.apply_scale(1e3)
    xc = ns["x_cut_mm"]

    def top(m, x, y):
        loc, _, _ = m.ray.intersects_location([[x, y, 100.0]], [[0, 0, -1.0]], multiple_hits=True)
        return float(loc[:, 2].max()) if len(loc) else np.nan

    def side(m, x, z):
        loc, _, _ = m.ray.intersects_location([[x, 100.0, z]], [[0, -1.0, 0]], multiple_hits=True)
        return float(loc[:, 1].max()) if len(loc) else np.nan
    for y in (0.05, 4.0):
        assert abs(top(n, xc - 0.15, y) - top(b, xc + 0.15, y)) < 0.5, y
    zs = np.arange(6.0, 22.0, 0.5)
    zw = zs[int(np.nanargmax([side(n, xc - 0.5, z) for z in zs]))]
    assert abs(side(n, xc - 0.15, zw) - side(b, xc + 0.15, zw)) < 0.5
