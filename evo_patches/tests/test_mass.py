"""Mass bookkeeping of the pinned evolution builder, with evo_patches/apply.py applied.

    EVO_CHECKOUT=<patched evolution checkout> [CAR_DIR=<built car>] [CAR_DIR_SHRUNK=<same layout, body shrunk>] \
        python -m pytest evo_patches/tests -q

CAR_DIR / CAR_DIR_SHRUNK: run_car / layout_car.build output folders (summary.json, layout_car.json, parts/,
body_half.stl; carhybrid/make_car.py puts them in <out>/evo). Without them the built-car tests are skipped and
only the cheap invariants run (CI).

(a) the referee's T3.6 weighs the car layout_car.mass weighs (supports included), within 0.05 g;
(b) the race mass falls when stage-2 numbers shrink the body of a car with no lead (and nothing resizes it);
(c) the invariants the fixes rely on: supports in the no-joints manufactured mass, a fixed skin never sized,
    printed parts weighed at the same density in the search path and the joints path.
"""
from __future__ import annotations

import importlib
import importlib.util
import json
import os
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pytest

EVO = os.environ.get("EVO_CHECKOUT")
CAR = os.environ.get("CAR_DIR")
CAR2 = os.environ.get("CAR_DIR_SHRUNK")
HERE = Path(__file__).resolve().parent
APPLY = HERE.parent / "apply.py"

pytestmark = pytest.mark.skipif(not EVO or not Path(EVO).is_dir(), reason="EVO_CHECKOUT not set")


def _use_evo():
    """The evolution checkout's modules, not this repo's (same flat names: car/build/run_car.py etc.)."""
    evo = Path(EVO).resolve()
    here_root = HERE.parents[1]
    for name, mod in list(sys.modules.items()):
        f = getattr(mod, "__file__", None) or ""
        if f and Path(f).resolve().is_relative_to(here_root) and not Path(f).resolve().is_relative_to(HERE):
            del sys.modules[name]
    sys.path[:] = [p for p in sys.path if not (p and Path(p).resolve().is_relative_to(here_root))]
    spec = importlib.util.spec_from_file_location("repo_paths", evo / "repo_paths.py")
    rp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rp)
    sys.modules["repo_paths"] = rp
    rp.add()
    sys.path.insert(0, str(evo))
    return evo


if EVO and Path(EVO).is_dir():
    _use_evo()

TOL_G = 0.05


def _car(path):
    p = Path(path)
    return p / "evo" if (p / "evo" / "summary.json").exists() else p


# ------------------------------------------------------------------ the patches are in
def test_patches_applied():
    r = subprocess.run([sys.executable, str(APPLY), EVO, "--check"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert all(line.endswith("already applied") for line in r.stdout.strip().splitlines()), r.stdout


# ------------------------------------------------------------------ (c) cheap invariants
def _fake_asm(sup_each_g=(2.393, 1.859), wheel_g=(1.178, 1.322), parts=None, halo_kg=0.00107):
    """An assembly as assembly.build returns it, its fixed_hardware_kwargs from assembly.fixed_hardware_kwargs."""
    parts = parts or {"nose": 5.6, "fwing": 1.5, "rwing": 3.7, "tethers": 0.15}
    wd = {"front": {"mass_g": wheel_g[0]}, "rear": {"mass_g": wheel_g[1]}}
    mf, mr = 2 * (wheel_g[0] + sup_each_g[0]) * 1e-3, 2 * (wheel_g[1] + sup_each_g[1]) * 1e-3
    fhk = {"rear_wing_mass_kg": sum(parts.values()) * 1e-3, "wheels_front_mass_kg": mf, "wheels_rear_mass_kg": mr,
           "wheels_axles_mass_kg": mf + mr, "halo_mass_kg": halo_kg}
    return {"parts_mass_g": dict(parts), "wheel_design": wd, "fixed_hardware_kwargs": fhk,
            "info": {"support_mass_g_each": list(sup_each_g)}}


def test_no_joints_manufactured_mass_includes_the_supports(monkeypatch, tmp_path):
    """P1: a --no-joints car as manufactured = body + printed parts + wheels + supports + halo, i.e. the race
    model's roll-up (layout_car.mass: foam + rear_wing_mass + wheels_front + wheels_rear + halo)."""
    import run_car as rc
    import unified_phi as up
    asm0 = _fake_asm()
    fake_p4 = types.SimpleNamespace(build=lambda *a, **k: json.loads(json.dumps(asm0)))
    monkeypatch.setitem(sys.modules, "assembly", fake_p4)
    monkeypatch.setitem(sys.modules, "joints", types.SimpleNamespace())
    body_g = 26.144
    monkeypatch.setattr(up, "compute_mass_com", lambda g: [types.SimpleNamespace(mass_kg=body_g * 1e-3)])
    a = types.SimpleNamespace(no_joints=True, ballast_g=0.0, W=138.0, x_front=55.4, d_halo=23.2)
    asm = rc.parts_and_joints(a, tmp_path, object(), tmp_path / "body_half.stl", "parts")
    f = asm["fixed_hardware_kwargs"]
    race_model_car_g = body_g + (f["rear_wing_mass_kg"] + f["wheels_front_mass_kg"] + f["wheels_rear_mass_kg"]
                                 + f["halo_mass_kg"]) * 1e3
    assert abs(rc.manufactured_mass_kg(asm) * 1e3 - race_model_car_g) < 1e-9
    assert abs(asm["joints"]["supports_g"] - 2 * sum(asm0["info"]["support_mass_g_each"])) < 1e-9


def test_fixed_skin_build_is_built_once_and_never_sized(monkeypatch, tmp_path):
    """(b)/P2: a fixed skin (every search build) builds the body once with no mass target, however far the car
    is from 48.2 g, and does not call it a failed sizing: what the design draws is what is weighed."""
    import run_car as rc
    import unified_phi as up
    targets = []

    def build_body(a, out):
        targets.append(a._body_target_kg)
        return types.SimpleNamespace(body=0.020), {}
    monkeypatch.setattr(rc, "build_body", build_body)
    monkeypatch.setattr(rc, "export_body", lambda g, p: {})
    monkeypatch.setattr(rc, "log", lambda *a, **k: None)
    monkeypatch.setattr(up, "compute_mass_com", lambda g: [types.SimpleNamespace(mass_kg=g.body)])
    monkeypatch.setattr(rc, "parts_and_joints", lambda a, out, g, *x, **k: {
        "body": g.body, "joints": {"mass_delta_kg": {}}, "info": {}})
    monkeypatch.setattr(rc, "manufactured_mass_kg", lambda asm: asm["body"] + 0.017)
    a = types.SimpleNamespace(skin_mm=0.0)
    geom, info, asm = rc.build_car(a, tmp_path)
    assert targets == [None]
    assert info["sized"] is False and info["mass_converged"] is None
    assert rc.mass_refusal(info["car_mass_error_g"] * 1e-3, True) is None


def test_param_body_fixed_skin_ignores_a_target():
    """param_body.build: a fixed skin returns make(skin) before any target is looked at."""
    import inspect
    import param_body as pb
    src = inspect.getsource(pb.build)
    i_skin = src.index("if skin_offset_mm is not None:")
    assert "return make(skin_offset_mm)" in src[i_skin:i_skin + 200]
    assert i_skin < src.index("if target_body_kg is None:")


def test_printed_parts_density_is_the_joints_density():
    """P3: the search path (assembly.py) weighs the nose, wings, struts and guides at the density the joints
    path (the finalist's manufacturing files) weighs them: ABS 1.04 g/cm3."""
    import assembly
    import component_contract as cc
    import joints
    rho = cc.DENSITY_G_CM3[assembly.PART_MATERIAL]
    for k in ("nose", "fwing", "rwing", "tethers"):
        assert abs(rho - joints.DENSITY_G_CM3[k]) < 1e-12, k


def test_race_mass_falls_with_the_foam_on_a_no_lead_car(tmp_path):
    """The race model itself: above the floor every gram of foam removed is a gram of race mass; below it the
    lead fills to 48.2 g (layout_car.mass), so the bookkeeping upstream is what has to be right."""
    import trimesh
    import layout_car as LC
    import layout_genes as LG
    d = LG.decode(LG.x_default())
    asm = _fake_asm()
    f = asm["fixed_hardware_kwargs"]
    f.update(rear_wing_com=[0.1, 0, 0.02], wheels_front_com=[0.055, 0, 0.014], wheels_rear_com=[0.19, 0, 0.014],
             halo_com=[0.09, 0, 0.034])
    (tmp_path / "parts").mkdir()
    (tmp_path / "parts" / "assembly.json").write_text(json.dumps(asm))
    out = []
    for h in (0.040, 0.036, 0.010):          # m: a foam box, shrunk, then far below the floor
        b = trimesh.creation.box(bounds=[[0.04, 0.0, 0.004], [0.24, 0.03, 0.004 + h]])
        b.export(str(tmp_path / "body_half.stl"))
        out.append(LC.mass(tmp_path, d))
    (m0, m1, m2) = out
    assert m0["lead_g"] == 0 and m1["lead_g"] == 0
    assert m1["race_mass_g"] < m0["race_mass_g"]
    assert abs((m0["race_mass_g"] - m1["race_mass_g"]) - (m0["foam_g"] - m1["foam_g"])) < 1e-9
    assert m2["lead_g"] > 0 and abs(m2["car_g"] + m2["lead_g"] - LC.FLOOR_G) < 1e-9


def test_stage2_shrinking_numbers_shrink_the_drawn_body_and_nose():
    """(b), cheap: the stage-2 numbers that pull the side in (u[30:45] < 0) take volume out of the drawn body
    behind Ref A and out of the nose ahead of it (freeform.field on a 1 mm grid, the team layout's seed)."""
    import freeform as fr
    import layout_genes as LG
    import simscale_shape as SS
    d = LG.decode(LG.x_default())
    seed = fr.Freeform.from_dict(d["freeform"])
    ref_a = d["x_front"] - 16.0
    h = 1.0
    X, Y, Z = np.meshgrid(np.arange(seed.x_tip_mm, seed.x_rear_mm + h, h), np.arange(0, 40, h) + h / 2,
                          np.arange(0, 55, h) + h / 2, indexing="ij")
    n = fr.n_genes(*SS.COARSE)

    def vol(u):
        ff = fr.from_unit(seed, u, SS.AMP_MM, *SS.COARSE)
        inside = fr.field(ff, X, Y, Z, ff.x_tip_mm, ff.x_rear_mm, axles=(d["x_front"], d["x_front"] + d["W"], None)) < 0
        return (inside & (X >= ref_a)).sum(), (inside & (X < ref_a)).sum()
    u = np.zeros(n); u[30:45] = -0.3
    (f0, n0), (f1, n1) = vol(np.zeros(n)), vol(u)
    assert f1 < 0.95 * f0 and n1 < n0


# ------------------------------------------------------------------ (a), (b), (c) on built cars
needs_car = pytest.mark.skipif(not CAR or not (_car(CAR or ".") / "summary.json").exists(), reason="CAR_DIR not set")


def _judged(car):
    import layout_car as LC
    S = json.loads((car / "summary.json").read_text())
    rec = json.loads((car / "layout_car.json").read_text())
    return S, rec, LC.mass(car, rec["decoded"])


@needs_car
def test_referee_T36_weighs_the_race_model_car():
    """(a): T3.6's mass (the referee) equals layout_car.mass's car_g (the race model), supports included."""
    car = _car(CAR)
    S, rec, m = _judged(car)
    t36 = S["legality"]["checks"]["T3.6_mass"]["margin"] + 48.0
    assert abs(t36 - m["car_g"]) <= TOL_G, (t36, m["car_g"])
    assert abs(S["mass"]["manufactured_mass_g"] - S["mass"]["competition_mass_g"]) <= TOL_G
    if m["car_g"] >= 48.0 + TOL_G:
        assert S["legality"]["checks"]["T3.6_mass"]["pass"]


@needs_car
def test_built_supports_are_in_both_rollups():
    """(c): the wheel masses the race model uses are the wheels plus the supports the support engine weighed,
    and those match the part files (part + hubcap at 1.04 g/cm3, less the slot's foam at 0.163)."""
    import trimesh
    car = _car(CAR)
    asm = json.loads((car / "parts" / "assembly.json").read_text())
    f, wd = asm["fixed_hardware_kwargs"], asm["wheel_design"]
    sup = asm["info"]["support_mass_g_each"]
    assert abs(f["wheels_front_mass_kg"] * 1e3 - 2 * (wd["front"]["mass_g"] + sup[0])) < 1e-6
    assert abs(f["wheels_rear_mass_kg"] * 1e3 - 2 * (wd["rear"]["mass_g"] + sup[1])) < 1e-6
    se = asm["info"]["support_engine"]
    rho = se["design"]["material_g_cm3"]
    whole = trimesh.load(str(car / "parts" / "supports.stl"), force="mesh")
    assert abs(abs(whole.volume) * 1e9 * rho * 1e-3
               - sum(se[t]["part_g_each"] - _cap_g(car, t, rho) for t in ("front", "rear"))) < 0.02
    for tag, g in zip(("front", "rear"), sup):
        slot = abs(trimesh.load(str(car / "parts" / f"support_{tag}_slot.stl"), force="mesh").volume) * 1e9
        assert abs(se[tag]["foam_removed_g_each"] - slot * 0.163e-3) < 0.005
        assert abs(g - se[tag]["net_g_each"]) < 1e-9


def _cap_g(car, tag, rho):
    import trimesh
    p = car / "parts" / f"hubcap_{tag}.stl"
    return abs(trimesh.load(str(p), force="mesh").volume) * 1e9 * rho * 1e-3 if p.exists() else 0.0


@needs_car
def test_built_printed_parts_are_weighed_at_abs():
    """(c)/P3: the printed nose (with the front wing's pylons) and the wings weigh their files' volume x 1.04."""
    import trimesh
    car = _car(CAR)
    S = json.loads((car / "summary.json").read_text())
    for k in ("nose", "fwing", "rwing"):
        v = abs(trimesh.load(str(car / "parts" / f"{k}.stl"), force="mesh").volume) * 1e9 * 2   # both halves
        assert abs(S["parts_mass_g"][k] - v * 1.04e-3) < 0.03 * S["parts_mass_g"][k], k


@pytest.mark.skipif(not CAR or not CAR2, reason="CAR_DIR and CAR_DIR_SHRUNK not set")
def test_shrunk_body_lowers_the_race_mass():
    """(b): the same layout with the body shrunk by stage-2 numbers: less foam, lower race mass, no lead, and the
    body built at its fixed skin (not resized back to a target)."""
    S0, r0, m0 = _judged(_car(CAR))
    S1, r1, m1 = _judged(_car(CAR2))
    assert r0["layout"] == r1["layout"]
    assert m0["lead_g"] == 0 and m1["lead_g"] == 0
    assert m1["foam_g"] < m0["foam_g"] - 0.5
    assert m1["race_mass_g"] < m0["race_mass_g"] - 0.5
    for S in (S0, S1):
        assert S["body"]["build_passes"] == 1 and S["body"]["build"]["skin_offset_mm"] == 0.0
