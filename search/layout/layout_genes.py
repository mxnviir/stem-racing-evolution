"""Genes of the whole-car search, read from layout_params.json (see layout_catalogue.py).

Every parameter whose mode is "evolve" is one gene in [0, 1] (CMA works on that cube, lo..hi);
"fixed" ones are held at their value. `decode` turns a gene vector into everything Parts 1-5
build from. The LAYOUT (W, rear, halo_rear) is decoded so every vector is a legal layout, and
`legal_layouts` lists every legal layout on a grid: the space the search moves in.

Layout, mm, x from the nose tip (car frame of Parts 1-5):
  W          wheelbase, front axle -> rear axle                      T7.3: 120.0-140.0
  rear       rear axle -> body rear face                             36 .. 207 - W (the 223 mm block)
  halo_rear  halo pocket FRONT edge -> body rear face
  x_front    = 16 + nose length (the nose tip is x = 0, Ref A is x_front - 16)

  d_halo (Part 1's halo parameter, Ref A -> pocket front) = W + rear + 16 - halo_rear.
  Halo limits (car/body/geometry_contract.py, project owner 2026-07-24):
    pocket front on or behind the front axle          d_halo >= 16          -> halo_rear <= W + rear
    pocket rear (front + 50) ahead of the rear axle   d_halo < W - 34       -> halo_rear > rear + 50
    pocket rear ahead of the cartridge bore (face - 52)  d_halo < W + rear - 86 -> halo_rear > 102
    16 + W + rear <= 223   T1.14: Ref A to the rear face fits the 223 mm block (never binds)
Cargo (T4.2): 60 long, top at 24 under the halo; cargo_x 0..1 slides it from its front on the
front axle line to its back on the rear axle line; cargo_flip >= 0.5 puts the wide end rearward.
"""
from __future__ import annotations

import sys
import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PARAMS = json.loads((HERE / "layout_params.json").read_text())["params"]
BY = {p["name"]: p for p in PARAMS}
GENES = [p for p in PARAMS if p["mode"] == "evolve"]
SHAPE_GENES = [p for p in PARAMS if p["mode"] == "shape"]    # held at value here; evolved by stage 2 (search/shape)
NAMES = [p["name"] for p in GENES]
LO = np.array([p["lo"] for p in GENES], float)
HI = np.array([p["hi"] for p in GENES], float)
DIM = len(GENES)
BLOCK_MM = 223.0
CARGO_LEN_MM, CARGO_Z_BASE_MM = 60.0, 14.0


def rear_range(W):
    """[lo, hi]: rear axle -> rear face for this wheelbase. The body from Ref A to its rear face must fit the
    223 mm block (T1.14); T9.4.2's 40 mm behind Ref B is the rear wing and its support only (team, 2026-11-22)."""
    return BY["rear"]["lo"], min(BY["rear"]["hi"], BLOCK_MM - 16.0 - W)


def halo_rear_range(W, rear):
    """(lo, hi]: legal halo pocket front -> rear face distances for this W and rear."""
    return max(rear + 50.0, 102.0), W + rear


def layout_legal(W, rear, halo_rear) -> bool:
    lo, hi = halo_rear_range(W, rear)
    rlo, rhi = rear_range(W)
    return (BY["W"]["lo"] <= W <= BY["W"]["hi"] and rlo <= rear <= rhi + 1e-9
            and lo < halo_rear <= hi + 1e-9 and 16.0 + W + rear <= BLOCK_MM + 1e-9)


def legal_layouts(step=1.0):
    """Every legal (W, rear, halo_rear) on a `step` mm grid."""
    out = []
    for W in np.arange(BY["W"]["lo"], BY["W"]["hi"] + 1e-9, step):
        for r in np.arange(rear_range(W)[0], rear_range(W)[1] + 1e-9, step):
            lo, hi = halo_rear_range(W, r)
            for h in np.arange(np.floor(lo / step) * step + step, hi + 1e-9, step):
                if layout_legal(W, r, h):
                    out.append((float(W), float(r), float(h)))
    return out


def _start(p):
    """A parameter's start value in [0, 1]. halo_rear's gene is the fraction of its legal interval."""
    if p["name"] == "rear":
        lo, hi = rear_range(BY["W"]["value"])
        return float(np.clip((p["value"] - lo) / (hi - lo), 0, 1))
    if p["name"] == "halo_rear":
        W, r = BY["W"]["value"], BY["rear"]["value"]
        lo, hi = halo_rear_range(W, r)
        return float(np.clip((p["value"] - lo - 0.5) / (hi - lo - 0.5), 0, 1))
    return float(np.clip((p["value"] - p["lo"]) / (p["hi"] - p["lo"]), 0, 1))


def x_default():
    """The start of the search: every evolved parameter at its value (the team-like car)."""
    return np.array([_start(p) for p in GENES])


def x_of(**values):
    """Gene vector with the named genes set (physical values; halo_rear as its 0..1 fraction of the
    legal interval), the rest at their defaults. For tests and hand-picked cars."""
    x = x_default()
    W = float(values.get("W", BY["W"]["value"]))
    for k, v in values.items():
        i = NAMES.index(k)
        if k == "halo_rear":
            x[i] = float(v)
        elif k == "rear":
            lo, hi = rear_range(W)
            x[i] = (float(v) - lo) / (hi - lo)
        else:
            x[i] = (float(v) - LO[i]) / (HI[i] - LO[i])
    if "W" in values and "rear" not in values and "rear" in NAMES:      # keep the default rear length, not its fraction
        lo, hi = rear_range(W)
        x[NAMES.index("rear")] = (BY["rear"]["value"] - lo) / (hi - lo)
    return np.clip(x, 0.0, 1.0)


def decode(x) -> dict:
    """[0, 1]^DIM (clipped) -> layout + Part 1 / Part 4 parameters. Always a legal layout."""
    x = np.clip(np.asarray(x, float), 0.0, 1.0)
    v = {p["name"]: p["value"] for p in PARAMS}
    t = dict(zip(NAMES, x))
    for p, xi in zip(GENES, x):
        if p["name"] not in ("halo_rear", "rear"):
            v[p["name"]] = p["lo"] + xi * (p["hi"] - p["lo"])
    W = round(float(v["W"]), 1)
    rlo, rhi = rear_range(W)
    rear = rlo + t["rear"] * (rhi - rlo) if "rear" in t else float(np.clip(v["rear"], rlo, rhi))
    rear = float(np.floor(rear * 10 + 0.5) / 10)
    if rear > rhi:
        rear = float(np.floor(rhi * 10) / 10)
    lo, hi = halo_rear_range(W, rear)
    if "halo_rear" in t:
        halo_rear = round(lo + 0.5 + t["halo_rear"] * (hi - lo - 0.5), 1)    # strictly above lo
    else:
        halo_rear = float(np.clip(v["halo_rear"], lo + 0.5, hi))
    nose_len = round(float(v["nose_len"]), 1)
    x_front = 16.0 + nose_len
    d_halo = W + rear + 16.0 - halo_rear
    s_x0u = float(v["s_x0u"]); s_x1u = min(1.0, s_x0u + float(v["s_lenu"]))
    body = dict(
        st_b=tuple(float(v[f"st_b{i}"]) for i in range(8)),
        st_zt=tuple(max(float(v[f"st_zt{i}"]), float(v[f"st_zb{i}"]) + 6.0) for i in range(8)),   # top above bottom
        st_zb=tuple(float(v[f"st_zb{i}"]) for i in range(8)),
        p=float(v["p"]), blend_mm=float(v["blend_mm"]),
        s_b=float(v["s_b"]), s_zt=max(float(v["s_zt"]), float(v["s_zb"]) + 6.0), s_zb=float(v["s_zb"]),
        s_taper=float(v["s_taper"]), s_p=float(v["s_p"]),
        ft_h=float(v["ft_h"]), ft_w=float(v["ft_w"]), ft_x0=float(v["ft_x0"]), ft_x1=float(v["ft_x1"]),
        ft_ramp=float(v["ft_ramp"]),
        modes=tuple(float(v[f"mode{i}"]) for i in range(56)))
    nose = dict(length_mm=nose_len, k=float(v["nose_k"]), p=float(v["nose_p"]), tip_z_mm=float(v["nose_tip_z"]),
                blend_after_ref_a_mm=float(v["nose_blend"]), root_scale=float(v["nose_root_scale"]),
                wall_mm=float(v["nose_wall"]))
    # the wings' rule limits on their placement (owner 2026-11-28: as far up, down, forward, back, in and out as the
    # rules let them): T8.2, the front wing within 40 mm ahead of Ref A; T9.4.2, the rear wing within 40 mm behind Ref B
    v["fw_te_gap"] = min(float(v["fw_te_gap"]), FRONT_REACH_MM - float(v["fw_chord"]) - RULE_MARGIN_MM)
    v["rw_x_le"] = min(float(v["rw_x_le"]), REAR_REACH_MM - float(v["rw_chord"]) - RULE_MARGIN_MM)
    fwing = front_wing_design(v)
    rwing = rear_wing_design(v)              # the searches' rear wing (rear_wing_engine): height, span, chord, place
    cargo_x0 = x_front + float(v["cargo_x"]) * (W - CARGO_LEN_MM)
    # a preference: Part 1 puts the cargo at the legal position nearest it (unified_phi, "prefer")
    cargo = {"x_start_m": cargo_x0 / 1e3, "z_base_m": CARGO_Z_BASE_MM / 1e3, "flip": bool(float(v["cargo_flip"]) >= 0.5),
             "prefer": True}
    tethers = (float(v["tether_front"]), float(v["tether_rear"]))
    support = support_design(v)
    out = dict(W=W, rear=rear, halo_rear=halo_rear, x_front=x_front, d_halo=round(d_halo, 2),
               s_x0u=s_x0u, s_x1u=s_x1u, body=body, nose=nose, fwing=fwing, rwing=rwing, cargo=cargo,
               tethers=tethers, support=support, skin_mm=float(v["skin_mm"]), legal_layout=layout_legal(W, rear, halo_rear))
    ff = organic_freeform(out)
    if ff is not None:
        out["freeform"] = ff
    return out


FRONT_REACH_MM = 40.0          # T8.2: nothing more than this ahead of Ref A
REAR_REACH_MM = 40.0           # T9.4.2: the rear wing within this behind Ref B
RULE_MARGIN_MM = 0.3
NOSE_LEAD_MM = 2.0             # the organic nose ends this far ahead of the front wing's leading edge
ORGANIC_SEED = HERE / "organic_seed.json"
_SEED = None


def organic_freeform(d) -> dict | None:
    """No cone (owner 2026-11-28): the car's body and nose as free-form surfaces (car/body/freeform.py): the team
    body's fitted surfaces (organic_seed.json) stretched from the nose tip to THIS car's rear face, and a nose
    grown from the body's own Ref A section to just ahead of THIS car's front wing (freeform.seed_nose)."""
    global _SEED
    if not ORGANIC_SEED.exists():
        return None
    if _SEED is None:
        _SEED = json.loads(ORGANIC_SEED.read_text())
    import freeform as fr
    ref_a = d["x_front"] - 16.0
    x_tip = ref_a - d["nose"]["length_mm"]
    ff = fr.Freeform.from_dict(dict(_SEED["freeform"], x_tip_mm=x_tip, x_rear_mm=d["x_front"] + d["W"] + d["rear"]))
    x_le = ref_a - d["fwing"]["te_gap_mm"] - d["fwing"]["chord_mm"]
    ff = fr.seed_nose_from_section(ff, _SEED["section"], ref_a, max(x_le - NOSE_LEAD_MM, x_tip))
    return ff.as_dict()


def support_design(v=None) -> dict:
    """support_engine.SupportDesign's numbers from parameter values (the catalogue's values by default)."""
    v = v or {p["name"]: p["value"] for p in PARAMS}
    return dict(chord_mm=float(v["sup_chord"]), t_mm=float(v["sup_t"]), taper=float(v["sup_taper"]),
                fillet_mm=float(v["sup_fillet"]),
                saddle_h_mm=float(v["sup_saddle_h"]), disc_r_mm=float(v["sup_disc_r"]),
                disc_recess_mm=float(v["sup_disc_recess"]))


def _strut(v, prefix):
    """wing_strut's numbers (front_wing_engine / rear_wing_engine .Support) from parameter values."""
    return dict(aim_out_deg=float(v[f"{prefix}_out"]), aim_fwd_deg=float(v[f"{prefix}_fwd"]),
                leave=float(v[f"{prefix}_leave"]), arrive=float(v[f"{prefix}_arrive"]),
                chord_mm=float(v[f"{prefix}_chord"]), thick_mm=float(v[f"{prefix}_thick"]),
                flare=float(v[f"{prefix}_flare"]), round=float(v[f"{prefix}_round"]),
                mid_x=float(v[f"{prefix}_mid_x"]), mid_y=float(v[f"{prefix}_mid_y"]), mid_z=float(v[f"{prefix}_mid_z"]),
                swell=float(v[f"{prefix}_swell"]),
                **{k: float(v[f"{prefix}_{k}"]) for k in ("root", "n") if f"{prefix}_{k}" in v})


def _supports(v, prefix, n):
    return [dict(width_mm=float(v[f"{prefix}{k}_width"]), y_frac=float(v[f"{prefix}{k}_y"]),
                 chord_frac=float(v[f"{prefix}{k}_chord"]), foot_scale=float(v[f"{prefix}{k}_foot"]),
                 waist=float(v[f"{prefix}{k}_waist"])) for k in range(n)]


def front_wing_design(v=None) -> dict:
    """front_wing_engine.FrontWingDesign's numbers (straight; its loft supports) from parameter values."""
    v = v or {p["name"]: p["value"] for p in PARAMS}
    return dict(pieces=2 if float(v["fw_pieces"]) >= 0.5 else 1, chord_mm=float(v["fw_chord"]),
                te_gap_mm=float(v["fw_te_gap"]), z_mm=float(v["fw_z"]), tip_y_mm=float(v["fw_tip_y"]),
                root_y_mm=float(v["fw_root_y"]),
                support=_strut(v, "fwsup"))


N_FW_STRUTS = 2


def rear_wing_design(v=None) -> dict:
    """rear_wing_engine.RearWingDesign's numbers (straight; its loft supports) from parameter values."""
    v = v or {p["name"]: p["value"] for p in PARAMS}
    return dict(chord_mm=float(v["rw_chord"]), t_mm=float(v["rw_t"]), half_span_mm=float(v["rw_half_span"]),
                x_le_mm=float(v["rw_x_le"]), lift_mm=float(v["rw_lift"]),
                support=_strut(v, "rwsup"))


N_RW_STRUTS = 2


def shape_from_unit(u, z_fw=None, rw_x_le=None, stage1=None) -> dict:
    """Stage 2: the shape genes (SHAPE_GENES, in order) from CMA's numbers: 0 is the catalogue's value,
    -1 its lo, +1 its hi (clipped), linear on each side. Returns {"support": ..., "fwing": ..., "rwing": ...};
    z_fw and rw_x_le keep stage 1's front wing height and rear wing position."""
    v = {p["name"]: p["value"] for p in PARAMS}
    for p, ui in zip(SHAPE_GENES, np.clip(np.asarray(u, float), -1.0, 1.0)):
        v[p["name"]] = p["value"] + ui * ((p["hi"] - p["value"]) if ui > 0 else (p["value"] - p["lo"]))
    if z_fw is not None:
        v["fw_z"] = z_fw
    if rw_x_le is not None:
        v["rw_x_le"] = rw_x_le
    if stage1 is not None:                    # the layout search's genes (layout x vector): their decoded values win
        dd = decode(stage1)
        fw, rw = dd["fwing"], dd["rwing"]
        v.update(fw_pieces=1.0 if fw["pieces"] == 2 else 0.0, fw_chord=fw["chord_mm"], fw_te_gap=fw["te_gap_mm"],
                 fw_z=fw["z_mm"], fw_tip_y=fw["tip_y_mm"], fw_root_y=fw["root_y_mm"], rw_chord=rw["chord_mm"],
                 rw_half_span=rw["half_span_mm"], rw_x_le=rw["x_le_mm"], rw_lift=rw["lift_mm"])
    return {"support": support_design(v), "fwing": front_wing_design(v), "rwing": rear_wing_design(v)}


def support_from_unit(u) -> dict:
    """Stage 2: the shape genes from CMA's numbers, one per SHAPE_GENES entry: 0 is the catalogue's value,
    -1 its lo, +1 its hi (clipped), linear on each side."""
    v = {p["name"]: p["value"] for p in PARAMS}
    for p, ui in zip(SHAPE_GENES, np.clip(np.asarray(u, float), -1.0, 1.0)):
        v[p["name"]] = p["value"] + ui * ((p["hi"] - p["value"]) if ui > 0 else (p["value"] - p["lo"]))
    return support_design(v)


def body_params(d):
    """BodyParams dict for Part 1: the sidepod's x span from its u fractions of Ref A -> rear face."""
    ref_a, rear_face = d["x_front"] - 16.0, d["x_front"] + d["W"] + d["rear"]
    L = rear_face - ref_a
    b = dict(d["body"]); b["s_x0"] = ref_a + d["s_x0u"] * L; b["s_x1"] = ref_a + d["s_x1u"] * L
    return b


if __name__ == "__main__":
    L = legal_layouts(1.0)
    print(f"{DIM} genes of {len(PARAMS)} parameters; {len(L)} legal layouts on a 1 mm grid "
          f"(W {BY['W']['lo']}-{BY['W']['hi']}, rear {BY['rear']['lo']}-{BY['rear']['hi']}); 2 mm grid {len(legal_layouts(2.0))}")
    d = decode(x_default())
    print({k: d[k] for k in ("W", "rear", "halo_rear", "x_front", "d_halo", "cargo", "legal_layout")})
    rng = np.random.default_rng(0)
    assert all(decode(rng.random(DIM))["legal_layout"] for _ in range(2000)), "decode left the legal set"
    print("2000 random vectors: all legal layouts")
