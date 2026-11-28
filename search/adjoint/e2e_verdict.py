"""E2E_TEST.md's marks from the downloaded results.
    python e2e_verdict.py <results dir> <out.json>
<results dir>/adjoint/result.json (the base job), <results dir>/fd_cases.json, <results dir>/fd/<id>/result.json."""
import json
import math
import sys
from pathlib import Path

NOISE_MN = 0.32          # coarse mesh-to-mesh sd (half car)


def main(root, outp):
    root = Path(root)
    base = json.loads((root / "adjoint" / "result.json").read_text())
    d = base.get("dafoam") or {}
    st = base.get("steady") or {}
    g = base.get("gradient") or {}
    pct = d.get("primal_vs_openfoam_pct")
    e1 = {"dafoam_ok": bool(d.get("ok")), "steady": bool(st.get("steady")), "scheme": st.get("scheme"),
          "primal_vs_openfoam_pct": pct, "primal_within_2pct": pct is not None and abs(pct) <= 2.0,
          "n_gradients": len(g.get("coarse_half_mN_per_unit") or [])}
    e1["pass"] = e1["dafoam_ok"] and e1["steady"] and e1["primal_within_2pct"] and e1["n_gradients"] == 45
    cases = json.loads((root / "fd_cases.json").read_text()) if (root / "fd_cases.json").exists() else []
    by = {}
    for c in cases:
        p = root / "fd" / c["id"] / "result.json"
        r = json.loads(p.read_text()) if p.exists() else {}
        D = (r.get("steady") or {}).get("D_half_mN") if (r.get("steady") or {}).get("steady") else None
        by.setdefault(c["gene"], {"delta": c["delta"], "pred": c["pred_half_mN_per_unit"]})[c["sign"]] = D
    rows = []
    for k, v in by.items():
        if v.get(1) is None or v.get(-1) is None:
            rows.append({"gene": k, "measured": None, "pred": v["pred"], "why": "a forward run failed or was not steady"})
            continue
        s = (v[1] - v[-1]) / (2 * v["delta"])
        se = NOISE_MN * math.sqrt(2) / (2 * v["delta"])
        rows.append({"gene": k, "delta": v["delta"], "D_plus": v[1], "D_minus": v[-1], "measured": s, "se": se,
                     "pred": v["pred"], "ratio": v["pred"] / s if s else None, "evidence": abs(s) >= 2 * se,
                     "same_sign": (s > 0) == (v["pred"] > 0)})
    ev = [r for r in rows if r.get("evidence")]
    if len(rows) < 3 or any(r.get("measured") is None for r in rows):
        e2 = {"verdict": "INCOMPLETE", "pass": False}
    elif len(ev) < 2:
        e2 = {"verdict": "INCONCLUSIVE", "pass": False}
    else:
        signs = all(r["same_sign"] for r in ev)
        inband = sum(1 for r in ev if r["ratio"] is not None and 0.5 <= r["ratio"] <= 4.0)
        e2 = {"signs_all_same": signs, "ratio_in_band": inband, "n_evidence": len(ev),
              "pass": signs and inband >= 2, "verdict": "PASS" if signs and inband >= 2 else "FAIL"}
    res = {"E1": e1, "E2": e2, "fd": rows, "PASS": e1["pass"] and e2["pass"],
           "reported": {"production": base.get("production"), "steady": st, "dafoam": d, "wall_s": base.get("wall_s"),
                        "n_active": g.get("n_active"), "n_points": g.get("n_points"),
                        "active_share_of_dD": g.get("active_share_of_dD"), "error": base.get("error")}}
    Path(outp).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({"E1": e1["pass"], "E2": e2["verdict"], "PASS": res["PASS"]}))
    for r in rows:
        print(json.dumps(r, default=float))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
