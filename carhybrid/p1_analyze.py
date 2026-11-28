"""P1 analysis against the marks in carhybrid/P1.md (fixed before any result).

    CARHYBRID_EVO=<evolution checkout> python carhybrid/p1_analyze.py carhybrid/results/p1 [more result folders]

Several folders are merged (a later folder's case replaces an earlier one's: re-runs such as p1r). Prints the
table and the marks; writes <first folder>/p1_analysis.json. Race time uses loop.py's frame_score (the evolution race
model at each car's own mass, frame-k drag and lift); paired comparisons over frames 0 and 1 with
SE = max(block SE, |Delta_0 - Delta_1| / 2) as P1.md fixes.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import loop  # noqa: E402

NAMES = ["B", "Tm", "Tp", "R1", "R2", "Dm", "Dp", "S1", "S2"]
K = (0, 1)
PRED = {"Tm-B_ms": -6.7, "Dm-B_mN_whole": -8.8}


def gather(folders, evo):
    out = {}
    for f in folders:
        for n in NAMES:
            for k in K:
                d = Path(f) / f"p1_{n}_k{k}"
                if not (d / "result.json").exists():
                    continue
                T, se, info = loop.frame_score(f, f"p1_{n}_k{k}", evo)
                if T is None:
                    continue
                r = json.load(open(d / "result.json"))
                out[(n, k)] = {"T": T, "seT": se, "D": 2 * float(r["drag_mN"]),
                               "seD": 2 * float(r.get("drag_block_se_mN") or 0), "L": 2 * float(r.get("lift_mN") or 0),
                               "m": info.get("mass_g"), "src": str(f)}
    return out


def pair(V, a, b, q):
    """Delta = a - b on q ('T' s or 'D' mN), paired by frame: (mean, SE) or None."""
    if not all((x, k) in V for x in (a, b) for k in K):
        return None
    d = [V[(a, k)][q] - V[(b, k)][q] for k in K]
    se_b = np.mean([np.hypot(V[(a, k)]["se" + q], V[(b, k)]["se" + q]) for k in K]) / np.sqrt(len(K))
    return float(np.mean(d)), float(max(se_b, abs(d[0] - d[1]) / 2)), d


def main(argv):
    folders = argv[1:] or [str(HERE / "results" / "p1")]
    evo = os.environ.get("CARHYBRID_EVO", "evo")
    V = gather(folders, evo)
    print("case  frame  D_whole_mN  L_whole_mN  mass_g  T_s")
    for n in NAMES:
        for k in K:
            v = V.get((n, k))
            print(f"{n:3s} k{k}  " + ("missing" if v is None else
                  f"{v['D']:8.2f}  {v['L']:7.2f}  {v['m']}  {v['T']:.5f}"))
    res = {"cases": {f"{n}_k{k}": v for (n, k), v in V.items()}}

    def cmp(a, b, q):
        p = pair(V, a, b, q)
        if p:
            s = 1e3 if q == "T" else 1.0
            print(f"{a} - {b} [{q}]: {p[0]*s:+.2f} +- {p[1]*s:.2f} ({'ms' if q == 'T' else 'mN whole'}; "
                  f"frames {', '.join(f'{x*s:+.2f}' for x in p[2])})")
            res[f"{a}-{b}_{q}"] = {"mean": p[0] * s, "se": p[1] * s, "frames": [x * s for x in p[2]]}
        return p

    tmtp = cmp("Tm", "Tp", "T"); tr = [cmp("Tm", r, "T") for r in ("R1", "R2")]
    dmdp = cmp("Dm", "Dp", "D"); dr = [cmp("Dm", s, "D") for s in ("S1", "S2")]
    tb = cmp("Tm", "B", "T"); db = cmp("Dm", "B", "D")
    cmp("Dp", "B", "D"); cmp("Tp", "B", "T"); cmp("Dm", "B", "T")

    def mark(main_, others):
        if main_ is None or any(o is None for o in others):
            return None
        return bool(main_[0] <= -3 * main_[1] and all(o[0] < 0 for o in others))
    res["P1-T"] = mark(tmtp, tr)
    res["P1-D"] = mark(dmdp, dr)
    if tb:
        res["slope_T"] = tb[0] * 1e3 / PRED["Tm-B_ms"]
    if db:
        res["slope_D"] = db[0] / PRED["Dm-B_mN_whole"]
    print(json.dumps({k: res.get(k) for k in ("P1-T", "P1-D", "slope_T", "slope_D")}))
    json.dump(res, open(Path(folders[0]) / "p1_analysis.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main(sys.argv)
