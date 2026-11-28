"""P1 (carhybrid/DESIGN.md): the pre-registered gradient check's cases, from a P0 gradient (grad.json).

    python carhybrid/p1_cases.py <grad.json> <layout> <batch name> [--h 0.25] [--hd 0.5] [--seed 20261128]
                                 [--frames 0,1] [--grid gu,gv,gw,gb] > cases.json

--frames: the frames every direction is scored on (default 0,1 = P1; P1F, carhybrid/FINE_ADJOINT.md: 0,1,2,3).
--grid: the free-form grid of the gradient's numbers (default stage 2's 5,2,3,4 = 45 numbers; "24,9,12,9" = the
  seed's own 720); written into every case as "grid" when given, so the build uses the same grid (it must match the
  length of grad.json's dT_du_s). Without these options the output is P1's, unchanged.
Directions in the gradient's numbers (unit length), each built at u = h * d and scored on every frame:
  T-   along the race-time descent      -dT/du / |dT/du|        (DESIGN.md's P1)
  T+   against it                       +dT/du / |dT/du|
  R1, R2  two random unit directions (numpy default_rng(seed), fixed before any result)
  D-   along the drag descent           -dD/du / |dD/du|        (added: the race-time direction is ~90 % mass on a
  D+   against it                       +dD/du / |dD/du|         car with no lead, so T+- tests the mass map; D+- tests
                                                                  the adjoint's drag)
  S1, S2  two more random unit directions, at the drag arm's step
  B    the base car (u = 0), the reference for the measured slopes
T and R at h (DESIGN.md: 0.25); D and S at hd (0.5: the drag arm's predicted half-car change at 0.25 is -2.2 mN, near
the frame-to-frame noise if the adjoint overstates size 2x as on the sidepod).
plus one DAFoam gradient on the base car as built now (the P0 gradient was computed before evo_patches' mass and
cap fixes): its cosine with the P0 direction is reported, so the directions' provenance is checked, not assumed.
"""
import json, sys

import numpy as np


def main(grad_p, layout, batch, h=0.25, hd=0.5, seed=20261128, frames=(0, 1), grid=None):
    g = json.load(open(grad_p))
    dT, dD = np.asarray(g["dT_du_s"]), np.asarray(g["dD_du_N"])
    n = len(dT)
    if grid is not None:
        gg = [int(v) for v in str(grid).split(",")]
        shapes = [(gg[0], gg[1]), (gg[0], gg[3] if len(gg) > 3 else gg[1]), (gg[0], gg[2])]   # freeform.coarse_grid_shapes
        if sum(a * b for a, b in shapes) != n:
            raise SystemExit(f"grid {grid} has {sum(a * b for a, b in shapes)} numbers, the gradient {n}")
    rng = np.random.default_rng(seed)
    r1, r2, r3, r4 = (rng.standard_normal(n) for _ in range(4))
    un = lambda v: v / np.linalg.norm(v)
    dirs = {"Tm": (h, -un(dT)), "Tp": (h, un(dT)), "R1": (h, un(r1)), "R2": (h, un(r2)),
            "Dm": (hd, -un(dD)), "Dp": (hd, un(dD)), "S1": (hd, un(r3)), "S2": (hd, un(r4)), "B": (0.0, np.zeros(n))}
    cases = []
    for name, (s, d) in dirs.items():
        u = json.dumps([round(float(v), 6) for v in s * d])
        for k in frames:
            c = {"name": f"{batch}_{name}_k{k}", "type": "score", "layout": layout, "k": k, "u": u}
            if grid is not None:
                c["grid"] = str(grid)
            cases.append(c)
    c = {"name": f"{batch}_B_dafoam", "type": "dafoam", "layout": layout}
    if grid is not None:
        c["grid"] = str(grid)
    cases.append(c)
    return cases


if __name__ == "__main__":
    a = sys.argv[1:]
    kw = {}
    if "--h" in a:
        kw["h"] = float(a[a.index("--h") + 1])
    if "--hd" in a:
        kw["hd"] = float(a[a.index("--hd") + 1])
    if "--seed" in a:
        kw["seed"] = int(a[a.index("--seed") + 1])
    if "--frames" in a:
        kw["frames"] = tuple(int(v) for v in a[a.index("--frames") + 1].split(","))
    if "--grid" in a:
        kw["grid"] = a[a.index("--grid") + 1]
    print(json.dumps(main(a[0], a[1], a[2], **kw), indent=1))
