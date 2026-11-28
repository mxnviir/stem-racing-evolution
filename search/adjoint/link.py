"""link.py -- the adjoint's drag sensitivity, per CFD surface point, turned into the gradient of drag with
respect to the free-form body's numbers (owner 2026-11-28: "Link the adjoint", fine-mesh DAFoam on AWS).

DAFoam (adjoint repo, research/sidepod_opt/dafoam/runScript.py) returns dD/dX_i for every point i of the
`car` patch: the exact derivative of its half-car drag with respect to moving that point, so for any small
change of the surface  dD = sum_i dD/dX_i . dX_i.  This module supplies dX_i for each body number.

The body (car/body/freeform.py, param_body.build) is a level set: the solid is phi(x; c) < level, where phi is
freeform.field (top, bottom, side control grids c, blended by a smooth max), cut at the milled body's front
(Ref A, or a foam nose's limit) and rear face, and level is the build's fixed skin (0 in the searches). A
point x on that surface moves, when the numbers change by dc, along the normal n = grad phi / |grad phi| by

    dX = -(dphi/dc . dc) / |grad phi| * n        so        dD/dc = sum_i  w_i * dphi_i/dc,
    w_i = -(dD/dX_i . grad phi_i) / |grad phi_i|^2

(only the normal part of DAFoam's dD/dX counts: a tangential slide of the points does not change the shape).
A point is ACTIVE when the free-form level set is what made the surface there: |phi(x) - level| < tol. Elsewhere
the surface is something no body number moves (the wings, nose mounts, supports, halo, the cut planes, the
cargo bore, the halo pocket, the T7.9 wheel zones), and dX = 0. The milling ball (3.175 mm) rounds concave
corners: points it moved by more than tol are inactive too (a small, conservative loss, reported).

Two gradients come out:
  * fine: one per control point of the seed's own grid (FIT_GRID 24 x 9 top, 24 x 9 bottom, 24 x 12 side:
    720 numbers) -- what an adjoint polish can use, since its cost does not grow with the number of numbers;
  * coarse: one per number of stage 2's CMA vector (freeform.from_unit: 45 numbers, mm per unit amp_mm) --
    the chain dD/du = dD/dc . dc/du, exact (from_unit is linear in u, apart from side >= 0).
Both are for the half car as DAFoam reports it and for the whole car (x 2, the shape is symmetric).

Units: DAFoam's file is SI (Xs in m, dD/dXs in N/m, half car). Gradients here: mN per mm of control point
(fine) and mN per unit of CMA's number (coarse).

    python link.py <dafoam_sens.npz> <job.json> <out.json>
job.json: {"freeform": <Freeform dict>, "coarse": [5,2,3,4], "amp_mm": [6,4,6], "v": [CMA vector, optional],
           "W": .., "x_front": .., "rear": .., "skin_mm": 0, "tol_mm": 1.0}

check_volume() is the exact test of this chain without CFD: with dD/dX_i replaced by the area vector of the
surface at point i, sum_i A_i n_i . dX_i is the volume change, so the "gradient" must equal the body's volume
derivative measured by rebuilding the body (test_link.py).
"""
from __future__ import annotations

import sys as _sys, pathlib as _pl  # the repo root, then every code folder (repo_paths.py)
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))
import repo_paths  # noqa: E402,F401
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

import freeform as fr

TOL_MM = 1.0             # a point within this of the free-form level set belongs to it (build grid 1 mm, milled)
EPS_MM = 1e-3            # finite step for the scalar partials of phi (phi is smooth in the surface values)


# ------------------------------------------------------------------ the body's field, as param_body draws it
@dataclasses.dataclass
class Frame:
    """Where the milled body starts and ends, and the axles (car frame, mm): what param_body.build adds to
    freeform.field."""
    x_front_cut: float              # Ref A, or a foam nose's front limit
    x_rear: float                   # rear face
    x_ref_a: float
    axles: tuple | None             # (x_front_axle, x_rear_axle, None): the 27.175 mm sidepod-top cap between them
    skin_mm: float = 0.0            # the build's fixed skin offset: the surface is phi = skin

    @classmethod
    def from_layout(cls, W, x_front, rear, nose_foam=False, skin_mm=0.0):
        ref_a = x_front - 16.0
        x_rear = x_front + W + rear
        cut = ref_a
        if nose_foam:
            import param_body as pb
            cut = max(ref_a - pb.FOAM_NOSE_OVERHANG_MM, x_rear - pb.BLOCK_LENGTH_MM + pb.BLOCK_MARGIN_MM)
        return cls(x_front_cut=cut, x_rear=x_rear, x_ref_a=ref_a, axles=(x_front, x_front + W, None), skin_mm=skin_mm)


def _phi_from(top, bot, side, P, ff, frame):
    """freeform.field from given surface values (so its partials in them can be taken), plus param_body's cuts."""
    X, Y, Z = P[:, 0], P[:, 1], P[:, 2]
    k = ff.blend_mm
    phi = fr._smax(fr._smax(Z - top, bot - Z, k), np.abs(Y) - side, k)
    phi = np.maximum(phi, np.maximum(ff.x_tip_mm - X, X - ff.x_rear_mm))
    if frame.axles is not None:
        xf, xr, _ = frame.axles
        between = (X > xf) & (X < xr)
        cap = np.minimum(Z - fr.SIDEPOD_TOP_MM, np.abs(Y) - fr.SPINE_HALF_MM)
        phi = np.where(between, np.maximum(phi, cap), phi)
    # the milled body's ends; ahead of Ref A (printed nose) the nose is the same field, so the front cut applies
    # only to the milled part: a point ahead of the cut is the printed nose's (level 0, no skin)
    milled = np.maximum(phi, np.maximum(frame.x_front_cut - X, X - frame.x_rear))
    return np.where(X >= frame.x_front_cut, milled, np.maximum(phi, X - frame.x_rear))


def level(P, frame):
    """The level the surface sits at: the skin on the milled body, 0 on a printed nose (built at phi = 0)."""
    return np.where(P[:, 0] >= frame.x_front_cut, frame.skin_mm, 0.0)


def phi_and_grad(ff, P, frame, h=0.05):
    top, bot, side = fr.surfaces(ff, P[:, 0], P[:, 1], P[:, 2], ff.x_tip_mm, ff.x_rear_mm)
    phi = _phi_from(top, bot, side, P, ff, frame)
    G = np.zeros_like(P)
    for a in range(3):
        d = np.zeros(3); d[a] = h
        Pp, Pm = P + d, P - d
        fp = _phi_from(*fr.surfaces(ff, Pp[:, 0], Pp[:, 1], Pp[:, 2], ff.x_tip_mm, ff.x_rear_mm), Pp, ff, frame)
        fm = _phi_from(*fr.surfaces(ff, Pm[:, 0], Pm[:, 1], Pm[:, 2], ff.x_tip_mm, ff.x_rear_mm), Pm, ff, frame)
        G[:, a] = (fp - fm) / (2 * h)
    return phi, G, (top, bot, side)


def surface_partials(ff, P, frame, surf):
    """dphi/dtop, dphi/dbottom, dphi/dside at each point (central differences in the surface values)."""
    top, bot, side = surf
    out = []
    for j in range(3):
        a = [top, bot, side]
        ap = list(a); am = list(a)
        ap[j] = a[j] + EPS_MM; am[j] = a[j] - EPS_MM
        out.append((_phi_from(*ap, P, ff, frame) - _phi_from(*am, P, ff, frame)) / (2 * EPS_MM))
    return out


def basis(ff, P):
    """The cubic-spline weight of every control point at every point: (N, nu*nv) for top and bottom (same
    coordinates), (N, nu*nw) for the side. _interp is linear in the grid, so the weight of control point j is
    _interp of the unit grid e_j."""
    u = (P[:, 0] - ff.x_tip_mm) / (ff.x_rear_mm - ff.x_tip_mm)
    v = np.abs(P[:, 1]) / fr.Y_MAX_MM
    w = P[:, 2] / fr.Z_MAX_MM
    def mat(n0, n1, a, b):
        M = np.zeros((len(P), n0 * n1))
        for j in range(n0 * n1):
            e = np.zeros(n0 * n1); e[j] = 1.0
            M[:, j] = fr._interp(e.reshape(n0, n1), a, b)
        return M
    return mat(ff.nu, ff.nv, u, v), mat(ff.nu, ff.nw, u, w)


# ------------------------------------------------------------------ the link
def fine_weights(ff, P_mm, frame, tol_mm=TOL_MM):
    """Per point: active mask, phi, grad phi and dphi/dc for every fine control point (top | bottom | side)."""
    phi, G, surf = phi_and_grad(ff, P_mm, frame)
    off = phi - level(P_mm, frame)
    g2 = (G * G).sum(1)
    act = (np.abs(off) < tol_mm) & (g2 > 0.04)
    dT, dB, dS = surface_partials(ff, P_mm[act], frame, tuple(s[act] for s in surf))
    Buv, Buw = basis(ff, P_mm[act])
    J = np.hstack([dT[:, None] * Buv, dB[:, None] * Buv, dS[:, None] * Buw])      # (n_active, 720) dphi/dc
    return {"active": act, "phi": phi, "off": off, "grad": G, "J": J}


def normal_weights(dD_dX, G):
    """w_i = -(dD/dX_i . grad phi_i) / |grad phi_i|^2 : the drag change per unit of phi's change at point i."""
    return -(dD_dX * G).sum(1) / np.maximum((G * G).sum(1), 1e-12)


def coarse_jacobian(seed, v, amp, coarse, h=1e-3):
    """dc/du (720 x 45): from_unit is linear in u (exact by central differences) except where the side is
    clamped at 0 (there the derivative is 0, as from_unit then is)."""
    n = fr.n_genes(*coarse)
    v = np.zeros(n) if v is None else np.asarray(v, float)[:n]
    def c(u):
        f = fr.from_unit(seed, u, amp, *coarse)
        return np.concatenate([np.asarray(f.top), np.asarray(f.bottom), np.asarray(f.side)])
    Jc = np.zeros((len(c(v)), n))
    for k in range(n):
        e = np.zeros(n); e[k] = h
        Jc[:, k] = (c(v + e) - c(v - e)) / (2 * h)
    return Jc


def gradient(dD_dX_SI, X_SI, ff, frame, seed=None, v=None, amp=None, coarse=None, tol_mm=TOL_MM):
    """DAFoam's (X_i [m], dD/dX_i [N/m], half car) -> gradients. ff: the body that was solved (Freeform);
    seed/v/amp/coarse: stage 2's parametrisation (seed + from_unit(v)), for the coarse gradient."""
    P = np.asarray(X_SI, float) * 1e3
    g = np.asarray(dD_dX_SI, float)                     # N/m = mN/mm
    fw = fine_weights(ff, P, frame, tol_mm)
    act = fw["active"]
    w = normal_weights(g[act], fw["grad"][act])          # mN per mm of phi
    fine_half = w @ fw["J"]                               # mN per mm of control point, half car
    out = {"n_points": int(len(P)), "n_active": int(act.sum()), "tol_mm": tol_mm,
           "active_share_of_dD": float(np.abs((g[act] * fw["grad"][act]).sum(1)).sum() /
                                       max(np.abs((g * np.nan_to_num(fw["grad"])).sum(1)).sum(), 1e-30)),
           "fine_half_mN_per_mm": fine_half.tolist(), "fine_car_mN_per_mm": (2 * fine_half).tolist(),
           "fine_layout": {"top": [ff.nu, ff.nv], "bottom": [ff.nu, ff.nv], "side": [ff.nu, ff.nw]}}
    if seed is not None:
        Jc = coarse_jacobian(seed, v, amp, coarse)
        cg = fine_half @ Jc
        out.update(coarse_half_mN_per_unit=cg.tolist(), coarse_car_mN_per_unit=(2 * cg).tolist(),
                   coarse=list(coarse), amp_mm=list(amp))
    return out


# ------------------------------------------------------------------ the exact check (no CFD)
def area_vectors(mesh):
    """Per vertex: the area vector A_i n_i (one third of each adjacent face's area vector). sum_i A_i.dX_i is
    the volume change of a closed mesh moved by dX (first order)."""
    F = mesh.faces
    fa = 0.5 * np.cross(mesh.vertices[F[:, 1]] - mesh.vertices[F[:, 0]], mesh.vertices[F[:, 2]] - mesh.vertices[F[:, 0]])
    A = np.zeros_like(mesh.vertices)
    for j in range(3):
        np.add.at(A, F[:, j], fa / 3.0)
    return A


def predicted_moment(mesh_mm, ff, frame, f=None, seed=None, v=None, amp=None, coarse=None, tol_mm=TOL_MM):
    """d/dc of  integral_V f dV  (f = 1: the volume), through the same chain as the drag: 'dD/dX_i' = f(x_i)
    A_i n_i (mm^3 per mm). Returns (fine, coarse or None, info)."""
    P = np.asarray(mesh_mm.vertices, float)
    A = area_vectors(mesh_mm) * (np.ones(len(P)) if f is None else f(P))[:, None]
    fw = fine_weights(ff, P, frame, tol_mm)
    act = fw["active"]
    w = normal_weights(A[act], fw["grad"][act])
    fine = w @ fw["J"]
    co = fine @ coarse_jacobian(seed, v, amp, coarse) if seed is not None else None
    return fine, co, {"n_active": int(act.sum()), "n_points": int(len(P)),
                      "active_area_share": float(np.linalg.norm(A[act], axis=1).sum() / np.linalg.norm(A, axis=1).sum())}


def main(argv):
    sens, job, outp = argv[1], json.loads(Path(argv[2]).read_text()), argv[3]
    z = np.load(sens)
    seed = fr.Freeform.from_dict(job["freeform_seed"]) if job.get("freeform_seed") else None
    ff = fr.Freeform.from_dict(job["freeform"])
    frame = Frame.from_layout(job["W"], job["x_front"], job["rear"], job["freeform"].get("nose_foam", False),
                              job.get("skin_mm", 0.0))
    res = gradient(z["dD_dXs"], z["Xs"], ff, frame, seed=seed, v=job.get("v"), amp=job.get("amp_mm"),
                   coarse=job.get("coarse"), tol_mm=job.get("tol_mm", TOL_MM))
    res["D_half_N"] = float(z["D_N"]) if "D_N" in z else None
    Path(outp).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("n_points", "n_active", "active_share_of_dD", "D_half_N")}))


if __name__ == "__main__":
    main(sys.argv)
