"""Runner-side patch of the pinned stem-racing-adjoint checkout's DAFoam scripts for the fine adjoint
(carhybrid/FINE_ADJOINT.md; job body carhybrid/fine_adjoint.sh).

    python carhybrid/patch_daf.py <adjoint repo checkout> [--check]

Idempotent, as evo_patches/apply.py: every patch asserts that either its exact anchor text or its already-patched
text is in the file, replaces the anchor once, and leaves an already-patched file alone. --check only reports (exit 1
if a patch is neither applicable nor applied). Defaults reproduce the original behaviour except that lift is computed
too (DAF_FUNCS="D" gives the original drag-only run).

research/sidepod_opt/dafoam/runScript.py
- D1: functions selectable by env DAF_FUNCS (default "D,L"): D = force along +x, L = force along +z (both on car +
  wheelF + wheelR, scale rho, half car, N). Same daOptions otherwise.
- D2: GMRES (or fixed-point) convergence recorded per adjoint solve: DAFoamSolver.solve_linear is wrapped and, after
  each solve, reads the KSP's iteration count, converged reason and last residual norm (petsc4py; DAFoam passes the
  same KSP object to its C++ solver). One solve per function (OpenMDAO rev mode), in the order of DAF_FUNCS.
  DAF_KEEP_UNCONVERGED=1 (default): a solve that misses its tolerance (DAFoam raises AnalysisError) is recorded as
  not converged and its gradient is still saved, flagged; 0 = DAFoam's behaviour (that gradient is lost).
- D3: compute_totals once per function (the preconditioner is built on the first and reused: DAFoam rebuilds it only
  when the solution folder is renamed, i.e. once per design), so a failed L adjoint still leaves D's gradient.
  dafoam_sens.npz gets Xs, dD_dXs, D_N and, with L, dL_dXs, L_N (compressed); dafoam_sens.json gets the function
  values, the per-function adjoint times and "gmres": {function: {iterations, reason, residual_norm, converged,
  seconds, error}}. The primal-only path (DAF_PRIMAL_ONLY=1) writes D and L to dafoam_primal.json.
research/sidepod_opt/dafoam/run_docker.sh
- R1: ranks from env DAF_NP when the np argument is absent; /dev/shm size from env DAF_SHM (default 16g: 2g is
  tight for OpenMPI's shared-memory transport at 64 ranks); DAF_FUNCS passed into the container.
"""
from __future__ import annotations

import sys
from pathlib import Path

RS = "research/sidepod_opt/dafoam/runScript.py"
RD = "research/sidepod_opt/dafoam/run_docker.sh"

PATCHES = [
    dict(
        id="D1_functions_D_and_L",
        file=RS,
        old='''    "function": {
        "D": {
            "type": "force",
            "source": "patchToFace",
            "patches": ["car", "wheelF", "wheelR"],
            "directionMode": "fixedDirection",
            "direction": [1.0, 0.0, 0.0],
            "scale": RHO,
        },
    },
''',
        new='''    # carhybrid/patch_daf.py D1: D (+x) and L (+z), selectable by DAF_FUNCS (default "D,L")
    "function": {f: {"type": "force", "source": "patchToFace", "patches": ["car", "wheelF", "wheelR"],
                     "directionMode": "fixedDirection", "direction": FUNC_DIRS[f], "scale": RHO} for f in FUNCS},
''',
    ),
    dict(
        id="D1b_function_list",
        file=RS,
        old='''daOptions = {
    "solverName": "DASimpleFoam",''',
        new='''# carhybrid/patch_daf.py D1: the functions, in adjoint-solve order (drag first: its gradient is the one never lost)
FUNC_DIRS = {"D": [1.0, 0.0, 0.0], "L": [0.0, 0.0, 1.0]}
FUNCS = [f.strip() for f in os.environ.get("DAF_FUNCS", "D,L").split(",") if f.strip()]
assert FUNCS and all(f in FUNC_DIRS for f in FUNCS), f"DAF_FUNCS must be a list of {sorted(FUNC_DIRS)}, got {FUNCS}"

daOptions = {
    "solverName": "DASimpleFoam",''',
    ),
    dict(
        id="D2_gmres_record",
        file=RS,
        old='''from dafoam.mphys import DAFoamBuilder
''',
        new='''from dafoam.mphys import DAFoamBuilder

# carhybrid/patch_daf.py D2: the adjoint solver's convergence, one record per solve_linear call (one per function)
from dafoam.mphys import mphys_dafoam as _MD

ADJ_SOLVES = []
_solve_linear0 = _MD.DAFoamSolver.solve_linear


# DAF_KEEP_UNCONVERGED=1 (default): an adjoint that misses its GMRES tolerance (DAFoam raises AnalysisError after it
# has stored psi) is recorded as not converged and its gradient is still computed and saved (gate F3 reads the
# record); 0 = DAFoam's own behaviour (the function's gradient is lost)
KEEP_UNCONVERGED = os.environ.get("DAF_KEEP_UNCONVERGED", "1") == "1"


def _solve_linear(self, d_outputs, d_residuals, mode):
    t, err = time.time(), None
    try:
        return _solve_linear0(self, d_outputs, d_residuals, mode)
    except Exception as e:
        err = repr(e)[:300]
        if not (KEEP_UNCONVERGED and "Adjoint solution failed" in str(e)):
            raise
    finally:
        r = {"seconds": round(time.time() - t, 1), "error": err}
        ksp = getattr(self.DASolver, "ksp", None)
        if ksp is not None and self.DASolver.getOption("adjEqnSolMethod") == "Krylov":
            try:
                r.update(iterations=int(ksp.getIterationNumber()), reason=int(ksp.getConvergedReason()),
                         residual_norm=float(ksp.getResidualNorm()))
                r["converged"] = r["reason"] > 0 and err is None
                r["gradient_kept_unconverged"] = err is not None and KEEP_UNCONVERGED
            except Exception as e2:
                r["ksp_error"] = repr(e2)[:200]
        else:
            r["converged"] = err is None
        ADJ_SOLVES.append(r)


_MD.DAFoamSolver.solve_linear = _solve_linear
''',
    ),
    dict(
        id="D3_totals_per_function",
        file=RS,
        old='''    D = float(prob.get_val("cruise.aero_post.D")[0])
    t1 = time.time()
    if gcomm.rank == 0:
        print(f"DAFOAM primal D_half = {D:.6f} N  ({t1 - t0:.0f} s)", flush=True)
    if os.environ.get("DAF_PRIMAL_ONLY") == "1":          # boundary-condition / steadiness check, no adjoint
        if gcomm.rank == 0:
            json.dump({"D_half_N": D, "primal_s": t1 - t0, "nproc": gcomm.size, "primal_only": True},
                      open("dafoam_primal.json", "w"), indent=1)
        return
    tot = prob.compute_totals(of=["cruise.aero_post.D"], wrt=["mesh.x_aero0"], get_remote=True)
    g = np.asarray(tot[("cruise.aero_post.D", "mesh.x_aero0")]).ravel()
    t2 = time.time()
    DASolver = prob.model.cruise.coupling.solver.DASolver
    Xs = np.asarray(DASolver.getSurfaceCoordinates(DASolver.designSurfacesGroup)).ravel()
    off = gcomm.exscan(len(Xs)) or 0
    local = g[off:off + len(Xs)]
    allX = gcomm.gather(Xs, root=0)
    allG = gcomm.gather(local, root=0)
    if gcomm.rank == 0:
        X = np.concatenate(allX).reshape(-1, 3)
        G = np.concatenate(allG).reshape(-1, 3)
        np.savez("dafoam_sens.npz", Xs=X, dD_dXs=G, D_N=D, primal_s=t1 - t0, adjoint_s=t2 - t1,
                 nproc=gcomm.size, gmres_tol=GMRES_TOL)
        json.dump({"D_half_N": D, "adjoint_method": ADJ, "n_surface_points": int(len(X)), "primal_s": t1 - t0, "adjoint_s": t2 - t1,
                   "nproc": gcomm.size, "sum_dD_dXs": G.sum(0).tolist()}, open("dafoam_sens.json", "w"), indent=1)
        print(f"DAFOAM adjoint done: {len(X)} surface points, {t2 - t1:.0f} s", flush=True)
''',
        new='''    # carhybrid/patch_daf.py D3: every function's value, then one compute_totals per function
    val = {f: float(prob.get_val(f"cruise.aero_post.{f}")[0]) for f in FUNCS}
    D = val.get("D")
    t1 = time.time()
    if gcomm.rank == 0:
        print("DAFOAM primal " + "  ".join(f"{f}_half = {v:.6f} N" for f, v in val.items()) + f"  ({t1 - t0:.0f} s)",
              flush=True)
    if os.environ.get("DAF_PRIMAL_ONLY") == "1":          # boundary-condition / steadiness check, no adjoint
        if gcomm.rank == 0:
            json.dump({**{f"{f}_half_N": v for f, v in val.items()}, "primal_s": t1 - t0, "nproc": gcomm.size,
                       "primal_only": True, "functions": FUNCS}, open("dafoam_primal.json", "w"), indent=1)
        return
    DASolver = prob.model.cruise.coupling.solver.DASolver
    Xs = np.asarray(DASolver.getSurfaceCoordinates(DASolver.designSurfacesGroup)).ravel()
    off = gcomm.exscan(len(Xs)) or 0
    allX = gcomm.gather(Xs, root=0)
    grads, gmres, secs = {}, {}, {}
    for f in FUNCS:
        ta, n0 = time.time(), len(ADJ_SOLVES)
        try:
            tot = prob.compute_totals(of=[f"cruise.aero_post.{f}"], wrt=["mesh.x_aero0"], get_remote=True)
            g = np.asarray(tot[(f"cruise.aero_post.{f}", "mesh.x_aero0")]).ravel()
            allG = gcomm.gather(g[off:off + len(Xs)], root=0)
            if gcomm.rank == 0:
                grads[f] = np.concatenate(allG).reshape(-1, 3)
        except Exception as e:                            # e.g. AnalysisError: this function's adjoint failed
            if gcomm.rank == 0:
                print(f"DAFOAM adjoint {f} FAILED: {e!r}"[:400], flush=True)
            if not ADJ_SOLVES[n0:]:
                ADJ_SOLVES.append({"error": repr(e)[:300], "converged": False})
        secs[f] = time.time() - ta
        rec = ADJ_SOLVES[n0:]
        gmres[f] = rec[-1] if len(rec) == 1 else {"solves": rec, "converged": bool(rec) and all(r.get("converged") for r in rec)}
        if gcomm.rank == 0:
            print(f"DAFOAM adjoint {f}: {json.dumps(gmres[f])} ({secs[f]:.0f} s)", flush=True)
    t2 = time.time()
    if gcomm.rank == 0:
        X = np.concatenate(allX).reshape(-1, 3)
        arr = {f"d{f}_dXs": G for f, G in grads.items()}
        arr.update({f"{f}_N": v for f, v in val.items()})
        np.savez_compressed("dafoam_sens.npz", Xs=X, **arr, primal_s=t1 - t0, adjoint_s=t2 - t1,
                            nproc=gcomm.size, gmres_tol=GMRES_TOL)
        json.dump({**{f"{f}_half_N": v for f, v in val.items()}, "functions": FUNCS, "gradients": sorted(grads),
                   "adjoint_method": ADJ, "n_surface_points": int(len(X)), "primal_s": t1 - t0, "adjoint_s": t2 - t1,
                   "adjoint_s_per_function": secs, "nproc": gcomm.size, "gmres_rel_tol": GMRES_TOL, "gmres": gmres,
                   **{f"sum_d{f}_dXs": G.sum(0).tolist() for f, G in grads.items()}},
                  open("dafoam_sens.json", "w"), indent=1)
        print(f"DAFOAM adjoint done: {len(X)} surface points, gradients {sorted(grads)}, {t2 - t1:.0f} s", flush=True)
''',
    ),
    dict(
        id="R1_ranks_shm_env",
        file=RD,
        old='''CASE=$(realpath "$1"); NP=${2:-4}; IMAGE=${3:-${DAFOAM_IMAGE:-dafoam/opt-packages:v5.1.1}}''',
        new='''CASE=$(realpath "$1"); NP=${2:-${DAF_NP:-4}}; IMAGE=${3:-${DAFOAM_IMAGE:-dafoam/opt-packages:v5.1.1}}  # carhybrid/patch_daf.py R1''',
    ),
    dict(
        id="R1b_shm_funcs",
        file=RD,
        old='''docker run --rm -u dafoamuser --shm-size=2g \\''',
        new='''docker run --rm -u dafoamuser --shm-size=${DAF_SHM:-16g} -e DAF_FUNCS \\''',
    ),
]


def main(argv):
    adj = Path(argv[0])
    check = "--check" in argv
    bad = 0
    for p in PATCHES:
        f = adj / p["file"]
        t = f.read_text()
        if p["new"] in t:
            print(f"{p['id']}: already applied")
        elif t.count(p["old"]) == 1:
            if check:
                print(f"{p['id']}: applicable")
            else:
                f.write_text(t.replace(p["old"], p["new"], 1))
                print(f"{p['id']}: applied")
        else:
            print(f"{p['id']}: NEITHER anchor nor patched text found in {f}")
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
