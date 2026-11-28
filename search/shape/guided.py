"""guided.py -- gradient-guided CMA-ES for stage 2 (search/shape/gradient_bench: arm G won the registered bench
against a trust-region polish and plain CMA, 2026-11-28, RESULTS.md there).

Every generation: an adjoint at CMA's mean gives g = dD/dv (the body numbers' part; search/adjoint); then the
generation is 16 of CMA's own samples plus 4 points on the adjoint's line

    m + sigma * t * d,   t in FAN = (0.5, 1, 2, 4),   d = -C g / sqrt(g' C g)   (CMA-preconditioned, unit length
                                                                                  in CMA's metric)

injected into pycma with force=True (kept as given; `tell` trims an injected step to Mahalanobis length
sqrt(n) + 2n/(n+2), Hansen 2011), all forward-run and told to CMA, so CMA's own selection checks the gradient
every generation and its 2.5x overstatement does not matter (no step length is taken from it). A generation
without a usable gradient (the adjoint failed, or the primal was not steady) is a plain one.

Drag only (owner 2026-11-28: drag scoring for now): the direction is the drag gradient's; CMA's ranking is the race
time as always.
"""
from __future__ import annotations

import numpy as np

FAN = (0.5, 1.0, 2.0, 4.0)


def covariance(es):
    """CMA's search covariance without sigma^2 (pycma: sigma_vec scaling x C x sigma_vec)."""
    C = np.asarray(es.sm.C, float)
    sv = getattr(getattr(es, "sigma_vec", None), "scaling", 1.0)
    sv = np.broadcast_to(np.asarray(sv, float), (C.shape[0],))
    return sv[:, None] * C * sv[None, :]


def direction(es, grad):
    """-C g, unit length in CMA's metric (|C^-1/2 d| = 1). grad: the full CMA vector's gradient (zeros where no
    adjoint reaches: the wing and support genes)."""
    g = np.asarray(grad, float)
    C = covariance(es)
    q = float(g @ C @ g)
    if not np.isfinite(q) or q <= 0:
        return None
    return -(C @ g) / np.sqrt(q)


def inject_line(es, grad, fan=FAN) -> dict:
    """Put the adjoint's line into the next `ask` (before ask). Returns what was injected (for the log)."""
    d = direction(es, grad)
    if d is None:
        return {"guided": False, "why": "zero gradient"}
    m = np.asarray(es.mean, float)
    pts = [m + es.sigma * t * d for t in fan]
    es.inject(pts, force=True)
    g = np.asarray(grad, float)
    return {"guided": True, "fan": list(fan), "sigma": float(es.sigma), "n_injected": len(pts),
            # first-order drag change predicted along the line per unit t (the adjoint's own, overstated ~2.5x)
            "dD_per_t_pred": float(es.sigma * (g @ d)), "grad_norm": float(np.linalg.norm(g))}


def mean_phenotype(es):
    """CMA's mean as a design vector (inside the bounds)."""
    return np.asarray(es.result.xfavorite, float)
