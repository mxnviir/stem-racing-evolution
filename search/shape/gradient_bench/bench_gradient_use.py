"""Gradient-use bench (PLAN.md in this folder): how should a fine-mesh adjoint gradient be used in stage 2
(search/shape/simscale_shape.py: CMA-ES, population 20, over the 45 coarse numbers of car/body/freeform.from_unit)?
No CFD: a known-answer 45-dimensional bench, numpy only.

    python bench_gradient_use.py run <out.json> --instances 41 42 43 44 45 --seeds 1 2 3 4 5 [--settings ...]
    python bench_gradient_use.py verdict <out.json>          # tables and the registered marks

THE BENCH (not a model of aerodynamics; a landscape with the structure the real one is believed to have).
x in [-2, 2]^45 = CMA's numbers (top 5x2, bottom 5x4, side 5x3 coarse displacement grids, as COARSE = (5, 2, 3, 4)).
x = 0 is the stage-1 car. Units are "mN".

  TRUE ("production") drag   f(x) = 1/2 (x - x*)' H (x - x*) + sum_j a_j (1 - cos(w_j . (x - x*)))   (f* = 0)
    * H = V diag(lam) V': V = the three grids' orthonormal 2-D DCT modes, mixed across grids by a random rotation
      near the identity (QR of I + 0.35 G / sqrt(45)); log10 lam = 2.5 (1 - k / k_max) + N(0, 0.25^2) - per-grid
      offset U(-0.3, 0.3), k = the mode's spatial frequency: smooth modes matter most (condition ~2000). lam is scaled
      so an isotropic sigma0 = 0.5 move costs SPREAD = 12 mN on average (1/2 sigma0^2 tr H);
    * ripple: 8 terms, w_j = 4 x a random unit direction in the 12 smoothest modes, a_j = RIPPLE D / 8 (RIPPLE 0.10):
      mildly nonconvex, smooth;
    * x* = V a, a_i ~ N(0, 1) lam_i^-0.25 (the seed is nearer optimal in the stiff modes), scaled so f(0) - f* = D =
      SCALE_D = 12 mN, then clipped to +-1.8 (inside the box, so f* = 0 is the optimum; f(0) a little under D then).
  STEADY MODEL f_s (what the adjoint differentiates; the steady continuation is not the production model):
      f_s(x) = 1/2 (x - x*_s)' H_s (x - x*_s) + the same ripple about x*_s; H_s = V_s diag(lam exp(N(0, 0.15^2))) V_s',
      V_s = V rotated by QR(I + 0.10 G / sqrt(45)); x*_s = x* + delta with 1/2 delta' H delta = 0.02 D.
  ADJOINT at x (one call), mode "exact" (primary): unit direction u = COS g^ + sqrt(1 - COS^2) (e_m + xi) / sqrt(2):
      g^ = grad f / |grad f|, e_m = the steady model's error direction at x (grad f_s's part orthogonal to g^, unit:
      systematic, the same for nearby calls), xi = a fresh random unit vector orthogonal to both (per call). So the
      cosine with the TRUE gradient is exactly COS (primary 0.75; rows 0.4, 0.9). Mode "model" (row `cosmodel`):
      grad f_s's direction plus random error up to cosine COS where the model is better than that, else grad f_s
      alone (its cosine falls near the optimum, where the model's error dominates). Magnitude (both): KAPPA x COS x
      |grad f|, i.e. the linear gain it predicts along its own direction is KAPPA (primary 2.5) times the true linear
      gain (the seed-evolution record: predicted -7.5 mN, measured -2.6).
  FORWARD run (one CFD): f(clip(x)) + N(0, sd), sd = NOISE x G_TYP, G_TYP = 5 mN (primary NOISE 0.40: 2.1 mN against
      gains of a few mN; reported 0.20).
  COST: a forward run = 1 unit; an adjoint = A units (primary 4; reported 2 and 8). ROUNDS (wall clock): a batch of
      up to 20 forward runs in parallel = 1 round; an adjoint = 1 round (forwards that do not need it may share it).

THE ARMS (all: CMA-ES, numpy implementation of the standard (mu/mu_w, 20) CMA-ES with Hansen's default constants,
x0 = 0, sigma0 = 0.5, box [-2, 2] handled as stage 2 does in effect: evaluated at clip(x), told f + 10 |x - clip(x)|^2).
Every arm stops its search phase when the next step would cut into the RESERVE (>= 10 % of the budget); what is left
is the hand-over re-run (same for all, below). "Direction" d of a gradient g at the CMA state (m, sigma, C) is the
CMA-preconditioned descent direction -C g, scaled to unit length in CMA's metric (|C^-1/2 d| = 1).
  C   plain CMA (control): generations of 20 forward runs.
  P   trust-region polish: CMA for its generations while units + 20 <= 0.65 b; then from CMA's mean x0 (2 runs + an
      adjoint, one round): step x1 = clip(x0 + Delta sigma d), Delta0 = 2, one forward run (one round). ACCEPT if
      y1 < ybar0 - 1.0 sd sqrt(1 + 1/k0) (k0 = runs of x0): then x1 is re-run once and a new adjoint is taken at it
      (one round), Delta doubles if measured/predicted gain >= 0.75, else stays. REJECT: Delta halves (same adjoint).
      After 4 rejections in a row the polish stops and the remaining units go back to CMA generations.
  G   gradient-guided CMA: every generation an adjoint at CMA's mean (one round, A units), then a batch of 20: 16 CMA
      samples + 4 injected points m + sigma t d on the adjoint's line, t in (0.5, 1, 2, 4) (one round); all 20
      forward-run and told to CMA (injected steps clipped to Mahalanobis length sqrt(n) + 2n/(n + 2), Hansen 2011
      "Injecting external solutions into CMA-ES"), so CMA's own selection verifies the gradient every generation.
  GP  "both" (reported): G while units + 20 + A <= 0.65 b, then P's polish from G's mean.
  G2  (reported): G with an adjoint only every other generation (guided and plain generations alternate).
  Glag (reported): G with the adjoint at the CURRENT mean taken in the same round as the batch and used for the NEXT
      generation (one round per generation instead of two; a one-generation-old gradient; generation 1 unguided).

HAND-OVER (all arms, the same rule): the reserve (the units the arm did not spend, >= 10 % of b) re-runs 5
candidates: the arm's current point (CMA's mean for C/G/Glag; the polished design for P/GP) and the 4 lowest-mean
other evaluated designs, round-robin in that order. The arm hands over the candidate with the lowest mean of all its
runs. SCORE = TRUE regret f(handed) - f*. Reported too: the regret of the current point, rounds, adjoints, cosines.
Each (arm, budget) is its own run; the random streams (CMA sampling, noise, adjoint error) are keyed by
(instance, seed), so P's CMA phase is exactly C's first generations (paired comparisons).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import math
import sys
import time
from multiprocessing import Pool

import numpy as np

GRIDS = ((5, 2), (5, 4), (5, 3))   # top, bottom, side: COARSE = (5, 2, 3, 4) of simscale_shape.py
N = sum(a * b for a, b in GRIDS)   # 45
SCALE_D = 12.0                     # f(0) - f* (mN); frozen in PLAN.md
GAMMA = 0.25                       # x*'s mode amplitudes ~ lam^-GAMMA
XS_MAX = 1.8                       # x* clipped to +-XS_MAX (inside the box, so f* = 0 is the optimum)
SPREAD = 12.0                      # mean drag cost (mN) of an isotropic sigma0 = 0.5 move (0.5 sigma0^2 tr H)
RIPPLE = 0.10                      # ripple amplitude in total, share of D
MODEL_OFF = 0.02                   # f(x*_s) - f* (quadratic part) as a share of D
G_TYP = 5.0                        # typical realised gain (mN) the noise ratio refers to
POP = 20
SIGMA0 = 0.5
BOX = 2.0
RESERVE = 0.10
POLISH_FROM = 0.65                 # P / GP: search phase while units + next step <= 0.65 b
FAN = (0.5, 1.0, 2.0, 4.0)

PRIMARY = {"cos": 0.75, "A": 4, "noise": 0.40, "kappa": 2.5}
SETTINGS = {
    "primary": PRIMARY,
    "cos0.4": dict(PRIMARY, cos=0.4),
    "cos0.9": dict(PRIMARY, cos=0.9),
    "A2": dict(PRIMARY, A=2),
    "A8": dict(PRIMARY, A=8),
    "noise0.2": dict(PRIMARY, noise=0.20),
    "kappa1": dict(PRIMARY, kappa=1.0),
    "cosmodel": dict(PRIMARY, mode="model"),
}
ARMS = ("C", "P", "G", "GP", "Glag", "G2")
BUDGETS = (100, 200, 400)


# ------------------------------------------------------------------ the landscape
def _dct(n):
    k = np.arange(n)[:, None]; i = np.arange(n)[None, :]
    M = np.cos(np.pi * (i + 0.5) * k / n) * np.sqrt(2.0 / n)
    M[0] /= np.sqrt(2.0)
    return M                       # rows = orthonormal modes


def _modes():
    """Columns: the orthonormal 2-D DCT modes of the three grids (block diagonal); freq: each mode's frequency."""
    V = np.zeros((N, N)); freq = np.zeros(N); grid = np.zeros(N, int); c = 0; r = 0
    for gi, (a, b) in enumerate(GRIDS):
        Da, Db = _dct(a), _dct(b)
        for ka in range(a):
            for kb in range(b):
                V[r:r + a * b, c] = np.outer(Da[ka], Db[kb]).ravel()
                freq[c] = math.hypot(ka / a, kb / b); grid[c] = gi; c += 1
        r += a * b
    return V, freq, grid


def _near_rot(rng, eps):
    Q, R = np.linalg.qr(np.eye(N) + eps * rng.standard_normal((N, N)) / np.sqrt(N))
    return Q * np.sign(np.diag(R))


class Problem:
    def __init__(self, inst):
        rng = np.random.default_rng([7919, inst])
        B, freq, grid = _modes()
        V = B @ _near_rot(rng, 0.35)
        off = rng.uniform(-0.3, 0.3, 3)
        lam = 10 ** (2.5 * (1 - freq / freq.max()) + rng.normal(0, 0.25, N) - off[grid])
        lam *= SPREAD / (0.5 * SIGMA0 ** 2 * lam.sum())        # a sigma0 sample costs SPREAD on average
        a0 = rng.normal(0, 1.0, N) * lam ** -GAMMA            # x* in modes: stiff modes nearer the seed
        xs = V @ a0
        # ripple along smooth directions
        smooth = np.argsort(freq)[:12]
        W = []
        for _ in range(8):
            u = B[:, smooth] @ rng.standard_normal(12); W.append(4.0 * u / np.linalg.norm(u))
        self.W = np.array(W)
        a = RIPPLE * SCALE_D / 8
        fq = lambda t: 0.5 * t * t * (xs @ (V * lam) @ V.T @ xs) + a * np.sum(1 - np.cos(self.W @ (t * xs)))
        lo, hi = 0.0, 50.0                                      # scale x* so that f(0) - f* = D
        for _ in range(100):
            mid = (lo + hi) / 2; lo, hi = (mid, hi) if fq(mid) < SCALE_D else (lo, mid)
        xs = np.clip(xs * (lo + hi) / 2, -XS_MAX, XS_MAX)
        self.lam = lam; self.V = V; self.xs = xs; self.a = a
        self.H = (V * self.lam) @ V.T
        self.fstar = 0.0
        # steady model: Hessian rotated and re-weighted, its own optimum MODEL_OFF D away (in truth), same ripple
        Vs = V @ _near_rot(rng, 0.10)
        lam_s = self.lam * np.exp(rng.normal(0, 0.15, N))
        self.Hs = (Vs * lam_s) @ Vs.T
        dl = rng.standard_normal(N); dl *= math.sqrt(2 * MODEL_OFF * SCALE_D / (dl @ self.H @ dl))
        self.xs_s = xs + dl

    def f(self, x):
        d = np.clip(x, -BOX, BOX) - self.xs
        return float(0.5 * d @ self.H @ d + self.a * np.sum(1 - np.cos(self.W @ d)))

    def grad(self, x):
        d = x - self.xs
        return self.H @ d + self.a * (self.W.T @ np.sin(self.W @ d))

    def grad_s(self, x):
        d = x - self.xs_s
        return self.Hs @ d + self.a * (self.W.T @ np.sin(self.W @ d))


def adjoint(prob, x, cos, kappa, rng, mode="exact"):
    """One adjoint call: (gradient as reported, its true cosine with grad f).
    exact: cosine with grad f exactly `cos`; the error (the part orthogonal to grad f) is half systematic (the steady
           model's own error direction at x) and half random per call (in variance).
    model: grad f_s plus random per-call error up to cosine `cos` when the model is better than that, else grad f_s
           alone (its cosine falls near the optimum, where the model's error dominates the true gradient)."""
    gt = prob.grad(x); gs = prob.grad_s(x)
    ut = gt / np.linalg.norm(gt); us = gs / np.linalg.norm(gs)
    if mode == "exact":
        em = us - (us @ ut) * ut; em /= np.linalg.norm(em)
        Q, _ = np.linalg.qr(np.stack([ut, em], 1))
        xi = rng.standard_normal(N); xi -= Q @ (Q.T @ xi); xi /= np.linalg.norm(xi)
        u = cos * ut + math.sqrt(1 - cos * cos) * (em + xi) / math.sqrt(2)
    else:
        cs = float(us @ ut)
        u = us.copy()
        if cs > cos:
            Q, _ = np.linalg.qr(np.stack([ut, us], 1))
            xi = rng.standard_normal(N); xi -= Q @ (Q.T @ xi); xi /= np.linalg.norm(xi)
            u = us + math.sqrt((cs / cos) ** 2 - 1) * xi
    u /= np.linalg.norm(u)
    return kappa * cos * np.linalg.norm(gt) * u, float(u @ ut)


# ------------------------------------------------------------------ CMA-ES (standard, Hansen's defaults)
class CMA:
    def __init__(self, n, sigma, rng, lam=POP):
        self.n, self.lam, self.rng = n, lam, rng
        self.mu = lam // 2
        w = np.log(self.mu + 0.5) - np.log(np.arange(1, self.mu + 1))
        self.w = w / w.sum(); self.mueff = 1 / np.sum(self.w ** 2)
        me = self.mueff
        self.cs = (me + 2) / (n + me + 5)
        self.ds = 1 + 2 * max(0, math.sqrt((me - 1) / (n + 1)) - 1) + self.cs
        self.cc = (4 + me / n) / (n + 4 + 2 * me / n)
        self.c1 = 2 / ((n + 1.3) ** 2 + me)
        self.cmu = min(1 - self.c1, 2 * (me - 2 + 1 / me) / ((n + 2) ** 2 + me))
        self.chin = math.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n * n))
        self.m = np.zeros(n); self.sigma = sigma
        self.C = np.eye(n); self.B = np.eye(n); self.D = np.ones(n)
        self.pc = np.zeros(n); self.ps = np.zeros(n); self.g = 0

    def sample(self, k):
        z = self.rng.standard_normal((k, self.n))
        return self.m + self.sigma * (z * self.D) @ self.B.T

    def direction(self, grad):
        """-C grad, unit length in CMA's metric."""
        d = -self.C @ grad
        return d / math.sqrt(grad @ self.C @ grad)

    def invsqrt(self, y):
        return self.B @ ((self.B.T @ y) / self.D)

    def tell(self, xs, fs, injected=None):
        xs = np.asarray(xs); n = self.n
        Y = (xs - self.m) / self.sigma
        if injected is not None:
            cap = math.sqrt(n) + 2 * n / (n + 2)
            for i in np.flatnonzero(injected):
                L = np.linalg.norm(self.invsqrt(Y[i]))
                if L > cap:
                    Y[i] *= cap / L
        order = np.argsort(fs)[:self.mu]
        Ys = Y[order]; yw = self.w @ Ys
        self.m = self.m + self.sigma * yw
        self.ps = (1 - self.cs) * self.ps + math.sqrt(self.cs * (2 - self.cs) * self.mueff) * self.invsqrt(yw)
        self.g += 1
        hs = np.linalg.norm(self.ps) / math.sqrt(1 - (1 - self.cs) ** (2 * self.g)) / self.chin < 1.4 + 2 / (n + 1)
        self.pc = (1 - self.cc) * self.pc + hs * math.sqrt(self.cc * (2 - self.cc) * self.mueff) * yw
        dh = (1 - hs) * self.cc * (2 - self.cc)
        self.C = ((1 - self.c1 - self.cmu + self.c1 * dh) * self.C + self.c1 * np.outer(self.pc, self.pc)
                  + self.cmu * (Ys.T * self.w) @ Ys)
        self.sigma *= math.exp((self.cs / self.ds) * (np.linalg.norm(self.ps) / self.chin - 1))
        self.C = (self.C + self.C.T) / 2
        D2, self.B = np.linalg.eigh(self.C)
        self.D = np.sqrt(np.maximum(D2, 1e-20))


# ------------------------------------------------------------------ one run of one arm
class Run:
    def __init__(self, prob, st, inst, seed, budget):
        self.p, self.st, self.b = prob, st, budget
        self.sd = st["noise"] * G_TYP
        self.r_cma = np.random.default_rng([inst, seed, 1])
        self.r_noise = np.random.default_rng([inst, seed, 2])
        self.r_adj = np.random.default_rng([inst, seed, 3])
        self.units = 0; self.rounds = 0; self.n_adj = 0
        self.designs = []           # [x (clipped), [measurements]]
        self.cos = []; self.accepts = 0; self.steps = 0
        self.limit = budget - math.ceil(RESERVE * budget)

    def forward(self, xs):
        """A batch of forward runs (rounds counted by the caller). Returns measurements and design indices."""
        ys, ids = [], []
        for x in xs:
            xc = np.clip(x, -BOX, BOX)
            y = self.p.f(xc) + self.sd * self.r_noise.standard_normal()
            self.designs.append([xc, [y]]); ids.append(len(self.designs) - 1); ys.append(y)
        self.units += len(xs)
        return np.array(ys), ids

    def rerun(self, i, k=1):
        for _ in range(k):
            self.designs[i][1].append(self.p.f(self.designs[i][0]) + self.sd * self.r_noise.standard_normal())
        self.units += k

    def adjoint(self, x):
        self.units += self.st["A"]; self.n_adj += 1
        g, c = adjoint(self.p, np.clip(x, -BOX, BOX), self.st["cos"], self.st["kappa"], self.r_adj,
                       self.st.get("mode", "exact"))
        self.cos.append(c)
        return g

    @staticmethod
    def told(x, ys):
        pen = 10.0 * np.sum((np.asarray(x) - np.clip(x, -BOX, BOX)) ** 2, axis=1)
        return ys + pen

    # -- phases
    def cma_gen(self, es):
        xs = es.sample(POP)
        ys, _ = self.forward(xs); self.rounds += 1
        es.tell(xs, self.told(xs, ys))

    def guided_gen(self, es, g):
        k = len(FAN)
        xs = es.sample(POP - k)
        d = es.direction(g)
        inj = np.array([es.m + es.sigma * t * d for t in FAN])
        X = np.vstack([xs, inj])
        ys, _ = self.forward(X); self.rounds += 1
        es.tell(X, self.told(X, ys), injected=np.r_[np.zeros(POP - k, bool), np.ones(k, bool)])

    def polish(self, es, stop):
        """P's trust-region polish from CMA's mean until `stop` units. Returns the current design index."""
        A = self.st["A"]
        if self.units + A + 3 > stop:
            return None
        x0 = np.clip(es.m, -BOX, BOX)
        _, (i0,) = self.forward([x0]); self.rerun(i0)
        g = self.adjoint(x0); self.rounds += 1
        d = es.direction(g); sig = es.sigma
        Delta, rej = 2.0, 0
        while self.units + 1 <= stop:
            x1 = np.clip(self.designs[i0][0] + Delta * sig * d, -BOX, BOX)
            (y1,), (i1,) = self.forward([x1]); self.rounds += 1; self.steps += 1
            runs0 = self.designs[i0][1]; y0 = float(np.mean(runs0))
            bar = 1.0 * self.sd * math.sqrt(1 + 1 / len(runs0))
            if y1 < y0 - bar:
                self.accepts += 1
                pred = float(g @ (self.designs[i0][0] - self.designs[i1][0]))
                if pred > 0 and (y0 - y1) / pred >= 0.75:
                    Delta *= 2
                i0, rej = i1, 0
                if self.units + A + 1 <= stop:
                    self.rerun(i1); g = self.adjoint(self.designs[i1][0]); self.rounds += 1
                    d = es.direction(g)
                else:
                    break
            else:
                Delta /= 2; rej += 1
                if rej == 4:
                    break
        return i0

    def handover(self, cur_x=None, cur_i=None):
        """The reserve re-runs 5 candidates; the lowest mean is handed over. Returns (handed x, current x)."""
        if cur_x is not None:
            _, (cur_i,) = self.forward([cur_x])      # the current point's first run is a reserve run
        means = [(np.mean(r), i) for i, (_, r) in enumerate(self.designs) if i != cur_i]
        cands = ([cur_i] if cur_i is not None else []) + [i for _, i in sorted(means)[:5]]
        cands = cands[:5]
        left = self.b - self.units
        j = 0
        while left > 0:
            self.rerun(cands[j % len(cands)]); left -= 1; j += 1
        self.rounds += math.ceil((j + (cur_x is not None)) / POP) if (j or cur_x is not None) else 0
        best = min(cands, key=lambda i: np.mean(self.designs[i][1]))
        return self.designs[best][0], (self.designs[cur_i][0] if cur_i is not None else None)


def run_arm(task):
    sname, arm, inst, seed, b = task
    st = SETTINGS[sname]
    p = Problem(inst)
    R = Run(p, st, inst, seed, b)
    es = CMA(N, SIGMA0, R.r_cma)
    A = st["A"]
    t0 = time.time()
    cur_x, cur_i = None, None
    if arm == "C":
        while R.units + POP <= R.limit:
            R.cma_gen(es)
        cur_x = es.m
    elif arm in ("G", "GP"):
        stop = R.limit if arm == "G" else POLISH_FROM * b
        while R.units + POP + A <= stop:
            R.guided_gen(es, R.adjoint(es.m)); R.rounds += 1
        if arm == "GP":
            cur_i = R.polish(es, R.limit)
            while R.units + POP <= R.limit:
                R.cma_gen(es)
            if cur_i is None:
                cur_x = es.m
        else:
            cur_x = es.m
    elif arm == "G2":
        k = 0
        while R.units + POP + (A if k % 2 == 0 else 0) <= R.limit:
            if k % 2 == 0:
                R.guided_gen(es, R.adjoint(es.m)); R.rounds += 1
            else:
                R.cma_gen(es)
            k += 1
        cur_x = es.m
    elif arm == "Glag":
        g = None
        while R.units + POP + A <= R.limit:
            # the adjoint at this mean runs in the same round as this generation's batch (skipped when no later
            # generation fits the budget: nothing could use it)
            gn = R.adjoint(es.m) if R.units + 2 * (POP + A) <= R.limit else None
            if g is None:
                R.cma_gen(es)
            else:
                R.guided_gen(es, g)
            g = gn
            if gn is None:
                while R.units + POP <= R.limit:
                    R.cma_gen(es)
                break
        cur_x = es.m
    elif arm == "P":
        while R.units + POP <= POLISH_FROM * b:
            R.cma_gen(es)
        cur_i = R.polish(es, R.limit)
        while R.units + POP <= R.limit:
            R.cma_gen(es)
        if cur_i is None:
            cur_x = es.m
    else:
        raise ValueError(arm)
    hx, cx = R.handover(np.clip(cur_x, -BOX, BOX) if cur_x is not None else None, cur_i)
    assert R.units == b, (arm, R.units, b)
    return {"setting": sname, "arm": arm, "inst": inst, "seed": seed, "budget": b,
            "regret": p.f(hx) - p.fstar, "regret_current": p.f(cx) - p.fstar, "f0": p.f(np.zeros(N)) - p.fstar,
            "measured_handed": float(min(np.mean(r) for x, r in R.designs if np.array_equal(x, hx))),
            "rounds": R.rounds, "adjoints": R.n_adj, "forwards": b - A * R.n_adj,
            "cos_mean": float(np.mean(R.cos)) if R.cos else None, "steps": R.steps, "accepts": R.accepts,
            "sigma_end": float(es.sigma), "secs": round(time.time() - t0, 3)}


# ------------------------------------------------------------------ driver
def cmd_run(a):
    tasks = []
    for s in a.settings:
        for arm in a.arms:
            if arm == "C" and SETTINGS[s]["noise"] == PRIMARY["noise"] and s != "primary" and "primary" in a.settings:
                continue                          # C does not depend on cos / A / kappa: reuse primary's
            for inst in a.instances:
                for seed in a.seeds:
                    for b in a.budgets:
                        tasks.append((s, arm, inst, seed, b))
    t0 = time.time()
    with Pool(a.workers) as pool:
        rows = pool.map(run_arm, tasks, chunksize=4)
    out = {"spec": {"SCALE_D": SCALE_D, "G_TYP": G_TYP, "settings": {s: SETTINGS[s] for s in a.settings},
                    "instances": a.instances, "seeds": a.seeds, "budgets": a.budgets, "arms": a.arms},
           "rows": rows, "secs": round(time.time() - t0, 1)}
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=0)
    print(f"{len(rows)} runs in {out['secs']} s -> {a.out}")


def _getter(rows, setting):
    """(arm, budget) -> {(inst, seed): row}; C is taken from 'primary' when the setting did not re-run it."""
    def get(arm, b):
        s = setting if any(r["setting"] == setting and r["arm"] == arm for r in rows) else "primary"
        return {(r["inst"], r["seed"]): r for r in rows if r["setting"] == s and r["arm"] == arm and r["budget"] == b}
    return get


def wins(ra, rb):
    keys = sorted(set(ra) & set(rb))
    return sum(ra[k]["regret"] < rb[k]["regret"] for k in keys), len(keys)


def wilcoxon_p(d):
    """Two-sided Wilcoxon signed-rank p (normal approximation)."""
    d = np.asarray([v for v in d if v != 0]); n = len(d)
    if n == 0:
        return 1.0
    r = np.argsort(np.argsort(np.abs(d))) + 1.0
    W = r[d > 0].sum(); mu = n * (n + 1) / 4; sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z = (W - mu) / sd
    return float(math.erfc(abs(z) / math.sqrt(2)))


MARK = 18


def verdict(rows, settings=None, mark=MARK):
    out = []
    budgets = sorted({r["budget"] for r in rows})
    settings = settings or [s for s in SETTINGS if any(r["setting"] == s for r in rows)]
    decision = None
    for s in settings:
        get = _getter(rows, s)
        out.append(f"\n### setting {s}: {SETTINGS[s]}")
        out.append("| arm | " + " | ".join(f"{b}" for b in budgets) + " | adjoints | rounds | mean cos |")
        out.append("|---|" + "---|" * (len(budgets) + 3))
        for arm in ARMS:
            cells, adj, rnd, cs = [], [], [], []
            for b in budgets:
                d = get(arm, b)
                if not d:
                    cells.append("-"); continue
                v = np.array([r["regret"] for r in d.values()])
                cells.append(f"{v.mean():.2f} ({v.std(ddof=1):.2f})")
                adj.append(f"{np.mean([r['adjoints'] for r in d.values()]):.1f}")
                rnd.append(f"{np.mean([r['rounds'] for r in d.values()]):.1f}")
                c = [r["cos_mean"] for r in d.values() if r["cos_mean"] is not None]
                cs.append(f"{np.mean(c):.2f}" if c else "-")
            if any(c != "-" for c in cells):
                out.append(f"| {arm} | " + " | ".join(cells) + f" | {'/'.join(adj)} | {'/'.join(rnd)} | {'/'.join(cs)} |")
        out.append("")
        out.append("| pair (X beats Y: X wins of n) | " + " | ".join(f"{b}" for b in budgets) + " |")
        out.append("|---|" + "---|" * len(budgets))
        W = {}
        for x, y in (("P", "C"), ("G", "C"), ("G", "P"), ("GP", "C"), ("GP", "G"), ("GP", "P"), ("Glag", "G"), ("Glag", "C"), ("G2", "G"), ("G2", "C")):
            cells = []
            for b in budgets:
                ra, rb = get(x, b), get(y, b)
                if not ra or not rb:
                    cells.append("-"); continue
                w, n = wins(ra, rb)
                keys = sorted(set(ra) & set(rb))
                dd = [ra[k]["regret"] - rb[k]["regret"] for k in keys]
                W[(x, y, b)] = w; W[(y, x, b)] = n - w - sum(v == 0 for v in dd)
                cells.append(f"{w}/{n} ({np.mean(dd):+.2f} +- {np.std(dd, ddof=1) / math.sqrt(n):.2f}, p {wilcoxon_p(dd):.1g})")
            out.append(f"| {x} vs {y} | " + " | ".join(cells) + " |")
        if s == "primary":
            decision = decide(W, budgets, mark)
            out.append(f"\n**Registered decision (primary): {decision}**")
    return "\n".join(out), decision


def decide(W, budgets, mark=MARK):
    beat = lambda x, y: [b for b in budgets if W.get((x, y, b), 0) >= mark]
    gp, pg = beat("G", "P"), beat("P", "G")
    gc, pc, cp = beat("G", "C"), beat("P", "C"), beat("C", "P")
    tag = f"[G>P at {gp}, P>G at {pg}, G>C at {gc}, P>C at {pc}, C>P at {cp}]"
    if len(gp) >= 2 and not pg and len(gc) >= 2:
        return "G " + tag
    if len(pg) >= 2 and not gp and len(pc) >= 1 and not cp:
        return "P " + tag
    if gp and pg:
        return "budget-dependent " + tag
    if not cp:
        return "P (tie rule) " + tag
    return "C (tie rule; P loses to C) " + tag


def cmd_verdict(a):
    rows = []
    for p in a.files:
        rows += json.load(open(p))["rows"]
    txt, _ = verdict(rows)
    print(txt)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("out")
    r.add_argument("--instances", type=int, nargs="+", default=[40])
    r.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    r.add_argument("--budgets", type=int, nargs="+", default=list(BUDGETS))
    r.add_argument("--settings", nargs="+", default=["primary"])
    r.add_argument("--arms", nargs="+", default=list(ARMS))
    r.add_argument("--workers", type=int, default=4)
    v = sub.add_parser("verdict"); v.add_argument("files", nargs="+")
    a = ap.parse_args()
    {"run": cmd_run, "verdict": cmd_verdict}[a.cmd](a)


if __name__ == "__main__":
    main()
