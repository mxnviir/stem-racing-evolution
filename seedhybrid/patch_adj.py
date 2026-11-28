"""Runner-side patch of the pinned stem-racing-adjoint checkout for seeded hybrid lineages (seedhybrid/README.md).

    python seedhybrid/patch_adj.py <adjoint repo checkout>

One change, opt-in by env STEP_REGION=1 (default off = the original code's behaviour):
- esi/make_step.py and esi/predict_candidates.py: the move taper is multiplied by a smoothed copy of the legal sidepod
  design region (gate.region_mask: the CarSpec box minus the T7.9.2 / T7.9.3 wheel clearances, the halo zone and the
  spine). The original hybrid's taper was a box x 70-200 mm, y > 12 mm, which reaches 20 mm past the region's rear end
  (x 179 mm). Cached per base surface in gate/original_rm/region_taper.npy.
"""
import sys
from pathlib import Path

adj = Path(sys.argv[1])
sp = adj / "research" / "sidepod_opt"

HELPER = '''

def region_taper(V, A):
    """STEP_REGION=1: smoothed legal-region weight per base vertex (1 inside, 0 outside, 2 mm Gaussian blend)."""
    if os.environ.get("STEP_REGION") != "1":
        return np.ones(len(V))
    cache = gate.ORIG / "region_taper.npy"
    if cache.exists():
        r = np.load(cache)
        if len(r) == len(V):
            return r
    s = gate.smooth_on_surface(V, gate.region_mask(V), A, 2.0e-3)
    r = smoothstep((s - 0.2) / 0.6)
    np.save(cache, r)
    return r
'''

ms = sp / "esi" / "make_step.py"
t = ms.read_text()
anchor = "def smooth_normals("
assert anchor in t and "def region_taper(" not in t
t = t.replace(anchor, HELPER.lstrip("\n") + "\n\n" + anchor, 1)
old = """             * smoothstep((V[:, 1] - Y0) / TAPER))
    # opt-in, larger steps"""
assert old in t, "make_step taper not found"
t = t.replace(old, """             * smoothstep((V[:, 1] - Y0) / TAPER)) * region_taper(V, A)
    # opt-in, larger steps""", 1)
ms.write_text(t)

pc = sp / "esi" / "predict_candidates.py"
t = pc.read_text()
old = "         * MS.smoothstep((V[:, 1] - MS.Y0) / MS.TAPER))\n"
assert old in t, "predict_candidates taper not found"
t = t.replace(old, "         * MS.smoothstep((V[:, 1] - MS.Y0) / MS.TAPER)) * MS.region_taper(V, A)\n", 1)
# the team pod's gradient files are only used for informational predictions; a lineage's base surface differs, so
# load them only when present and skip predictions that cannot be mapped
t = t.replace("G = {k: np.load(p) for k, p in GRADS.items()}", "G = {k: np.load(p) for k, p in GRADS.items() if p.exists()}", 1)
pc.write_text(t)
# opt-in hooks used by the confirm stage (same as the drag audit's workflow; defaults reproduce the original exactly):
# P2_BG_DIV = background-cell divisor (12 = the scoring mesh), PRODCHECK_TURB = turbulence model
oc = adj / "part2-simulation" / "openfoam_case.py"
t = oc.read_text()
for a in ("max(lx / 12.0, 1e-4)", "max(ref_len / 12.0, 1e-4)"):
    assert t.count(a) == 1, a
    t = t.replace(a, a.replace("12.0", 'float(os.environ.get("P2_BG_DIV", "12"))'))
oc.write_text(t)
pk = sp / "prodcheck.py"
t = pk.read_text()
a = 'extra["stage_timeout_s"] = int(os.environ.get("PRODCHECK_TIMEOUT", 18000))\n'
assert a in t
t = t.replace(a, a + 'if os.environ.get("PRODCHECK_TURB", "kOmegaSST") != "kOmegaSST":\n    extra["turbulence_model"] = os.environ["PRODCHECK_TURB"]\n', 1)
pk.write_text(t)
# opt-in: GATE_DOMAIN_JSON = [[x0, y0, z0], [x1, y1, z1]] (m) replaces gate.py's fixed meshing frame (the trial car's,
# 257 mm long), for whole cars that do not fit inside it (carhybrid); unset = the original frame
gp = sp / "gate.py"
t = gp.read_text()
a = "DOMAIN = ((-0.005, 0.0, 0.0), (0.265, 0.045, 0.066))"
assert t.count(a) == 1, "gate DOMAIN"
t = t.replace(a, a + "\nif os.environ.get(\"GATE_DOMAIN_JSON\"):\n    DOMAIN = tuple(tuple(map(float, c)) for c in json.loads(os.environ[\"GATE_DOMAIN_JSON\"]))", 1)
gp.write_text(t)
print("patched", ms, pc, oc, pk, gp)
