"""Runner-side geometry of one lineage's seed (seedhybrid/README.md).

    python seedhybrid/prep.py geom <adjoint repo checkout> <lineage json> <out dir>

Builds the seed as a legal car folder, puts it at gate/original (the meshed surface's source) and its even-triangle
remesh at gate/original_rm (the surface every move is defined on and every mesh is morphed from), and packs both into
<out dir>/geom.tgz. Seed kinds:
- {"kind": "seedevo", "campaign": "seedevo/campaign_poc", "id": 0}: a design of a seed-evolution campaign, rebuilt
  from its stored field through the rules (3-axis mill, 3 mm walls, one piece, cargo, clearance): legal by construction.
- {"kind": "cad"}: the team's exact CAD (not legalised: the original hybrid's start).
"""
import json, os, shutil, subprocess, sys, tarfile
from pathlib import Path


def geom(adj, lin, out):
    sp = Path(adj).resolve() / "research" / "sidepod_opt"
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    L = json.loads(lin) if lin.strip().startswith("{") else json.loads(Path(lin).read_text())
    seed = L["seed"]
    g = sp / "gate"; g.mkdir(exist_ok=True)
    for d in ("original", "original_rm"):
        shutil.rmtree(g / d, ignore_errors=True)
    py = sys.executable
    if seed["kind"] == "seedevo":
        tmp = out / "build"
        r = subprocess.run([py, "-W", "ignore", "-m", "seedevo.campaign", "build", seed["campaign"], str(seed["id"]), str(tmp)],
                           cwd=sp, capture_output=True, text=True)
        (out / "build.log").write_text(r.stdout[-5000:] + r.stderr[-5000:])
        if r.returncode != 0:
            raise SystemExit("seed build failed: " + (r.stdout + r.stderr)[-500:])
        shutil.copytree(tmp / "car", g / "original")
        shutil.rmtree(tmp)
    elif seed["kind"] == "cad":
        r = subprocess.run([py, "-W", "ignore", "-c", "from pathlib import Path; from geometry import Car\n"
                            "c = Car(); c.export_car_folder(c.body_mesh, Path('gate/original'), 'original')"],
                           cwd=sp, capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit("CAD export failed: " + (r.stdout + r.stderr)[-500:])
    else:
        raise SystemExit(f"unknown seed kind {seed['kind']}")
    # gate.cmd_remesh reuses gate_results/adjoint_*.npy (the TEAM pod's surface) when present: move them aside so a
    # seed is remeshed from its own surface (pymeshlab isotropic remesh, ~1 mm triangles, 30 deg edges kept)
    hide = sp / "gate_results" / "_hidden_adjoint"
    if seed["kind"] != "cad":
        hide.mkdir(exist_ok=True)
        for f in (sp / "gate_results").glob("adjoint_*.npy"):
            f.rename(hide / f.name)
    r = subprocess.run([py, "-W", "ignore", "gate.py", "remesh"], cwd=sp, capture_output=True, text=True)
    (out / "remesh.log").write_text(r.stdout[-3000:] + r.stderr[-3000:])
    if r.returncode != 0:
        raise SystemExit("remesh failed: " + (r.stdout + r.stderr)[-500:])
    # the symmetry face must sit exactly at y = 0 (make_step holds |y| < 1e-9 there)
    import numpy as np, trimesh
    sys.path.insert(0, str(sp))
    import cadio as CI
    m = trimesh.load(g / "original_rm" / "car_half.stl", process=False)
    m.merge_vertices()
    v = m.vertices.copy()
    near = np.abs(v[:, 1]) < 2e-6
    v[near, 1] = 0.0
    m.vertices = v
    CI.write_ascii_named(m, g / "original_rm" / "car_half.stl", "car", m.volume, 0.002)
    info = {"seed": seed, "vertices": int(len(m.vertices)), "faces": int(len(m.faces)), "snapped_to_y0": int(near.sum()),
            "volume_mm3_half": float(m.volume * 1e9), "watertight": bool(m.is_watertight)}
    # the region taper once (cached in gate/original_rm/region_taper.npy, read by the patched make_step)
    r = subprocess.run([py, "-W", "ignore", "-c", "import os, sys; os.environ['STEP_REGION']='1'; sys.path.insert(0,'esi');"
                        "import numpy as np, trimesh, gate, make_step as MS; b = trimesh.load(gate.ORIG/'car_half.stl');"
                        "r = MS.region_taper(b.vertices, MS.vertex_area(b)); print(float((r>0.5).mean()))"],
                       cwd=sp, capture_output=True, text=True)
    info["region_share_of_vertices"] = r.stdout.strip()[-20:]
    if r.returncode != 0:
        raise SystemExit("region taper failed: " + r.stderr[-500:])
    (out / "geom.json").write_text(json.dumps(info, indent=1))
    with tarfile.open(out / "geom.tgz", "w:gz") as tf:
        tf.add(g / "original", arcname="gate/original")
        tf.add(g / "original_rm", arcname="gate/original_rm")
    print(json.dumps(info))


if __name__ == "__main__":
    if sys.argv[1] == "geom":
        geom(*sys.argv[2:5])
