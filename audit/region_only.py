"""The original hybrid's final design with every move outside the legal sidepod region removed (audit batch b3).

    python audit/region_only.py <adjoint checkout> <out car_half.stl>

hybrid/final_design/car_half.stl is gate/original_rm's surface (same triangles) plus the loop's moves. Here each vertex's
move is multiplied by gate.region_mask at the BASE vertex (1 inside the region, 0 outside): exactly the moves a
region-confined loop could have made, no smoothing, so no move inside the region is weakened.
"""
import sys
from pathlib import Path
import numpy as np


def main(adj, out):
    sp = Path(adj).resolve() / "research" / "sidepod_opt"
    sys.path.insert(0, str(sp)); sys.path.insert(0, str(sp.parents[1] / "part2-simulation"))
    import trimesh, gate, cadio as CI
    v = np.load(sp / "gate_results" / "adjoint_verts.npy"); f = np.load(sp / "gate_results" / "adjoint_faces.npy")
    fin = trimesh.load(sp / "hybrid" / "final_design" / "car_half.stl", process=False)
    V = fin.vertices.reshape(-1, 3, 3)
    assert V.shape[0] == len(f), "final design is not original_rm's triangulation"
    base = v[f]
    D = V - base                                             # per face-corner moves
    w = gate.region_mask(base.reshape(-1, 3)).reshape(-1, 3)
    new = base + D * w[..., None]
    m = trimesh.Trimesh(new.reshape(-1, 3), np.arange(len(f) * 3).reshape(-1, 3), process=True)
    kept = float(np.abs(D * w[..., None]).sum() / max(np.abs(D).sum(), 1e-12))
    CI.write_ascii_named(m, Path(out), "car", m.volume, 0.01)
    print(f"region-only design written; share of the move magnitude kept {kept:.3f}; volume {m.volume * 1e9:.0f} mm3")


if __name__ == "__main__":
    main(*sys.argv[1:3])
