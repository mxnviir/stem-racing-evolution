"""Every code folder of the repo on sys.path, so the modules import each other by name.

Modules keep flat names (the saved CMA states are pickles that refer to e.g. `methods.CMA`), and the
folders are put at the FRONT of sys.path: car/parts/wheel.py must win over pip's `wheel` package.

    import repo_paths  # noqa: F401      (after putting the repo root on sys.path)
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FOLDERS = {
    "body": ROOT / "car" / "body",          # the milled body: layout, rule masks, loft, milling, mass
    "race": ROOT / "car" / "race",          # race-time model, thrust data, mass/COM roll-up
    "parts": ROOT / "car" / "parts",        # nose, wings, wheels, halo, tether guides, joints
    "build": ROOT / "car" / "build",        # one whole car: build + referee
    "cfd": ROOT / "cfd",                    # SimScale route: air volume, surface, evaluator
    "sidepod": ROOT / "search" / "sidepod", # CMA on the team car's sidepods
    "layout": ROOT / "search" / "layout",   # stage 1: the layout numbers (Bayesian)
    "shape": ROOT / "search" / "shape",     # stage 2: free-form body (CMA)
    "adjoint": ROOT / "search" / "adjoint", # the adjoint link and the fine (AWS) adjoint job
}


def add() -> None:
    for p in reversed(list(FOLDERS.values())):
        s = str(p)
        if s in sys.path:
            sys.path.remove(s)
        sys.path.insert(0, s)
    os.environ.setdefault("PART2_PATH", str(FOLDERS["race"]))     # Part 1 looks for the race model there


add()
