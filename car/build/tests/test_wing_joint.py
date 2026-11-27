"""T9.6's fillet allowance is measured from where a support joins the wing, buried part included."""
import trimesh

import legality


def _box(x0, x1, y0, y1, z0, z1):
    return trimesh.creation.box(bounds=[[x0, y0, z0], [x1, y1, z1]])


def test_buried_pylon_joint_spans_its_width():
    wing = _box(0.0, 15.0, 0.0, 25.5, 49.2, 51.84)          # right half of a wing, mm
    pylon = _box(3.0, 13.0, 0.0, 1.0, 43.0, 50.5)            # a 2 mm centre pylon, its top 1.3 mm inside the wing
    lo, hi = legality.support_joint_y(wing, pylon)
    assert lo == 0.0 and abs(hi - 1.0) < 1e-9                # the face at y = 1 counts, not only the centre plane


def test_support_off_the_wing_falls_back_to_its_bounds():
    wing = _box(0.0, 15.0, 0.0, 25.5, 49.2, 51.84)
    strut = _box(3.0, 13.0, 30.0, 32.0, 40.0, 45.0)          # not touching the wing
    assert legality.support_joint_y(wing, strut) == (30.0, 32.0)
