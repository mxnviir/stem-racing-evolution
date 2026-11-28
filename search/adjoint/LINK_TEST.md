# The adjoint link's exact test: results (2026-11-28)

Spec and marks: `test_link.py`'s docstring. The marks were written before the first run, with one dated amendment
after gene 7's result. There is no CFD in this test.

DAFoam's dD/dX is replaced by the surface's area vector (times f). The link must then return the derivative of the
body's volume (f = 1) and of its first moments (f = x, f = z). That is measured by rebuilding the body with
`param_body.build`, the searches' own build, at u_k ± 0.25.

**The car:** the team layout, with its own free-form body as the seed. The CMA numbers tested:
- 7: a top number (station 4 of 5, outer);
- 19: a bottom number;
- 33: a side number.

## Verdict as registered: FAIL, by one of nine checks (gene 33's x-moment at 0.887 against the 0.90 bar)

The table gives the predicted change divided by the rebuilt one.

| Gene | Moment | Before the amendment | With the rule masks (the amendment) | Mark (within 10 %) |
|---|---|---|---|---|
| 7 (top; set the threshold) | V | 1.304 | 1.014 | ok |
| 7 | x | 1.364 | 1.024 | ok |
| 7 | z | 1.458 | 1.051 | ok |
| **19 (bottom; held out)** | V | 1.012 | 0.959 | ok |
| 19 | x | 1.009 | 0.951 | ok |
| 19 | z | 0.975 | 0.918 | ok |
| **33 (side; held out)** | V | 0.934 | 0.930 | ok |
| 33 | x | 0.888 | **0.887** | **fail** |
| 33 | z | 0.911 | 0.908 | ok |

### What the test found

1. **A real defect, fixed.** Surface that the build's hard rule masks hold doesn't move with the numbers, and the
   link counted it.
   - About a third of gene 7's prediction sat on the spine (|y| < 13 mm, z 24–43 mm), on the chamber wall and the
     halo pocket.
   - The seed field was fitted to the built body, so it lies within 1 mm of those surfaces too. The "|phi − level| <
     1 mm" test can't tell them apart.
   - Per-point displacements showed those points held in both directions.
   - `link.Masks` / `build_masks` now exclude points within 0.6 build cells of a mask cell. Building the masks takes
     12 s per layout. The threshold was set on gene 7 alone.
   - The carhybrid `grad.py` (the shared core) has the same test and the same defect.
2. **Gene 33 is curved over the test step.** Its one-sided slopes differ by 20 % on the x-moment: +0.25 gives 245k,
   −0.25 gives 293k, both in mm⁴ per unit.
   - The central difference (269k) is therefore not the derivative at u = 0, to within the mark's 10 %.
   - The prediction (238k) sits at the near end of the one-sided range.
   - Likely sources: the side's clamp at zero half-width in `from_unit`, and the milling at the side.
   - The mask fix doesn't touch it (0.888 → 0.887).
3. **No chain errors.** All nine ratios lie in 0.89–1.05 with the right signs. A wrong frame, unit, sign, normal or
   chain rule would show up as factors of 10³, sign flips, or ratios far from 1.

### What it means for use

The link gives each number's movement of the body to within about 5–11 %. For the guided search that is ample:
- CMA uses only the gradient's direction;
- its own measurements choose the step;
- the adjoint itself overstates gains about 2.5×.

A finer check of gene 33, with a smaller step and a finer build grid, can be run if wanted. It wasn't, since the
result would not change how the gradient is used.

### Files

- `results/test_link_team_unmasked.json`
- `results/test_link_team_masked.json`
- `results/xcheck_p0_team.json`: the link against `grad.py` on the p0 team-car adjoint, cosine 1.0000.
