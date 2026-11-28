# Builder test suite runs

`.github/workflows/evo_tests.yml` runs the pinned evolution builder's own tests (`claude/loving-keller-q72mvv` at
`e9c896b`), split 20 ways with 320 minutes per shard (the repo's own CI stops two parts shards at 90 minutes), with
this branch's runtime fixes (`evo_patches/apply.py`) applied when present, plus `evo_patches/tests`.
Each run's per-test results and `SUMMARY.md` land in `evo_tests/results/<run id>/`. Launch: edit `request.json`.

| Run | What | Result |
|---|---|---|
| 1 | baseline, no patches | running |
