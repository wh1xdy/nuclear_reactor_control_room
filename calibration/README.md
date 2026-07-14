# OpenMC calibration pipeline

Anchors ReactorSim's reduced-order core constants to continuous-energy Monte
Carlo transport (OpenMC + ENDF/B-VII.1, NNDC HDF5 library).

- `pincell.py` — 3.1 w/o PWR pin cell: boron worth [pcm/ppm], Doppler (FTC)
  and moderator (MTC) coefficients via finite differences (bias cancels).
- `rodcurve.py` — 17×17 B4C-rodded 3D assembly, uniform-fuel box source,
  60 inactive / 150 batches (tall-core source convergence!): k vs insertion
  + 15-bin axial fission profiles.
- `aggregate.py` — merges results into `calibration.json`. The integral
  rod-worth SHAPE comes from φ²-weighting of the converged ARO axial profile
  (perturbation theory) — the raw deep-insertion k-worth is an
  infinite-lattice artifact and is kept only as reference.

The app bundles `calibration.json` (Sources/ReactorSim/Resources/); at
startup `PlantParams.applyCalibration()` overrides boron worth, FTC, MTC and
the rod-shape table for PWR/SMR. The TOTAL bank worth stays at the design
−5000 pcm (a single assembly cannot measure all banks). BWR keeps its own
void-dominated profile.

Rebuild: OpenMC from source (arm64: libomp flags + vendored fmt), data from
openmc.org; run the two models, then `aggregate.py`. Full recipe in the
project memory / session logs.

Result 2026-07-13: boron −8.78 pcm/ppm (trainer had −8), FTC −1.47 (had −2),
MTC −18.0 (had −30), rod S-curve within a few % of the analytic
integrated-cosine — the hand calibration held up honourably.
