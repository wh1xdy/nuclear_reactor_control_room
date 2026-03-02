# Model Limits and Verification/Validation Notes

This simulator is educational and not certified for operator licensing or safety-case use.

## Scope of validity (intended)
- Qualitative transient behavior for PWR/BWR/RBMK operation around normal and upset conditions.
- Relative trend exploration for control actions, alarm/trip latching, and simplified balance-of-plant coupling.

## Not represented in high fidelity
- Full 3D neutronics (nodal diffusion/transport), detailed reactor vessel CFD, or full plant network hydraulics.
- Certified instrumentation uncertainty analysis and complete safety-system actuation diversity.
- Regulatory acceptance criteria for licensing analyses.

## Current V&V approach in repo
- Deterministic transient checks in `validation.py`.
- Benchmark-envelope acceptance in `data/benchmark_envelopes.json`.
- Protection workflow checks: trip latching, acknowledge, shelving, reset permissives, and channel disagreement diagnostics.

## Recommended next V&V steps
- Add scenario matrix with pass/fail rationale and traceability IDs.
- Add external benchmark references with cited sources and uncertainty decomposition.
- Add Monte Carlo perturbation tests for sensor bias/failure combinations.
