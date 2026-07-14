#!/usr/bin/env python3
"""Merge pincell + rodcurve results into calibration.json for ReactorSim."""
import json, subprocess, sys, os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
pin = json.load(open("results/pincell.json"))
rod = json.load(open("results/rodcurve.json"))

curve = rod["worth_curve"]                 # [[insertion, pcm-vs-ARO]...] (raw)
total = curve[-1][1]
# The k-based deep-end is an infinite-lattice artifact (reflective assembly →
# fully-rodded k collapses). The SOUND shape comes from perturbation theory:
# differential worth ∝ φ²(z) at the rod tip, so the integral worth for a
# TOP-entry rod inserted to fraction x is ∫ over the covered top span of the
# measured ARO axial flux squared.
phi = rod["axial"]["0.0"]                  # ARO axial profile, 15 bins bottom→top
n = len(phi)
phi2 = [p * p for p in phi]
tot2 = sum(phi2)
xs, shape = [], []
for i in range(n + 1):
    x = i / n                              # insertion fraction (from the top)
    covered = phi2[n - i:]                 # top i bins are rodded
    xs.append(round(x, 4))
    shape.append(round(sum(covered) / tot2, 5))

commit = subprocess.run(["git", "-C", "src", "log", "-1", "--format=%h"],
                        capture_output=True, text=True).stdout.strip()
out = {
    "source": f"OpenMC {commit} + ENDF/B-VII.1 (NNDC HDF5)",
    "models": "3.1 w/o pin cell (coefficients) + 17x17 B4C-rodded 3D assembly (rod curve, axial)",
    "boron_pcm_per_ppm": round(pin["boron_pcm_per_ppm"], 3),
    "ftc_pcm_per_K": round(pin["ftc_pcm_per_K"], 3),
    "mtc_pcm_per_K": round(pin["mtc_pcm_per_K"], 3),
    "rod_worth_total_pcm_raw_inflattice": round(total, 1),
    "rod_shape_x": xs,
    "rod_shape_w": shape,
    "axial_shapes": rod["axial"],
}
json.dump(out, open("calibration.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "axial_shapes"}, indent=1))
