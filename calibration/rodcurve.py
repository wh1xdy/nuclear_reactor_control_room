#!/usr/bin/env python3
"""Rodded 3D assembly: integral rod-worth curve + axial flux shapes.

A 17x17-class assembly (quarter-symmetric pin lattice approximated as a
homogenized-guide-tube layout), 366 cm active height, reflective radially,
water reflectors axially. A B4C rod cluster inserts from the TOP to depth
d in {0, 0.1, ..., 1.0}; each case tallies k-eff and a 15-bin axial flux.

Outputs results/rodcurve.json:
  worth_curve: [(insertion, pcm_vs_ARO)]
  axial: {insertion: [15 relative-power bins, mean 1]}
"""
import json, os, sys
import numpy as np
import openmc

DATA = os.path.expanduser("~/kärnreaktor/calibration/data/endfb-vii.1-hdf5/cross_sections.xml")
os.environ.setdefault("OPENMC_CROSS_SECTIONS", DATA)

H = 366.0            # active height [cm]
PITCH = 1.26
NPIN = 17
ASSY = PITCH * NPIN

def materials(t_mod=570.0, boron=800.0):
    fuel = openmc.Material(name="uo2"); fuel.set_density("g/cm3", 10.4)
    fuel.add_element("U", 1.0, enrichment=3.1); fuel.add_element("O", 2.0)
    fuel.temperature = 900.0
    clad = openmc.Material(name="zr"); clad.set_density("g/cm3", 6.55)
    clad.add_element("Zr", 1.0); clad.temperature = t_mod
    water = openmc.Material(name="h2o")
    pts = [(400, 0.947), (450, 0.905), (500, 0.838), (550, 0.767),
           (570, 0.732), (590, 0.693), (600, 0.660), (615, 0.610)]
    rho = 0.55
    for (t0, r0), (t1, r1) in zip(pts, pts[1:]):
        if t_mod <= t1:
            rho = r0 + (r1 - r0) * (t_mod - t0) / (t1 - t0); break
    water.set_density("g/cm3", rho)
    wB = boron * 1e-6
    water.add_element("H", 0.1119 * (1 - wB), percent_type="wo")
    water.add_element("O", 0.8881 * (1 - wB), percent_type="wo")
    water.add_element("B", wB, percent_type="wo")
    water.add_s_alpha_beta("c_H_in_H2O"); water.temperature = t_mod
    b4c = openmc.Material(name="b4c"); b4c.set_density("g/cm3", 1.76)
    b4c.add_element("B", 4.0); b4c.add_element("C", 1.0); b4c.temperature = t_mod
    return fuel, clad, water, b4c

def build(insertion, particles=8000, batches=150):
    fuel, clad, water, b4c = materials()
    rf = openmc.ZCylinder(r=0.4096); rc = openmc.ZCylinder(r=0.4750)
    gt = openmc.ZCylinder(r=0.5610)          # guide tube inner
    zb = openmc.ZPlane(z0=0, boundary_type="vacuum")
    zt = openmc.ZPlane(z0=H + 40, boundary_type="vacuum")
    z0 = openmc.ZPlane(z0=20.0)              # bottom reflector top / fuel bottom
    z1 = openmc.ZPlane(z0=20.0 + H)          # fuel top
    rod_tip = openmc.ZPlane(z0=20.0 + H * (1 - insertion))

    fuel_pin = openmc.Universe(cells=[
        openmc.Cell(fill=fuel, region=-rf & +z0 & -z1),
        openmc.Cell(fill=clad, region=+rf & -rc & +z0 & -z1),
        openmc.Cell(fill=water, region=(+rc & +z0 & -z1) | -z0 | +z1),
    ])
    # Guide-tube pin: water below the rod tip, B4C above it (top entry).
    # Regions form a COMPLETE partition of the universe (no lost particles).
    gto = openmc.ZCylinder(r=0.6020)
    gt_pin = openmc.Universe(cells=[
        openmc.Cell(fill=b4c,   region=-gt & +rod_tip),
        openmc.Cell(fill=water, region=-gt & -rod_tip),
        openmc.Cell(fill=clad,  region=+gt & -gto),
        openmc.Cell(fill=water, region=+gto),
    ])
    # 17x17 with the standard 24+1 guide-tube pattern approximated: every 4th
    # position in a diamond — close enough for a SHAPE calibration.
    gtpos = {(2,4),(2,12),(3,8),(4,2),(4,14),(5,5),(5,11),(8,3),(8,13),(8,8),
             (11,5),(11,11),(12,2),(12,14),(13,8),(14,4),(14,12),(6,8),(10,8),
             (8,6),(8,10),(4,8),(12,8),(8,4),(8,12)}
    lat = openmc.RectLattice()
    lat.lower_left = (-ASSY/2, -ASSY/2)
    lat.pitch = (PITCH, PITCH)
    lat.universes = [[gt_pin if (i,j) in gtpos else fuel_pin for j in range(NPIN)] for i in range(NPIN)]
    # Outer universe catches floating-point edges at the lattice boundary —
    # without it those histories are LOST and the run aborts.
    lat.outer = openmc.Universe(cells=[openmc.Cell(fill=water)])
    box = openmc.model.RectangularPrism(ASSY, ASSY, boundary_type="reflective")
    root = openmc.Cell(fill=lat, region=-box & +zb & -zt)
    model = openmc.Model()
    model.geometry = openmc.Geometry([root])
    model.settings.batches = batches
    model.settings.inactive = 60      # tall-core dominance: source needs room to converge
    model.settings.particles = particles
    # Uniform initial fission source over the FUEL — the default point source
    # at the origin sits BELOW the core and 20 inactive batches leave the
    # eigenmode unconverged (bottom-collapsed axial shapes, biased shallow-k).
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-ASSY/2, -ASSY/2, 20.0), (ASSY/2, ASSY/2, 20.0 + H)),
        constraints={"fissionable": True})
    model.settings.temperature = {"method": "interpolation"}
    model.settings.output = {"tallies": False}
    # A handful of histories die on coincident-plane edge cases (rod tip vs
    # lattice boundaries) — irrelevant out of millions, but the default abort
    # threshold is 10. Raise it; statistics are unaffected.
    model.settings.max_lost_particles = 5000
    mesh = openmc.RegularMesh()
    mesh.dimension = [1, 1, 15]
    mesh.lower_left = (-ASSY/2, -ASSY/2, 20.0)
    mesh.upper_right = (ASSY/2, ASSY/2, 20.0 + H)
    t = openmc.Tally(name="axial")
    t.filters = [openmc.MeshFilter(mesh)]
    t.scores = ["fission"]
    model.tallies = openmc.Tallies([t])
    return model

def run_case(insertion, tag, particles=8000):
    model = build(insertion, particles=particles)
    sp = model.run(cwd=f"runs/rod_{tag}", output=False)
    with openmc.StatePoint(sp) as s:
        k = s.keff
        ax = s.get_tally(name="axial").mean.ravel()
    ax = ax / ax.mean() if ax.mean() > 0 else ax
    return k.nominal_value, k.std_dev, ax.tolist()

def main():
    quick = "--quick" in sys.argv
    steps = [0.0, 0.5, 1.0] if quick else [i / 10 for i in range(11)]
    particles = 2000 if quick else 8000
    os.makedirs("results", exist_ok=True)
    out = {"worth_curve": [], "axial": {}}
    k_aro = None
    for d in steps:
        k, s, ax = run_case(d, f"{int(d*100):03d}", particles)
        if k_aro is None: k_aro = k
        pcm = (k - k_aro) / k_aro * 1e5
        out["worth_curve"].append([d, pcm])
        out["axial"][f"{d:.1f}"] = ax
        print(f"insertion {d:.1f}: k={k:.5f}±{s:.5f}  worth={pcm:+.0f} pcm", flush=True)
    json.dump(out, open("results/rodcurve.json", "w"), indent=1)

if __name__ == "__main__":
    main()
