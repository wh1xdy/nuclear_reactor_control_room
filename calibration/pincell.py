#!/usr/bin/env python3
"""Pin-cell coefficient sweeps for ReactorSim calibration.

Builds a standard 17x17-class PWR pin cell (UO2 3.1 w/o, Zr clad, borated
light water) and sweeps:
  - boron concentration  -> boron worth   [pcm/ppm]
  - fuel temperature     -> Doppler (FTC) [pcm/K]
  - moderator temp+dens  -> MTC           [pcm/K]

Writes results to results/pincell.json. Each k-eff point runs quickly at
modest statistics; worth coefficients come from finite differences, so the
absolute k bias cancels to first order.
"""
import json, os, sys
import openmc

DATA = os.path.expanduser("~/kärnreaktor/calibration/data/endfb-vii.1-hdf5/cross_sections.xml")
os.environ.setdefault("OPENMC_CROSS_SECTIONS", DATA)

PITCH = 1.26            # cm
R_FUEL, R_CLAD_IN, R_CLAD_OUT = 0.4096, 0.4180, 0.4750

def water_density(t_k):
    """Liquid density at 15.5 MPa — IAPWS-anchored piecewise-linear table."""
    pts = [(400, 0.947), (450, 0.905), (500, 0.838), (550, 0.767),
           (570, 0.732), (590, 0.693), (600, 0.660), (615, 0.610)]
    for (t0, r0), (t1, r1) in zip(pts, pts[1:]):
        if t_k <= t1:
            return r0 + (r1 - r0) * (t_k - t0) / (t1 - t0)
    return 0.55

def build(boron_ppm, t_fuel, t_mod):
    fuel = openmc.Material(name="uo2")
    fuel.set_density("g/cm3", 10.4)
    fuel.add_element("U", 1.0, enrichment=3.1)
    fuel.add_element("O", 2.0)
    fuel.temperature = t_fuel

    clad = openmc.Material(name="zirc")
    clad.set_density("g/cm3", 6.55)
    clad.add_element("Zr", 1.0)
    clad.temperature = t_mod

    water = openmc.Material(name="water")
    water.set_density("g/cm3", water_density(t_mod))
    wB = boron_ppm * 1e-6                     # ppm BY WEIGHT, done properly
    water.add_element("H", 0.1119 * (1 - wB), percent_type="wo")
    water.add_element("O", 0.8881 * (1 - wB), percent_type="wo")
    if wB > 0:
        water.add_element("B", wB, percent_type="wo")
    water.add_s_alpha_beta("c_H_in_H2O")
    water.temperature = t_mod

    f = openmc.ZCylinder(r=R_FUEL)
    ci = openmc.ZCylinder(r=R_CLAD_IN)
    co = openmc.ZCylinder(r=R_CLAD_OUT)
    box = openmc.model.RectangularPrism(PITCH, PITCH, boundary_type="reflective")
    cells = [
        openmc.Cell(fill=fuel, region=-f),
        openmc.Cell(region=+f & -ci),                 # gap
        openmc.Cell(fill=clad, region=+ci & -co),
        openmc.Cell(fill=water, region=+co & -box),
    ]
    model = openmc.Model()
    model.geometry = openmc.Geometry(cells)
    model.settings.batches = 60
    model.settings.inactive = 15
    model.settings.particles = 5000
    model.settings.temperature = {"method": "interpolation"}
    model.settings.output = {"tallies": False}
    return model

def keff(boron, tf, tm, tag):
    model = build(boron, tf, tm)
    sp = model.run(cwd=f"runs/{tag}", output=False)
    with openmc.StatePoint(sp) as s:
        k = s.keff
    return k.nominal_value, k.std_dev

def main():
    os.makedirs("results", exist_ok=True)
    res = {}
    base = (800.0, 900.0, 570.0)
    k0, s0 = keff(*base, "base")
    res["k_base"] = [k0, s0]
    # Boron worth: +200 ppm
    k1, s1 = keff(1000.0, 900.0, 570.0, "boron")
    res["boron_pcm_per_ppm"] = (k1 - k0) / k0 * 1e5 / 200.0
    # Doppler: +300 K fuel
    k2, s2 = keff(800.0, 1200.0, 570.0, "doppler")
    res["ftc_pcm_per_K"] = (k2 - k0) / k0 * 1e5 / 300.0
    # MTC: +20 K moderator (density follows)
    k3, s3 = keff(800.0, 900.0, 590.0, "mtc")
    res["mtc_pcm_per_K"] = (k3 - k0) / k0 * 1e5 / 20.0
    json.dump(res, open("results/pincell.json", "w"), indent=1)
    print(json.dumps(res, indent=1))

if __name__ == "__main__":
    main()
