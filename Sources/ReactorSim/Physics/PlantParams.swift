// PlantParams.swift — Tunable physical constants for a PWR.
// All values taken from open literature; see physics_engine.py for references.
// Designed for extension: swap a different PlantParams to tune or add a new reactor.

import Foundation

struct PlantParams: Sendable {
    // MARK: — Six-group delayed neutron data (U-235 thermal fission)
    let beta: [Double] = [2.15e-4, 1.424e-3, 1.274e-3, 2.568e-3, 7.48e-4, 2.73e-4]
    let lambdaD: [Double] = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]

    // Prompt neutron generation time Λ [s] — ~20 µs for PWR thermal spectrum
    var lambdaPrompt: Double = 2.0e-5

    // MARK: — Reactor power
    var nominalPower: Double = 3.0e9          // W (3000 MWt, typical PWR)

    // MARK: — Thermal capacities [J/K]
    var fuelHeatCapacity: Double     = 1.0e8
    var coolantHeatCapacity: Double  = 4.0e7
    var sgHeatCapacity: Double       = 1.0e8

    // MARK: — Heat transfer coefficients [W/K]
    // Calibrated so that at n=1.0, flow=1.0, turbine=1.0 the steady-state
    // temperatures match the nominal values (T_fuel=900K, T_cool=550K, T_sg=490K):
    //   h_FC = P_nom / (T_fuel − T_cool) = 3e9 / 350 = 8.57e6
    //   h_CS = P_nom / (T_cool − T_sg)   = 3e9 / 60  = 5.0e7
    //   h_ST = P_nom / (T_sg  − T_cond)  = 3e9 / 180 = 1.67e7
    var hFuelToCoolant: Double   = 8.57e6
    var hCoolantToSG: Double     = 5.0e7
    var hSGToTurbine: Double     = 1.67e7
    // Condenser cold-side temperature (turbine exhaust sink) [K]
    var condenserTempK: Double   = 310.0

    // MARK: — Turbine
    var turbineEfficiency: Double = 0.33      // ~33% thermal efficiency

    // MARK: — Control rod reactivity
    var rodWorth: Double        = -0.05       // Δk/k fully inserted
    var scramExtraWorth: Double = -0.12       // additional emergency boration+drop
    var rodDropTau: Double      = 2.0         // s — gravity drop time
    // CRDM stepping speed [fraction of full travel per s].
    // Realistic PWR: ~72 steps/min of a 228-step stroke ≈ 0.0053/s (~3 min full travel).
    var rodSpeed: Double        = 0.0053

    // MARK: — Temperature feedback coefficients [Δk/k per K]
    var fuelTempCoeff: Double    = -2.0e-5    // Doppler, −2 pcm/K
    var coolantTempCoeff: Double = -3.0e-4    // moderator, −20 to −50 pcm/K

    // MARK: — Xenon / Iodine dynamics
    var gammaXe: Double         = 0.065       // effective I-135 fission yield
    var lambdaI: Double         = 0.6931 / (6.57 * 3600)   // I-135 decay [1/s]
    var lambdaXe: Double        = 0.6931 / (9.17 * 3600)   // Xe-135 decay [1/s]
    var xenonBurnCoeff: Double  = 2.1e-5      // σ_a(Xe)·φ at full power [1/s]
    var xenonReactivityCoeff: Double = 1.6e-5 // maps Xe inventory → Δk/k

    // MARK: — Nominal operating temperatures [K]
    var nominalFuelTemp: Double    = 900.0
    var nominalCoolantTemp: Double = 550.0
    var nominalSGTemp: Double      = 490.0   // T_sg = T_cool − P/(h_CS) = 550−60 = 490

    // MARK: — Integration
    // Thermal/xenon substep. Fastest thermal time constant is the coolant node
    // (~0.7 s), so 0.05 s keeps explicit Euler 14× under its stability limit
    // while making 600× time compression ~5× cheaper per frame. Kinetics has
    // its own adaptive sub-stepping below this.
    var internalDt: Double = 0.05
    var externalReactivity: Double = 0.0      // boron, instructor override

    // MARK: — Derived
    var betaTotal: Double { beta.reduce(0, +) }
}
