// PlantParams.swift — Tunable physical constants for a PWR.
// All values taken from open literature; see physics_engine.py for references.
// Designed for extension: swap a different PlantParams to tune or add a new reactor.

import Foundation

/// Which reactor the parameter set describes. The shared sub-models (point
/// kinetics, decay heat, xenon) are identical across kinds; only the thermal
/// profile, the void-feedback term, and the balance-of-plant differ.
enum ReactorKind: String, Sendable, CaseIterable { case pwr, bwr, smr }

struct PlantParams: Sendable {
    // MARK: — Reactor identity / profile
    // PWR defaults below keep `PlantParams()` bit-identical to the validated
    // model; `bwr()` / `smr()` factories override only what genuinely differs.
    var kind: ReactorKind = .pwr
    var hasPressurizer: Bool = true        // PWR/SMR: PZR holds primary pressure
    var hasBoron: Bool = true              // chemical shim (PWR/SMR); BWR has none in normal ops
    var hasSteamGenerator: Bool = true     // indirect cycle; BWR is direct (steam from core)
    var naturalCirculation: Bool = false   // SMR: buoyancy-driven flow, no RCPs
    var nominalPressureMPa: Double = 15.5  // primary (PWR/SMR) or steam-dome (BWR)
    /// Void-fraction reactivity coefficient [Δk/k per unit void], referenced to
    /// the nominal void below so the steady operating point is unchanged. ≈0 for
    /// subcooled PWR/SMR; strongly negative for a BWR (its defining feedback).
    var voidCoeff: Double = 0.0
    var nominalVoidFraction: Double = 0.0  // core-average void at full power & flow

    // MARK: — Six-group delayed neutron data (U-235 thermal fission)
    let beta: [Double] = [2.15e-4, 1.424e-3, 1.274e-3, 2.568e-3, 7.48e-4, 2.73e-4]
    let lambdaD: [Double] = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]

    // Prompt neutron generation time Λ [s] — ~20 µs for PWR thermal spectrum
    var lambdaPrompt: Double = 2.0e-5

    // Intrinsic/startup neutron source [relative n · s⁻¹]. Holds a source-range
    // floor when deeply subcritical so the approach to criticality gives real
    // 1/M behaviour (subcritical multiplication of a fixed source). Utterly
    // negligible at power — at n≈1 the precursor source is ~3×10² s⁻¹ — so the
    // validated at-power physics is unchanged.
    var neutronSource: Double = 1.0e-5

    // MARK: — Reactor power
    var nominalPower: Double = 3.0e9          // W (3000 MWt, typical PWR)

    // MARK: — Thermal capacities [J/K]
    // Coolant is split into hot-leg (core outlet) and cold-leg (SG outlet)
    // lumps so the loop carries a real ΔT and a transport delay.
    var fuelHeatCapacity: Double     = 1.0e8
    var hotLegCapacity: Double       = 2.0e7
    var coldLegCapacity: Double      = 2.0e7
    var sgHeatCapacity: Double       = 1.0e8

    // MARK: — Heat transfer coefficients [W/K]
    // Calibrated so that at n=1.0, flow=1.0, turbine=1.0 the steady state matches
    // nominal (T_fuel=900K, T_avg=550K, T_sg=553K — a realistic ~6.5 MPa S/G):
    //   h_FC = P_nom / (T_fuel  − T_avg) = 3e9 / 350 = 8.57e6
    //   h_CS = P_nom / (T_hot   − T_sg)  = 3e9 / 12  = 2.5e8   (evaporator pinch off the HOT leg)
    //   h_ST = P_nom / (T_sg    − T_cond)= 3e9 / 243 = 1.235e7
    // The S/G boils at 553K (~280°C); the primary HOT leg (565K) drives it with a
    // 12 K pinch — tight but physical, and the ceiling for this T-avg. These are
    // CONVECTIVE coefficients; they scale with flow^0.8 (Dittus–Boelter).
    var hFuelToCoolant: Double   = 8.57e6
    var hCoolantToSG: Double     = 2.5e8
    var hSGToTurbine: Double     = 1.235e7
    // Condenser cold-side temperature (turbine exhaust sink) [K]
    var condenserTempK: Double   = 310.0

    // Primary mass-advection conductance at flow=1 [W/K]. Sets the core ΔT:
    // P_nom / advection = ΔT_core → 3e9 / 1e8 = 30 K. Advection ∝ ṁ (LINEAR),
    // unlike the convective HTCs above.
    var primaryAdvection: Double     = 1.0e8
    var nominalCoreDeltaT: Double    = 30.0   // T_hot − T_cold at full power/flow [K]
    // Loop transport delay at full flow [s] — coolant transit core→SG and back.
    // Scales as 1/flow (slower coolant = longer transit).
    var tauHotLeg: Double            = 4.0
    var tauColdLeg: Double           = 4.0

    // MARK: — Turbine
    var turbineEfficiency: Double = 0.33      // ~33% thermal efficiency (legacy seed)
    // Gross efficiency is Carnot-scaled: η = fraction·(1 − T_cond/T_sg). At the
    // realistic S/G (T_sg=553, T_cond=310) the Carnot ceiling is 1−310/553 ≈ 0.44,
    // so 0.751·0.44 ≈ 0.33 net — the ~10-pt loss to irreversibilities/aux load that
    // makes 33% credible at 6.5 MPa steam (NOT the 0.9·0.37 of the old 1.9 MPa state).
    var turbineCarnotFraction: Double = 0.751

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
    // Boron differential worth [pcm/ppm] (negative). Overridden by the OpenMC
    // calibration when calibration.json is bundled.
    var boronWorthPcmPerPpm: Double = -8.0

    // Optional OpenMC-derived integral rod-worth SHAPE (normalized 0…1 over
    // insertion 0…1). The TOTAL bank worth stays `rodWorth` — a single rodded
    // assembly can't measure all banks; transport gives us the curve's shape.
    var rodShapeX: [Double]? = nil
    var rodShapeW: [Double]? = nil

    /// Integral rod-worth shape w(x) ∈ 0…1 at insertion x ∈ 0…1: the OpenMC
    /// table when calibrated, else the analytic integrated-cosine S-curve.
    func rodShape(_ x: Double) -> Double {
        let xc = max(0, min(1, x))
        guard let xs = rodShapeX, let ws = rodShapeW, xs.count == ws.count, xs.count >= 2 else {
            return 3 * xc * xc - 2 * xc * xc * xc
        }
        for i in 1..<xs.count where xc <= xs[i] {
            let f = (xc - xs[i - 1]) / max(1e-9, xs[i] - xs[i - 1])
            return ws[i - 1] + (ws[i] - ws[i - 1]) * f
        }
        return ws[ws.count - 1]
    }

    /// Apply the OpenMC calibration file (bundled as calibration.json) on top
    /// of the defaults. Absent/malformed file → defaults, silently.
    mutating func applyCalibration() {
        struct Cal: Decodable {
            var boron_pcm_per_ppm: Double?
            var ftc_pcm_per_K: Double?
            var mtc_pcm_per_K: Double?
            var rod_shape_x: [Double]?
            var rod_shape_w: [Double]?
        }
        guard let url = Bundle.module.url(forResource: "calibration", withExtension: "json"),
              let data = try? Data(contentsOf: url),
              let cal = try? JSONDecoder().decode(Cal.self, from: data) else { return }
        if let v = cal.boron_pcm_per_ppm, v < 0 { boronWorthPcmPerPpm = v }
        if let v = cal.ftc_pcm_per_K, v < 0 { fuelTempCoeff = v * 1e-5 }
        if let v = cal.mtc_pcm_per_K, v < 0 { coolantTempCoeff = v * 1e-5 }
        if let xs = cal.rod_shape_x, let ws = cal.rod_shape_w, xs.count == ws.count, xs.count >= 2 {
            rodShapeX = xs; rodShapeW = ws
        }
    }

    // MARK: — Xenon / Iodine dynamics
    var gammaXe: Double         = 0.065       // effective I-135 fission yield
    var lambdaI: Double         = 0.6931 / (6.57 * 3600)   // I-135 decay [1/s]
    var lambdaXe: Double        = 0.6931 / (9.17 * 3600)   // Xe-135 decay [1/s]
    var xenonBurnCoeff: Double  = 2.1e-5      // σ_a(Xe)·φ at full power [1/s]
    var xenonReactivityCoeff: Double = 1.6e-5 // maps Xe inventory → Δk/k

    // MARK: — Nominal operating temperatures [K]
    var nominalFuelTemp: Double    = 900.0
    var nominalCoolantTemp: Double = 550.0
    var nominalSGTemp: Double      = 553.0   // ~280°C saturation → ~6.5 MPa secondary

    // MARK: — Integration
    // Thermal/xenon substep. Fastest thermal time constant is the coolant node
    // (~0.7 s), so 0.05 s keeps explicit Euler 14× under its stability limit
    // while making 600× time compression ~5× cheaper per frame. Kinetics has
    // its own adaptive sub-stepping below this.
    var internalDt: Double = 0.05
    var externalReactivity: Double = 0.0      // boron, instructor override

    // MARK: — Derived
    var betaTotal: Double { beta.reduce(0, +) }

    // MARK: — Reactor-kind factories
    // PWR is the unmodified baseline. BWR adds void feedback + a direct cycle.
    // SMR is an integral PWR scaled to ~200 MWt with natural circulation.

    static func pwr() -> PlantParams { PlantParams() }

    /// Boiling Water Reactor: direct cycle, ~7 MPa steam dome, no pressurizer or
    /// chemical shim. Power is shaped by recirculation flow through the void
    /// coefficient — raise flow → collapse voids → +reactivity; raise power →
    /// more boiling → −reactivity (a strong, stabilizing negative power coeff).
    static func bwr() -> PlantParams {
        var p = PlantParams()
        p.kind                = .bwr
        p.hasPressurizer      = false
        p.hasBoron            = false      // SLC is emergency-only; not a control lever here
        p.hasSteamGenerator   = false      // steam is raised in the core itself
        p.nominalPressureMPa  = 7.0        // reactor steam dome
        p.nominalVoidFraction = 0.40       // core-average void at full power
        p.voidCoeff           = -0.04      // Δk/k per unit void deviation — the defining feedback
        p.coolantTempCoeff    = -1.5e-4    // moderator-temp coeff weaker; void carries moderation
        return p
    }

    /// Small Modular Reactor: integral PWR (pressurizer + SG inside the vessel)
    /// at ~200 MWt with passive natural-circulation flow. Heat capacities and
    /// conductances scale with power so temperatures and time-constants match a
    /// full PWR — only the absolute power (and MWe) shrink.
    static func smr() -> PlantParams {
        var p = PlantParams()
        p.kind               = .smr
        p.naturalCirculation = true
        p.nominalPressureMPa = 13.8
        let r = 200.0e6 / p.nominalPower   // power ratio 200 MWt : 3000 MWt
        p.nominalPower      *= r
        p.fuelHeatCapacity  *= r
        p.hotLegCapacity    *= r
        p.coldLegCapacity   *= r
        p.sgHeatCapacity    *= r
        p.hFuelToCoolant    *= r
        p.hCoolantToSG      *= r
        p.hSGToTurbine      *= r
        p.primaryAdvection  *= r
        return p
    }
}
