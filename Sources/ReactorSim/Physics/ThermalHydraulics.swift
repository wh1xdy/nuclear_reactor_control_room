// ThermalHydraulics.swift — Lumped-parameter primary + secondary thermal model.
//
// Primary loop is a TWO-node coolant circuit with a transport delay:
//   fuel → hot-leg (core outlet) →[τ]→ steam generator → cold-leg (SG outlet) →[τ]→ core
// so the loop carries a real ΔT (~30 K) and power changes reach the SG (and
// return to the core) with a realistic ~4 s-per-leg lag instead of instantly.
//
// Refinements over the old single-node model:
//  • hot-leg / cold-leg split → distinct T_hot, T_cold; T_avg drives feedback
//  • Dittus–Boelter: convective HTC scales with flow^0.8 (not linearly)
//  • mass advection scales LINEARLY with flow (∝ ṁ) — the physically correct split
//  • saturation steam-pressure model: SG secondary pressure = P_sat(T_sg)
//  • Carnot-scaled gross efficiency tied to the SG/condenser temperatures

import Foundation

final class ThermalHydraulics {
    let params: PlantParams

    var tFuel: Double       // fuel temperature [K]
    var tHot:  Double       // core-outlet / hot-leg coolant [K]
    var tCold: Double       // SG-outlet / cold-leg coolant (= core inlet) [K]
    var tHL:   Double       // coolant delivered to the SG (hot leg, transport-lagged)
    var tCL:   Double       // coolant delivered to the core (cold leg, transport-lagged)
    var tSG:   Double       // steam-generator secondary (saturation) node [K]
    var electricPower: Double
    private(set) var steamPressureMPa: Double = 0
    /// Core-average void fraction. Inert for subcooled PWR/SMR (stays 0); for a
    /// BWR it tracks boiling intensity (∝ power) against flow, and the plant
    /// reads it for the void-reactivity feedback. Kept OUT of the temperature
    /// ODEs so PWR thermal behavior is unchanged.
    private(set) var voidFraction: Double

    /// Reactor coolant average temperature — drives moderator feedback & display.
    var tAvg: Double { 0.5 * (tHot + tCold) }
    /// Back-compatible alias (old single-node name).
    var tCool: Double { tAvg }

    init(_ params: PlantParams) {
        self.params = params
        tFuel = params.nominalFuelTemp                  // 900 K
        let avg    = params.nominalCoolantTemp          // 550 K
        let halfDT = 0.5 * params.nominalCoreDeltaT     // 15 K
        tHot  = avg + halfDT                            // 565 K
        tCold = avg - halfDT                            // 535 K
        tHL   = tHot
        tCL   = tCold
        tSG   = params.nominalSGTemp                    // 490 K
        electricPower = params.nominalPower * params.turbineEfficiency
        steamPressureMPa = ThermalHydraulics.satPressureMPa(tSG)
        voidFraction = params.nominalVoidFraction       // start at equilibrium void
    }

    /// Advance by dt seconds given thermal power P_th [W] and control inputs.
    /// Returns gross electric power [W].
    @discardableResult
    func step(dt: Double, thermalPower: Double, flow: Double, turbineValve: Double) -> Double {
        let p = params
        let f = max(flow, 0)

        // Dittus–Boelter: convective heat-transfer coefficient ∝ Re^0.8 ∝ flow^0.8.
        let fHT = pow(f, 0.8)
        let hFC = p.hFuelToCoolant * fHT
        let hCS = p.hCoolantToSG   * fHT
        // Mass advection is LINEAR in mass flow (carries enthalpy ∝ ṁ).
        let adv = p.primaryAdvection * f

        let tAvgCore = 0.5 * (tHot + tCold)
        let qFC = hFC * (tFuel - tAvgCore)                 // fuel → coolant
        // SG evaporator boils at ~553K (6.5 MPa); the hot leg (≈565K) drives it
        // with a ~12 K pinch — the ceiling for this primary T-avg.
        let qCS = hCS * (tHL - tSG)                        // coolant → SG secondary
        // Steam removal to the condenser, modulated by the turbine/dump valve.
        let qST = p.hSGToTurbine * max(turbineValve, 0) * max(0, tSG - p.condenserTempK)

        // ── Primary nodes ───────────────────────────────────────────────────
        // Fuel: heat generated minus heat to coolant
        tFuel += (thermalPower - qFC) / p.fuelHeatCapacity * dt
        // Hot leg (core coolant): fuel heat in, advect the cold-leg inlet through
        tHot  += (qFC + adv * (tCL - tHot)) / p.hotLegCapacity * dt
        // Cold leg (SG outlet): advect the hot-leg inlet in, SG removes heat
        tCold += (adv * (tHL - tCold) - qCS) / p.coldLegCapacity * dt

        // Transport delay (first-order lag ≈ plug-flow transit). Transit time
        // lengthens as flow drops (τ ∝ 1/flow), so a coastdown stretches the lag.
        let tauH = p.tauHotLeg  / max(f, 0.05)
        let tauC = p.tauColdLeg / max(f, 0.05)
        tHL += (tHot  - tHL) / tauH * dt
        tCL += (tCold - tCL) / tauC * dt

        // ── Secondary node ──────────────────────────────────────────────────
        tSG += (qCS - qST) / p.sgHeatCapacity * dt

        // Saturation steam pressure from the SG secondary temperature.
        steamPressureMPa = ThermalHydraulics.satPressureMPa(tSG)

        // ── Core void fraction (BWR feedback driver; inert when nominal=0) ────
        // Boiling scales with power; flow sweeps voids out / subcools the core.
        // First-order lag (~1.5 s) — voids respond fast but not instantly.
        if p.nominalVoidFraction > 0 {
            let pf = max(0, thermalPower / max(p.nominalPower, 1))
            let target = min(0.9, p.nominalVoidFraction * pf / max(f, 0.1))
            voidFraction += (target - voidFraction) / 1.5 * dt
            voidFraction = max(0, min(0.95, voidFraction))
        }

        // Carnot-scaled gross efficiency: warmer steam / colder sink → more work.
        let eta = p.turbineCarnotFraction * max(0, 1 - p.condenserTempK / max(tSG, 1))
        electricPower = max(qST, 0) * eta
        return electricPower
    }

    /// Antoine saturation pressure for water, T[K] → P[MPa]. Constants re-anchored
    /// to steam-table points across the plant's range — condenser (310K→6.2 kPa),
    /// mid (490K→2.1 MPa) and PWR S/G (553K→6.4 MPa) — since the old 379–573 K set
    /// read ~19% low near 553 K and would have shown a sub-realistic S/G pressure.
    static func satPressureMPa(_ tK: Double) -> Double {
        let A = 5.444, B = 1951.8, C = -16.5
        let logPbar = A - B / (max(tK, 50) + C)   // P in bar
        return pow(10, logPbar) / 10.0            // bar → MPa
    }

    /// Exact inverse of satPressureMPa: saturation temperature [K] at P [MPa].
    /// (15.5 MPa → ~616 K, 7 MPa → ~559 K, 6.4 MPa → ~553 K.)
    static func satTempK(_ pMPa: Double) -> Double {
        let A = 5.444, B = 1951.8, C = -16.5
        return B / (A - log10(max(pMPa, 0.001) * 10)) - C
    }
}
