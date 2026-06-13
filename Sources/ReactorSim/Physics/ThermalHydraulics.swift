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
        // SG primary side sees the delivered hot leg averaged with its outlet.
        let tSGprim = 0.5 * (tHL + tCold)
        let qCS = hCS * (tSGprim - tSG)                    // coolant → SG secondary
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

        // Carnot-scaled gross efficiency: warmer steam / colder sink → more work.
        let eta = p.turbineCarnotFraction * max(0, 1 - p.condenserTempK / max(tSG, 1))
        electricPower = max(qST, 0) * eta
        return electricPower
    }

    /// Antoine saturation pressure for water, T[K] → P[MPa]. NIST coefficients
    /// for the 379–573 K range (covers SG and condenser operating points).
    static func satPressureMPa(_ tK: Double) -> Double {
        let A = 4.6543, B = 1435.264, C = -64.848
        let logPbar = A - B / (max(tK, 50) + C)   // P in bar
        return pow(10, logPbar) / 10.0            // bar → MPa
    }
}
