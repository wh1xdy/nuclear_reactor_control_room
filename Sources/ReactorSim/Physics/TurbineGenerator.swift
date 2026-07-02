// TurbineGenerator.swift — turbine-generator electrical + mechanical supervision.
//
// Upgrades the mimic's TURBINE-GENERATOR / ELECTRICAL readouts from "static
// numbers scaled by load" to small real models:
//
//  • Synchronous machine on an infinite bus (V = 1.0 pu, Xs = 1.8 pu on the
//    1100 MVA machine base). Active power comes from the turbine; the AVR is a
//    first-order excitation loop driving reactive power to its programmed
//    setpoint (≈0.92 pf at rated load). From the machine equations
//        P = E·V·sinδ / Xs        Q = (E·V·cosδ − V²) / Xs
//    we get real MVAr, load angle, stator current, and field V/A that respond
//    to load changes and de-excite on a trip.
//
//  • First-order thermal lags for bearing-metal / stator / lube-oil / H₂-gas
//    temperatures (τ = minutes) — gauges now trail load moves realistically
//    instead of teleporting.
//
//  • Shaft speed with a real coastdown on trip (τ ≈ 8 min) and a vibration
//    bump while coasting through the rotor criticals.
//
// Still deliberately simple — no rotor dynamics, no saturation — but every
// displayed number is now a STATE with correct qualitative behaviour.

import Foundation

final class TurbineGenerator {
    // ── Machine constants (per-unit on the machine base) ────────────────────
    // The MVA base scales with the plant rating (990 MWe PWR → 1100 MVA;
    // a 66 MWe SMR gets a proportionally small machine, not an idling giant).
    private var sBaseMVA  = 1100.0   // machine rating [MVA]
    private let vBus      = 1.0      // infinite-bus / terminal voltage [pu]
    private let xs        = 1.8      // synchronous reactance [pu]
    private let pfSet     = 0.92     // AVR reactive program: Q = P·tan(acos(pf))
    private let tauExc    = 1.5      // excitation (AVR) time constant [s]
    private let fieldAperE = 1620.0  // field amps per unit internal emf
    private let fieldOhms  = 0.105   // effective field resistance [V/A]

    // ── States ──────────────────────────────────────────────────────────────
    private(set) var emf: Double = 1.0        // internal emf E [pu]
    private(set) var rpm: Double = 3000
    private(set) var brgMetalC: Double = 92   // bearing metal [°C]
    private(set) var statorC:  Double = 97    // stator winding [°C]
    private(set) var lubeC:    Double = 52    // lube oil [°C]
    private(set) var h2GasC:   Double = 46    // H₂ cold-gas [°C]
    private(set) var vibMil:   Double = 1.8   // shaft vibration [mil]

    // ── Derived electrical outputs (updated each step) ──────────────────────
    private(set) var mvar: Double = 426       // reactive power [MVAr]
    private(set) var statorKA: Double = 29.6  // stator current at 21 kV [kA]
    private(set) var loadAngleDeg: Double = 40
    private(set) var powerFactor: Double = 0.92
    var fieldA: Double { emf * fieldAperE }
    var fieldV: Double { fieldA * fieldOhms }

    func step(dt: Double, grossMWe: Double, tripped: Bool, ratedMWe: Double = 990) {
        sBaseMVA = max(1, ratedMWe) / pfSet   // rating at the design power factor
        let rated = max(1, ratedMWe)
        // ── Shaft speed: governed at 3000 while on the grid; exponential
        //    coastdown (τ ≈ 8 min) once tripped.
        if tripped {
            rpm += (0 - rpm) / 480.0 * dt
            if rpm < 1 { rpm = 0 }
        } else {
            rpm += (3000 - rpm) / 20.0 * dt   // resync/run-up is governor-fast
        }

        // ── Excitation / AVR ────────────────────────────────────────────────
        let pPU = tripped ? 0 : max(0, grossMWe) / sBaseMVA
        if tripped {
            // Field breaker opens → de-excite toward no-load emf.
            emf += (1.0 - emf) / (tauExc * 4) * dt
            mvar = 0; statorKA = 0; loadAngleDeg = 0; powerFactor = 0
        } else {
            // AVR demand: the emf that delivers Q_ref at the present P.
            let qRef = pPU * tan(acos(pfSet))
            let eD1 = qRef * xs / vBus + vBus       // E·cosδ demanded
            let eD2 = pPU * xs / vBus               // E·sinδ demanded
            let eDemand = (eD1 * eD1 + eD2 * eD2).squareRoot()
            emf += (eDemand - emf) / tauExc * dt

            // Machine equations with the ACTUAL emf (lags the demand).
            let sinD = min(0.99, pPU * xs / (emf * vBus))
            let cosD = (1 - sinD * sinD).squareRoot()
            loadAngleDeg = asin(sinD) * 180 / .pi
            let qPU = (emf * vBus * cosD - vBus * vBus) / xs
            mvar = qPU * sBaseMVA
            let sPU = (pPU * pPU + qPU * qPU).squareRoot()
            statorKA = sPU * sBaseMVA / (1.732 * 21.0)      // at the 21 kV terminals
            powerFactor = sPU > 1e-6 ? pPU / sPU : 1.0
        }

        // ── Thermal lags (first-order toward load-dependent targets) ────────
        let load = tripped ? 0 : min(1.2, max(0, grossMWe / rated))
        brgMetalC += ((60 + 32 * load) - brgMetalC) / 300 * dt
        statorC   += ((55 + 42 * load) - statorC)  / 600 * dt
        lubeC     += ((44 +  8 * load) - lubeC)    / 400 * dt
        h2GasC    += ((40 +  6 * load) - h2GasC)   / 900 * dt

        // ── Vibration: load-dependent baseline; bump through the rotor
        //    criticals (~1400–1900 rpm) during a coastdown, then quiet.
        let inCriticals = rpm > 1200 && rpm < 2000 && tripped
        let vibTarget = rpm < 10 ? 0.05 : (inCriticals ? 3.4 : 1.1 + 0.7 * load)
        vibMil += (vibTarget - vibMil) / 25 * dt
    }
}
