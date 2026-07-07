// Containment.swift — lumped containment atmosphere model.
//
// Two-state model (atmosphere temperature + free-steam inventory) with the
// pieces an operator actually watches: pressure, temperature, sump level,
// and the mitigation systems — fan coolers (AC-powered, always trying),
// passive wall heat sink, and the containment SPRAY that actuates on HI-HI
// pressure and knocks the steam down. Mass/energy sources: a primary-side
// break (the malfunction menu's LOCA), and PORV/relief discharge (the relief
// tank bursts its disc on sustained lift and vents to containment).
//
// Calibration is design-basis-flavoured rather than exact: atmospheric
// 101 kPa / 310 K normally; a 50 kg/s small LOCA reaches the 135 kPa HI-HI
// spray setpoint in a few minutes and the spray caps it well under the
// 400 kPa design pressure.

import Foundation

final class Containment {
    // ── States ──────────────────────────────────────────────────────────────
    private(set) var tempK: Double = 310
    private(set) var steamKg: Double = 0          // free steam in the atmosphere
    private(set) var sumpM3: Double = 0           // collected water
    private(set) var sprayOn = false

    // ── Constants ───────────────────────────────────────────────────────────
    private let heatCapJK   = 4.0e8               // lumped atmosphere+structures near-field
    private let flashFrac   = 0.45                // hot primary water flashing to steam
    private let hFlash      = 1.6e6               // J/kg energy carried per kg released
    private let sprayCond   = 25.0                // kg/s steam condensed by spray
    private let sprayCoolW  = 3.0e7               // spray heat removal [W]

    /// Absolute pressure [kPa]: air partial pressure (ideal-gas with T) plus
    /// the steam partial pressure from the free-steam inventory.
    var pressureKPa: Double { 101.0 * (tempK / 310.0) + steamKg * 0.012 }

    /// Advance by dt. `releaseKgs` = primary mass entering containment,
    /// `acPowered` gates the fan coolers, `esfAvailable` gates the spray pumps.
    func step(dt: Double, releaseKgs: Double, acPowered: Bool, esfAvailable: Bool) {
        // Break/relief inflow: part flashes to steam, the rest rains to the sump.
        let flashed = releaseKgs * flashFrac * dt
        steamKg += flashed
        sumpM3  += releaseKgs * (1 - flashFrac) * dt / 1000 * 1.0   // m³ (ρ≈1000)
        let eIn = releaseKgs * hFlash * dt

        // Spray: actuates at HI-HI, resets with margin (latched behaviour is
        // the supervisor's alarm board job; the pumps themselves cycle).
        if esfAvailable, pressureKPa > 135 { sprayOn = true }
        if pressureKPa < 110 || !esfAvailable { sprayOn = false }

        // Heat removal: fan coolers (AC), passive walls, spray.
        let dT = max(0, tempK - 310)
        var qOut = 2.0e6 * dT / 50.0                       // walls
        if acPowered { qOut += 4.0e6 * dT / 50.0 }         // fan coolers
        if sprayOn   { qOut += sprayCoolW }

        tempK += (eIn - qOut * dt) / heatCapJK
        tempK = max(305, min(450, tempK))

        // Steam condensation: coolers/walls take some; spray takes a lot.
        var cond = (0.4 + (acPowered ? 0.8 : 0)) * dT / 50.0 * dt
        if sprayOn { cond += sprayCond * dt }
        cond = min(cond, steamKg)
        steamKg -= cond
        sumpM3  += cond / 1000
    }
}
