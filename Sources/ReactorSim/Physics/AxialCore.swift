// AxialCore.swift — 1-D axial nodal flux + xenon model.
//
// Sits ON TOP of the (calibrated, tested) point-kinetics. The TOTAL core power
// still comes from PointKinetics; this module solves the axial SHAPE and the
// per-node xenon so that axial offset (ΔI), peaking (Fz / Fq), peak clad
// temperature and DNBR are REAL — read off the flux distribution — and axial
// xenon oscillations emerge from the flux↔xenon coupling. (Verified: a rod
// kick at equilibrium xenon produces a damped ΔI oscillation with a ~32 h
// period — textbook BOL-core behaviour.)
//
// It is DRIVEN by the global state (power, rod insertion, inlet/outlet temps,
// flow, pressure) and does NOT feed back into the global reactivity, so the
// validated 0-D behaviour (and its tests) are unchanged; only the
// shape-dependent limits are new.
//
// Reactor kinds: PWR/SMR rods enter from the TOP; BWR control blades enter
// from the BOTTOM (mirrored coverage), and the BWR additionally carries a
// bottom-peaking void-profile tilt (voids build up the channel, so moderation
// is better low in the core). voidCoeff = 0 for PWR/SMR keeps that term inert.
//
// Shape = fundamental mode of a 1-D one-group diffusion operator on N nodes
// with Dirichlet-ish boundaries (→ chopped cosine, Fz≈1.3 at nominal), tilted
// by the per-node reactivity DEVIATION in real Δk/k units. Solved by a
// warm-started power iteration each substep — the flux is quasi-static
// relative to the xenon / thermal timescales. All per-substep scratch arrays
// are preallocated (PointKinetics pattern): this runs 200×/frame at 600×.

import Foundation

final class AxialCore {
    let params: PlantParams
    let N: Int

    private(set) var phi: [Double]   // relative axial power, mean = 1 (chopped cosine at BOL)
    private var iod: [Double]         // per-node I-135
    private var xen: [Double]         // per-node Xe-135

    // Preallocated per-substep scratch — no heap churn on the 60 Hz hot path.
    private var tm:   [Double]        // moderator temperature profile
    private var tf:   [Double]        // fuel temperature profile
    private var tilt: [Double]        // per-node reactivity deviation → shape tilt
    private var phiN: [Double]        // next-iterate buffer for the power iteration

    // ── Shape solver ────────────────────────────────────────────────────────
    private let coupling: Double = 1.0    // inter-node diffusion coupling D
    private let relax:    Double = 0.08   // power-iteration step (relax·4D = 0.32 < 0.5 → stable)
    private let sweeps           = 6      // warm-started iterations per substep

    // ── Shape tilt ──────────────────────────────────────────────────────────
    // Per-node reactivity DEVIATION (Δk/k) using the REAL coefficients, so all
    // effects share one unit system (the per-node xenon inventory is O(1e3), so
    // it MUST enter as reactivity, not raw inventory — else it collapses the
    // shape). A single gain converts local reactivity → flux tilt; xeMult sets
    // the axial-xenon stability (1.6 → damped ~32 h oscillation; raise toward
    // ~2.5 for the marginal/divergent EOL-core behaviour).
    private let shapeGain: Double = 8.0   // Δk/k local reactivity → flux tilt
    private let xeMult:    Double = 1.6   // axial-xenon feedback amplification (>1 = less stable)

    // ── DNBR / clad / peaking calibration ──────────────────────────────────
    private let chfRef:    Double = 2.66  // → nominal min-DNBR ≈ 2.0
    private let filmDT:    Double = 78.0  // clad film ΔT at unit local heat flux, FULL flow [K]
    private let radialNom: Double = 1.40  // nominal radial peaking (this model is axial-only)

    init(_ params: PlantParams, nodes: Int = 15) {
        self.params = params
        self.N = nodes
        // Dirichlet fundamental (chopped cosine), normalised to mean 1.
        var p = (0..<nodes).map { sin(Double.pi * (Double($0) + 1) / Double(nodes + 1)) }
        let m = p.reduce(0, +) / Double(nodes)
        for i in 0..<nodes { p[i] /= m }
        phi  = p
        iod  = [Double](repeating: 0, count: nodes)
        xen  = [Double](repeating: 0, count: nodes)
        tm   = [Double](repeating: 0, count: nodes)
        tf   = [Double](repeating: 0, count: nodes)
        tilt = [Double](repeating: 0, count: nodes)
        phiN = [Double](repeating: 0, count: nodes)
    }

    /// Fraction of node i (0 = bottom) covered by control rods inserted `rod`
    /// (0…1). PWR/SMR CRDMs enter from the TOP; BWR blades from the BOTTOM.
    private func rodded(_ i: Int, _ rod: Double) -> Double {
        let j = params.kind == .bwr ? (N - 1 - i) : i   // mirror for bottom-entry blades
        let top = 1.0 - Double(j + 1) / Double(N)
        let bot = 1.0 - Double(j)     / Double(N)
        let covered = (min(rod, bot) - min(rod, top)) / (bot - top)
        return max(0, min(1, covered))
    }

    /// Local moderator temperature profile for the current shape (inlet at the
    /// bottom rising to outlet at the top by cumulative power), into `tm`.
    private func fillModeratorTemps(inlet: Double, outlet: Double) {
        let dT = outlet - inlet
        var cum = 0.0
        for i in 0..<N {
            let frac = (cum + 0.5 * phi[i]) / Double(N)
            cum += phi[i]
            tm[i] = inlet + dT * frac
        }
    }

    func step(dt: Double, power n: Double, rodPos: Double,
              inlet: Double, outlet: Double, flow: Double) {
        let p = params
        let nPos = max(n, 0)
        fillModeratorTemps(inlet: inlet, outlet: outlet)
        for i in 0..<N {
            tf[i] = tm[i] + (p.nominalFuelTemp - p.nominalCoolantTemp) * phi[i] * nPos
        }
        var tmBar = 0.0, tfBar = 0.0, xeBar = 0.0
        for i in 0..<N { tmBar += tm[i]; tfBar += tf[i]; xeBar += xen[i] }
        tmBar /= Double(N); tfBar /= Double(N); xeBar /= Double(N)

        // BWR: void fraction builds up the channel (cumulative boiling), so the
        // local void deviation is ≈ nominalVoid·(2·cumFrac − 1)·(power/flow).
        // With voidCoeff < 0 that suppresses the top and bottom-peaks the shape
        // — the defining BWR axial driver. Inert for PWR/SMR (voidCoeff = 0).
        // The 0.30 attenuation accounts for the flux-shape solver seeing the
        // FULL static void span otherwise (ΔI would pin at −50%); calibrated to
        // a realistic BWR bottom-peak of ΔI ≈ −15% and Fz ≈ 1.5 at rated.
        let voidScale = p.voidCoeff != 0
            ? 0.30 * p.nominalVoidFraction * nPos / max(flow, 0.1)
            : 0

        // Per-node reactivity DEVIATION (Δk/k) → shape tilt (mean removed so it
        // only redistributes flux, leaving global criticality to point-kinetics).
        var cum = 0.0
        for i in 0..<N {
            let cumFrac = (cum + 0.5 * phi[i]) / Double(N)
            cum += phi[i]
            let localVoidDev = voidScale * (2 * cumFrac - 1)
            let rho = p.rodWorth           * rodded(i, rodPos)
                    + p.coolantTempCoeff   * (tm[i]  - tmBar)
                    + p.fuelTempCoeff      * (tf[i]  - tfBar)
                    + p.voidCoeff          * localVoidDev
                    - xeMult * p.xenonReactivityCoeff * (xen[i] - xeBar)
            tilt[i] = shapeGain * rho
        }
        var tiltBar = 0.0
        for i in 0..<N { tiltBar += tilt[i] }
        tiltBar /= Double(N)
        for i in 0..<N { tilt[i] -= tiltBar }

        // Quasi-static shape: warm-started power iteration on D·∇² + tilt.
        for _ in 0..<sweeps {
            for i in 0..<N {
                let left  = i > 0     ? phi[i - 1] : 0.0   // Dirichlet ghost = 0
                let right = i < N - 1 ? phi[i + 1] : 0.0
                let lap = left - 2 * phi[i] + right
                phiN[i] = max(1e-4, phi[i] + relax * (coupling * lap + tilt[i] * phi[i]))
            }
            var m = 0.0
            for i in 0..<N { m += phiN[i] }
            m /= Double(N)
            for i in 0..<N { phi[i] = phiN[i] / m }
        }

        // Per-node xenon / iodine, driven by the LOCAL flux (phi·n).
        for i in 0..<N {
            let flux = phi[i] * nPos
            let dI = p.gammaXe * flux - p.lambdaI * iod[i]
            let dX = p.lambdaI * iod[i] - (p.lambdaXe + p.xenonBurnCoeff * flux) * xen[i]
            iod[i] += dI * dt
            xen[i] += dX * dt
            if iod[i] < 0 { iod[i] = 0 }
            if xen[i] < 0 { xen[i] = 0 }
        }
    }

    // ── Diagnostics ─────────────────────────────────────────────────────────

    /// Axial offset ΔI = (P_top − P_bottom)/(P_top + P_bottom) [%].
    var axialOffsetPct: Double {
        let half = N / 2
        var top = 0.0, bot = 0.0
        for i in 0..<N {
            if i >= N - half { top += phi[i] }
            else if i < half { bot += phi[i] }
        }
        let s = top + bot
        return s > 0 ? (top - bot) / s * 100 : 0
    }

    /// Axial peaking Fz = peak local / average.
    var fz: Double { phi.max() ?? 1 }
    /// Total heat-flux hot-channel factor (axial × nominal radial).
    var fq: Double { fz * radialNom }
    /// Enthalpy-rise hot-channel factor (heuristic mapping from Fz).
    var fdh: Double { (1.0 + (fz - 1.0) * 0.62) * radialNom }

    /// Peak clad surface temperature [K] — hottest of (local coolant + film ΔT).
    /// The film ΔT rises as flow drops (Dittus–Boelter, h ∝ flow^0.8), so a
    /// loss-of-flow now heats the clad instead of leaving it pinned to nominal.
    func peakCladTempK(inlet: Double, outlet: Double, power n: Double, flow: Double) -> Double {
        fillModeratorTemps(inlet: inlet, outlet: outlet)
        let nPos = max(n, 0)
        let film = filmDT / pow(max(flow, 0.08), 0.8)
        var peak = 0.0
        for i in 0..<N { peak = max(peak, tm[i] + film * phi[i] * nPos) }
        return peak
    }

    /// Minimum DNBR across the core — simplified W-3-style CHF ratio on local
    /// conditions (pressure, mass flux, subcooling) vs local heat flux. The
    /// saturation anchor uses the real inverse-Antoine T_sat at the LIVE
    /// primary pressure, so depressurisation erodes the margin correctly (and
    /// the anchor is right for BWR 7 MPa / SMR 13.8 MPa, not just PWR 15.5).
    func minDNBR(power n: Double, flow: Double, pressureMPa: Double,
                 inlet: Double, outlet: Double) -> Double {
        fillModeratorTemps(inlet: inlet, outlet: outlet)
        let nPos = max(n, 0.001)
        let pRat = max(pressureMPa, 1) / params.nominalPressureMPa
        let tSat = ThermalHydraulics.satTempK(pressureMPa)
        let gFac = pow(max(flow, 0.05), 0.5) * pow(pRat, 0.3)
        var minR = 9.99
        for i in 0..<N {
            let qLoc = phi[i] * nPos
            let subcool = max(0, tSat - tm[i])
            let chf = chfRef * gFac * (0.6 + 0.008 * subcool)
            minR = min(minR, chf / max(qLoc, 0.001))
        }
        return min(9.99, minR)
    }
}
