// AxialCore.swift — 1-D axial nodal flux + xenon model.
//
// Sits ON TOP of the (calibrated, tested) point-kinetics. The TOTAL core power
// still comes from PointKinetics; this module solves the axial SHAPE and the
// per-node xenon so that axial offset (ΔI), peaking (Fq / FΔH), peak clad
// temperature and DNBR are REAL — read off the flux distribution — and axial
// xenon oscillations emerge from the flux↔xenon coupling.
//
// It is DRIVEN by the global state (power, rod insertion, inlet/outlet temps,
// flow, pressure) and does NOT feed back into the global reactivity, so the
// validated 0-D behaviour (and the 20 physics tests) are unchanged; only the
// shape-dependent limits are new.
//
// Shape = fundamental mode of a 1-D one-group diffusion operator on N nodes with
// Dirichlet-ish boundaries (→ chopped cosine, Fz≈1.5 at BOL), tilted by the
// per-node reactivity DEVIATION: rods insert from the top; local moderator +
// Doppler feedback; per-node xenon. Solved by warm-started power iteration each
// substep — the flux is quasi-static relative to the xenon / thermal timescales.

import Foundation

final class AxialCore {
    let params: PlantParams
    let N: Int

    private(set) var phi: [Double]   // relative axial power, mean = 1 (chopped cosine at BOL)
    private var iod: [Double]         // per-node I-135
    private var xen: [Double]         // per-node Xe-135

    // ── Shape solver ────────────────────────────────────────────────────────
    private let coupling: Double = 1.0    // inter-node diffusion coupling D
    private let relax:    Double = 0.08   // power-iteration step (relax·4D = 0.32 < 0.5 → stable)
    private let sweeps           = 6      // warm-started iterations per substep

    // ── Shape tilt ──────────────────────────────────────────────────────────
    // Per-node reactivity DEVIATION (Δk/k) using the REAL coefficients, so all
    // effects share one unit system (the per-node xenon inventory is O(1e3), so
    // it MUST enter as reactivity, not raw inventory — else it collapses the
    // shape). A single gain converts local reactivity → flux tilt; xeMult lets
    // the axial-xenon feedback be amplified toward the realistic marginal
    // instability of a tall PWR core (a single-mode model under-predicts it).
    private let shapeGain: Double = 8.0   // Δk/k local reactivity → flux tilt
    private let xeMult:    Double = 1.6   // axial-xenon feedback amplification (>1 = less stable)

    // ── DNBR / clad / peaking calibration ──────────────────────────────────
    private let chfRef:    Double = 2.66  // → nominal min-DNBR ≈ 2.0
    private let filmDT:    Double = 78.0  // clad film ΔT at unit local heat flux [K]
    private let radialNom: Double = 1.40  // nominal radial peaking (this model is axial-only)

    init(_ params: PlantParams, nodes: Int = 15) {
        self.params = params
        self.N = nodes
        // Dirichlet fundamental (chopped cosine), normalised to mean 1.
        var p = (0..<nodes).map { sin(Double.pi * (Double($0) + 1) / Double(nodes + 1)) }
        let m = p.reduce(0, +) / Double(nodes)
        for i in 0..<nodes { p[i] /= m }
        phi = p
        iod = [Double](repeating: 0, count: nodes)
        xen = [Double](repeating: 0, count: nodes)
    }

    /// Fraction of node i (0 = bottom) covered by control rods inserted `rod`
    /// (0…1) from the TOP of the core.
    private func rodded(_ i: Int, _ rod: Double) -> Double {
        let top = 1.0 - Double(i + 1) / Double(N)   // distance-from-top of node's lower edge
        let bot = 1.0 - Double(i)     / Double(N)    // distance-from-top of node's upper edge
        let covered = (min(rod, bot) - min(rod, top)) / (bot - top)
        return max(0, min(1, covered))
    }

    /// Local moderator temperature profile for the current shape (inlet at the
    /// bottom rising to outlet at the top by cumulative power).
    private func moderatorTemps(inlet: Double, outlet: Double) -> [Double] {
        let dT = outlet - inlet
        var t = [Double](repeating: 0, count: N)
        var cum = 0.0
        for i in 0..<N {
            let frac = (cum + 0.5 * phi[i]) / Double(N)
            cum += phi[i]
            t[i] = inlet + dT * frac
        }
        return t
    }

    func step(dt: Double, power n: Double, rodPos: Double,
              inlet: Double, outlet: Double, flow: Double) {
        let p = params
        let nPos = max(n, 0)
        let Tm = moderatorTemps(inlet: inlet, outlet: outlet)
        var Tf = [Double](repeating: 0, count: N)
        for i in 0..<N {
            Tf[i] = Tm[i] + (p.nominalFuelTemp - p.nominalCoolantTemp) * phi[i] * nPos
        }
        let TmBar = Tm.reduce(0, +) / Double(N)
        let TfBar = Tf.reduce(0, +) / Double(N)
        let XBar  = xen.reduce(0, +) / Double(N)

        // Per-node reactivity DEVIATION (Δk/k) → shape tilt (mean removed so it
        // only redistributes flux, leaving global criticality to point-kinetics).
        var tilt = [Double](repeating: 0, count: N)
        for i in 0..<N {
            let rho = p.rodWorth           * rodded(i, rodPos)
                    + p.coolantTempCoeff   * (Tm[i]  - TmBar)
                    + p.fuelTempCoeff      * (Tf[i]  - TfBar)
                    - xeMult * p.xenonReactivityCoeff * (xen[i] - XBar)
            tilt[i] = shapeGain * rho
        }
        let tiltBar = tilt.reduce(0, +) / Double(N)
        for i in 0..<N { tilt[i] -= tiltBar }

        // Quasi-static shape: warm-started power iteration on D·∇² + tilt.
        for _ in 0..<sweeps {
            var np = phi
            for i in 0..<N {
                let left  = i > 0     ? phi[i - 1] : 0.0   // Dirichlet ghost = 0
                let right = i < N - 1 ? phi[i + 1] : 0.0
                let lap = left - 2 * phi[i] + right
                np[i] = phi[i] + relax * (coupling * lap + tilt[i] * phi[i])
                if np[i] < 1e-4 { np[i] = 1e-4 }
            }
            let m = np.reduce(0, +) / Double(N)
            for i in 0..<N { phi[i] = np[i] / m }
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

    /// Axial offset ΔI = (P_top − P_bottom)/(P_top + P_bottom) [%]. Rods (from
    /// the top) drive it negative (bottom-peaked).
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
    /// Enthalpy-rise hot-channel factor.
    var fdh: Double { (1.0 + (fz - 1.0) * 0.62) * radialNom }

    /// Peak clad surface temperature [K] — hottest of (local coolant + film ΔT).
    func peakCladTempK(inlet: Double, outlet: Double, power n: Double) -> Double {
        let Tm = moderatorTemps(inlet: inlet, outlet: outlet)
        let nPos = max(n, 0)
        var peak = 0.0
        for i in 0..<N { peak = max(peak, Tm[i] + filmDT * phi[i] * nPos) }
        return peak
    }

    /// Minimum DNBR across the core — simplified W-3-style CHF ratio on local
    /// conditions (pressure, mass flux, subcooling) vs local heat flux.
    func minDNBR(power n: Double, flow: Double, pressureMPa: Double,
                 inlet: Double, outlet: Double) -> Double {
        let Tm = moderatorTemps(inlet: inlet, outlet: outlet)
        let nPos = max(n, 0.001)
        let pRat = max(pressureMPa, 1) / 15.5
        let tSat = 618.0 * pow(pRat, 0.05)                 // ≈ Tsat at the primary pressure
        let gFac = pow(max(flow, 0.05), 0.5) * pow(pRat, 0.3)
        var minR = 9.99
        for i in 0..<N {
            let qLoc = phi[i] * nPos
            let subcool = max(0, tSat - Tm[i])
            let chf = chfRef * gFac * (0.6 + 0.008 * subcool)
            minR = min(minR, chf / max(qLoc, 0.001))
        }
        return min(9.99, minR)
    }
}
