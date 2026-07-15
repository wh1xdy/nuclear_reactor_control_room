// PointKinetics.swift — Six-group point kinetics, RK2 (midpoint) integrator.
// Faithful Swift port of PointKinetics in physics_engine.py.
// Hot path: adaptive sub-stepping, all temporaries pre-allocated.

import Foundation

final class PointKinetics {
    let params: PlantParams

    var n: Double = 1.0                       // relative neutron population
    var C: [Double]                           // delayed precursor concentrations

    // Pre-allocated working arrays — avoid heap allocation in the hot loop
    private var C0:   [Double]
    private var Cmid: [Double]
    private var dC1:  [Double]
    private var dC2:  [Double]
    private let bL:   [Double]   // β_i / Λ — constant, precomputed once
    private let maxLam: Double   // fastest precursor group decay constant

    init(_ params: PlantParams) {
        self.params = params
        // Steady-state precursor concentrations: C_i = β_i / (λ_i · Λ)
        C = zip(params.beta, params.lambdaD).map { b, l in b / (l * params.lambdaPrompt) }
        let n6 = params.beta.count
        C0   = [Double](repeating: 0, count: n6)
        Cmid = [Double](repeating: 0, count: n6)
        dC1  = [Double](repeating: 0, count: n6)
        dC2  = [Double](repeating: 0, count: n6)
        bL   = params.beta.map { $0 / params.lambdaPrompt }
        maxLam = params.lambdaD.max() ?? 3.01
    }

    /// Force the flux to a given relative level with precursors in equilibrium
    /// for it (used to boot a cold, subcritical hot-standby state).
    func setLevel(_ level: Double) {
        n = max(0, level)
        for i in 0..<C.count {
            C[i] = params.beta[i] / (params.lambdaD[i] * params.lambdaPrompt) * n
        }
    }

    func step(dt: Double, rho: Double) {
        guard dt > 0 else { return }
        let p    = params
        let lam  = p.lambdaD
        let Lam  = p.lambdaPrompt
        let bT   = p.betaTotal
        let rbl  = (rho - bT) / Lam      // (ρ − β) / Λ — loop invariant
        let S    = p.neutronSource       // fixed source (source-range floor)

        // Deep-shutdown fast path: prompt-jump approximation.
        // With |ρ−β|/Λ ≈ 9000/s the prompt term is ultra-stiff, but n is then
        // slaved to the precursor source (n ≈ Σλ_iC_i · Λ/(β−ρ)) on a µs
        // timescale. Integrating it explicitly costs ~900 substeps per call;
        // the quasi-static solution is exact to O(µs) and O(1) per step.
        if rbl < 0, n < 1e-3, -rbl * dt > 20 {
            let nQS    = max(1, Int(ceil(dt * maxLam / 0.3)))
            let qsDt   = dt / Double(nQS)
            for _ in 0..<nQS {
                let src = (0..<6).reduce(0.0) { $0 + lam[$1] * C[$1] }
                n = (src + S) / -rbl
                for i in 0..<6 {
                    let e = exp(-lam[i] * qsDt)
                    C[i] = C[i] * e + bL[i] * n / lam[i] * (1.0 - e)
                }
            }
            return
        }

        // Adaptive sub-stepping.
        // dtStab: RK2 stability bound for the stiff prompt term — |rbl|·subDt
        // must stay ≲ 0.5 or a deep-scram step (rbl ≈ −9000/s) diverges to NaN.
        // dtSafe: accuracy bound, 30% change in n per substep.
        let nSafe  = max(n, 1e-6)
        let ps     = (0..<6).reduce(0.0) { $0 + lam[$1] * C[$1] }
        let dm     = abs(rbl * nSafe) + abs(ps)
        let dtStab = 0.5 / max(abs(rbl), 1e-10)
        let dtSafe: Double
        if n < 1e-4 {
            dtSafe = dtStab                       // power negligible: stability only
        } else {
            dtSafe = min(0.3 * nSafe / max(dm, 1e-10), dtStab)
        }
        // Divide by min(dt, dtSafe) — when dtSafe < dt we need MORE sub-steps, not fewer
        let nSub  = min(5000, max(1, Int(ceil(dt / min(dt, dtSafe)))))
        let subDt = dt / Double(nSub)
        let hdt   = 0.5 * subDt

        for _ in 0..<nSub {
            let n0 = n
            for i in 0..<6 { C0[i] = C[i] }

            // Stage 1
            let dn1 = rbl * n0 + (0..<6).reduce(0.0) { $0 + lam[$1] * C0[$1] } + S
            for i in 0..<6 { dC1[i] = bL[i] * n0 - lam[i] * C0[i] }
            let nMid = n0 + hdt * dn1
            for i in 0..<6 { Cmid[i] = C0[i] + hdt * dC1[i] }

            // Stage 2
            let dn2 = rbl * nMid + (0..<6).reduce(0.0) { $0 + lam[$1] * Cmid[$1] } + S
            for i in 0..<6 { dC2[i] = bL[i] * nMid - lam[i] * Cmid[i] }
            n += subDt * dn2
            for i in 0..<6 { C[i] += subDt * dC2[i] }
            if n < 0 { n = 0 }   // clamp inside the loop — a negative excursion
                                 // must not feed the next substep
        }
    }
}
