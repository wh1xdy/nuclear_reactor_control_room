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

    init(_ params: PlantParams) {
        self.params = params
        // Steady-state precursor concentrations: C_i = β_i / (λ_i · Λ)
        C = zip(params.beta, params.lambdaD).map { b, l in b / (l * params.lambdaPrompt) }
        let n6 = params.beta.count
        C0   = [Double](repeating: 0, count: n6)
        Cmid = [Double](repeating: 0, count: n6)
        dC1  = [Double](repeating: 0, count: n6)
        dC2  = [Double](repeating: 0, count: n6)
    }

    func step(dt: Double, rho: Double) {
        guard dt > 0 else { return }
        let p    = params
        let lam  = p.lambdaD
        let bet  = p.beta
        let Lam  = p.lambdaPrompt
        let bT   = p.betaTotal
        let rbl  = (rho - bT) / Lam      // (ρ − β) / Λ — loop invariant
        let iLam = 1.0 / Lam

        // β_i / Λ — precomputed per-group
        let bL = bet.map { $0 * iLam }

        // Adaptive sub-stepping: dt_safe = 0.3 · n / |dn/dt|
        let nSafe = max(n, 1e-6)
        let ps    = (0..<6).reduce(0.0) { $0 + lam[$1] * C[$1] }
        let dm    = abs(rbl * nSafe) + abs(ps)
        let nSub: Int
        if n < 1e-4 {
            nSub = 1
        } else {
            let dtSafe = 0.3 * nSafe / max(dm, 1e-10)
            // Divide by min(dt, dtSafe) — when dtSafe < dt we need MORE sub-steps, not fewer
            nSub = min(200, max(1, Int(ceil(dt / min(dt, dtSafe)))))
        }
        let subDt = dt / Double(nSub)
        let hdt   = 0.5 * subDt

        for _ in 0..<nSub {
            let n0 = n
            for i in 0..<6 { C0[i] = C[i] }

            // Stage 1
            let dn1 = rbl * n0 + (0..<6).reduce(0.0) { $0 + lam[$1] * C0[$1] }
            for i in 0..<6 { dC1[i] = bL[i] * n0 - lam[i] * C0[i] }
            let nMid = n0 + hdt * dn1
            for i in 0..<6 { Cmid[i] = C0[i] + hdt * dC1[i] }

            // Stage 2
            let dn2 = rbl * nMid + (0..<6).reduce(0.0) { $0 + lam[$1] * Cmid[$1] }
            for i in 0..<6 { dC2[i] = bL[i] * nMid - lam[i] * Cmid[i] }
            n += subDt * dn2
            for i in 0..<6 { C[i] += subDt * dC2[i] }
        }

        if n < 0 { n = 0 }
    }
}
