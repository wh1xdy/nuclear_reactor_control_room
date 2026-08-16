// XenonIodine.swift — Fission-product poison model (Xe-135 / I-135 + Sm-149 / Pm-149).
// Port of XenonIodine in physics_engine.py, extended with the samarium chain.
// Euler integration (both chains are slow; Euler is accurate at dt ≤ 0.05 s).

import Foundation

final class XenonIodine {
    let params: PlantParams

    var X: Double   // Xe-135 inventory (dimensionless, relative to equilibrium)
    var I: Double   // I-135 inventory
    var Sm: Double  // Sm-149 inventory (dimensionless)
    var Pm: Double  // Pm-149 inventory

    init(_ params: PlantParams) {
        self.params = params
        // Start poison-free (fresh reload). Both chains build over sim-hours;
        // learning to withdraw rods / dilute as they accumulate is the training
        // objective. Xenon peaks in ~6–8 h; samarium climbs monotonically toward
        // its (flux-independent) equilibrium and then a permanent post-shutdown level.
        X = 0
        I = 0
        Sm = 0
        Pm = 0
    }

    func step(dt: Double, n: Double) {
        let p  = params
        // Iodine / xenon
        let dI = p.gammaXe * n - p.lambdaI * I
        let dX = p.lambdaI * I - (p.lambdaXe + p.xenonBurnCoeff * n) * X
        I += dI * dt
        X += dX * dt
        if I < 0 { I = 0 }
        if X < 0 { X = 0 }

        // Promethium / samarium. Pm is produced by fission and β-decays to Sm;
        // Sm is stable and disappears ONLY by neutron capture (∝ n). At shutdown
        // (n=0) the burn term vanishes and the Pm bank drains entirely into Sm —
        // the permanent "dead samarium" that never decays back out.
        let dPm = p.pmYield * n - p.lambdaPm * Pm
        let dSm = p.lambdaPm * Pm - p.smBurnCoeff * n * Sm
        Pm += dPm * dt
        Sm += dSm * dt
        if Pm < 0 { Pm = 0 }
        if Sm < 0 { Sm = 0 }
    }

    /// Total fission-product poison reactivity [Δk/k] — always negative.
    var reactivity: Double {
        -params.xenonReactivityCoeff * X - params.smReactivityCoeff * Sm
    }
}
