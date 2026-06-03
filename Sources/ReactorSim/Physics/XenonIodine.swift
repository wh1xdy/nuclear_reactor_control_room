// XenonIodine.swift — Xe-135 / I-135 poison model.
// Port of XenonIodine in physics_engine.py.
// Euler integration (xenon dynamics are slow; Euler is accurate enough at dt ≤ 0.01 s).

import Foundation

final class XenonIodine {
    let params: PlantParams

    var X: Double   // Xe-135 inventory (dimensionless, relative to equilibrium)
    var I: Double   // I-135 inventory

    init(_ params: PlantParams) {
        self.params = params
        // Start xenon-free (fresh reload). Xenon builds naturally over 6-8 sim-hours.
        // Operator learns to withdraw rods as xenon accumulates — the training objective.
        X = 0
        I = 0
    }

    func step(dt: Double, n: Double) {
        let p  = params
        let dI = p.gammaXe * n - p.lambdaI * I
        let dX = p.lambdaI * I - (p.lambdaXe + p.xenonBurnCoeff * n) * X
        I += dI * dt
        X += dX * dt
        if I < 0 { I = 0 }
        if X < 0 { X = 0 }
    }

    /// Xenon reactivity contribution [Δk/k] — always negative.
    var reactivity: Double {
        -params.xenonReactivityCoeff * X
    }
}
