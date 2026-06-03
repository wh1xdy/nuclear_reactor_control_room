// DecayHeat.swift — Fission product decay heat.
// Step 1: simplified ANS 5.1 power-law approximation.
// Step 4 upgrade path: replace with full 23-group ANS 5.1-2014 model.

import Foundation

struct DecayHeat {
    private var operatingSeconds: Double = 0   // cumulative n-weighted time
    private(set) var fraction: Double = 0      // decay heat / nominal power

    mutating func step(dt: Double, n: Double) -> Double {
        operatingSeconds += dt * max(n, 0)
        // ANS 5.1 simplified: P_d/P_0 ≈ 0.066 · t^{−0.2}  for t > 10 s
        // At t = 0 ramp linearly to avoid discontinuity
        let t = operatingSeconds
        let f: Double
        if t < 10 {
            f = 0.066 * pow(max(t, 0.1), -0.2) * (t / 10.0)
        } else {
            f = 0.066 * pow(t, -0.2)
        }
        fraction = min(f, 0.08)       // physical ceiling ~7% at shutdown
        return fraction
    }
}
