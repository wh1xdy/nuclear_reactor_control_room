// DecayHeat.swift — Fission-product decay heat, full ANS/ANSI-5.1 23-group model.
//
// Each group i obeys  dh_i/dt = α_i·F(t) − λ_i·h_i  where F(t) is the fission
// rate normalized to rated power (here: kinetics.n). Decay power fraction is
// Σ h_i / Q with Q = 200 MeV/fission (ANS normalization convention).
//
// Groups are advanced with the EXACT exponential solution for a constant F over
// the substep — unconditionally stable at any dt (incl. 600× time compression):
//   h_i(t+dt) = h_i e^{−λ_i dt} + (α_i/λ_i)·F·(1 − e^{−λ_i dt})
//
// Behavior this restores vs. the old power-law stub:
//  • decay heat BUILDS UP toward ~6.5% of operating power during operation
//  • after shutdown it DECAYS along the ANS curve (the stub froze it forever)

import Foundation

struct DecayHeat {
    // ANS-5.1 U-235 thermal-fission 23-group fit: α [MeV/(fission·s)], λ [1/s]
    static let alpha: [Double] = [
        6.5057e-01, 5.1264e-01, 2.4384e-01, 1.3850e-01, 5.5440e-02,
        2.2225e-02, 3.3088e-03, 9.3015e-04, 8.0943e-04, 1.9567e-04,
        3.2535e-05, 7.5595e-06, 2.5232e-06, 4.9948e-07, 1.8531e-07,
        2.6608e-08, 2.2398e-09, 8.1641e-12, 8.7797e-11, 2.5131e-14,
        3.2176e-16, 4.5038e-17, 7.4791e-17,
    ]
    static let lambda: [Double] = [
        2.2138e+01, 5.1587e-01, 1.9594e-01, 1.0314e-01, 3.3656e-02,
        1.1681e-02, 3.5870e-03, 1.3930e-03, 6.2630e-04, 1.8906e-04,
        5.4988e-05, 2.0958e-05, 1.0010e-05, 2.5438e-06, 6.6361e-07,
        1.2290e-07, 2.7213e-08, 4.3714e-09, 7.5780e-10, 2.4786e-10,
        2.2384e-13, 2.4600e-14, 1.5699e-14,
    ]

    /// Recoverable energy per fission [MeV] used for normalization.
    static let qFission: Double = 200.0

    /// Decay-heat fraction at infinite-operation equilibrium: Σ(α_i/λ_i)/Q ≈ 0.0654.
    /// Used by the plant to split rated power into fission + decay shares.
    static let equilibriumFraction: Double =
        zip(alpha, lambda).reduce(0.0) { $0 + $1.0 / $1.1 } / qFission

    // Group energies h_i [MeV/fission-normalized]. Fresh core: all zero —
    // decay heat builds naturally over the first hours of operation.
    private var h: [Double] = Array(repeating: 0, count: 23)
    private(set) var fraction: Double = 0

    /// Advance by dt seconds at normalized fission power n; returns the
    /// decay-heat fraction of rated power.
    mutating func step(dt: Double, n: Double) -> Double {
        let f = max(n, 0)
        var total = 0.0
        for i in 0..<23 {
            let lam = Self.lambda[i]
            let e   = exp(-lam * dt)
            h[i] = h[i] * e + Self.alpha[i] / lam * f * (1.0 - e)
            total += h[i]
        }
        fraction = total / Self.qFission
        return fraction
    }
}
