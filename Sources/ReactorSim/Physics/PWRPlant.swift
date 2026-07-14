// PWRPlant.swift — Top-level reactor plant model.
// Composes PointKinetics + XenonIodine + ThermalHydraulics + DecayHeat.
// Reactor-kind-agnostic: the same machinery models PWR, BWR (void feedback +
// direct cycle) and SMR (scaled, natural circulation); the differences live in
// the injected PlantParams (see PlantParams.bwr()/.smr()).

import Foundation

struct ControlInputs {
    var rodPosition: Double  = 0.0
    var primaryFlow: Double  = 1.0
    var turbineValve: Double = 1.0   // effective SG heat-removal valve (turbine or steam dump)
    var turbineTripped: Bool = false // generator off the grid → 0 MWe even if dump flows
    var scram: Bool          = false
    /// Live primary pressure from the supervisor's pressurizer model [MPa].
    /// Feeds the DNBR saturation anchor so depressurisation erodes the margin;
    /// nil (tests, standalone plant) falls back to the nominal pressure.
    var primaryPressureMPa: Double? = nil
}

struct PlantSnapshot {
    var time:               Double
    var powerFraction:      Double
    var thermalPowerW:      Double
    var electricPowerW:     Double
    var fuelTempK:          Double
    var coolantTempK:       Double          // RCS T-avg = (T_hot + T_cold)/2
    var sgTempK:            Double
    var reactivity:         Double
    var xenonInventory:     Double
    var iodineInventory:    Double
    var rodPosition:        Double
    var scrammed:           Bool
    var decayHeatFraction:  Double
    // Two-node primary detail + saturation steam pressure (defaults keep older
    // PlantSnapshot(...) constructions valid).
    var hotLegTempK:        Double = 565
    var coldLegTempK:       Double = 535
    var steamPressureMPa:   Double = 1.9
    /// Core-average void fraction (BWR feedback driver; 0 for PWR/SMR).
    var voidFraction:       Double = 0
    // Axial-nodal outputs (real, from AxialCore — no longer mimic proxies).
    var axialOffsetPct:     Double = 0
    var fz:                 Double = 1.5     // axial peaking
    var fq:                 Double = 2.1     // total heat-flux hot-channel factor
    var fdh:                Double = 1.55    // enthalpy-rise hot-channel factor
    var peakCladTempK:      Double = 620
    var minDNBR:            Double = 2.0
    var axialProfile:       [Double] = []    // per-node relative power (for an axial flux display)
    // Radial model outputs (real, from AxialCore's ring + azimuthal states).
    var radialRings:        [Double] = []    // centre / mid / periphery relative power
    var tiltX:              Double = 0       // azimuthal first-harmonic components
    var tiltY:              Double = 0
    var qptr:               Double = 1.0     // quadrant power tilt ratio (real)
}

/// Back-compatible name — the plant now models every reactor kind via params.
typealias PWRPlant = ReactorPlant

final class ReactorPlant {
    var params: PlantParams      // var so supervisor can tune boron/external reactivity
    private let kinetics:  PointKinetics
    private let xenon:     XenonIodine
    private let thermal:   ThermalHydraulics
    private var decayHeat: DecayHeat
    private let axial:     AxialCore

    private(set) var time: Double = 0
    private(set) var scrammed: Bool = false
    // Tracks gradual rod drop; private(set) so resetScram can write via the method
    private var rodPosEffective: Double = 0.0

    init(params: PlantParams = PlantParams()) {
        self.params = params
        kinetics  = PointKinetics(params)
        xenon     = XenonIodine(params)
        thermal   = ThermalHydraulics(params)
        decayHeat = DecayHeat()
        axial     = AxialCore(params)
    }

    private func computeReactivity(_ ctrl: ControlInputs) -> Double {
        let p = params
        let x = rodPosEffective
        let w = p.rodShape(x)            // calibrated table or integrated-cosine S-curve
        var rho = w * p.rodWorth
        if scrammed { rho += p.scramExtraWorth }
        rho += p.fuelTempCoeff    * (thermal.tFuel - p.nominalFuelTemp)
        rho += p.coolantTempCoeff * (thermal.tAvg  - p.nominalCoolantTemp)
        // Void feedback (BWR). Referenced to the nominal void so the steady
        // operating point is unchanged; only deviations contribute. Zero for
        // subcooled PWR/SMR (voidCoeff = 0).
        rho += p.voidCoeff        * (thermal.voidFraction - p.nominalVoidFraction)
        rho += xenon.reactivity
        rho += p.externalReactivity
        return max(-0.9, min(p.betaTotal * 1.5, rho))
    }

    func step(dt: Double, ctrl: ControlInputs) -> PlantSnapshot {
        let p = params
        if ctrl.scram { scrammed = true }

        let nSteps = max(1, Int(ceil(dt / p.internalDt)))
        let dtSub  = dt / Double(nSteps)

        // Rated power splits into fission + decay-heat shares so that at
        // equilibrium (n=1, decay≈6.5%) total thermal power equals nominal.
        let fissionShare = 1.0 - DecayHeat.equilibriumFraction

        for _ in 0..<nSteps {
            if scrammed {
                rodPosEffective = min(1.0, rodPosEffective + dtSub / p.rodDropTau)
            } else {
                // CRDM rate limit: rods WALK toward demand, never teleport.
                // (Also prevents a reactivity step when a scram is reset with
                // the demand lever left withdrawn.)
                let demand  = max(0, min(1, ctrl.rodPosition))
                let maxStep = p.rodSpeed * dtSub
                rodPosEffective += max(-maxStep, min(maxStep, demand - rodPosEffective))
            }
            let rho   = computeReactivity(ctrl)
            kinetics.step(dt: dtSub, rho: rho)
            xenon.step(dt: dtSub, n: kinetics.n)
            let decay = decayHeat.step(dt: dtSub, n: kinetics.n)
            let pTh   = (kinetics.n * fissionShare + decay) * p.nominalPower
            thermal.step(dt: dtSub, thermalPower: pTh, flow: ctrl.primaryFlow, turbineValve: ctrl.turbineValve)
            // Axial shape + per-node xenon (driven by the just-updated global
            // state; does not feed back into the 0-D reactivity above).
            axial.step(dt: dtSub, power: kinetics.n, rodPos: rodPosEffective,
                       inlet: thermal.tCold, outlet: thermal.tHot, flow: ctrl.primaryFlow)
            time += dtSub
        }

        return PlantSnapshot(
            time:              time,
            // Fraction of RATED THERMAL power (fission share + decay heat) —
            // reads 100% at steady full power regardless of decay-heat buildup.
            powerFraction:     kinetics.n * fissionShare + decayHeat.fraction,
            thermalPowerW:     (kinetics.n * fissionShare + decayHeat.fraction) * p.nominalPower,
            electricPowerW:    ctrl.turbineTripped ? 0 : thermal.electricPower,
            fuelTempK:         thermal.tFuel,
            coolantTempK:      thermal.tAvg,
            sgTempK:           thermal.tSG,
            reactivity:        computeReactivity(ctrl),
            xenonInventory:    xenon.X,
            iodineInventory:   xenon.I,
            rodPosition:       rodPosEffective,
            scrammed:          scrammed,
            decayHeatFraction: decayHeat.fraction,
            hotLegTempK:       thermal.tHot,
            coldLegTempK:      thermal.tCold,
            steamPressureMPa:  thermal.steamPressureMPa,
            voidFraction:      thermal.voidFraction,
            axialOffsetPct:    axial.axialOffsetPct,
            fz:                axial.fz,
            fq:                axial.fq,
            fdh:               axial.fdh,
            peakCladTempK:     axial.peakCladTempK(inlet: thermal.tCold, outlet: thermal.tHot,
                                                   power: kinetics.n * fissionShare + decayHeat.fraction,
                                                   flow: ctrl.primaryFlow),
            minDNBR:           axial.minDNBR(power: kinetics.n * fissionShare + decayHeat.fraction,
                                             flow: ctrl.primaryFlow,
                                             pressureMPa: ctrl.primaryPressureMPa ?? p.nominalPressureMPa,
                                             inlet: thermal.tCold, outlet: thermal.tHot),
            axialProfile:      axial.phi,
            radialRings:       axial.ring,
            tiltX:             axial.tiltX,
            tiltY:             axial.tiltY,
            qptr:              axial.qptr
        )
    }

    func resetScram() {
        scrammed = false
        // Leave rodPosEffective where it is (rods stay fully inserted until operator withdraws)
    }

    /// End-of-cycle core mode: strengthens the axial-xenon feedback toward the
    /// marginal/divergent EOL oscillation (the hard flyspeck drill).
    func setEndOfCycle(_ on: Bool) { axial.eolFactor = on ? 1.55 : 1.0 }

    /// Excite the azimuthal flux tilt (stuck / dropped rod malfunctions).
    func kickTilt(x: Double, y: Double) { axial.kickTilt(x: x, y: y) }
}
