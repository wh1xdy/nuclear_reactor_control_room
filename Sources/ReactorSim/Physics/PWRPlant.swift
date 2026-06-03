// PWRPlant.swift — Top-level PWR plant model.
// Composes PointKinetics + XenonIodine + ThermalHydraulics + DecayHeat.

import Foundation

struct ControlInputs {
    var rodPosition: Double  = 0.0
    var primaryFlow: Double  = 1.0
    var turbineValve: Double = 1.0
    var scram: Bool          = false
}

struct PlantSnapshot {
    var time:               Double
    var powerFraction:      Double
    var thermalPowerW:      Double
    var electricPowerW:     Double
    var fuelTempK:          Double
    var coolantTempK:       Double
    var sgTempK:            Double
    var reactivity:         Double
    var xenonInventory:     Double
    var iodineInventory:    Double
    var rodPosition:        Double
    var scrammed:           Bool
    var decayHeatFraction:  Double
}

final class PWRPlant {
    var params: PlantParams      // var so supervisor can tune boron/external reactivity
    private let kinetics:  PointKinetics
    private let xenon:     XenonIodine
    private let thermal:   ThermalHydraulics
    private var decayHeat: DecayHeat

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
    }

    private func computeReactivity(_ ctrl: ControlInputs) -> Double {
        let p = params
        let x = rodPosEffective
        let w = 3*x*x - 2*x*x*x          // S-curve: integrated cosine profile
        var rho = w * p.rodWorth
        if scrammed { rho += p.scramExtraWorth }
        rho += p.fuelTempCoeff    * (thermal.tFuel - p.nominalFuelTemp)
        rho += p.coolantTempCoeff * (thermal.tCool - p.nominalCoolantTemp)
        rho += xenon.reactivity
        rho += p.externalReactivity
        return max(-0.9, min(p.betaTotal * 1.5, rho))
    }

    func step(dt: Double, ctrl: ControlInputs) -> PlantSnapshot {
        let p = params
        if ctrl.scram { scrammed = true }

        let nSteps = max(1, Int(ceil(dt / p.internalDt)))
        let dtSub  = dt / Double(nSteps)

        for _ in 0..<nSteps {
            if scrammed {
                rodPosEffective = min(1.0, rodPosEffective + dtSub / p.rodDropTau)
            } else {
                rodPosEffective = max(0, min(1, ctrl.rodPosition))
            }
            let rho   = computeReactivity(ctrl)
            kinetics.step(dt: dtSub, rho: rho)
            xenon.step(dt: dtSub, n: kinetics.n)
            let decay = decayHeat.step(dt: dtSub, n: kinetics.n)
            let pTh   = (kinetics.n + decay) * p.nominalPower
            thermal.step(dt: dtSub, thermalPower: pTh, flow: ctrl.primaryFlow, turbineValve: ctrl.turbineValve)
            time += dtSub
        }

        return PlantSnapshot(
            time:              time,
            powerFraction:     kinetics.n,
            thermalPowerW:     (kinetics.n + decayHeat.fraction) * p.nominalPower,
            electricPowerW:    thermal.electricPower,
            fuelTempK:         thermal.tFuel,
            coolantTempK:      thermal.tCool,
            sgTempK:           thermal.tSG,
            reactivity:        computeReactivity(ctrl),
            xenonInventory:    xenon.X,
            iodineInventory:   xenon.I,
            rodPosition:       rodPosEffective,
            scrammed:          scrammed,
            decayHeatFraction: decayHeat.fraction
        )
    }

    func resetScram() {
        scrammed = false
        // Leave rodPosEffective where it is (rods stay fully inserted until operator withdraws)
    }
}
