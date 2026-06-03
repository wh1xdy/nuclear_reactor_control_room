// ThermalHydraulics.swift — Lumped-parameter primary + secondary thermal model.
// Two primary nodes (fuel, coolant) + one secondary node (steam generator).
// Port of ThermalHydraulics in physics_engine.py.

import Foundation

final class ThermalHydraulics {
    let params: PlantParams

    var tFuel: Double       // fuel temperature [K]
    var tCool: Double       // primary coolant temperature [K]
    var tSG:   Double       // steam generator temperature [K]
    var electricPower: Double

    init(_ params: PlantParams) {
        self.params = params
        tFuel  = params.nominalFuelTemp     // 900 K
        tCool  = params.nominalCoolantTemp  // 550 K
        tSG    = params.nominalSGTemp       // 490 K
        electricPower = params.nominalPower * params.turbineEfficiency
    }

    /// Advance by dt seconds given thermal power P_th [W] and control inputs.
    /// Returns electric power [W].
    @discardableResult
    func step(dt: Double, thermalPower: Double, flow: Double, turbineValve: Double) -> Double {
        let p    = params
        let hFC  = p.hFuelToCoolant * max(flow, 0)
        let hCS  = p.hCoolantToSG   * max(flow, 0)
        let hST  = p.hSGToTurbine   * max(turbineValve, 0)

        let qFC  = hFC * (tFuel - tCool)
        let qCS  = hCS * (tCool - tSG)
        let qST  = hST * (tSG - p.condenserTempK)         // condenser cold-side sink

        tFuel += (thermalPower - qFC) / p.fuelHeatCapacity * dt
        tCool += (qFC - qCS)          / p.coolantHeatCapacity * dt
        tSG   += (qCS - qST)          / p.sgHeatCapacity * dt

        let elec = max(qST, 0) * p.turbineEfficiency
        electricPower = elec
        return elec
    }
}
