// PlantSupervisor.swift — @Observable wrapper around PWRPlant.
// Adds protection system, BOP simulation, alarm management, and historian.
// Architecture note: runs on @MainActor for Step 1.
// Upgrade path (Step 3): move plant.step() to a background actor; publish
// snapshots via AsyncStream to keep the main thread free.

import Foundation
import Observation

struct ReactorAlarm: Identifiable {
    let id: String
    let message: String
    let priority: Int           // 1 = trip/emergency, 2 = warning, 3 = advisory
    var state: String           // "unack" | "ack"
    var isTrip: Bool
    var isFirstOut: Bool
    let simTime: Double
}

@Observable @MainActor
final class PlantSupervisor {

    // MARK: — Operator controls (UI writes)
    var rodPosition:    Double = 0.0
    var primaryFlow:    Double = 1.0
    var turbineValve:   Double = 1.0
    var feedwaterValve: Double = 0.7
    var borationRate:    Double = 0.0
    var dilutionRate:    Double = 0.0
    var scrammed:        Bool   = false
    // Fault / permit switches
    var startupPermit:   Bool   = true
    var turbineTrip:     Bool   = false
    var pumpDegraded:    Bool   = false
    var feedwaterFault:  Bool   = false

    // MARK: — Readable plant state (UI reads)
    private(set) var snapshot:     PlantSnapshot
    private(set) var pressureMPa:  Double = 15.5
    private(set) var boronPPM:     Double = 800.0
    private(set) var feedwaterInv: Double = 1.0
    private(set) var steamInv:     Double = 1.0
    private(set) var omegaRCP:     Double = 1.0
    private(set) var porvOpen:     Bool   = false
    private(set) var eccsActuated: Bool   = false
    private(set) var condTempK:    Double = 305.0
    private(set) var alarms:       [ReactorAlarm] = []
    private(set) var trips:        [String] = []
    private(set) var scramMessage: String = ""

    // Fast historian — 600 points (~10 s at 60 Hz)
    private(set) var histPower:  [Double]   // core power %
    private(set) var histReact:  [Double]   // reactivity Δk/k
    private(set) var histFuelT:  [Double]   // fuel temp K
    private(set) var histDecay:  [Double]   // decay heat %
    private(set) var histPress:  [Double]   // primary pressure MPa
    private(set) var histElec:   [Double]   // electric power MWe
    private(set) var histSteamT: [Double]   // SG temperature K
    private(set) var histCoolT:  [Double]   // coolant temperature K

    // MARK: — Private
    private let plant: PWRPlant
    private var _alarmMap:    [String: ReactorAlarm] = [:]
    private var _firstOutSet: Bool = false
    private var _histIdx:     Int  = 0
    private let _histLen:     Int  = 600
    // Xenon peak tracking for transient detection
    private var _xeMax:       Double = 0.0

    init() {
        let p = PlantParams()
        plant = PWRPlant(params: p)
        snapshot = PlantSnapshot(
            time: 0, powerFraction: 1, thermalPowerW: 3e9, electricPowerW: 990e6,
            fuelTempK: 900, coolantTempK: 550, sgTempK: 550,
            reactivity: 0, xenonInventory: 0, iodineInventory: 0,
            rodPosition: 0, scrammed: false, decayHeatFraction: 0
        )
        let n = 600
        histPower  = Array(repeating: 100.0,  count: n)
        histReact  = Array(repeating: 0.0,    count: n)
        histFuelT  = Array(repeating: 900.0,  count: n)
        histDecay  = Array(repeating: 0.0,    count: n)
        histPress  = Array(repeating: 15.5,   count: n)
        histElec   = Array(repeating: 990.0,  count: n)
        histSteamT = Array(repeating: 550.0,  count: n)
        histCoolT  = Array(repeating: 550.0,  count: n)
    }

    // MARK: — Step

    func step(dt: Double) {
        let ctrl = ControlInputs(
            rodPosition:  rodPosition,
            primaryFlow:  primaryFlow * omegaRCP,
            turbineValve: turbineValve,
            scram:        scrammed
        )
        snapshot = plant.step(dt: dt, ctrl: ctrl)
        scrammed = snapshot.scrammed

        updateBOP(dt: dt)
        updateBoron(dt: dt)
        updateAlarms()
        appendHistory()
    }

    // MARK: — Operator actions

    func triggerScram() {
        scrammed = true
    }

    func resetScram() {
        guard scrammed else { scramMessage = "Not scrammed"; return }
        guard snapshot.rodPosition > 0.94 else {
            scramMessage = "RESET BLOCKED — rods not fully inserted"; return
        }
        guard trips.isEmpty else {
            scramMessage = "RESET BLOCKED — active trips present"; return
        }
        scrammed = false
        plant.resetScram()
        scramMessage = "SCRAM RESET APPROVED"
    }

    func acknowledgeAllAlarms() {
        for key in _alarmMap.keys { _alarmMap[key]?.state = "ack" }
        _rebuildAlarmLists()
    }

    // MARK: — BOP

    private func updateBOP(dt: Double) {
        // RCP coast-down on pump fault
        if pumpDegraded {
            omegaRCP = max(0.04, omegaRCP - dt / 12.0)
        } else {
            omegaRCP = min(1.0, omegaRCP + dt / 30.0)
        }

        // Steam inventory: SG produces steam proportional to heat input, turbine removes it
        let steamProduction = snapshot.powerFraction * turbineValve * 0.06 * dt
        let steamConsumption = turbineTrip ? 0.0 : turbineValve * 0.08 * dt
        steamInv = max(0, min(2.0, steamInv + steamProduction - steamConsumption))

        // Condenser temp: rises with steam flow, drops toward cooling water temp
        let steamHeat  = turbineTrip ? 0.0 : turbineValve * snapshot.powerFraction * 5.0
        condTempK += (steamHeat - (condTempK - 305.0) * 0.8) * dt
        condTempK  = max(300.0, min(370.0, condTempK))

        // Feedwater balance (fault or turbine trip stops feedwater)
        let fwIn  = (feedwaterFault || turbineTrip) ? 0.0 : feedwaterValve * 0.08 * dt
        let fwOut = max(0, snapshot.powerFraction) * 0.07 * dt
        feedwaterInv = max(0, min(1.2, feedwaterInv + fwIn - fwOut))

        // Pressurizer: slow drift toward 15.5 MPa + small coolant-temp coupling
        // Coefficient kept small so normal ops don't trip HIGH_PRESS
        let pTarget = 15.5 + (snapshot.coolantTempK - 550.0) * 0.02
        pressureMPa += (pTarget - pressureMPa) * 0.02 * dt
        pressureMPa = max(10.0, min(17.5, pressureMPa))

        // PORV
        porvOpen = pressureMPa > 16.55

        // ECCS arms on low pressure
        if pressureMPa < 11.5 { eccsActuated = true }
        if pressureMPa > 13.0 && feedwaterInv > 0.3 { eccsActuated = false }
    }

    private func updateBoron(dt: Double) {
        let vPrimary = 300.0   // m³ approximate primary coolant volume
        let dB = (borationRate * 150.0 - dilutionRate * boronPPM) / vPrimary * dt
        boronPPM = max(0, boronPPM + dB)
        // Reactivity: −8 pcm/ppm, nominal 800 ppm → 0 reactivity
        plant.params.externalReactivity = -8e-5 * (boronPPM - 800.0)
    }

    // MARK: — Alarms

    private func updateAlarms() {
        let s = snapshot
        let nomP = 15.5

        let pTrip = nomP + 1.5   // 17.0 MPa trip (15.5 + 1.5), realistic PWR HIHI
        _check("HIGH_FLUX",   s.powerFraction > 1.20, 1, "HIGH NEUTRON FLUX — TRIP",       isTrip: true)
        _check("HIGH_FUEL_T", s.fuelTempK > 1500,     1, "HIGH FUEL TEMPERATURE — TRIP",   isTrip: true)
        _check("HIGH_PRESS",  pressureMPa > pTrip,    1, "HIGH REACTOR COOLANT PRESS TRIP", isTrip: true)
        _check("HIGH_COOL_T", s.coolantTempK > 620,   1, "HIGH COOLANT TEMPERATURE TRIP",   isTrip: true)
        _check("ECCS_ACT",    eccsActuated,            1, "ECCS ACTUATION",                  isTrip: true)
        _check("WARN_PRESS",  pressureMPa > nomP+0.8, 2, "HIGH PRESSURE WARNING",            isTrip: false)
        _check("PORV_OPEN",   porvOpen,                2, "PRESSURIZER PORV OPEN",            isTrip: false)
        _check("LOW_FEED",    feedwaterInv < 0.1,      2, "LOW FEEDWATER INVENTORY",          isTrip: false)
        // Xenon transient: Xe well above equilibrium AND power below 50%
        _xeMax = max(_xeMax, s.xenonInventory)
        let xeEq = PlantParams().gammaXe / (PlantParams().lambdaXe + PlantParams().xenonBurnCoeff)
        _check("XE_TRANSIENT", s.xenonInventory > xeEq * 1.4 && s.powerFraction < 0.5,
               3, "XENON TRANSIENT IN PROGRESS", isTrip: false)

        // Auto-trip on any unacknowledged trip
        if _alarmMap.values.contains(where: { $0.isTrip }) {
            scrammed = true
        }

        _rebuildAlarmLists()
    }

    private func _check(_ id: String, _ cond: Bool, _ priority: Int,
                        _ msg: String, isTrip: Bool) {
        if cond {
            if _alarmMap[id] == nil {
                let fo = !_firstOutSet
                if fo { _firstOutSet = true }
                _alarmMap[id] = ReactorAlarm(id: id, message: msg, priority: priority,
                                             state: "unack", isTrip: isTrip,
                                             isFirstOut: fo, simTime: snapshot.time)
            }
        } else {
            _alarmMap.removeValue(forKey: id)
        }
    }

    private func _rebuildAlarmLists() {
        alarms = _alarmMap.values.sorted { $0.priority < $1.priority }
        trips  = alarms.filter(\.isTrip).map(\.message)
        if _alarmMap.isEmpty { _firstOutSet = false }
    }

    // MARK: — Historian

    private func appendHistory() {
        let i = _histIdx % _histLen
        histPower[i]  = snapshot.powerFraction * 100
        histReact[i]  = snapshot.reactivity
        histFuelT[i]  = snapshot.fuelTempK
        histDecay[i]  = snapshot.decayHeatFraction * 100
        histPress[i]  = pressureMPa
        histElec[i]   = snapshot.electricPowerW / 1e6
        histSteamT[i] = snapshot.sgTempK
        histCoolT[i]  = snapshot.coolantTempK
        _histIdx += 1
    }

    // Returns the historian as a time-ordered array starting from oldest
    func orderedHistory<T>(_ arr: [T]) -> [T] {
        let n = arr.count
        let start = _histIdx % n
        return Array(arr[start...]) + Array(arr[..<start])
    }
}
