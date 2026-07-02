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

    // MARK: — Switchyard breakers (operator-toggleable on the mimic one-line).
    // 52G = generator breaker; line1/line2 = the two 400 kV outgoing circuits.
    var genBreakerOpen:   Bool  = false
    var line1BreakerOpen: Bool  = false
    var line2BreakerOpen: Bool  = false

    // Sim-clock pause (Esc). Lives here (a reference type) so the 60 Hz timer,
    // which captures the ContentView struct, reads it live rather than stale.
    var simPaused:        Bool  = false
    // Core-map popup (opened by tapping the core on the mimic; Esc closes).
    var coreMapOpen:      Bool  = false
    /// End-of-cycle core: axial-xenon feedback pushed toward divergence.
    var eolCore: Bool = false { didSet { plant.setEndOfCycle(eolCore) } }
    /// Control-room audio (procedural: turbine hum, annunciator chime, breaker clunk).
    let sound = SoundEngine()
    var soundEnabled: Bool = true { didSet { sound.enabled = soundEnabled } }
    private var _lastAlarmCount = 0

    // MARK: — Automation (all default OFF — manual operation is the trainer's
    // baseline; auto control is the realistic option, not the default)
    var rodAutoEnabled:  Bool   = false   // holds RCS T-avg at 550 K
    var pzrAutoEnabled:  Bool   = false   // heaters + spray hold 15.5 MPa
    var fwAutoEnabled:   Bool   = false   // holds FW inventory at 1.0
    var autoStartup:     Bool   = false   // sequencer: shutdown → full power
    private(set) var startupPhase: String = "STANDBY"
    private var _seqPhase: Int = 0        // 0 standby, 1 rod withdrawal, 2 ascension

    // MARK: — Readable plant state (UI reads)
    private(set) var snapshot:     PlantSnapshot
    private(set) var pressureMPa:  Double = 15.5
    private(set) var boronPPM:     Double = 800.0
    private(set) var feedwaterInv: Double = 1.0
    private(set) var steamInv:     Double = 1.0
    private(set) var omegaRCP:     Double = 1.0
    private(set) var porvOpen:     Bool   = false
    private(set) var eccsActuated: Bool   = false
    private(set) var steamDumpValve: Double = 0   // condenser dump (post-trip heat sink)
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
    private(set) var histTime:   [Double]   // SIM time [s] — normalises rate reads (SUR) at any speed
    private(set) var histAO:     [Double]   // axial offset ΔI [%] — flyspeck trail

    // Slow historian for post-run analysis / CSV export: one sample every 5 sim-s,
    // 4320 points ≈ 6 sim-hours. Columns are fixed (see exportCSV).
    private(set) var slowHist: [[Double]] = []
    private var _slowLastT: Double = -1e9

    // MARK: — Private
    private let plant: PWRPlant
    /// Turbine-generator electrical (AVR/excitation → MVAr) + mechanical
    /// supervision states — the mimic reads these instead of static gauges.
    let turbineGen = TurbineGenerator()
    private var _alarmMap:    [String: ReactorAlarm] = [:]
    private var _firstOutSet: Bool = false
    private var _histIdx:     Int  = 0
    private let _histLen:     Int  = 600
    // Xenon peak tracking for transient detection
    private var _xeMax:       Double = 0.0
    // Equilibrium Xe inventory at full power (computed once, not per frame)
    private let _xeEquilibrium: Double = {
        let p = PlantParams()
        return p.gammaXe / (p.lambdaXe + p.xenonBurnCoeff)
    }()

    /// The reactor kind in service. UI reads it to label kind-specific values.
    var reactorKind: ReactorKind { plant.params.kind }
    // Kind profile flags — the mimic gates pressurizer/SG/boron/pump drawing on
    // these instead of hard-coding the PWR plant picture.
    var hasPressurizer:  Bool { plant.params.hasPressurizer }
    var hasSteamGenerator: Bool { plant.params.hasSteamGenerator }
    var hasBoron:        Bool { plant.params.hasBoron }
    var isNaturalCirc:   Bool { plant.params.naturalCirculation }
    var nominalPressureMPa: Double { plant.params.nominalPressureMPa }
    var nominalMWe:      Double { plant.params.nominalPower * plant.params.turbineEfficiency / 1e6 }
    var nominalMWt:      Double { plant.params.nominalPower / 1e6 }

    init(kind: ReactorKind = .pwr) {
        let p: PlantParams
        switch kind {
        case .pwr: p = .pwr()
        case .bwr: p = .bwr()
        case .smr: p = .smr()
        }
        plant = ReactorPlant(params: p)
        let elecMW = p.nominalPower * p.turbineEfficiency / 1e6
        pressureMPa = p.nominalPressureMPa
        snapshot = PlantSnapshot(
            time: 0, powerFraction: 1, thermalPowerW: p.nominalPower,
            electricPowerW: p.nominalPower * p.turbineEfficiency,
            fuelTempK: 900, coolantTempK: 550, sgTempK: 553,
            reactivity: 0, xenonInventory: 0, iodineInventory: 0,
            rodPosition: 0, scrammed: false, decayHeatFraction: 0
        )
        let n = 600
        histPower  = Array(repeating: 100.0,                 count: n)
        histReact  = Array(repeating: 0.0,                   count: n)
        histFuelT  = Array(repeating: 900.0,                 count: n)
        histDecay  = Array(repeating: 0.0,                   count: n)
        histPress  = Array(repeating: p.nominalPressureMPa,  count: n)
        histElec   = Array(repeating: elecMW,                count: n)
        histSteamT = Array(repeating: 553.0,                 count: n)
        histCoolT  = Array(repeating: 550.0,                 count: n)
        histTime   = Array(repeating: 0.0,                   count: n)
        histAO     = Array(repeating: 0.0,                   count: n)
        sound.start()
    }

    // MARK: — Step

    func step(dt: Double) {
        updateAutomation(dt: dt)
        // Turbine trips automatically on reactor scram (standard interlock).
        if scrammed { turbineTrip = true }
        // With the turbine tripped, condenser STEAM DUMP becomes the SG heat
        // sink: a proportional valve holds no-load T_avg ≈ 550 K. Without it,
        // stored + decay heat bottles up the SG and the plant overheats with
        // no recovery path.
        let effectiveValve = turbineTrip
            ? max(0, min(0.2, 0.02 * (snapshot.coolantTempK - 550.0)))
            : turbineValve
        steamDumpValve = turbineTrip ? effectiveValve : 0
        // Natural circulation (SMR): coolant flow is buoyancy-driven and rises
        // with power — the operator can't pump harder. Otherwise flow is the
        // pump/recirc demand scaled by RCP speed.
        let effFlow = plant.params.naturalCirculation
            ? min(1.0, 0.25 + 0.75 * max(0, snapshot.powerFraction).squareRoot())
            : primaryFlow * omegaRCP
        let ctrl = ControlInputs(
            rodPosition:    rodPosition,
            primaryFlow:    effFlow,
            turbineValve:   effectiveValve,
            turbineTripped: turbineTrip,
            scram:          scrammed,
            primaryPressureMPa: pressureMPa   // live PZR pressure → DNBR anchor
        )
        snapshot = plant.step(dt: dt, ctrl: ctrl)
        scrammed = snapshot.scrammed

        turbineGen.step(dt: dt, grossMWe: snapshot.electricPowerW / 1e6,
                        tripped: turbineTrip || genBreakerOpen || (line1BreakerOpen && line2BreakerOpen),
                        ratedMWe: nominalMWe)
        updateBOP(dt: dt)
        updateBoron(dt: dt)
        updateAlarms()
        appendHistory()

        // Audio: hum tracks shaft speed; chime once per newly-raised alarm.
        sound.update(rpmFraction: turbineGen.rpm / 3000.0)
        if alarms.count > _lastAlarmCount { sound.chime() }
        _lastAlarmCount = alarms.count
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
            scramMessage = "RESET BLOCKED — active trips present (acknowledge first)"; return
        }
        scrammed = false
        plant.resetScram()
        // Demand follows actual position after a trip reset: rods stay inserted
        // until the operator deliberately withdraws (real CRDM behavior).
        rodPosition = snapshot.rodPosition
        // Re-close the switchyard so the unit can resync after a reset.
        genBreakerOpen = false; line1BreakerOpen = false; line2BreakerOpen = false
        scramMessage = "SCRAM RESET APPROVED"
    }

    // MARK: — Switchyard breaker operation
    /// Opening the generator breaker while synchronised is a load rejection —
    /// the turbine trips (no electrical load → overspeed protection).
    func toggleGenBreaker() {
        genBreakerOpen.toggle()
        sound.clunk()
        if genBreakerOpen && !turbineTrip { turbineTrip = true }
    }
    /// One outgoing circuit is redundant (the other carries full load). Opening
    /// BOTH with the generator breaker closed is a full load rejection → trip.
    func toggleLineBreaker(_ i: Int) {
        if i == 0 { line1BreakerOpen.toggle() } else { line2BreakerOpen.toggle() }
        sound.clunk()
        if line1BreakerOpen && line2BreakerOpen && !genBreakerOpen && !turbineTrip {
            turbineTrip = true
        }
    }

    func acknowledgeAllAlarms() {
        for key in _alarmMap.keys { _alarmMap[key]?.state = "ack" }
        _rebuildAlarmLists()
    }

    // MARK: — Automation

    private func updateAutomation(dt: Double) {
        // ── Startup sequencer (semi-automatic power ascension, APR1400-style) ──
        if autoStartup {
            if scrammed {
                startupPhase = "HOLD — SCRAM (reset first)"
                _seqPhase = 0
            } else if !startupPermit {
                startupPhase = "HOLD — NO STARTUP PERMIT"
            } else {
                let pf = snapshot.powerFraction
                if _seqPhase == 0 { _seqPhase = 1 }
                if _seqPhase == 1 {
                    // Turbine stays TRIPPED here: criticality is approached
                    // against the steam dump holding T-avg ≈ 550 K. Re-latching
                    // a wide-open turbine at 3% power would be a cold-water
                    // reactivity excursion.
                    startupPhase = "ROD WITHDRAWAL — DUMP HOLDS T-AVG"
                    turbineTrip = true
                    rodPosition = max(0, rodPosition - 0.0053 * dt)   // CRDM max rate
                    feedwaterValve = max(feedwaterValve, 0.30)
                    fwAutoEnabled = true
                    if pf >= 0.15 {
                        _seqPhase = 2
                        turbineTrip  = false
                        turbineValve = 0.05           // roll turbine from hot standby
                        rodAutoEnabled = true         // T-avg program takes the rods
                    }
                }
                if _seqPhase == 2 {
                    startupPhase = "POWER ASCENSION — RAMPING TURBINE"
                    turbineValve = min(1.0, turbineValve + 0.005 * dt)
                    if pf >= 0.98 && turbineValve >= 0.999 {
                        startupPhase = "AT POWER — SEQ COMPLETE"
                        autoStartup = false           // hands over to rod auto
                        _seqPhase = 0
                    }
                }
            }
        } else {
            if !scrammed { startupPhase = "STANDBY" }
            _seqPhase = 0
        }

        // ── Rod auto-control: hold RCS T-avg at 550 K (0.5 K deadband).
        // Drives the DEMAND at CRDM speed; actual rods are rate-limited anyway.
        if rodAutoEnabled && !scrammed && snapshot.powerFraction > 0.02 {
            let err = snapshot.coolantTempK - 550.0
            if abs(err) > 0.5 {
                // hot → insert (demand up), cold → withdraw
                let rate = max(-1.0, min(1.0, err * 0.2)) * 0.0053
                rodPosition = max(0, min(1, rodPosition + rate * dt))
            }
        }

        // ── Feedwater auto: power feedforward + inventory correction, slewed.
        if fwAutoEnabled && !feedwaterFault {
            let pf = max(0, snapshot.powerFraction)
            let target = 0.7 * pf + 0.8 * (1.0 - feedwaterInv)
            let dv = max(-0.05 * dt, min(0.05 * dt, target - feedwaterValve))
            feedwaterValve = max(0, min(1, feedwaterValve + dv))
        }
    }

    // MARK: — BOP

    private func updateBOP(dt: Double) {
        // RCP coast-down on pump fault. Natural-circulation plants (SMR) have
        // no reactor coolant pumps — flow is set in step(), so pin omega to 1.
        if plant.params.naturalCirculation {
            omegaRCP = 1.0
        } else if pumpDegraded {
            omegaRCP = max(0.04, omegaRCP - dt / 12.0)
        } else {
            omegaRCP = min(1.0, omegaRCP + dt / 30.0)
        }

        // Steam inventory: SG production follows core heat; the turbine removes it.
        // Matched coefficients (0.08/0.08) so it balances at steady state, and
        // consumption scales with available inventory so the balance is
        // self-regulating (settles at inv ≈ power/valve) instead of drifting
        // to a clamp on any tiny power/valve mismatch.
        let steamProduction  = max(0, snapshot.powerFraction) * 0.08 * dt
        let steamConsumption = turbineTrip
            ? 0.0
            : turbineValve * 0.08 * min(1.0, max(0.0, steamInv)) * dt
        steamInv = max(0, min(2.0, steamInv + steamProduction - steamConsumption))

        // Condenser temp: first-order ODE solved EXACTLY (stable at 600× speed;
        // explicit Euler diverges for dt > 2.5 s with this 0.8/s rate constant).
        let steamHeat  = turbineTrip ? 0.0 : turbineValve * snapshot.powerFraction * 5.0
        let condTarget = 305.0 + steamHeat / 0.8
        condTempK = condTarget + (condTempK - condTarget) * exp(-0.8 * dt)
        condTempK = max(300.0, min(370.0, condTempK))

        // Feedwater balance (fault or turbine trip stops feedwater).
        // In = valve·0.10, out = power·0.07 → balances at valve 0.70 for 100%
        // power, matching the default valve position (no spurious LOW_FEED).
        let fwIn  = (feedwaterFault || turbineTrip) ? 0.0 : feedwaterValve * 0.10 * dt
        let fwOut = max(0, snapshot.powerFraction) * 0.07 * dt
        feedwaterInv = max(0, min(1.2, feedwaterInv + fwIn - fwOut))

        // Primary/dome pressure. PWR & SMR hold pressure with a pressurizer
        // (heaters + spray); the target tracks the nominal set pressure with an
        // asymmetric thermal coupling so an overheat can still reach the HIGH
        // trip while a cooldown holds above the ECCS setpoint. A BWR has no
        // pressurizer: the steam dome swells with steam inventory and is relieved
        // by the turbine. All thresholds are referenced to the nominal pressure.
        let nomP = plant.params.nominalPressureMPa
        if plant.params.hasPressurizer {
            let dT = snapshot.coolantTempK - 550.0
            let pTarget = pzrAutoEnabled
                ? nomP + (dT > 0 ? dT * 0.03 : dT * 0.004)
                : nomP + (dT > 0 ? dT * 0.05 : dT * 0.01)
            let pRate = pzrAutoEnabled ? 0.05 : 0.02
            pressureMPa = pTarget + (pressureMPa - pTarget) * exp(-pRate * dt)
            pressureMPa = max(10.0, min(17.5, pressureMPa))
        } else {
            // BWR steam dome: pressure rises with stored steam, turbine relieves it.
            let pTarget = nomP + (steamInv - 1.0) * 1.5
            pressureMPa = pTarget + (pressureMPa - pTarget) * exp(-0.1 * dt)
            pressureMPa = max(2.0, min(nomP + 2.0, pressureMPa))
        }

        // Relief valve (PORV / BWR SRV) lifts above the set pressure.
        porvOpen = pressureMPa > nomP + 1.05

        // ECCS arms on low pressure (loss of inventory), resets when restored.
        if pressureMPa < nomP - 4.0 { eccsActuated = true }
        if pressureMPa > nomP - 2.5 && feedwaterInv > 0.3 { eccsActuated = false }
    }

    private func updateBoron(dt: Double) {
        // BWRs have no soluble-boron chemical shim in normal operation.
        guard plant.params.hasBoron else { return }
        let vPrimary = 300.0   // m³ approximate primary coolant volume
        let dB = (borationRate * 150.0 - dilutionRate * boronPPM) / vPrimary * dt
        boronPPM = max(0, boronPPM + dB)
        // Reactivity: −8 pcm/ppm, nominal 800 ppm → 0 reactivity
        plant.params.externalReactivity = -8e-5 * (boronPPM - 800.0)
    }

    // MARK: — Alarms

    private func updateAlarms() {
        let s = snapshot
        let nomP = plant.params.nominalPressureMPa

        let pTrip = nomP + 1.5   // HIHI trip 1.5 MPa over set pressure (17.0 PWR / 8.5 BWR)
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
        _check("XE_TRANSIENT", s.xenonInventory > _xeEquilibrium * 1.4 && s.powerFraction < 0.5,
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
            // Annunciator latching: an alarm whose condition has cleared stays
            // on the board until ACKNOWLEDGED (so first-out evidence survives,
            // and a trip can't silently reset itself).
            if let a = _alarmMap[id], a.state == "ack" {
                _alarmMap.removeValue(forKey: id)
            }
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
        histTime[i]   = snapshot.time
        histAO[i]     = snapshot.axialOffsetPct
        _histIdx += 1

        // Slow historian: one row per 5 sim-seconds, ~6 sim-hours retained.
        if snapshot.time - _slowLastT >= 5.0 {
            _slowLastT = snapshot.time
            slowHist.append([snapshot.time,
                             snapshot.powerFraction * 100,
                             snapshot.coolantTempK,
                             pressureMPa,
                             snapshot.steamPressureMPa,
                             snapshot.electricPowerW / 1e6,
                             snapshot.axialOffsetPct,
                             snapshot.minDNBR,
                             -1.6 * snapshot.xenonInventory,
                             boronPPM])
            if slowHist.count > 4320 { slowHist.removeFirst(slowHist.count - 4320) }
        }
    }

    /// Write the slow historian to ~/Downloads as CSV. Returns the file path
    /// (surfaced via scramMessage-style toast by the caller) or nil on failure.
    func exportCSV() -> String? {
        let header = "sim_time_s,power_pct,t_avg_K,rcs_p_MPa,steam_p_MPa,gross_MWe,axial_offset_pct,min_dnbr,xenon_pcm,boron_ppm"
        var rows = [header]
        rows.reserveCapacity(slowHist.count + 1)
        for r in slowHist {
            rows.append(r.map { String(format: "%.4g", $0) }.joined(separator: ","))
        }
        let dir = FileManager.default.urls(for: .downloadsDirectory, in: .userDomainMask).first
            ?? URL(fileURLWithPath: NSHomeDirectory())
        let stamp = Int(snapshot.time)
        let url = dir.appendingPathComponent("ReactorSim-\(reactorKind.rawValue)-t\(stamp)s.csv")
        do {
            try rows.joined(separator: "\n").write(to: url, atomically: true, encoding: .utf8)
            return url.path
        } catch { return nil }
    }

    // Returns the historian as a time-ordered array starting from oldest
    func orderedHistory<T>(_ arr: [T]) -> [T] {
        let n = arr.count
        let start = _histIdx % n
        return Array(arr[start...]) + Array(arr[..<start])
    }
}
