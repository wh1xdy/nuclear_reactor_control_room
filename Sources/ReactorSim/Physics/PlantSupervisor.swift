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
    /// Generator paralleled to the grid. True at power; false during a cold
    /// startup until the operator syncs (rolls the turbine to speed, catches the
    /// synchroscope at 12 o'clock, closes 52G). Unsynced ⇒ 0 grid MWe.
    var genSynced: Bool = true

    // MARK: — Startup state
    enum StartMode { case atPower, coldStartup }
    /// Set while a cold startup is in progress (subcritical → on the grid); the
    /// hot-standby temperature hold and sync gating key off it. Cleared once the
    /// plant is generating on the grid above ~10 %.
    private(set) var coldStartActive = false

    // MARK: — Auxiliary (station) power — the LOOP/SBO state machine.
    // GRID: offsite or the unit's own generator feeds the aux buses.
    // HOUSE LOAD: grid lost but the generator islands onto its own auxiliaries
    //   (turbine runback to ~8% — Swedish plants drill exactly this).
    // DIESELS: all AC lost → EDGs crank ~10 s, then carry the ESF buses only
    //   (AFW motor pump, chargers — NOT the RCPs or circ-water pumps).
    // BLACKOUT: diesels failed/cranking → batteries only, clock ticking.
    enum AuxPowerState: String { case grid = "GRID", houseLoad = "HOUSE LOAD",
                                 diesel = "DIESELS", sbo = "BLACKOUT" }
    private(set) var auxPower: AuxPowerState = .grid
    private(set) var dieselsRunning = false
    private(set) var dieselTimer: Double = 0          // crank progress [s]
    private(set) var batteryMin: Double = 120         // SBO battery clock [min]
    private(set) var afwRunning = false               // turbine-driven aux feed
    private(set) var srvOpen = false                  // SG safety valves venting
    var dieselFault = false                           // malfunction hook (menu)
    private var _prevOffsite = true
    private var _prevDiesels = false
    private var _prevSBO = false
    /// Station AC available for the big non-ESF loads (RCPs, circ water, main FW).
    var stationACPowered: Bool { (!turbineTrip && !genBreakerOpen) || !(line1BreakerOpen && line2BreakerOpen) }

    // Sim-clock pause (Esc). Lives here (a reference type) so the 60 Hz timer,
    // which captures the ContentView struct, reads it live rather than stale.
    // Pausing also gates the audio (hum/horn/speech go quiet with the physics).
    var simPaused:        Bool  = false { didSet { sound.setPaused(simPaused) } }
    // Core-map popup (opened by tapping the core on the mimic; Esc closes).
    var coreMapOpen:      Bool  = false
    // Startup-physics popup (tap the NEUTRONICS dock): 1/M, ECP, rod worth.
    var startupPanelOpen: Bool  = false
    /// 1/M approach-to-critical points: (withdrawn fraction, CR₀/CR).
    private(set) var invMPoints: [(x: Double, invM: Double)] = []
    private var _invMBaseN: Double = 0
    private var _invMLastRod: Double = 1
    // Settings sheet visibility — the mimic pauses its 120 Hz canvas behind
    // modals so sheet scrolling stays smooth (Canvas redraw fights the sheet).
    var settingsOpen:     Bool  = false
    /// End-of-cycle core: axial-xenon feedback pushed toward divergence.
    var eolCore: Bool = false { didSet { plant.setEndOfCycle(eolCore) } }
    /// Control-room audio (procedural: turbine hum, annunciator chime, breaker clunk).
    let sound = SoundEngine()
    var soundEnabled: Bool = true { didSet { sound.enabled = soundEnabled } }
    /// Spoken annunciator: major events (trips, safety injection) always speak;
    /// this toggle additionally speaks EVERY newly-raised alarm window.
    var voiceAllAlarms: Bool = false
    private var _lastAlarmCount = 0
    private var _newAlarmMsgs: [(id: String, msg: String)] = []   // raised this step (voice queue)
    private var _prevScrammed = false
    private var _prevTbnTrip  = false
    private var _prevECCS     = false
    // Rate-of-change tracking (sampled every ~2 sim-s, exponentially smoothed).
    private var _rateT: Double = -1
    private var _ratePrev: (tAvg: Double, p: Double, pw: Double) = (550, 15.5, 100)
    private var _rate: (dT: Double, dP: Double, dPw: Double) = (0, 0, 0)

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
    // Per-pump states (4 RCPs on a PWR, 2 recirc pumps on a BWR, 0 on the
    // nat-circ SMR). The thermal model consumes the AGGREGATE (omegaRCP =
    // mean over installed pumps), so the validated physics is untouched;
    // individual pumps give the trainer single-pump-loss scenarios.
    var rcpRunning: [Bool] = [true, true, true, true]
    private(set) var rcpOmega: [Double] = [1, 1, 1, 1]
    /// Installed reactor-coolant/recirc pump count for the active kind.
    var rcpCount: Int { plant.params.naturalCirculation ? 0 : (plant.params.kind == .bwr ? 2 : 4) }
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
    /// Containment atmosphere (pressure/temp/sump + spray). Sources: PORV
    /// relief-tank venting and the malfunction menu's primary leak.
    let containment = Containment()
    /// Primary-side break flow into containment [kg/s] — set by the LOCA
    /// malfunction; 0 in normal operation.
    var primaryLeakKgs: Double = 0

    // MARK: — Malfunctions (instructor menu). Each is an independent latched
    // fault with real physics; CLEAR ALL restores the plant systems (the
    // transient they caused still has to be recovered by the operator).
    var malfMenuOpen = false                 // instructor panel visibility
    var sgtrLeakKgs: Double = 0              // SG tube rupture: primary → secondary
    var stuckRod = false                     // one RCCA frozen: tilt + degraded scram worth
    var droppedRod = false                   // one RCCA on the bottom: −350 pcm + tilt
    var atwsFault = false                    // reactor protection fails to scram
    var msivClosed = false                   // spurious main-steam isolation
    private var _rodMalfApplied = false
    /// Any instructor fault engaged (tints the MALFUNCTIONS chip).
    var anyMalfunctionActive: Bool {
        primaryLeakKgs > 0 || sgtrLeakKgs > 0 || stuckRod || droppedRod
            || atwsFault || msivClosed || dieselFault || pumpDegraded || feedwaterFault
            || rcpRunning.prefix(rcpCount).contains(false)
    }
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

    init(kind: ReactorKind = .pwr, mode: StartMode = .atPower) {
        var p: PlantParams
        switch kind {
        case .pwr: p = .pwr()
        case .bwr: p = .bwr()
        case .smr: p = .smr()
        }
        // OpenMC calibration applies to the PWR-type lattices; BWR physics is
        // void-dominated and keeps its own profile.
        if kind != .bwr { p.applyCalibration() }
        plant = ReactorPlant(params: p)
        let elecMW = p.nominalPower * p.turbineEfficiency / 1e6
        pressureMPa = p.nominalPressureMPa
        let cold = (mode == .coldStartup)
        snapshot = PlantSnapshot(
            time: 0, powerFraction: cold ? 0 : 1,
            thermalPowerW: cold ? 0 : p.nominalPower,
            electricPowerW: cold ? 0 : p.nominalPower * p.turbineEfficiency,
            fuelTempK: cold ? 550 : 900, coolantTempK: 550, sgTempK: cold ? 550 : 553,
            reactivity: 0, xenonInventory: 0, iodineInventory: 0,
            rodPosition: cold ? 1 : 0, scrammed: false, decayHeatFraction: 0
        )
        let n = 600
        histPower  = Array(repeating: cold ? 0 : 100.0,      count: n)
        histReact  = Array(repeating: 0.0,                   count: n)
        histFuelT  = Array(repeating: cold ? 550 : 900.0,    count: n)
        histDecay  = Array(repeating: 0.0,                   count: n)
        histPress  = Array(repeating: p.nominalPressureMPa,  count: n)
        histElec   = Array(repeating: cold ? 0 : elecMW,     count: n)
        histSteamT = Array(repeating: cold ? 550 : 553.0,    count: n)
        histCoolT  = Array(repeating: 550.0,                 count: n)
        histTime   = Array(repeating: 0.0,                   count: n)
        histAO     = Array(repeating: 0.0,                   count: n)

        if cold {
            // Hot standby: rods in, source-range flux, turbine offline, generator
            // off the bus. The operator drives the approach to criticality.
            plant.coldStandby()
            coldStartActive = true
            rodPosition   = 1.0
            turbineValve  = 0.0
            turbineTrip   = true
            feedwaterValve = 0.0
            genBreakerOpen = true
            genSynced      = false
            turbineGen.park()
        }
        sound.start()
    }

    // MARK: — Step

    func step(dt: Double) {
        updateAutomation(dt: dt)
        updateAuxPower(dt: dt)
        // Turbine trips automatically on reactor scram (standard interlock).
        if scrammed { turbineTrip = true }
        // House-load operation: islanded on our own auxiliaries — the governor
        // runs the turbine back to just carry the station load (~8%).
        if auxPower == .houseLoad { turbineValve = min(turbineValve, 0.08) }
        // With the turbine tripped, condenser STEAM DUMP becomes the SG heat
        // sink IF the condenser still has vacuum (circ-water pumps are big
        // non-ESF loads — a LOOP takes them, the vacuum collapses, and the SG
        // safety valves become the only heat path, venting to atmosphere).
        let vacuumOK = stationACPowered && condTempK < 340
        let dumpValve = vacuumOK ? max(0, min(0.35, 0.02 * (snapshot.coolantTempK - 550.0))) : 0
        // SRVs: lift/reseat with hysteresis on SG (steam) pressure. Physical
        // valves — they backstop EVERY configuration (trip, house load, SBO).
        if plant.params.hasSteamGenerator {
            if snapshot.steamPressureMPa > 7.9 { srvOpen = true }
            if snapshot.steamPressureMPa < 7.4 { srvOpen = false }
        } else { srvOpen = false }
        let srvValve = srvOpen ? min(0.15, 0.08 * (snapshot.steamPressureMPa - 7.4)) : 0
        // SG heat removal: turbine when latched; on a house-load runback the
        // DUMP carries the shed load (that's how a real runback survives —
        // turbine to ~8 %, dump takes the rest, rods walk power down); after a
        // trip the dump holds no-load T-avg while vacuum lasts.
        // Spurious MSIV closure: turbine AND dump are downstream of the MSIVs —
        // only the SG safety valves can still relieve. The turbine trips on
        // loss of steam.
        if msivClosed && !turbineTrip { turbineTrip = true }
        let effectiveValve: Double
        if msivClosed {
            effectiveValve = srvValve
        } else if turbineTrip {
            effectiveValve = dumpValve + srvValve
        } else if auxPower == .houseLoad {
            effectiveValve = turbineValve + dumpValve + srvValve
        } else {
            effectiveValve = turbineValve + srvValve
        }
        steamDumpValve = dumpValve > 0 && !msivClosed && (turbineTrip || auxPower == .houseLoad) ? dumpValve : 0

        // Rod malfunctions: apply the flux tilt (and degraded scram worth for a
        // stuck RCCA) once on engagement; restore the worth on clear.
        if (stuckRod || droppedRod) && !_rodMalfApplied {
            _rodMalfApplied = true
            plant.kickTilt(x: droppedRod ? -0.05 : 0.04, y: 0.025)
            if stuckRod { plant.params.scramExtraWorth = -0.102 }   // one RCCA short
        }
        if !(stuckRod || droppedRod) && _rodMalfApplied {
            _rodMalfApplied = false
            plant.params.scramExtraWorth = -0.12
        }
        // Natural circulation (SMR): coolant flow is buoyancy-driven and rises
        // with power — the operator can't pump harder. Otherwise flow is the
        // pump/recirc demand scaled by RCP speed, with a natural-circulation
        // FLOOR once the pumps coast down (a dead-pump PWR still convects —
        // this is what carries decay heat through a station blackout).
        var effFlow = plant.params.naturalCirculation
            ? min(1.0, 0.25 + 0.75 * max(0, snapshot.powerFraction).squareRoot())
            : primaryFlow * omegaRCP
        if !plant.params.naturalCirculation {
            let natFloor = 0.06 + 0.10 * max(0, snapshot.powerFraction).squareRoot()
            effFlow = max(effFlow, natFloor)
        }
        let ctrl = ControlInputs(
            rodPosition:    rodPosition,
            primaryFlow:    effFlow,
            turbineValve:   effectiveValve,
            turbineTripped: turbineTrip,
            gridConnected:  genSynced,
            scram:          scrammed,
            primaryPressureMPa: pressureMPa   // live PZR pressure → DNBR anchor
        )
        snapshot = plant.step(dt: dt, ctrl: ctrl)
        scrammed = snapshot.scrammed

        // Turbine-generator. Off the bus (pre-sync) the machine free-runs on the
        // throttle; the synchroscope shows the slip.
        let steaming = !turbineTrip && effectiveValve > 0.03 && snapshot.powerFraction > 0.02
        turbineGen.step(dt: dt, grossMWe: snapshot.electricPowerW / 1e6,
                        tripped: turbineTrip || (genSynced && (genBreakerOpen || (line1BreakerOpen && line2BreakerOpen))),
                        synced: genSynced, steaming: steaming,
                        ratedMWe: nominalMWe)
        // Cold startup completes once the unit is paralleled and carrying load.
        if coldStartActive && genSynced && snapshot.powerFraction > 0.10 { coldStartActive = false }
        // Containment atmosphere: fed by any primary break (LOCA malfunction)
        // plus PORV discharge once the relief tank vents; coolers need AC,
        // spray needs the ESF buses (diesels suffice).
        containment.step(dt: dt,
                         releaseKgs: primaryLeakKgs + (porvOpen ? 6 : 0),
                         acPowered: stationACPowered || dieselsRunning,
                         esfAvailable: stationACPowered || dieselsRunning)
        updateBOP(dt: dt)
        updateBoron(dt: dt)
        updateAlarms()
        appendHistory()
        recordInvM()

        // ── Audio + spoken annunciator ────────────────────────────────────
        sound.update(rpmFraction: turbineGen.rpm / 3000.0)
        // Major events: urgent horn on a reactor trip; calm voice callouts.
        // ── 52G follows the unit state (fixes "clunk but nothing happens"):
        // the reverse-power relay OPENS the gen breaker on a turbine trip, and
        // auto-sync recloses it when the trip clears — so the drawn breaker,
        // the stored state, and the sound are always the same event.
        // Only while already paralleled — during a cold startup the operator
        // syncs manually, so the reverse-power/auto-sync relay stays out of it.
        if genSynced {
            if turbineTrip && !_prevTbnTrip && !genBreakerOpen {
                genBreakerOpen = true
                sound.clunk(closing: false)
            }
            if !turbineTrip && _prevTbnTrip && genBreakerOpen {
                genBreakerOpen = false
                sound.clunk(closing: true)
            }
        }

        let majorEvent = (scrammed && !_prevScrammed) || (eccsActuated && !_prevECCS)
        if scrammed && !_prevScrammed {
            sound.horn()
            sound.announce(key: "reactor-trip", text: "Reactor trip. Reactor trip.", priority: true)
        }
        if turbineTrip && !_prevTbnTrip && !scrammed {
            sound.announce(key: "turbine-trip", text: "Turbine trip.", priority: true)
        }
        if eccsActuated && !_prevECCS {
            sound.announce(key: "safety-injection", text: "Safety injection initiated.", priority: true)
        }
        _prevScrammed = scrammed; _prevTbnTrip = turbineTrip; _prevECCS = eccsActuated
        // LOOP-ladder callouts (grid gone / diesels up / blackout) + diesel sound.
        let offsiteNow = !(line1BreakerOpen && line2BreakerOpen)
        if !offsiteNow && _prevOffsite {
            sound.announce(key: "loop", text: "Loss of offsite power.", priority: true)
        }
        _prevOffsite = offsiteNow
        if dieselsRunning && !_prevDiesels {
            sound.announce(key: "dg_run", text: "Diesel generators running.")
        }
        _prevDiesels = dieselsRunning
        let sboNow = auxPower == .sbo && dieselTimer > 12
        if sboNow && !_prevSBO {
            sound.announce(key: "sbo", text: "Station blackout. Station blackout.", priority: true)
        }
        _prevSBO = sboNow
        sound.setDiesel(dieselsRunning ? 2 : (auxPower == .sbo && dieselTimer > 0 && dieselTimer < 10 && !dieselFault ? 1 : 0))
        // Ordinary new alarms: chime, plus per-alarm voice if enabled — capped
        // per step, and window messages that duplicate a just-spoken major
        // callout (the trip windows) are skipped rather than double-announced.
        // The alarm ID doubles as the pre-rendered callout key.
        if alarms.count > _lastAlarmCount { sound.chime() }
        _lastAlarmCount = alarms.count
        if voiceAllAlarms {
            for item in _newAlarmMsgs.prefix(3) where !(majorEvent && item.msg.contains("TRIP")) {
                sound.announce(key: item.id.lowercased(), text: item.msg)
            }
        }
        _newAlarmMsgs.removeAll()
    }

    // MARK: — Operator actions

    func triggerScram() {
        // ATWS malfunction: the reactor protection system fails to act — the
        // drill is emergency boration, not the scram button.
        guard !atwsFault else { scramMessage = "RPS FAILURE — RODS DID NOT INSERT"; return }
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
        // Re-close the transmission circuits; 52G stays with the turbine —
        // auto-sync recloses it only when the turbine trip actually clears.
        line1BreakerOpen = false; line2BreakerOpen = false
        scramMessage = "SCRAM RESET APPROVED"
    }

    // MARK: — Switchyard breaker operation
    /// Opening the generator breaker while synchronised is a load rejection —
    /// the turbine trips (no electrical load → overspeed protection).
    func toggleGenBreaker() {
        // Pre-sync (cold startup): closing 52G PARALLELS the generator — allowed
        // only at speed and near phase (catch the synchroscope at 12 o'clock).
        if !genSynced {
            guard genBreakerOpen else { return }        // already off the bus
            guard turbineGen.readyToSync else {
                scramMessage = "52G CLOSE BLOCKED — turbine not at speed (roll to 3000 rpm)"; return
            }
            let a = abs(turbineGen.syncAngleDeg)
            guard a < 20 || a > 340 else {
                scramMessage = "52G — OUT OF PHASE, close as the synchroscope passes 12"; return
            }
            genSynced = true
            genBreakerOpen = false
            sound.clunk(closing: true)
            scramMessage = "GENERATOR SYNCHRONISED — on the grid"
            return
        }
        // Sync check: closing 52G onto a tripped/dead turbine is REJECTED (no
        // state change, no sound) — you can't parallel a machine that isn't
        // at speed. Clear the turbine trip first; auto-sync then recloses.
        if genBreakerOpen && turbineTrip {
            scramMessage = "52G CLOSE BLOCKED — SYNC CHECK (turbine tripped)"
            return
        }
        genBreakerOpen.toggle()
        sound.clunk(closing: !genBreakerOpen)
        if genBreakerOpen && !turbineTrip { turbineTrip = true }
    }
    /// One outgoing circuit is redundant (the other carries full load).
    /// Opening BOTH while generating no longer force-trips the unit — the
    /// generator ISLANDS onto house load (governor runback); the aux-power
    /// state machine takes it from there. Losing the gen too → diesels/SBO.
    func toggleLineBreaker(_ i: Int) {
        if i == 0 { line1BreakerOpen.toggle() } else { line2BreakerOpen.toggle() }
        sound.clunk(closing: !(i == 0 ? line1BreakerOpen : line2BreakerOpen))
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
                    // Startup-rate hold: withdraw only while the power rise is
                    // controlled (a real operator paces the approach to
                    // criticality on the startup-rate meter — continuous
                    // max-rate withdrawal rides the period straight into the
                    // HIGH FLUX trip after a scram).
                    if _rate.dPw < 0.25 {
                        rodPosition = max(0, rodPosition - 0.0053 * dt)
                    } else {
                        startupPhase = "ROD HOLD — STARTUP RATE"
                    }
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
                    // Auto-sync the generator once it rolls up to speed (a cold
                    // startup boots off the bus; the sequencer parallels it).
                    if !genSynced {
                        startupPhase = "ROLLING TURBINE — AUTO SYNC"
                        turbineValve = max(turbineValve, 0.08)   // roll on modest steam
                        if turbineGen.readyToSync {
                            genSynced = true; genBreakerOpen = false
                            sound.clunk(closing: true)
                        }
                    } else {
                        // Power-limited ascension: open the throttle only while
                        // below 100 % and while power isn't already climbing hard,
                        // so the reactor follows without riding the period into a
                        // high-flux trip (a wide-open ramp from a low-power sync
                        // overshoots to ~120 %).
                        startupPhase = "POWER ASCENSION — RAMPING TURBINE"
                        // Ramp the throttle fully open, pausing whenever power is
                        // climbing fast (rate brake) so the reactor tracks the
                        // steam demand without riding the period into a trip.
                        if turbineValve < 0.999 && _rate.dPw < 0.20 {
                            turbineValve = min(1.0, turbineValve + 0.0025 * dt)
                        }
                        if turbineValve >= 0.999 && pf >= 0.95 {
                            startupPhase = "AT POWER — SEQ COMPLETE"
                            autoStartup = false           // hands over to rod auto
                            _seqPhase = 0
                        }
                    }
                }
            }
        } else {
            if !scrammed { startupPhase = "STANDBY" }
            _seqPhase = 0
        }

        // ── Rod auto-control: hold RCS T-avg at 550 K (0.5 K deadband), with
        // a FLUX-LIMITER channel on top — above 103% RTP the controller drives
        // rods IN regardless of temperature (every real rod controller has a
        // power-limiting channel; T-avg alone reacts too late on an overshoot
        // because the dump/turbine hold the temperature while flux runs).
        if rodAutoEnabled && !scrammed && snapshot.powerFraction > 0.02 {
            let pf = snapshot.powerFraction
            if pf > 1.03 {
                // Flux limiter: drive rods IN regardless of temperature.
                let drive = min(1.0, (pf - 1.03) * 25)          // full rate by ~107%
                rodPosition = max(0, min(1, rodPosition + 0.0053 * drive * dt))
            } else {
                // Rod speed proportional to T-avg error (0.5 K deadband). The rod
                // is an integrator, so proportional velocity control already holds
                // T-avg offset-free — no explicit integral needed. The plant
                // self-regulates power to the turbine demand via moderator
                // feedback (open-loop stable); the rods only TRIM T-avg back to
                // 550, so the gain is deliberately LOW — an aggressive controller
                // fights an already-stable plant and limit-cycles. The withdrawal
                // block keeps the T-avg channel from pushing power into the flux
                // limiter (the two channels bang-banged at high gain).
                let tErr = snapshot.coolantTempK - 550.0
                if abs(tErr) > 0.5 {
                    var drive = max(-1.0, min(1.0, tErr * 0.06))
                    if drive < 0 && pf > 1.005 { drive = 0 }
                    rodPosition = max(0, min(1, rodPosition + drive * 0.0053 * dt))
                }
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

    /// The LOOP/SBO ladder. Evaluated every step; transitions announce below.
    private func updateAuxPower(dt: Double) {
        let genOnline = !turbineTrip && !genBreakerOpen
        let offsite   = !(line1BreakerOpen && line2BreakerOpen)
        if offsite || genOnline {
            auxPower = offsite && !genOnline ? .grid : (offsite ? .grid : .houseLoad)
            dieselsRunning = false
            dieselTimer = 0
            batteryMin = min(120, batteryMin + dt / 120)   // chargers restore slowly
        } else if dieselsRunning {
            auxPower = .diesel
        } else if dieselFault {
            auxPower = .sbo
            batteryMin = max(0, batteryMin - dt / 60)
        } else {
            // EDGs cranking (~10 s) — on batteries until they pick up.
            auxPower = .sbo
            dieselTimer += dt
            batteryMin = max(0, batteryMin - dt / 60)
            if dieselTimer >= 10 { dieselsRunning = true }
        }
    }

    private func updateBOP(dt: Double) {
        // Per-pump coast/run: a pump spins only if IT is selected running, the
        // station has AC (RCPs are huge non-ESF loads — diesels don't carry
        // them), and the legacy degradation fault is clear. The thermal model
        // sees the aggregate. Natural-circulation plants have no pumps.
        if plant.params.naturalCirculation {
            omegaRCP = 1.0
        } else {
            let n = rcpCount
            for i in 0..<n {
                if rcpRunning[i] && stationACPowered && !pumpDegraded {
                    rcpOmega[i] = min(1.0, rcpOmega[i] + dt / 30.0)
                } else {
                    rcpOmega[i] = max(0.02, rcpOmega[i] - dt / 12.0)
                }
            }
            omegaRCP = rcpOmega.prefix(n).reduce(0, +) / Double(n)
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
        // Without station AC the circ-water pumps stop and the vacuum collapses.
        let steamHeat  = turbineTrip ? 0.0 : turbineValve * snapshot.powerFraction * 5.0
        let condTarget = stationACPowered ? 305.0 + steamHeat / 0.8 : 368.0
        let condRate   = stationACPowered ? 0.8 : 0.05
        condTempK = condTarget + (condTempK - condTarget) * exp(-condRate * dt)
        condTempK = max(300.0, min(370.0, condTempK))

        // Feedwater balance. Main FW pumps are station-AC loads; when they die
        // the TURBINE-DRIVEN aux-feed pump carries decay-heat feed as long as
        // the SG can supply steam (works even in a blackout — that's its job),
        // helped by the motor AFW pump once the diesels are up.
        let fwPowered = stationACPowered && !feedwaterFault && !turbineTrip
        var fwIn = fwPowered ? feedwaterValve * 0.10 * dt : 0.0
        let sgSteaming = snapshot.steamPressureMPa > 0.8
        afwRunning = !fwPowered && !feedwaterFault && sgSteaming && feedwaterInv < 0.6
        if afwRunning {
            fwIn += 0.030 * dt                                   // TD-AFW
            if dieselsRunning || stationACPowered { fwIn += 0.015 * dt }   // motor AFW
        }
        // SGTR: ruptured tubes feed PRIMARY water into the SG — the level
        // rises without feedwater (the classic identification cue).
        fwIn += sgtrLeakKgs * 4e-4 * dt
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
            // Inventory-loss malfunctions bleed pressure faster than the PZR
            // heaters can hold; safety injection fights back once actuated.
            pressureMPa -= (primaryLeakKgs * 4e-4 + sgtrLeakKgs * 3e-4) * dt
            if eccsActuated { pressureMPa += 0.004 * dt }
            pressureMPa = max(6.0, min(17.5, pressureMPa))
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
        // Rod-malfunction reactivity rides on the external term with boron:
        // a dropped RCCA is ≈ −350 pcm sitting on the bottom of the core.
        let malfRho = droppedRod ? -350e-5 : 0
        // BWRs have no soluble-boron chemical shim in normal operation.
        guard plant.params.hasBoron else {
            plant.params.externalReactivity = malfRho
            return
        }
        let vPrimary = 300.0   // m³ approximate primary coolant volume
        let dB = (borationRate * 150.0 - dilutionRate * boronPPM) / vPrimary * dt
        boronPPM = max(0, boronPPM + dB)
        // Differential boron worth from PlantParams (OpenMC-calibrated when
        // calibration.json is bundled); nominal 800 ppm → 0 reactivity.
        plant.params.externalReactivity = plant.params.boronWorthPcmPerPpm * 1e-5 * (boronPPM - 800.0) + malfRho
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

        // Rate-of-change annunciators (sampled every ~2 sim-s, smoothed so a
        // frame blip can't chatter the board). Power rate alarms only on the
        // way UP — a scram's plunge is already announced by the trip windows.
        if _rateT < 0 { _rateT = s.time; _ratePrev = (s.coolantTempK, pressureMPa, s.powerFraction * 100) }
        let rdt = s.time - _rateT
        if rdt >= 2.0 {
            let raw = ((s.coolantTempK - _ratePrev.tAvg) / rdt,
                       (pressureMPa - _ratePrev.p) / rdt,
                       (s.powerFraction * 100 - _ratePrev.pw) / rdt)
            // Stride-aware smoothing (τ ≈ 4 sim-s regardless of the stride):
            // at ×600 a step IS the stride (10 s), so a fixed alpha would slow
            // the annunciator 5× — the board must not depend on playback speed.
            let a = 1 - exp(-rdt / 4.0)
            _rate = ((1 - a) * _rate.dT + a * raw.0,
                     (1 - a) * _rate.dP + a * raw.1,
                     (1 - a) * _rate.dPw + a * raw.2)
            _rateT = s.time
            _ratePrev = (s.coolantTempK, pressureMPa, s.powerFraction * 100)
        }
        _check("T_RATE",  abs(_rate.dT) > 0.5,  2, "RCS T-AVG HIGH RATE OF CHANGE", isTrip: false)
        _check("P_RATE",  abs(_rate.dP) > 0.05, 2, "PZR PRESSURE HIGH RATE",        isTrip: false)
        _check("PW_RATE", _rate.dPw > 1.5,       2, "REACTOR POWER HIGH RATE",       isTrip: false)

        // LOOP / SBO ladder. Low RCS flow at power is a REACTOR TRIP (RCP bus
        // undervoltage / low-flow trip — this is what actually scrams the
        // plant when the grid goes away at power).
        _check("LOW_FLOW", !plant.params.naturalCirculation && omegaRCP < 0.85 && s.powerFraction > 0.3,
               1, "LOW RCS FLOW — TRIP", isTrip: true)
        _check("LOOP",     !(line1BreakerOpen == false || line2BreakerOpen == false),
               2, "LOSS OF OFFSITE POWER", isTrip: false)
        _check("DG_RUN",   dieselsRunning,      3, "DIESEL GENERATORS RUNNING",     isTrip: false)
        _check("SBO",      auxPower == .sbo && dieselTimer > 12, 1, "STATION BLACKOUT", isTrip: false)
        _check("SRV_LIFT", srvOpen,             2, "SG SAFETY VALVES LIFTING",      isTrip: false)
        _check("AFW_RUN",  afwRunning,          3, "AUX FEEDWATER RUNNING",         isTrip: false)
        _check("BATT_LOW", batteryMin < 30 && auxPower == .sbo, 1, "STATION BATTERIES LOW", isTrip: false)

        // Containment.
        _check("CTMT_HI",   containment.pressureKPa > 115, 2, "CONTAINMENT PRESSURE HIGH",              isTrip: false)
        _check("CTMT_HIHI", containment.pressureKPa > 135, 1, "CTMT PRESS HI-HI — SPRAY ACTUATED",      isTrip: false)
        _check("CTMT_SUMP", containment.sumpM3 > 20,       2, "CONTAINMENT SUMP LEVEL HIGH",            isTrip: false)

        // Malfunction windows.
        _check("SGTR_RAD",  sgtrLeakKgs > 0,   1, "AIR EJECTOR RADIATION HIGH — SGTR", isTrip: false)
        _check("ROD_DROP",  droppedRod,        1, "CONTROL ROD DROP",                  isTrip: false)
        _check("ROD_DEV",   stuckRod,          2, "ROD POSITION DEVIATION",            isTrip: false)
        _check("MSIV_CLSD", msivClosed,        1, "MAIN STEAM ISOLATION",              isTrip: false)
        _check("ATWS",      atwsFault && _alarmMap.values.contains(where: { $0.isTrip }) && !scrammed,
               1, "RPS FAILURE TO SCRAM — ATWS", isTrip: false)

        // Auto-trip on any unacknowledged trip — unless the RPS itself has
        // failed (the ATWS malfunction): then the board demands boration.
        if !atwsFault, _alarmMap.values.contains(where: { $0.isTrip }) {
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
                _newAlarmMsgs.append((id, msg))    // spoken-annunciator queue
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

    // MARK: — Startup physics (1/M · ECP · rod worth)

    /// Record 1/M points during a subcritical rod withdrawal: count rate ∝ n,
    /// baseline taken at the first sample; a fresh scram clears the plot.
    private func recordInvM() {
        if scrammed { invMPoints.removeAll(); _invMBaseN = 0; _invMLastRod = 1; return }
        let rod = snapshot.rodPosition
        guard snapshot.powerFraction < 0.10, rod < _invMLastRod - 0.01 else { return }
        _invMLastRod = rod
        let n = max(snapshot.powerFraction, 1e-9)
        if _invMBaseN == 0 { _invMBaseN = n }
        invMPoints.append((x: 1 - rod, invM: min(1.2, _invMBaseN / n)))
        if invMPoints.count > 80 { invMPoints.removeFirst() }
    }

    /// Integral rod-worth shape (calibrated table or analytic S-curve).
    func rodWorthShape(_ x: Double) -> Double { plant.params.rodShape(x) }

    /// Reactivity components [Δk/k] for the ECP block.
    var rhoComponents: (boron: Double, xenon: Double, mod: Double, dop: Double) {
        let p = plant.params
        return (p.externalReactivity,
                -p.xenonReactivityCoeff * snapshot.xenonInventory,
                p.coolantTempCoeff * (snapshot.coolantTempK - p.nominalCoolantTemp),
                p.fuelTempCoeff * (snapshot.fuelTempK - p.nominalFuelTemp))
    }

    /// Estimated critical position [SWD], from inverting the rod S-curve
    /// against the non-rod reactivity balance. nil = criticality precluded at
    /// any rod position (e.g. deep in the xenon peak) — wait or dilute.
    var ecpSWD: Int? {
        let p = plant.params
        let c = rhoComponents
        let others = c.boron + c.xenon + c.mod + c.dop
        let wNeeded = others / -p.rodWorth        // rodWorth < 0
        guard wNeeded >= 0 else { return 228 }    // excess without rods → out
        guard wNeeded <= 1 else { return nil }    // even fully inserted can't... precluded end
        var lo = 0.0, hi = 1.0
        for _ in 0..<40 {
            let mid = (lo + hi) / 2
            if plant.params.rodShape(mid) < wNeeded { lo = mid } else { hi = mid }
        }
        return Int((228 * (1 - (lo + hi) / 2)).rounded())
    }

    // Returns the historian as a time-ordered array starting from oldest
    func orderedHistory<T>(_ arr: [T]) -> [T] {
        let n = arr.count
        let start = _histIdx % n
        return Array(arr[start...]) + Array(arr[..<start])
    }
}
