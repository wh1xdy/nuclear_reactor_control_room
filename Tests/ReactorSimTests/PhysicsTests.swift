// PhysicsTests.swift — headless verification of the reactor physics.
// Run with `swift test`. These are the checks RECAP.md refers to.

import XCTest
@testable import ReactorSim

final class DecayHeatTests: XCTestCase {

    func testEquilibriumFractionIsANS() {
        // Σ(α_i/λ_i)/200 MeV ≈ 6.5% for U-235 thermal
        XCTAssertEqual(DecayHeat.equilibriumFraction, 0.0654, accuracy: 0.003)
    }

    func testBuildsUpDuringOperation() {
        var dh = DecayHeat()
        XCTAssertEqual(dh.fraction, 0)
        // 1 day at full power (exact integrator → big dt is fine)
        for _ in 0..<864 { _ = dh.step(dt: 100, n: 1.0) }
        XCTAssertGreaterThan(dh.fraction, 0.04)            // most of equilibrium
        XCTAssertLessThanOrEqual(dh.fraction, DecayHeat.equilibriumFraction)
    }

    func testDecaysAfterShutdown() {
        // The old power-law stub froze decay heat forever after scram.
        var dh = DecayHeat()
        for _ in 0..<864 { _ = dh.step(dt: 100, n: 1.0) }   // 1 day operation
        let atShutdown = dh.fraction
        _ = dh.step(dt: 60, n: 0)                           // 1 min after
        let after1min = dh.fraction
        _ = dh.step(dt: 3540, n: 0)                         // 1 h after
        let after1h = dh.fraction
        XCTAssertLessThan(after1min, atShutdown)
        XCTAssertLessThan(after1h, after1min)
        XCTAssertGreaterThan(after1h, 0.005)                // but never zero quickly
    }
}

final class PlantTests: XCTestCase {

    /// Steady full power: 10 sim-minutes at nominal controls.
    func testSteadyStateFullPower() {
        let plant = PWRPlant()
        let ctrl  = ControlInputs(rodPosition: 0, primaryFlow: 1, turbineValve: 1, scram: false)
        var snap  = plant.step(dt: 1, ctrl: ctrl)
        for _ in 0..<600 { snap = plant.step(dt: 1, ctrl: ctrl) }

        XCTAssertTrue(snap.powerFraction.isFinite)
        XCTAssertEqual(snap.powerFraction, 1.0, accuracy: 0.03)   // ~100% rated thermal
        XCTAssertEqual(snap.fuelTempK,    900, accuracy: 15)
        XCTAssertEqual(snap.coolantTempK, 550, accuracy: 10)
        XCTAssertEqual(snap.electricPowerW / 1e6, 990, accuracy: 40)
        XCTAssertFalse(snap.scrammed)
    }

    /// Scram: power collapses to decay-heat levels and decay heat then FALLS.
    func testScramShutdownAndDecay() {
        let plant = PWRPlant()
        var ctrl  = ControlInputs(rodPosition: 0, primaryFlow: 1, turbineValve: 1, scram: false)
        for _ in 0..<600 { _ = plant.step(dt: 1, ctrl: ctrl) }    // 10 min at power

        ctrl.scram = true
        var snap = plant.step(dt: 10, ctrl: ctrl)
        XCTAssertLessThan(snap.powerFraction, 0.10)               // fission gone in 10 s
        let decayAt10s = snap.decayHeatFraction
        XCTAssertGreaterThan(decayAt10s, 0.01)                    // decay heat present

        for _ in 0..<60 { snap = plant.step(dt: 10, ctrl: ctrl) } // +10 min
        XCTAssertLessThan(snap.decayHeatFraction, decayAt10s)     // ANS curve falls
        XCTAssertGreaterThan(snap.decayHeatFraction, 0)
    }

    /// CRDM rate limit: rods walk toward demand, never teleport.
    func testRodMotionIsRateLimited() {
        let plant = PWRPlant()
        let ctrl  = ControlInputs(rodPosition: 1.0, primaryFlow: 1, turbineValve: 1, scram: false)
        let snap  = plant.step(dt: 10, ctrl: ctrl)                // demand a full insertion
        // 10 s × 0.0053/s ≈ 0.053 — nowhere near the demanded 1.0
        XCTAssertEqual(snap.rodPosition, 0.053, accuracy: 0.01)
    }
}

@MainActor
final class SupervisorTests: XCTestCase {

    /// BOP must stay stable and balanced at 600× time compression
    /// (dt = 10 s/frame — the old explicit-Euler condenser diverged here).
    func testBOPStableAt600x() {
        let sup = PlantSupervisor()
        for _ in 0..<600 { sup.step(dt: 10.0) }                   // 100 sim-minutes

        XCTAssertTrue(sup.condTempK.isFinite)
        XCTAssertEqual(sup.condTempK, 311, accuracy: 15)          // settled, not banging clamps
        XCTAssertGreaterThan(sup.feedwaterInv, 0.4)               // no steady-state drain
        XCTAssertGreaterThan(sup.steamInv, 0.5)                   // self-regulating balance
        XCTAssertLessThan(sup.steamInv, 1.5)
        // No spurious alarms during normal full-power operation
        XCTAssertTrue(sup.trips.isEmpty, "unexpected trips: \(sup.trips)")
    }

    /// Annunciator latching: a cleared alarm stays on the board until acknowledged.
    func testAlarmLatchesUntilAcknowledged() {
        let sup = PlantSupervisor()
        // Drain feedwater via fault until LOW_FEED comes in
        sup.feedwaterFault = true
        for _ in 0..<20 { sup.step(dt: 1.0) }
        XCTAssertTrue(sup.alarms.contains { $0.id == "LOW_FEED" })

        // Clear the fault and refill above the setpoint
        sup.feedwaterFault = false
        sup.feedwaterValve = 1.0
        for _ in 0..<60 { sup.step(dt: 1.0) }
        XCTAssertGreaterThan(sup.feedwaterInv, 0.1)               // condition cleared
        XCTAssertTrue(sup.alarms.contains { $0.id == "LOW_FEED" },
                      "alarm must latch until acknowledged")

        sup.acknowledgeAllAlarms()
        sup.step(dt: 1.0)
        XCTAssertFalse(sup.alarms.contains { $0.id == "LOW_FEED" })
    }

    /// Full cycle: trip → reset → AUTO STARTUP sequencer returns the plant to
    /// power with no protection trips along the way.
    func testAutoStartupSequencer() {
        let sup = PlantSupervisor()
        for _ in 0..<30 { sup.step(dt: 1.0) }
        sup.triggerScram()
        for _ in 0..<40 { sup.step(dt: 1.0) }                     // rods in, dump holding
        sup.acknowledgeAllAlarms()
        sup.step(dt: 1.0)
        sup.resetScram()
        XCTAssertFalse(sup.scrammed, sup.scramMessage)

        sup.autoStartup = true
        // 1000 s sim at 10 s frames: withdrawal (~190 s) + ascension (~200 s)
        for _ in 0..<100 { sup.step(dt: 10.0) }
        XCTAssertFalse(sup.scrammed, "sequencer must not trip the plant: \(sup.trips)")
        XCTAssertGreaterThan(sup.snapshot.powerFraction, 0.90,
                             "phase: \(sup.startupPhase)")
        XCTAssertTrue(sup.rodAutoEnabled)                          // handed over
        XCTAssertEqual(sup.snapshot.coolantTempK, 550, accuracy: 8)
    }

    /// Scram reset syncs the rod DEMAND to the inserted position — no
    /// reactivity step from a stale withdrawn lever.
    func testScramResetSyncsRodDemand() {
        let sup = PlantSupervisor()
        for _ in 0..<10 { sup.step(dt: 1.0) }
        sup.triggerScram()
        for _ in 0..<30 { sup.step(dt: 1.0) }                     // rods drop (τ = 2 s)
        XCTAssertGreaterThan(sup.snapshot.rodPosition, 0.94)

        sup.acknowledgeAllAlarms()                                // clear any latched trips
        sup.step(dt: 1.0)
        sup.rodPosition = 0.0                                     // stale withdrawn demand
        sup.resetScram()
        XCTAssertFalse(sup.scrammed, sup.scramMessage)
        XCTAssertGreaterThan(sup.rodPosition, 0.94, "demand must follow actual position")
    }
}
