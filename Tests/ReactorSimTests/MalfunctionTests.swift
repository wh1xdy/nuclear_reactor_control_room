import XCTest
@testable import ReactorSim

/// Instructor-menu faults: each must produce its signature transient.
final class MalfunctionTests: XCTestCase {

    @MainActor
    func testATWSBlocksScramUntilBoration() {
        let sup = PlantSupervisor()
        for _ in 0..<300 { sup.step(dt: 1) }
        sup.atwsFault = true
        sup.triggerScram()
        sup.step(dt: 1)
        XCTAssertFalse(sup.snapshot.scrammed, "ATWS must defeat the scram button")
        // Emergency boration drives power down without rods (against the
        // still-open turbine the descent is demand-limited — assert progress).
        sup.borationRate = 1.0
        for _ in 0..<900 { sup.step(dt: 1) }
        XCTAssertLessThan(sup.snapshot.powerFraction, 0.6,
                          "boration must be shutting the core down (got \(sup.snapshot.powerFraction))")
    }

    @MainActor
    func testDroppedRodTiltsAndDepressesPower() {
        let sup = PlantSupervisor()
        for _ in 0..<300 { sup.step(dt: 1) }
        let q0 = sup.snapshot.qptr
        sup.droppedRod = true
        // Prompt dip, then reactor-follows-turbine recovery — assert the DIP,
        // the persistent QPTR tilt, and no trip.
        var minPf = 1.0
        for _ in 0..<120 { sup.step(dt: 1); minPf = min(minPf, sup.snapshot.powerFraction) }
        XCTAssertLessThan(minPf, 0.93, "−350 pcm must visibly dip power (min \(minPf))")
        XCTAssertGreaterThan(sup.snapshot.qptr, q0 + 0.01, "a dropped rod must persist in QPTR")
        XCTAssertFalse(sup.snapshot.scrammed)
    }

    @MainActor
    func testSGTRRaisesSGLevelAndAlarm() {
        let sup = PlantSupervisor()
        for _ in 0..<300 { sup.step(dt: 1) }
        let inv0 = sup.feedwaterInv
        sup.sgtrLeakKgs = 15
        sup.fwAutoEnabled = false
        for _ in 0..<300 { sup.step(dt: 1) }
        XCTAssertGreaterThan(sup.feedwaterInv, inv0 + 0.01, "SG fills from the primary side")
        XCTAssertTrue(sup.alarms.contains { $0.id == "SGTR_RAD" }, "radiation window must raise")
    }

    @MainActor
    func testMSIVClosureTripsTurbineAndLiftsSRVs() {
        let sup = PlantSupervisor()
        for _ in 0..<300 { sup.step(dt: 1) }
        sup.msivClosed = true
        var sawSRV = false
        for _ in 0..<600 { sup.step(dt: 1); sawSRV = sawSRV || sup.srvOpen }
        XCTAssertTrue(sup.turbineTrip, "no steam → turbine trips")
        XCTAssertTrue(sawSRV, "SRVs are the only relief path with the MSIVs shut")
    }
}

/// Four-loop primary: single-pump loss behaviour.
final class FourLoopTests: XCTestCase {

    /// Losing ONE of four RCPs at full power must reach the low-flow trip.
    @MainActor
    func testSinglePumpTripAtPowerScrams() {
        let sup = PlantSupervisor()
        for _ in 0..<300 { sup.step(dt: 1) }
        XCTAssertEqual(sup.rcpCount, 4)
        sup.rcpRunning[2] = false
        var tripped = false
        for _ in 0..<60 { sup.step(dt: 1); if sup.snapshot.scrammed { tripped = true; break } }
        XCTAssertTrue(tripped, "one-of-four pump loss at 100% must trip on low flow")
    }

    /// The same loss at LOW power is ride-through: three pumps carry the plant.
    @MainActor
    func testSinglePumpTripAtLowPowerRidesThrough() {
        let sup = PlantSupervisor()
        for _ in 0..<300 { sup.step(dt: 1) }
        sup.rodAutoEnabled = false
        sup.turbineValve = 0.2                        // run back to ~20% steam demand
        for _ in 0..<600 { sup.step(dt: 1) }
        guard sup.snapshot.powerFraction < 0.3 else {
            return XCTFail("setup: power should be below the low-flow trip gate (got \(sup.snapshot.powerFraction))")
        }
        sup.rcpRunning[0] = false
        for _ in 0..<300 { sup.step(dt: 1) }
        XCTAssertFalse(sup.snapshot.scrammed, "3/4 pumps at low power is a ride-through")
        XCTAssertLessThan(sup.omegaRCP, 0.85)
        XCTAssertGreaterThan(sup.omegaRCP, 0.7)
    }
}

/// Startup kit: ECP sanity.
final class StartupKitTests: XCTestCase {
    @MainActor
    func testECPExistsAtSteadyState() {
        let sup = PlantSupervisor()
        for _ in 0..<600 { sup.step(dt: 1) }
        let ecp = sup.ecpSWD
        XCTAssertNotNil(ecp, "steady plant must have a critical position")
        XCTAssertGreaterThan(ecp ?? 0, 150, "near-rods-out core → ECP high in SWD")
        let c = sup.rhoComponents
        XCTAssertLessThan(abs(c.boron + c.xenon + c.mod + c.dop) * 1e5, 300,
                          "components ≈ balance at steady state")
    }
}
