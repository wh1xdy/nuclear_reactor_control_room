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
