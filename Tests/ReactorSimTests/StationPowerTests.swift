import XCTest
@testable import ReactorSim

/// LOOP / station-blackout ladder: grid loss islands to house load; losing the
/// generator too cranks the diesels; RCPs coast to natural circulation; the
/// low-flow trip scrams the plant at power; AFW + SG safety valves carry the
/// decay heat with the condenser vacuum gone.
final class StationPowerTests: XCTestCase {

    /// Opening both lines while generating must NOT trip the unit — it islands
    /// onto house load (governor runback) with the RCPs still powered.
    @MainActor
    func testDualLineOpenIslandsToHouseLoad() {
        let sup = PlantSupervisor()
        for _ in 0..<600 { sup.step(dt: 1) }
        sup.rodAutoEnabled = true                     // graceful runback
        sup.toggleLineBreaker(0)
        sup.toggleLineBreaker(1)
        for _ in 0..<300 { sup.step(dt: 1) }
        XCTAssertEqual(sup.auxPower, .houseLoad)
        XCTAssertFalse(sup.snapshot.scrammed, "islanding must not scram the plant")
        XCTAssertGreaterThan(sup.omegaRCP, 0.95, "RCPs stay powered on house load")
        XCTAssertLessThan(sup.snapshot.powerFraction, 0.45, "governor ran the unit back")
    }

    /// Full-power loss of ALL station AC: low-flow trip, diesel start, natural
    /// circulation + AFW + SRVs hold the plant safe for 30 minutes.
    @MainActor
    func testStationLossAtPowerTripsAndRidesThrough() {
        let sup = PlantSupervisor()
        for _ in 0..<600 { sup.step(dt: 1) }
        sup.toggleLineBreaker(0)
        sup.toggleLineBreaker(1)
        sup.turbineTrip = true                        // gen gone too → no AC
        for _ in 0..<20 { sup.step(dt: 1) }
        XCTAssertTrue(sup.snapshot.scrammed, "low-flow trip must scram at power")
        XCTAssertTrue(sup.trips.contains { $0.contains("LOW RCS FLOW") })
        XCTAssertTrue(sup.dieselsRunning, "EDGs pick up within ~10 s")
        XCTAssertEqual(sup.auxPower, .diesel)

        var sawSRV = false, sawAFW = false
        for _ in 0..<1800 {
            sup.step(dt: 1)
            sawSRV = sawSRV || sup.srvOpen
            sawAFW = sawAFW || sup.afwRunning
        }
        XCTAssertLessThan(sup.omegaRCP, 0.1, "RCPs coasted down")
        XCTAssertTrue(sawAFW, "aux feedwater must carry the SG")
        XCTAssertTrue(sawSRV, "SG safety valves must lift with the condenser gone")
        XCTAssertGreaterThan(sup.feedwaterInv, 0.05, "AFW keeps SG inventory")
        XCTAssertLessThan(sup.snapshot.fuelTempK, 1500, "no fuel damage on natural circulation")
        XCTAssertLessThan(sup.snapshot.sgTempK, 600, "SRVs bound the SG temperature")
    }
}
