import XCTest
@testable import ReactorSim

/// Containment response to a primary-side break: pressurization to the HI-HI
/// spray setpoint, spray knockdown well under design pressure, sump fill, and
/// recovery once the break is isolated.
final class ContainmentTests: XCTestCase {
    @MainActor
    func testSmallLOCAPressurizesSprayCapsAndRecovers() {
        let sup = PlantSupervisor()
        for _ in 0..<300 { sup.step(dt: 1) }
        XCTAssertLessThan(sup.containment.pressureKPa, 110, "normal ops ≈ atmospheric")

        sup.primaryLeakKgs = 50                          // small LOCA
        var sawHiHi = false, sawSpray = false, peak = 0.0
        for _ in 0..<1200 {
            sup.step(dt: 1)
            peak = max(peak, sup.containment.pressureKPa)
            sawHiHi = sawHiHi || sup.containment.pressureKPa > 135
            sawSpray = sawSpray || sup.containment.sprayOn
        }
        XCTAssertTrue(sawHiHi, "50 kg/s must reach the HI-HI setpoint (peak \(peak) kPa)")
        XCTAssertTrue(sawSpray, "spray must actuate on HI-HI")
        XCTAssertLessThan(peak, 400, "spray caps pressure under design")
        XCTAssertGreaterThan(sup.containment.sumpM3, 5, "break inventory collects in the sump")

        sup.primaryLeakKgs = 0                           // break isolated
        for _ in 0..<1800 { sup.step(dt: 1) }
        XCTAssertLessThan(sup.containment.pressureKPa, 120, "coolers recover the atmosphere")
    }
}
