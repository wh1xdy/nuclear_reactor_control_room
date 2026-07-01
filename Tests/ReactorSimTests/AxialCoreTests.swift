import XCTest
@testable import ReactorSim

/// Fast sanity checks for the 1-D axial nodal core. (The full multi-day
/// axial-xenon-oscillation study lives outside CI; these guard the static
/// behaviour and the 0-D invariants.)
final class AxialCoreTests: XCTestCase {

    private func settle() -> (ReactorPlant, PlantSnapshot) {
        let plant = ReactorPlant()
        let ctrl = ControlInputs(rodPosition: 0, primaryFlow: 1, turbineValve: 1, scram: false)
        var snap = plant.step(dt: 1, ctrl: ctrl)
        for _ in 0..<600 { snap = plant.step(dt: 1, ctrl: ctrl) }   // 10 min
        return (plant, snap)
    }

    /// Nominal full-power shape is sane (no collapse): a chopped-cosine peaking,
    /// a positive DNBR margin, a credible peak clad temperature, near-zero ΔI.
    func testNominalShapeSane() {
        let (_, s) = settle()
        XCTAssertEqual(s.powerFraction, 1.0, accuracy: 0.02, "0-D power unchanged")
        XCTAssertTrue((1.05...1.9).contains(s.fz),   "Fz \(s.fz) out of range")
        XCTAssertTrue((1.6...3.2).contains(s.minDNBR), "DNBR \(s.minDNBR) out of range")
        XCTAssertTrue((600...720).contains(s.peakCladTempK), "clad \(s.peakCladTempK) out of range")
        XCTAssertLessThan(abs(s.axialOffsetPct), 25, "|ΔI| \(s.axialOffsetPct) too large at nominal")
        XCTAssertEqual(s.axialProfile.reduce(0,+) / Double(s.axialProfile.count), 1.0, accuracy: 0.01,
                       "axial profile must be mean-normalised")
    }

    /// Inserting rods from the top drives the flux peak DOWN (ΔI more negative)
    /// and sharpens the axial peak (Fz up) — the real shape coupling, not a proxy.
    /// (DNBR actually IMPROVES here: the peak moves to the cooler, more-subcooled
    /// bottom and total power drops — so we only require it to stay sane.)
    func testRodInsertionShiftsFluxDown() {
        let (plant, s0) = settle()
        let ctrl = ControlInputs(rodPosition: 0.30, primaryFlow: 1, turbineValve: 1, scram: false)
        var s = s0
        for _ in 0..<40 { s = plant.step(dt: 10, ctrl: ctrl) }      // step rods to 30%
        XCTAssertLessThan(s.axialOffsetPct, s0.axialOffsetPct - 3, "ΔI should swing negative on top rod insertion")
        XCTAssertGreaterThan(s.fz, s0.fz, "top rodding should sharpen the axial peak")
        XCTAssertGreaterThan(s.minDNBR, 1.3, "DNBR must remain a sane positive margin")
    }
}
