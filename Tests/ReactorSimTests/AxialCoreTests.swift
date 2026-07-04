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

/// Regression: the startup sequencer after a scram+reset must ascend to power
/// WITHOUT riding the period into the HIGH FLUX trip (user-reported: "auto
/// can't handle startup after scram — it scrams due to too high power").
/// Guards: startup-rate hold on rod withdrawal + the rod-auto flux limiter.
final class StartupAfterScramTests: XCTestCase {
    @MainActor
    func testAutoStartupAfterScramNoRetrip() {
        let sup = PlantSupervisor()
        for _ in 0..<600 { sup.step(dt: 1) }                 // settle at power
        sup.triggerScram()
        for _ in 0..<180 { sup.step(dt: 1) }                 // rods in, plant cooling
        sup.acknowledgeAllAlarms()
        for _ in 0..<10 { sup.step(dt: 1) }                  // acked+cleared windows drop
        sup.resetScram()
        sup.step(dt: 1)     // snapshot refreshes on the next step
        XCTAssertFalse(sup.snapshot.scrammed, "scram reset should be approved: \(sup.scramMessage)")

        sup.startupPermit = true
        sup.autoStartup = true
        var maxPf = 0.0
        var retripped = false
        for _ in 0..<5400 {
            sup.step(dt: 1)
            maxPf = max(maxPf, sup.snapshot.powerFraction)
            if sup.snapshot.scrammed { retripped = true; break }
        }
        XCTAssertFalse(retripped, "auto startup re-scrammed (peak \(maxPf * 100)% RTP)")
        XCTAssertLessThan(maxPf, 1.19, "power overshoot must stay clear of the 120% trip")
        XCTAssertGreaterThan(sup.snapshot.powerFraction, 0.85,
                             "ascension should reach power (got \(sup.snapshot.powerFraction * 100)%)")
    }
}
