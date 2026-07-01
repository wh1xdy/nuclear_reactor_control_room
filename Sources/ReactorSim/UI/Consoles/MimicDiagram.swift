// MimicDiagram.swift — detailed full-screen PWR plant mimic.
// A purpose-built schematic (not the small panel P&ID) with proper component
// shapes, level indication, valves, flow animation, and live values embedded
// at every system — sized to fill the whole console.
//
// Layout philosophy: the process flows left→right and uses the WHOLE canvas.
//   • Primary island (left): pressure vessel + pressurizer + RCP, full height.
//   • Steam generator (centre-left): tall U-tube vessel, primary↔secondary.
//   • Turbine hall (upper right): main steam → HP → MSR → LP → generator → grid.
//   • Heat sink + feed train (lower band): condenser → CP → heaters → FW pump → SG.
// Every node carries a live ISA-style readout so the empty field is gone.

import SwiftUI

struct MimicDiagram: View {
    let supervisor: PlantSupervisor

    var body: some View {
        GeometryReader { geo in
            TimelineView(.animation) { context in
                let t = context.date.timeIntervalSinceReferenceDate
                Canvas { ctx, size in
                    draw(ctx, size, snap: supervisor.snapshot, sup: supervisor, t: t)
                }
            }
            .contentShape(Rectangle())
            .gesture(SpatialTapGesture().onEnded { ev in
                handleTap(ev.location, size: geo.size)
            })
        }
        .background(Theme.schematicBg)
    }

    /// Operator taps on the mimic. Currently the switchyard breakers are live —
    /// hit boxes must track the positions drawn in drawSwitchyard (52G / feeders).
    private func handleTap(_ p: CGPoint, size: CGSize) {
        let w = size.width, h = size.height
        let hitR: CGFloat = 15
        let g52 = CGPoint(x: 0.852 * w, y: 0.148 * h)   // 52G generator breaker
        let l1  = CGPoint(x: 0.873 * w, y: 0.085 * h)   // LINE 1 breaker
        let l2  = CGPoint(x: 0.938 * w, y: 0.085 * h)   // LINE 2 breaker
        if      hypot(p.x - g52.x, p.y - g52.y) < hitR { supervisor.toggleGenBreaker() }
        else if hypot(p.x - l1.x,  p.y - l1.y)  < hitR { supervisor.toggleLineBreaker(0) }
        else if hypot(p.x - l2.x,  p.y - l2.y)  < hitR { supervisor.toggleLineBreaker(1) }
    }

    private func draw(_ ctx: GraphicsContext, _ size: CGSize,
                      snap: PlantSnapshot, sup: PlantSupervisor, t: Double) {
        let W = size.width, H = size.height
        func fx(_ v: CGFloat) -> CGFloat { v * W }
        func fy(_ v: CGFloat) -> CGFloat { v * H }
        func P(_ x: CGFloat, _ y: CGFloat) -> CGPoint { CGPoint(x: fx(x), y: fy(y)) }

        let scram   = snap.scrammed
        let flowF   = scram ? 0 : max(0, Double(sup.primaryFlow) * Double(sup.omegaRCP))
        let steamF  = scram ? 0 : Double(sup.turbineValve)
        let pf      = max(0, snap.powerFraction)

        // ── Derived plant numbers (calibrated display values) ────────────────
        let rcsKgs   = 18_800.0 * flowF                       // total RCS flow
        let steamKgs = scram ? 0 : 1_650.0 * pf               // main steam flow ∝ power
        // Feed delivery: the FW REG valve meters rated flow at its 70% balanced
        // point, so 0.70 → 1650 kg/s ≈ steam. Under/overfeeding now diverges from
        // steam AND is integrated into the S/G level below — they're coupled.
        let fwKgs    = sup.feedwaterFault ? 0 : 1_650.0 * Double(sup.feedwaterValve) / 0.70
        let deltaT   = snap.hotLegTempK - snap.coldLegTempK
        // Pressurizer level follows the T-avg program (~25% no-load → ~60% at the
        // 550 K full-power band), NOT pressure — level and pressure are distinct.
        let pzrLvl   = max(0.05, min(0.95, 0.25 + 0.35 * (snap.coolantTempK - 530) / 20))
        // S/G narrow-range level tracks the feed↔steam mass balance (feedwaterInv
        // integrates FW in − steam out): underfeed → level falls → LOW FEED. No
        // longer decoupled from the flows the operator sees.
        let sgLvl    = max(0.08, min(0.92, sup.feedwaterInv * 0.62))
        let rpm      = (scram || sup.turbineTrip) ? 0 : 3000
        let condKPa  = magnusKPa(sup.condTempK)               // condenser backpressure

        // ════════════════════════════════════════════════════════════════════
        // COMPONENT FRAMES
        // ════════════════════════════════════════════════════════════════════
        let vessel = CGRect(x: fx(0.050), y: fy(0.34), width: fx(0.070), height: fy(0.44))
        let hotY   = vessel.minY + fy(0.07)
        let coldY  = vessel.maxY - fy(0.10)
        let pzr    = CGRect(x: fx(0.158), y: fy(0.085), width: fx(0.040), height: fy(0.185))
        let rcp    = CGRect(x: fx(0.170), y: coldY - fy(0.045), width: fx(0.052), height: fy(0.090))
        let sg     = CGRect(x: fx(0.268), y: fy(0.155), width: fx(0.086), height: fy(0.585))
        let sgHotIn   = CGPoint(x: sg.minX, y: sg.maxY - fy(0.10))
        let sgColdOut = CGPoint(x: sg.minX, y: coldY)                       // straight wall, at the RCP centerline (no dog-leg)
        let sgSteam   = CGPoint(x: sg.maxX, y: sg.minY + sg.height * 0.26)  // below taperY (0.20·h) → lands on the straight wall, not the dome taper
        let sgFeed    = CGPoint(x: sg.maxX, y: sg.minY + fy(0.20))

        // Turbine casings share one shaft centerline (hpT.midY = fy(0.205)) so the
        // HP→MSR→LP crossovers and the steam header run dead level.
        let hpT  = CGRect(x: fx(0.470), y: fy(0.150),  width: fx(0.072), height: fy(0.110))
        let msr  = CGRect(x: fx(0.556), y: fy(0.1675), width: fx(0.028), height: fy(0.075))
        let lpT  = CGRect(x: fx(0.598), y: fy(0.115),  width: fx(0.110), height: fy(0.180))
        let msiv = CGPoint(x: fx(0.405), y: hpT.midY)                       // on the level steam run
        let gen  = CGRect(x: fx(0.730), y: fy(0.160), width: fx(0.075), height: fy(0.100))

        let cond  = CGRect(x: fx(0.620), y: fy(0.470), width: fx(0.158), height: fy(0.130))
        let feedY = fy(0.745)
        let cp    = CGPoint(x: fx(0.635), y: feedY)        // condensate pump
        let lpFwh = CGRect(x: fx(0.553), y: feedY - fy(0.030), width: fx(0.054), height: fy(0.060))
        let dae   = CGRect(x: fx(0.492), y: feedY - fy(0.034), width: fx(0.056), height: fy(0.068))
        let fp    = CGPoint(x: fx(0.452), y: feedY)        // main feed pump
        let hpFwh = CGRect(x: fx(0.388), y: feedY - fy(0.030), width: fx(0.054), height: fy(0.060))

        // ════════════════════════════════════════════════════════════════════
        // PIPING  (drawn first so equipment sits on top)
        // ════════════════════════════════════════════════════════════════════
        // Primary hot leg: vessel hot nozzle → over → down into SG inlet.
        let hot = [(vessel.maxX, hotY), (fx(0.232), hotY), (fx(0.232), sgHotIn.y), (sgHotIn.x, sgHotIn.y)]
        // Primary cold leg: SG outlet → across through RCP → vessel cold nozzle.
        let coldA = [(sgColdOut.x, sgColdOut.y), (rcp.maxX, sgColdOut.y), (rcp.maxX, coldY)]
        let coldB = [(rcp.minX, coldY), (vessel.maxX, coldY)]
        // Pipe color = the medium it carries / its own temperature. The hot leg
        // renders calm teal at ~565 K and only ramps warm if it heats abnormally;
        // the cold leg reads blue. Color is the instrument, not decoration.
        let cHot   = Theme.colorFor(snap.hotLegTempK, trip: 620)
        let cCold  = Theme.colorFor(snap.coldLegTempK, trip: 620)
        let cSteam = steamF > 0.05 ? Theme.fluidTwoPhase : Theme.fluidExhaust
        pipe(ctx, hot,   cHot,  4)
        pipe(ctx, coldA, cCold, 4)
        pipe(ctx, coldB, cCold, 4)
        // Surge line PZR → hot leg.
        pipe(ctx, [(pzr.midX, pzr.maxY), (pzr.midX, hotY)], Theme.ink.opacity(0.45), 2)

        // Main steam: short stub off the nozzle (so the riser clears the tapered
        // dome instead of running up alongside it), then up to the shaft line, level to HP.
        let steamRiserX = sgSteam.x + fx(0.018)
        let steamHdr = [(sgSteam.x, sgSteam.y), (steamRiserX, sgSteam.y), (steamRiserX, hpT.midY), (msiv.x, hpT.midY), (hpT.minX, hpT.midY)]
        pipe(ctx, steamHdr, cSteam, 3)
        // HP exhaust → MSR → LP (lower-energy exhaust steam) — all on the shaft line.
        pipe(ctx, [(hpT.maxX, hpT.midY), (msr.minX, msr.midY)], Theme.fluidExhaust, 2)
        pipe(ctx, [(msr.maxX, msr.midY), (lpT.minX, lpT.midY)], Theme.fluidExhaust, 2)
        // LP exhaust → condenser (start on the sloped casing bottom at midX, not below it).
        pipe(ctx, [(lpT.midX, lpT.maxY - lpT.height * 0.30 * 0.5), (lpT.midX, cond.minY)], Theme.fluidExhaust.opacity(0.85), 3)
        // Shaft LP → generator (the gen→switchyard lead is drawn with the switchyard).
        pipe(ctx, [(lpT.maxX, lpT.midY), (gen.minX, gen.midY)], Theme.ink.opacity(0.55), 3)

        // Feed / condensate train: subcooled condensate (blue) becomes pumped
        // feedwater (cyan-blue) after the FW pump — three readable secondary media.
        // Drain taps the hotwell at the condenser centre.
        pipe(ctx, [(cond.midX, cond.maxY - fy(0.02)), (cond.midX, feedY), (cp.x, feedY)], Theme.fluidSubcooled, 3)
        pipe(ctx, [(cp.x, feedY), (lpFwh.maxX, feedY)], Theme.fluidSubcooled, 3)
        pipe(ctx, [(lpFwh.minX, feedY), (dae.maxX, feedY)], Theme.fluidSubcooled, 3)
        pipe(ctx, [(dae.minX, feedY), (fp.x, feedY)], Theme.fluidSubcooled, 3)
        pipe(ctx, [(fp.x, feedY), (hpFwh.maxX, feedY)], Theme.fluidFeedwater, 3)
        // FW riser: off the HP-heater LEFT face → out clear of the box → up → S/G feed nozzle.
        let fwRiserX = hpFwh.minX - fx(0.014)
        pipe(ctx, [(hpFwh.minX, feedY), (fwRiserX, feedY), (fwRiserX, sgFeed.y), (sgFeed.x, sgFeed.y)], Theme.fluidFeedwater, 3)

        // Steam dump (post-trip heat sink). Taps the header at the MSIV, drops down
        // LEFT of the STEAM CYCLE dock, crosses ABOVE it, then into the condenser via
        // the clear gap (dumpX) between the dock and the condenser — never over the panel.
        if sup.steamDumpValve > 0.001 {
            // Clean 3-bend: down off the header (left of the dock), level above the
            // dock, then straight down into the condenser LEFT of the LP-exhaust drop
            // — no stair-step, no crossing the LP exhaust, no manifold.
            let dumpX = cond.minX + fx(0.012)
            let dumpTopY = fy(0.40)
            pipe(ctx, [(msiv.x, hpT.midY), (msiv.x, dumpTopY), (dumpX, dumpTopY), (dumpX, cond.minY)],
                 Theme.steam.opacity(0.7), 2)
            txt(ctx, "STM DUMP", CGPoint(x: msiv.x + 5, y: (hpT.midY + dumpTopY)/2), .leading, Theme.caution, 8)
        }

        // ── Flow animation: dash speed AND spacing encode mass flow (q). At q≈0
        // the dashes are sparse and nearly frozen, so loss-of-flow reads instantly.
        flowDash(ctx, hot,   cHot,  flowF, t)
        flowDash(ctx, coldA, cCold, flowF, t)
        flowDash(ctx, coldB, cCold, flowF, t)
        flowDash(ctx, steamHdr, cSteam, steamF, t)
        flowDash(ctx, [(cond.midX, cond.maxY - fy(0.02)), (cond.midX, feedY), (fwRiserX, feedY),
                       (fwRiserX, sgFeed.y), (sgFeed.x, sgFeed.y)], Theme.fluidFeedwater,
                 sup.feedwaterFault ? 0 : Double(sup.feedwaterValve), t)

        // ── Valves ────────────────────────────────────────────────────────────
        valve(ctx, msiv, open: steamF > 0.02)

        // ════════════════════════════════════════════════════════════════════
        // PRIMARY ISLAND
        // ════════════════════════════════════════════════════════════════════
        let fuelAlarm = snap.fuelTempK > 1200
        vesselRPV(ctx, vessel, snap: snap, alarm: fuelAlarm ? Theme.alarm : nil)
        // nozzle stubs
        nozzle(ctx, CGPoint(x: vessel.maxX, y: hotY),  Theme.hotLeg)
        nozzle(ctx, CGPoint(x: vessel.maxX, y: coldY), Theme.water)

        // Reactor data block under the vessel (bottom is clear).
        let pcx = vessel.midX
        reading(ctx, CGPoint(x: pcx, y: vessel.maxY + fy(0.050)),
                String(format: "%.1f%%", pf * 100),
                pf > 1.1 ? Theme.alarm : Theme.ink, 15, .center)
        txt(ctx, "RTP", CGPoint(x: pcx, y: vessel.maxY + fy(0.078)), .center, Theme.textDim, 8)
        txt(ctx, String(format: "FUEL  %.0f K", snap.fuelTempK),
            CGPoint(x: pcx, y: vessel.maxY + fy(0.105)), .center,
            Theme.colorForFuel(snap.fuelTempK), 9)
        txt(ctx, "\(Int((228 * (1 - snap.rodPosition)).rounded())) SWD",
            CGPoint(x: pcx, y: vessel.maxY + fy(0.128)), .center, Theme.textDim, 9)

        // Pressurizer with level + heaters + PORV.
        let pressAlarm = sup.pressureMPa > 17
        pressurizer(ctx, pzr, level: pzrLvl, heatersOn: sup.pzrAutoEnabled || sup.pressureMPa < 15.4,
                    alarm: pressAlarm ? Theme.alarm : nil)
        reading(ctx, CGPoint(x: pzr.maxX + fx(0.010), y: pzr.minY + fy(0.045)),
                String(format: "%.2f", sup.pressureMPa), pressAlarm ? Theme.alarm : Theme.ink, 12, .leading)
        txt(ctx, "MPa", CGPoint(x: pzr.maxX + fx(0.010), y: pzr.minY + fy(0.073)), .leading, Theme.textDim, 8)
        txt(ctx, String(format: "LVL %.0f%%", pzrLvl * 100),
            CGPoint(x: pzr.maxX + fx(0.010), y: pzr.minY + fy(0.103)), .leading, Theme.textDim, 9)
        if sup.porvOpen {
            txt(ctx, "● PORV OPEN", CGPoint(x: pzr.maxX + fx(0.010), y: pzr.minY + fy(0.130)), .leading, Theme.caution, 9)
        }

        // RCP on the cold leg.
        equipBox(ctx, rcp, "", sup.pumpDegraded ? Theme.alarm : nil)
        pump(ctx, CGPoint(x: rcp.midX, y: rcp.midY), r: min(rcp.width, rcp.height) * 0.30,
             running: flowF > 0.1, tint: sup.pumpDegraded ? Theme.alarm : Theme.accent, mirrored: true)
        txt(ctx, "RCP-1", CGPoint(x: rcp.midX, y: rcp.minY - fy(0.020)), .center, Theme.textHdr, 9)
        txt(ctx, String(format: "%.0f%%  %.0f", sup.omegaRCP * 100, rcsKgs),
            CGPoint(x: rcp.midX, y: rcp.maxY + fy(0.022)), .center, Theme.textDim, 8)
        txt(ctx, "kg/s", CGPoint(x: rcp.midX, y: rcp.maxY + fy(0.042)), .center, Theme.textDim, 7)

        // Hot / cold leg temperature tags on clear pipe runs.
        reading(ctx, CGPoint(x: vessel.maxX + fx(0.030), y: hotY - fy(0.022)),
                String(format: "T-HOT  %.0f K", snap.hotLegTempK), cHot, 10, .leading)
        reading(ctx, CGPoint(x: vessel.maxX + fx(0.006), y: coldY + fy(0.026)),
                String(format: "T-COLD  %.0f K", snap.coldLegTempK), cCold, 9, .leading)
        txt(ctx, String(format: "ΔT %.0f K", deltaT),
            CGPoint(x: fx(0.232) + fx(0.006), y: (hotY + sgHotIn.y)/2), .leading, Theme.textDim, 9)
        // T-avg deviation from the 550 K program (±0.5 K controller deadband).
        reading(ctx, CGPoint(x: vessel.maxX + fx(0.006), y: coldY + fy(0.046)),
                String(format: "Δ %+.1f K", snap.coolantTempK - 550),
                Theme.setpointDev(snap.coolantTempK, 550, 0.5), 8, .leading)

        // ════════════════════════════════════════════════════════════════════
        // STEAM GENERATOR
        // ════════════════════════════════════════════════════════════════════
        steamGen(ctx, sg, level: sgLvl, alarm: sup.steamInv < 0.25 ? Theme.alarm : nil)
        reading(ctx, CGPoint(x: sg.midX, y: sg.maxY + fy(0.038)),
                String(format: "%.2f MPa", snap.steamPressureMPa), Theme.ink, 11, .center)
        txt(ctx, String(format: "%.0f K   LVL %.0f%%", snap.sgTempK, sgLvl * 100),
            CGPoint(x: sg.midX, y: sg.maxY + fy(0.065)), .center, Theme.textDim, 9)
        // steam / feed flow tags — placed to the RIGHT of the risers so neither
        // pipe runs through its label.
        reading(ctx, CGPoint(x: steamRiserX + fx(0.010), y: sgSteam.y - fy(0.030)),
                String(format: "STM %.0f kg/s", steamKgs), Theme.ink, 9, .leading)
        reading(ctx, CGPoint(x: steamRiserX + fx(0.010), y: sgFeed.y - fy(0.004)),
                String(format: "FW %.0f kg/s", fwKgs),
                sup.feedwaterFault ? Theme.alarm : Theme.textDim, 9, .leading)

        // ════════════════════════════════════════════════════════════════════
        // TURBINE HALL
        // ════════════════════════════════════════════════════════════════════
        let trip = sup.turbineTrip
        turbine(ctx, hpT, tripped: trip, stages: 4)
        equipBox(ctx, msr, "", nil)
        txt(ctx, "MSR", CGPoint(x: msr.midX, y: msr.midY), .center, Theme.textHdr, 8)
        turbine(ctx, lpT, tripped: trip, stages: 6)
        txt(ctx, "HP", CGPoint(x: hpT.midX, y: hpT.minY - fy(0.024)), .center, Theme.textHdr, 9)
        txt(ctx, "LP TURBINE", CGPoint(x: lpT.midX, y: lpT.minY - fy(0.010)), .center, Theme.textHdr, 9)
        reading(ctx, CGPoint(x: lpT.midX, y: lpT.maxY + fy(0.030)),
                trip ? "TRIPPED" : "\(rpm) rpm", trip ? Theme.alarm : Theme.textDim, 9, .center)

        // Generator (driven off the LP shaft).
        equipBox(ctx, gen, "", trip ? Theme.alarm : nil)
        let gc = CGPoint(x: gen.midX, y: gen.midY)
        ctx.stroke(Path(ellipseIn: CGRect(x: gc.x - fx(0.016), y: gc.y - fx(0.016),
                                          width: fx(0.032), height: fx(0.032))),
                   with: .color(Theme.ink.opacity(0.55)), lineWidth: 1)
        txt(ctx, "G", gc, .center, Theme.ink, 12)
        txt(ctx, "GEN", CGPoint(x: gen.midX, y: gen.minY - fy(0.024)), .center, Theme.textHdr, 9)
        reading(ctx, CGPoint(x: gen.midX, y: gen.maxY + fy(0.032)),
                trip ? "0 MWe" : String(format: "%.0f MWe", snap.electricPowerW / 1e6),
                trip ? Theme.alarm : Theme.ink, 13, .center)
        // 400 kV switchyard one-line rising into the top-right corner.
        drawSwitchyard(ctx, w: W, h: H, gen: gen, snap: snap, sup: sup)

        // ════════════════════════════════════════════════════════════════════
        // CONDENSER + FEED TRAIN
        // ════════════════════════════════════════════════════════════════════
        condenser(ctx, cond, vacOK: condKPa < 12)
        reading(ctx, CGPoint(x: cond.midX, y: cond.minY + fy(0.050)),
                String(format: "%.1f kPa", condKPa), condKPa < 12 ? Theme.ink : Theme.caution, 11, .center)
        txt(ctx, String(format: "%.0f K", sup.condTempK),
            CGPoint(x: cond.midX, y: cond.minY + fy(0.075)), .center, Theme.textDim, 9)

        // Feed-train components.
        pump(ctx, cp, r: fy(0.022), running: !sup.feedwaterFault && sup.feedwaterValve > 0.02, tint: Theme.accent, mirrored: true)
        txt(ctx, "CP", CGPoint(x: cp.x, y: cp.y + fy(0.040)), .center, Theme.textDim, 8)
        heaterBox(ctx, lpFwh, "LP HTR")
        heaterBox(ctx, dae, "DAERATOR")
        pump(ctx, fp, r: fy(0.024), running: !sup.feedwaterFault && sup.feedwaterValve > 0.05, tint: Theme.accent, mirrored: true)
        txt(ctx, "FW PUMP", CGPoint(x: fp.x, y: fp.y + fy(0.042)), .center,
            sup.feedwaterFault ? Theme.alarm : Theme.textDim, 8)
        heaterBox(ctx, hpFwh, "HP HTR")
        reading(ctx, CGPoint(x: fp.x, y: fp.y - fy(0.040)),
                String(format: "FW %.0f%%", sup.feedwaterValve * 100),
                sup.feedwaterInv < 0.1 ? Theme.alarm : Theme.textDim, 9, .center)

        // ════════════════════════════════════════════════════════════════════
        // INSTRUMENT DOCKS + MARGINS STRIP — fill the dead zones with the data
        // operators actually scan. Every value derives from the live snapshot.
        // ════════════════════════════════════════════════════════════════════
        drawMarginsStrip(ctx, CGRect(x: fx(0.420), y: fy(0.034), width: fx(0.395), height: fy(0.060)), snap: snap, sup: sup)
        drawNeutronicsDock(ctx, CGRect(x: fx(0.012), y: fy(0.045), width: fx(0.140), height: fy(0.205)), snap: snap, sup: sup)
        // Electrical picture tucked 4 px under the GSU transformer circles (the
        // lower circle bottom = gen.midY + xfR), beside the generator. Computed in
        // px, not a fixed fraction, so the 4 px gap holds at any window size.
        let elecTop  = gen.midY + fx(0.013) + 4
        let elecRect = CGRect(x: fx(0.815), y: elecTop, width: fx(0.170), height: fy(0.250))
        drawElectricalDock(ctx, elecRect, snap: snap, sup: sup, t: t)
        // Under it: turbine-generator mechanical + excitation, filling down to the data page.
        let tgTop = elecRect.maxY + fy(0.012)
        drawTurbineGen(ctx, CGRect(x: fx(0.815), y: tgTop, width: fx(0.170), height: fy(0.824) - tgTop), snap: snap, sup: sup)
        // Fill the empty bottom-left (primary loop) and the centre void (steam cycle).
        drawPrimaryDock(ctx, CGRect(x: fx(0.150), y: fy(0.815), width: fx(0.175), height: fy(0.165)), snap: snap, sup: sup)
        drawSteamCycle(ctx, CGRect(x: fx(0.442), y: fy(0.430), width: fx(0.150), height: fy(0.150)), snap: snap, sup: sup)
        // Fill the lower-right with a dense engineering data page (raw numbers).
        // Sits low + clear of the feed-train labels above it.
        drawCoreLimits(ctx, CGRect(x: fx(0.420), y: fy(0.838), width: fx(0.565), height: fy(0.152)), snap: snap, sup: sup)
    }

    // MARK: — instrument docks -------------------------------------------------

    /// Flat engineering bezel (NOT a glass card): faint structural tint, square
    /// corners, hairline border, title top-left.
    private func dockField(_ ctx: GraphicsContext, _ r: CGRect, _ title: String) {
        ctx.fill(Path(r), with: .color(Theme.dockTint))
        ctx.stroke(Path(r), with: .color(Theme.ink.opacity(0.18)), lineWidth: 1)
        ctx.draw(Text(title).font(.system(size: 9, weight: .semibold, design: .monospaced))
                    .foregroundColor(Theme.textHdr),
                 at: CGPoint(x: r.minX + 8, y: r.minY + 5), anchor: .topLeading)
    }

    private func tagRow(_ ctx: GraphicsContext, _ r: CGRect, _ y: CGFloat,
                        _ label: String, _ value: String, _ color: Color) {
        ctx.draw(Text(label).font(.system(size: 8, design: .monospaced)).foregroundColor(Theme.textDim),
                 at: CGPoint(x: r.minX + 8, y: y), anchor: .leading)
        ctx.draw(Text(value).font(.system(size: 10, weight: .semibold, design: .monospaced)).foregroundColor(color),
                 at: CGPoint(x: r.maxX - 8, y: y), anchor: .trailing)
    }

    /// Zero-centered signed bar (e.g. net reactivity, deviation).
    private func signedBar(_ ctx: GraphicsContext, _ r: CGRect, value: Double, span: Double, _ color: Color) {
        ctx.fill(Path(roundedRect: r, cornerRadius: 2, style: .continuous), with: .color(Theme.ink.opacity(0.08)))
        let cx = r.midX
        let f = max(-1, min(1, value / max(1e-9, span)))
        let w = (r.width / 2 - 1) * CGFloat(abs(f))
        let fr = f >= 0 ? CGRect(x: cx, y: r.minY + 1, width: w, height: r.height - 2)
                        : CGRect(x: cx - w, y: r.minY + 1, width: w, height: r.height - 2)
        ctx.fill(Path(fr), with: .color(color))
        ctx.stroke(Path { p in p.move(to: CGPoint(x: cx, y: r.minY)); p.addLine(to: CGPoint(x: cx, y: r.maxY)) },
                   with: .color(Theme.ink.opacity(0.5)), lineWidth: 0.5)
    }

    /// Rail/standing cue: green when comfortable, warming toward the trip.
    private func marginColor(_ pct: Double) -> Color {
        if pct <= 0  { return Theme.alarm }
        if pct < 15  { return Theme.alarm }
        if pct < 30  { return Theme.warning }
        if pct < 45  { return Theme.caution }
        return Theme.isFlat ? Theme.ink.opacity(0.25) : Theme.statusNormal   // comfortable: neutral on ISA skins
    }

    /// Five whole-plant margin-to-trip tiles where the eye lands. The big value
    /// stays CALM ink while comfortable and only warms as the parameter nears its
    /// real trip; a thin left rail carries the standing OK/approach cue, so a
    /// healthy plant is quiet (no wall of green) even on the ISA skins.
    private func drawMarginsStrip(_ ctx: GraphicsContext, _ r: CGRect, snap: PlantSnapshot, sup: PlantSupervisor) {
        let tiles: [(String, Double)] = [
            ("FLUX",   (1.20 - snap.powerFraction) / 0.20),
            ("FUEL-T", (1500 - snap.fuelTempK) / 600),
            ("RCS-P",  (17.0 - sup.pressureMPa) / 1.5),
            ("COOL-T", (620 - snap.coolantTempK) / 70),
            ("SG-LVL", (sup.steamInv - 0.25) / 0.75),
        ]
        let tw = r.width / CGFloat(tiles.count)
        ctx.draw(Text("MARGIN TO TRIP").font(.system(size: 7, design: .monospaced)).foregroundColor(Theme.textDim),
                 at: CGPoint(x: r.minX, y: r.minY - 4), anchor: .bottomLeading)
        for (i, t) in tiles.enumerated() {
            let pct = max(0, min(100, t.1 * 100))
            let cell = CGRect(x: r.minX + tw * CGFloat(i) + 2, y: r.minY, width: tw - 4, height: r.height)
            ctx.fill(Path(roundedRect: cell, cornerRadius: 2, style: .continuous), with: .color(Theme.ink.opacity(0.06)))
            let rail = marginColor(pct)
            ctx.fill(Path(CGRect(x: cell.minX, y: cell.minY, width: 2, height: cell.height)), with: .color(rail))
            let valC = pct >= 45 ? Theme.ink : rail   // calm unless approaching trip
            ctx.draw(Text(t.0).font(.system(size: 7, design: .monospaced)).foregroundColor(Theme.textDim),
                     at: CGPoint(x: cell.midX, y: cell.minY + 7), anchor: .center)
            ctx.draw(Text("\(Int(pct))%").font(.system(size: 12, weight: .bold, design: .monospaced)).foregroundColor(valC),
                     at: CGPoint(x: cell.midX, y: cell.maxY - 5), anchor: .bottom)
        }
    }

    /// Left dock: the neutronics an operator watches on the approach to power.
    private func drawNeutronicsDock(_ ctx: GraphicsContext, _ r: CGRect, snap: PlantSnapshot, sup: PlantSupervisor) {
        dockField(ctx, r, "NEUTRONICS")
        let pcm = snap.reactivity * 1e5
        let xeWorth = -1.6 * snap.xenonInventory            // model coeff 1.6e-5 → pcm
        let decayPct = snap.decayHeatFraction * 100
        // Period from the historian's recent power slope (qualitative — honest).
        let h = sup.orderedHistory(sup.histPower)
        let p1 = h.last ?? 100, p0 = h.count > 40 ? h[h.count - 40] : (h.first ?? 100)
        let period: (String, Color) = (p1 > 0.2 && p0 > 0.2)
            ? (p1 > p0 * 1.01 ? ("RISING", Theme.caution)
               : p1 < p0 * 0.99 ? ("FALLING", Theme.fluidSubcooled) : ("STABLE", Theme.statusNormal))
            : ("—", Theme.textDim)
        var y = r.minY + 24
        ctx.draw(Text("NET REACTIVITY").font(.system(size: 8, design: .monospaced)).foregroundColor(Theme.textDim),
                 at: CGPoint(x: r.minX + 8, y: y), anchor: .leading)
        ctx.draw(Text(fmtPcm(pcm)).font(.system(size: 9, weight: .semibold, design: .monospaced))
                    .foregroundColor(abs(pcm) > 500 ? Theme.alarm : abs(pcm) > 100 ? Theme.caution : Theme.statusNormal),
                 at: CGPoint(x: r.maxX - 8, y: y), anchor: .trailing)
        y += 10
        signedBar(ctx, CGRect(x: r.minX + 8, y: y, width: r.width - 16, height: 7),
                  value: pcm, span: 500, abs(pcm) > 100 ? Theme.caution : Theme.statusNormal)
        y += 20
        let nPitch = max(14, (r.height - 62) / 3)   // 4 rows / 3 gaps, derived from height
        tagRow(ctx, r, y, "PWR TREND", period.0, period.1); y += nPitch
        tagRow(ctx, r, y, "BORON", String(format: "%.0f ppm", sup.boronPPM), Theme.ink); y += nPitch
        tagRow(ctx, r, y, "XENON", fmtPcm(xeWorth), Theme.ink); y += nPitch
        tagRow(ctx, r, y, "DECAY HT", String(format: "%.1f %%", decayPct), Theme.ink)
    }

    /// Right dock: the electrical / grid picture, with the iconic synchroscope.
    private func drawElectricalDock(_ ctx: GraphicsContext, _ r: CGRect, snap: PlantSnapshot, sup: PlantSupervisor, t: Double) {
        dockField(ctx, r, "ELECTRICAL")
        let gross = snap.electricPowerW / 1e6
        let aux = gross * 0.05
        let synced = !sup.turbineTrip && !snap.scrammed && gross > 1
        let gcol = synced ? Theme.elecGold : Theme.deEnergized
        // Values band (top) — kept clear of the synchroscope band below.
        var y = r.minY + 20
        tagRow(ctx, r, y, "GROSS", String(format: "%.0f MWe", gross), gcol); y += 13
        tagRow(ctx, r, y, "NET", String(format: "%.0f MWe", gross - aux), gcol); y += 13
        tagRow(ctx, r, y, "MVAr*", String(format: "%.0f", gross * 0.43), gcol); y += 13
        let hzCol = synced ? (Theme.isFlat ? Theme.ink : Theme.statusNormal) : Theme.deEnergized
        tagRow(ctx, r, y, "GRID Hz", synced ? "50.00" : "--.--", hzCol); y += 11
        ctx.draw(Text("*0.92 pf").font(.system(size: 6, design: .monospaced)).foregroundColor(Theme.textDim),
                 at: CGPoint(x: r.minX + 8, y: y), anchor: .leading)
        // Synchroscope sized/placed within the band BELOW the values, so it never
        // climbs into the value rows even as the window shrinks.
        let valuesBottom = r.minY + 78
        let dr = max(12, min(r.width * 0.18, (r.maxY - valuesBottom - 16) / 2))
        synchroscope(ctx, CGPoint(x: r.midX, y: valuesBottom + dr + 6), dr, synced: synced, t: t)
    }

    /// Generator-paralleling synchroscope: pointer dead-still at 12 o'clock with a
    /// green lock pip when synchronized; rotates at slip rate otherwise.
    private func synchroscope(_ ctx: GraphicsContext, _ c: CGPoint, _ rad: CGFloat, synced: Bool, t: Double) {
        ctx.stroke(Path(ellipseIn: CGRect(x: c.x - rad, y: c.y - rad, width: rad * 2, height: rad * 2)),
                   with: .color(Theme.ink.opacity(0.3)), lineWidth: 1)
        for i in 0..<12 {
            let a = Double(i) / 12 * 2 * .pi
            let p0 = CGPoint(x: c.x + cos(a) * rad * 0.82, y: c.y + sin(a) * rad * 0.82)
            let p1 = CGPoint(x: c.x + cos(a) * rad * 0.95, y: c.y + sin(a) * rad * 0.95)
            ctx.stroke(Path { p in p.move(to: p0); p.addLine(to: p1) }, with: .color(Theme.ink.opacity(0.2)), lineWidth: 0.5)
        }
        let ang = synced ? -Double.pi / 2 : (-Double.pi / 2 + t * 0.4 * 2 * .pi).truncatingRemainder(dividingBy: 2 * .pi)
        let tip = CGPoint(x: c.x + cos(ang) * rad * 0.7, y: c.y + sin(ang) * rad * 0.7)
        ctx.stroke(Path { p in p.move(to: c); p.addLine(to: tip) },
                   with: .color(synced ? (Theme.isFlat ? Theme.ink : Theme.statusNormal) : Theme.caution), lineWidth: 1.5)
        if synced {
            ctx.fill(Path(ellipseIn: CGRect(x: c.x - 2, y: c.y - rad - 1, width: 4, height: 4)),
                     with: .color(Theme.isFlat ? Theme.textHdr : Theme.statusNormal))
        }
        ctx.draw(Text("SYNC").font(.system(size: 6, design: .monospaced)).foregroundColor(Theme.textDim),
                 at: CGPoint(x: c.x, y: c.y + rad + 6), anchor: .center)
    }

    /// Saturation temperature [K] at pressure [MPa] — inverse of satPressureMPa.
    private func tsatK(_ pMPa: Double) -> Double {
        let pbar = max(0.01, pMPa * 10)
        return 1951.8 / (5.444 - log10(pbar)) + 16.5
    }

    /// Bottom-left dock: the primary-loop picture that complements the neutronics.
    private func drawPrimaryDock(_ ctx: GraphicsContext, _ r: CGRect, snap: PlantSnapshot, sup: PlantSupervisor) {
        dockField(ctx, r, "RCS PRIMARY")
        let flowF = snap.scrammed ? 0 : max(0, Double(sup.primaryFlow) * sup.omegaRCP)
        let subcool = tsatK(sup.pressureMPa) - snap.hotLegTempK     // margin to saturation
        let pitch = max(14, (r.height - 32) / 4)                    // 5 rows / 4 gaps, derived from height
        var y = r.minY + 24
        tagRow(ctx, r, y, "RCS FLOW", String(format: "%.0f kg/s", 18_800 * flowF), Theme.ink); y += pitch
        tagRow(ctx, r, y, "RCP-1 ΔP", String(format: "%.2f MPa", 0.62 * sup.omegaRCP * sup.omegaRCP), Theme.ink); y += pitch
        tagRow(ctx, r, y, "RCP-1 AMPS", String(format: "%.0f A", 6_000 * sup.omegaRCP), Theme.ink); y += pitch
        tagRow(ctx, r, y, "SUBCOOL", String(format: "%.0f K", subcool),
               subcool < 15 ? Theme.alarm : subcool < 30 ? Theme.caution : Theme.statusNormal); y += pitch
        tagRow(ctx, r, y, "RC PUMPS", sup.pumpDegraded ? "DEGRADED" : "1 / 1 RUN",
               sup.pumpDegraded ? Theme.alarm : Theme.textDim)
    }

    /// Centre dock: the secondary thermodynamic picture, sited in the steam path.
    private func drawSteamCycle(_ ctx: GraphicsContext, _ r: CGRect, snap: PlantSnapshot, sup: PlantSupervisor) {
        dockField(ctx, r, "STEAM CYCLE")
        let eta = snap.thermalPowerW > 1 ? snap.electricPowerW / snap.thermalPowerW * 100 : 0
        let steamKgs = snap.scrammed ? 0 : 1_650.0 * max(0, snap.powerFraction)
        let pitch = max(14, (r.height - 32) / 4)                    // 5 rows / 4 gaps, derived from height
        var y = r.minY + 24
        tagRow(ctx, r, y, "HDR PRESS", String(format: "%.2f MPa", snap.steamPressureMPa), Theme.ink); y += pitch
        tagRow(ctx, r, y, "SAT TEMP", String(format: "%.0f K", snap.sgTempK), Theme.ink); y += pitch
        tagRow(ctx, r, y, "STM FLOW", String(format: "%.0f kg/s", steamKgs), Theme.fluidTwoPhase); y += pitch
        tagRow(ctx, r, y, "CYCLE η", String(format: "%.1f %%", eta), Theme.elecGold); y += pitch
        tagRow(ctx, r, y, "TBN", sup.turbineTrip ? "TRIPPED" : "3000 rpm",
               sup.turbineTrip ? Theme.alarm : Theme.textDim)
    }

    /// Lower-right dock: an engineering DATA PAGE of the advanced reactor / thermal
    /// numbers — fewer, BIG, readable. Derived values directionally correct + labeled.
    private func drawCoreLimits(_ ctx: GraphicsContext, _ r: CGRect, snap: PlantSnapshot, sup: PlantSupervisor) {
        dockField(ctx, r, "CORE · THERMAL DATA")
        let pf   = max(0.001, snap.powerFraction)
        let flow = max(0.05, Double(sup.primaryFlow) * sup.omegaRCP)
        let pRat = max(0.3, sup.pressureMPa / 15.5)
        let tFac = max(0.5, 1 - (snap.coolantTempK - 550) / 130)
        let dnbr = min(9.99, (flow.squareRoot() / pf.squareRoot()) * pRat * tFac * 1.93)
        let ao   = -snap.rodPosition * 32 + (1 - flow) * 6
        let lhr  = 13.0 * pf * (1 + abs(ao) / 110)
        let qptr = 1.00 + abs(ao) / 1800
        let pct  = snap.coolantTempK + 80 * pf + max(0, snap.fuelTempK - 900) * 0.4
        let subcool = tsatK(sup.pressureMPa) - snap.hotLegTempK
        let pcm  = snap.reactivity * 1e5
        let sdm  = max(0, 5.2 - snap.rodPosition * 1.5)
        let sgL  = max(8, min(92, sup.feedwaterInv * 62))
        let mwe  = snap.electricPowerW / 1e6
        let eta  = snap.thermalPowerW > 1 ? mwe / (snap.thermalPowerW / 1e6) * 100 : 0
        let hp = sup.orderedHistory(sup.histPower)
        let p1 = hp.last ?? 100, p0 = hp.count > 50 ? hp[hp.count - 50] : p1
        let rate = (p1 > 0.1 && p0 > 0.1) ? (p1 - p0) / max(0.1, p0) : 0
        let period = abs(rate) < 0.0015 ? 999.0 : min(999, 0.8 / abs(rate))
        let cA = { (c: Bool) in c ? Theme.alarm : Theme.ink }
        let dnbrC = dnbr < 1.30 ? Theme.alarm : dnbr < 1.55 ? Theme.caution : Theme.statusNormal

        let items: [(String, String, Color)] = [
            ("DNBR",     String(format: "%.2f", dnbr), dnbrC),
            ("PK LHR",   String(format: "%.1f", lhr) + " kW/ft", cA(lhr > 20)),
            ("AXIAL ΔI", String(format: "%+.0f%%", ao), abs(ao) > 15 ? Theme.caution : Theme.ink),
            ("PK CLAD",  String(format: "%.0f K", pct), cA(pct > 1200)),
            ("QPTR",     String(format: "%.3f", qptr), cA(qptr > 1.02)),
            ("SUR",      String(format: "%+.1f dpm", rate * 60), Theme.ink),
            ("PERIOD",   period >= 999 ? "∞ s" : String(format: "%.0f s", period), Theme.ink),
            ("ρ NET",    fmtPcm(pcm), abs(pcm) > 100 ? Theme.caution : Theme.ink),
            ("SDM",      String(format: "%.1f%%Δk", sdm), sdm < 1.3 ? Theme.alarm : Theme.ink),
            ("BORON",    String(format: "%.0f ppm", sup.boronPPM), Theme.ink),
            ("T-HOT",    String(format: "%.0f K", snap.hotLegTempK), Theme.ink),
            ("T-COLD",   String(format: "%.0f K", snap.coldLegTempK), Theme.ink),
            ("ΔT",       String(format: "%.0f K", snap.hotLegTempK - snap.coldLegTempK), Theme.ink),
            ("SUBCOOL",  String(format: "%.0f K", subcool), subcool < 15 ? Theme.alarm : subcool < 30 ? Theme.caution : Theme.ink),
            ("RCS-P",    String(format: "%.1f MPa", sup.pressureMPa), Theme.ink),
            ("RCS FLOW", String(format: "%.0f", 18_800 * flow), Theme.ink),
            ("SG-P",     String(format: "%.1f MPa", snap.steamPressureMPa), Theme.ink),
            ("SG LVL",   String(format: "%.0f%%", sgL), sgL < 20 ? Theme.alarm : Theme.ink),
            ("GROSS",    String(format: "%.0f MWe", mwe), Theme.elecGold),
            ("CYCLE η",  String(format: "%.1f%%", eta), Theme.ink),
        ]

        let cols = 4
        let rpc = (items.count + cols - 1) / cols       // 5 rows
        let gx: CGFloat = 10
        let colW = (r.width - gx * 2) / CGFloat(cols)
        let top = r.minY + 18
        let rowH = (r.maxY - top - 5) / CGFloat(rpc)
        for (i, it) in items.enumerated() {
            let col = i / rpc, row = i % rpc
            let cx = r.minX + gx + CGFloat(col) * colW
            let cy = top + CGFloat(row) * rowH + rowH / 2
            ctx.draw(Text(it.0).font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                     at: CGPoint(x: cx, y: cy), anchor: .leading)
            ctx.draw(Text(it.1).font(.system(size: 11, weight: .semibold, design: .monospaced)).foregroundColor(it.2),
                     at: CGPoint(x: cx + colW - 24, y: cy), anchor: .trailing)
        }
    }

    // MARK: — derived helpers ---------------------------------------------------

    /// Saturation pressure (kPa) via the Magnus formula — condenser backpressure.
    /// (Magnus already yields kPa: ~6.6 kPa at 311 K — healthy turbine vacuum.)
    private func magnusKPa(_ tK: Double) -> Double {
        let tC = tK - 273.15
        return 0.61094 * exp(17.625 * tC / (tC + 243.04))
    }

    /// pcm string without the ugly "−0" that "%+.0f" prints near zero.
    private func fmtPcm(_ v: Double) -> String { abs(v) < 0.5 ? "0 pcm" : String(format: "%+.0f pcm", v) }

    // MARK: — primitives --------------------------------------------------------

    /// Bevel every interior corner of a polyline into a 45° miter (PCB-style /
    /// schematic convention) — auto-shrinking the bevel on short segments so it
    /// never overshoots. Endpoints are preserved, so component connections hold.
    private func chamfered(_ pts: [(CGFloat, CGFloat)], _ bevel: CGFloat = 10) -> [CGPoint] {
        let p = pts.map { CGPoint(x: $0.0, y: $0.1) }
        guard p.count >= 3 else { return p }
        var out: [CGPoint] = [p[0]]
        for i in 1..<(p.count - 1) {
            let a = p[i - 1], b = p[i], c = p[i + 1]
            let dIn = hypot(b.x - a.x, b.y - a.y), dOut = hypot(c.x - b.x, c.y - b.y)
            guard dIn > 1, dOut > 1 else { out.append(b); continue }
            let ci = min(bevel, dIn / 2), co = min(bevel, dOut / 2)
            out.append(CGPoint(x: b.x - (b.x - a.x) / dIn * ci, y: b.y - (b.y - a.y) / dIn * ci))
            out.append(CGPoint(x: b.x + (c.x - b.x) / dOut * co, y: b.y + (c.y - b.y) / dOut * co))
        }
        out.append(p[p.count - 1])
        return out
    }

    private func polyPath(_ pts: [CGPoint]) -> Path {
        var path = Path()
        guard let first = pts.first else { return path }
        path.move(to: first)
        for q in pts.dropFirst() { path.addLine(to: q) }
        return path
    }

    private func pipe(_ ctx: GraphicsContext, _ pts: [(CGFloat, CGFloat)], _ color: Color, _ w: CGFloat) {
        guard pts.count >= 2 else { return }
        ctx.stroke(polyPath(chamfered(pts)), with: .color(color),
                   style: StrokeStyle(lineWidth: w, lineJoin: .round))
    }

    private func flowDash(_ ctx: GraphicsContext, _ pts: [(CGFloat, CGFloat)], _ tint: Color, _ q: Double, _ t: Double) {
        guard pts.count >= 2, q > 0.01 else { return }
        let qq = max(0, min(1, q))
        ctx.stroke(polyPath(chamfered(pts)), with: .color(tint.opacity(0.95)),
                   style: StrokeStyle(lineWidth: 1 + qq, lineCap: .butt,
                                      dash: [3 + 1.5 * qq, 16 - 9 * qq],
                                      dashPhase: CGFloat(-(t * (26 + 150 * qq)))))
    }

    private func equipBox(_ ctx: GraphicsContext, _ r: CGRect, _ label: String, _ alarm: Color?) {
        let s = Path(roundedRect: r, cornerRadius: 4, style: .continuous)
        ctx.fill(s, with: .color(Theme.equipFill))
        ctx.stroke(s, with: .color(alarm ?? Theme.ink.opacity(0.30)), lineWidth: alarm == nil ? 1 : 2)
        if !label.isEmpty {
            ctx.draw(Text(label).font(.system(size: 10, design: .monospaced)).foregroundColor(alarm ?? Theme.textHdr),
                     at: CGPoint(x: r.midX, y: r.minY + 9), anchor: .top)
        }
    }

    /// Reactor pressure vessel: domed heads, CRDM stubs, core lattice, inserted
    /// control rods, downcomer annulus, power-dependent core glow.
    private func vesselRPV(_ ctx: GraphicsContext, _ body: CGRect, snap: PlantSnapshot, alarm: Color?) {
        let domeH = body.width * 0.42
        var v = Path()
        v.move(to: CGPoint(x: body.minX, y: body.minY))
        v.addQuadCurve(to: CGPoint(x: body.maxX, y: body.minY),
                       control: CGPoint(x: body.midX, y: body.minY - domeH))
        v.addLine(to: CGPoint(x: body.maxX, y: body.maxY))
        v.addQuadCurve(to: CGPoint(x: body.minX, y: body.maxY),
                       control: CGPoint(x: body.midX, y: body.maxY + domeH * 0.85))
        v.closeSubpath()
        ctx.fill(v, with: .color(Theme.equipFill))
        ctx.stroke(v, with: .color(alarm ?? Theme.ink.opacity(0.32)), lineWidth: alarm == nil ? 1 : 2)

        // CRDM stubs above the head — short, uniform drive housings.
        let stubBase = body.minY - domeH * 0.38
        for i in 0..<5 {
            let x = body.minX + body.width * (0.20 + 0.15 * CGFloat(i))
            ctx.stroke(Path { p in
                p.move(to: CGPoint(x: x, y: stubBase)); p.addLine(to: CGPoint(x: x, y: stubBase - body.width * 0.20))
            }, with: .color(Theme.ink.opacity(0.32)), lineWidth: 1.5)
        }

        // Downcomer annulus.
        for dx in [body.minX + 2.5, body.maxX - 2.5] {
            ctx.stroke(Path { p in
                p.move(to: CGPoint(x: dx, y: body.minY + 4)); p.addLine(to: CGPoint(x: dx, y: body.maxY - 4))
            }, with: .color(Theme.ink.opacity(0.14)), lineWidth: 1)
        }

        // Core region with glow + assembly lattice.
        let core = CGRect(x: body.minX + body.width * 0.16, y: body.minY + body.height * 0.42,
                          width: body.width * 0.68, height: body.height * 0.36)
        // Core incandescence by power: radial flux peak, subcritical-dark → amber
        // at full power → white-hot, with an overpower red rim before the trip.
        let g = max(0, min(1, snap.powerFraction / 1.10))
        let glow = Theme.fluxGlow(g)
        ctx.fill(Path(core), with: .radialGradient(
            Gradient(colors: [glow.center, glow.edge]),
            center: CGPoint(x: core.midX, y: core.midY),
            startRadius: 0, endRadius: core.width * 0.7))
        if snap.powerFraction > 1.10 {
            ctx.stroke(Path(core), with: .color(Theme.alarm.opacity(min(1, (snap.powerFraction - 1.10) / 0.10))), lineWidth: 2)
        } else {
            ctx.stroke(Path(core), with: .color(Theme.ink.opacity(0.20)), lineWidth: 1)
        }
        for gx in 1..<4 {
            let x = core.minX + core.width * CGFloat(gx) / 4
            ctx.stroke(Path { p in p.move(to: CGPoint(x: x, y: core.minY)); p.addLine(to: CGPoint(x: x, y: core.maxY)) },
                       with: .color(Theme.ink.opacity(0.12)), lineWidth: 0.75)
        }
        for gy in 1..<3 {
            let y = core.minY + core.height * CGFloat(gy) / 3
            ctx.stroke(Path { p in p.move(to: CGPoint(x: core.minX, y: y)); p.addLine(to: CGPoint(x: core.maxX, y: y)) },
                       with: .color(Theme.ink.opacity(0.12)), lineWidth: 0.75)
        }

        // Control rods: enter from the top head, inserted depth ∝ rodPosition.
        let rodPos = max(0, min(1, snap.rodPosition))
        let rodTop = body.minY + 4
        let rodSpan = core.maxY - rodTop
        for i in 0..<4 {
            let x = core.minX + core.width * (0.12 + 0.25 * CGFloat(i))
            let depth = rodSpan * CGFloat(rodPos)
            if depth > 1 {
                ctx.stroke(Path { p in
                    p.move(to: CGPoint(x: x, y: rodTop)); p.addLine(to: CGPoint(x: x, y: rodTop + depth))
                }, with: .color(Theme.ink.opacity(0.55)), lineWidth: 2)
            }
        }

        ctx.draw(Text("RPV").font(.system(size: 10, weight: .semibold, design: .monospaced))
                    .foregroundColor(alarm ?? Theme.textHdr),
                 at: CGPoint(x: body.midX, y: body.minY + 12), anchor: .top)
    }

    private func nozzle(_ ctx: GraphicsContext, _ c: CGPoint, _ color: Color) {
        let r = CGRect(x: c.x - 1, y: c.y - 3, width: 6, height: 6)
        ctx.fill(Path(r), with: .color(color.opacity(0.8)))
    }

    /// Pressurizer: vertical cylinder, level fill, heater coil hint at the base.
    private func pressurizer(_ ctx: GraphicsContext, _ r: CGRect, level: Double, heatersOn: Bool, alarm: Color?) {
        let domeH = r.width * 0.4
        var s = Path()
        s.move(to: CGPoint(x: r.minX, y: r.minY))
        s.addQuadCurve(to: CGPoint(x: r.maxX, y: r.minY), control: CGPoint(x: r.midX, y: r.minY - domeH))
        s.addLine(to: CGPoint(x: r.maxX, y: r.maxY)); s.addLine(to: CGPoint(x: r.minX, y: r.maxY)); s.closeSubpath()
        ctx.fill(s, with: .color(Theme.equipFill))
        let lv = max(0, min(1, level))
        let liq = CGRect(x: r.minX + 1.5, y: r.maxY - r.height * CGFloat(lv) + 1,
                         width: r.width - 3, height: r.height * CGFloat(lv) - 2)
        ctx.fill(Path(roundedRect: liq, cornerRadius: 2, style: .continuous), with: .color(Theme.water.opacity(0.32)))
        ctx.stroke(Path(CGRect(x: liq.minX, y: liq.minY, width: liq.width, height: 1)),
                   with: .color(Theme.water.opacity(0.7)), lineWidth: 1)
        ctx.stroke(s, with: .color(alarm ?? Theme.ink.opacity(0.32)), lineWidth: alarm == nil ? 1 : 2)
        // heater hint
        ctx.fill(Path(CGRect(x: r.minX + 3, y: r.maxY - 4, width: r.width - 6, height: 2)),
                 with: .color(heatersOn ? Theme.caution : Theme.ink.opacity(0.2)))
        ctx.draw(Text("PZR").font(.system(size: 9, weight: .semibold, design: .monospaced)).foregroundColor(alarm ?? Theme.textHdr),
                 at: CGPoint(x: r.midX, y: r.minY + 8), anchor: .top)
    }

    /// Steam generator: fat evaporator shell tapering to a slimmer steam drum,
    /// with a secondary-side level fill.
    private func steamGen(_ ctx: GraphicsContext, _ r: CGRect, level: Double, alarm: Color?) {
        let neck = r.width * 0.10            // gentle shoulder, not a bullet taper
        let taperY = r.minY + r.height * 0.20
        let domeH = r.width * 0.22
        var s = Path()
        s.move(to: CGPoint(x: r.minX, y: r.maxY))
        s.addLine(to: CGPoint(x: r.minX, y: taperY))
        s.addLine(to: CGPoint(x: r.minX + neck, y: r.minY + domeH))
        s.addQuadCurve(to: CGPoint(x: r.maxX - neck, y: r.minY + domeH),
                       control: CGPoint(x: r.midX, y: r.minY - domeH * 0.15))
        s.addLine(to: CGPoint(x: r.maxX, y: taperY))
        s.addLine(to: CGPoint(x: r.maxX, y: r.maxY))
        s.closeSubpath()
        ctx.fill(s, with: .color(Theme.equipFill))

        // Secondary level fill in the lower shell.
        let lv = max(0, min(1, level))
        let h = (r.maxY - taperY) * CGFloat(lv)
        let liq = CGRect(x: r.minX + 1.5, y: r.maxY - h + 1, width: r.width - 3, height: h - 2)
        if liq.height > 0 {
            ctx.fill(Path(liq), with: .color(Theme.water.opacity(0.30)))
            ctx.stroke(Path(CGRect(x: liq.minX, y: liq.minY, width: liq.width, height: 1)),
                       with: .color(Theme.water.opacity(0.7)), lineWidth: 1)
        }
        // U-tube bundle hint in the evaporator section.
        let bx = r.midX
        for off in [r.width * 0.14, r.width * 0.27] {
            ctx.stroke(Path { p in
                p.move(to: CGPoint(x: bx - off, y: r.maxY - r.height * 0.12))
                p.addLine(to: CGPoint(x: bx - off, y: r.minY + r.height * 0.34))
                p.addQuadCurve(to: CGPoint(x: bx + off, y: r.minY + r.height * 0.34),
                               control: CGPoint(x: bx, y: r.minY + r.height * 0.26))
                p.addLine(to: CGPoint(x: bx + off, y: r.maxY - r.height * 0.12))
            }, with: .color(Theme.ink.opacity(0.13)), lineWidth: 1)
        }

        ctx.stroke(s, with: .color(alarm ?? Theme.ink.opacity(0.32)), lineWidth: alarm == nil ? 1 : 2)
        ctx.draw(Text("S/G").font(.system(size: 11, weight: .semibold, design: .monospaced)).foregroundColor(alarm ?? Theme.textHdr),
                 at: CGPoint(x: r.midX, y: r.minY + r.height * 0.08), anchor: .top)
    }

    private func pump(_ ctx: GraphicsContext, _ c: CGPoint, r: CGFloat, running: Bool, tint: Color, mirrored: Bool = false) {
        let rect = CGRect(x: c.x - r, y: c.y - r, width: r*2, height: r*2)
        ctx.fill(Path(ellipseIn: rect), with: .color(Theme.equipFill))
        ctx.stroke(Path(ellipseIn: rect), with: .color(Theme.ink.opacity(0.32)), lineWidth: 1)
        // Impeller arrow points in the flow direction (mirror = right→left).
        let d: CGFloat = mirrored ? -1 : 1
        var tri = Path()
        tri.move(to: CGPoint(x: c.x - d*r*0.4, y: c.y - r*0.5))
        tri.addLine(to: CGPoint(x: c.x - d*r*0.4, y: c.y + r*0.5))
        tri.addLine(to: CGPoint(x: c.x + d*r*0.55, y: c.y))
        tri.closeSubpath()
        ctx.fill(tri, with: .color(running ? tint : Theme.deEnergized))
    }

    private func valve(_ ctx: GraphicsContext, _ c: CGPoint, open: Bool) {
        let s: CGFloat = 6
        var p = Path()
        p.move(to: CGPoint(x: c.x - s, y: c.y - s)); p.addLine(to: CGPoint(x: c.x + s, y: c.y + s))
        p.addLine(to: CGPoint(x: c.x + s, y: c.y - s)); p.addLine(to: CGPoint(x: c.x - s, y: c.y + s)); p.closeSubpath()
        ctx.fill(p, with: .color(open ? Theme.statusNormal.opacity(0.9) : Theme.deEnergized.opacity(0.8)))
        ctx.stroke(p, with: .color(Theme.ink.opacity(0.5)), lineWidth: 1)
    }

    /// Turbine casing: trapezoid widening to the exhaust end, with stage blades.
    private func turbine(_ ctx: GraphicsContext, _ r: CGRect, tripped: Bool, stages: Int) {
        var p = Path()
        let inset = r.height * 0.30
        p.move(to: CGPoint(x: r.minX, y: r.minY + inset))
        p.addLine(to: CGPoint(x: r.maxX, y: r.minY))
        p.addLine(to: CGPoint(x: r.maxX, y: r.maxY))
        p.addLine(to: CGPoint(x: r.minX, y: r.maxY - inset))
        p.closeSubpath()
        ctx.fill(p, with: .color(Theme.equipFill))
        ctx.stroke(p, with: .color(tripped ? Theme.alarm : Theme.ink.opacity(0.32)), lineWidth: tripped ? 2 : 1)
        for i in 1..<stages {
            let f = CGFloat(i) / CGFloat(stages)
            let x = r.minX + r.width * f
            let topInset = inset * (1 - f)
            ctx.stroke(Path { pa in
                pa.move(to: CGPoint(x: x, y: r.minY + topInset)); pa.addLine(to: CGPoint(x: x, y: r.maxY - topInset))
            }, with: .color(Theme.ink.opacity(0.18)), lineWidth: 1)
        }
    }

    private func transformer(_ ctx: GraphicsContext, _ c: CGPoint, r: CGFloat, live: Bool) {
        let col = (live ? Theme.elecGold : Theme.deEnergized).opacity(0.85)
        ctx.stroke(Path(ellipseIn: CGRect(x: c.x - r, y: c.y - r, width: r*2, height: r*2)), with: .color(col), lineWidth: 1.3)
        ctx.stroke(Path(ellipseIn: CGRect(x: c.x - r, y: c.y - r*0.2, width: r*2, height: r*2)), with: .color(col), lineWidth: 1.3)
    }

    /// 400 kV switchyard one-line filling the top-right corner: generator breaker
    /// 52G → GSU transformer → HV bus → two transmission feeders (LINE 1 in
    /// service, LINE 2 on standby). Shows component + breaker states plus the
    /// headline grid numbers (kV / Hz / kA). Brass when energized; the gen side
    /// de-energizes on a unit trip while the 400 kV grid bus stays live. All
    /// colour collapses to graphite in the AUTHENTIC skins.
    private func drawSwitchyard(_ ctx: GraphicsContext, w: CGFloat, h: CGFloat,
                                gen: CGRect, snap: PlantSnapshot, sup: PlantSupervisor) {
        func fx(_ v: CGFloat) -> CGFloat { v * w }
        func fy(_ v: CGFloat) -> CGFloat { v * h }
        let trip    = sup.turbineTrip
        let g52Open = trip || sup.genBreakerOpen                  // gen breaker: open on trip or by hand
        let genCol  = g52Open ? Theme.deEnergized : Theme.elecGold
        let bus     = Theme.elecGold                              // 400 kV grid bus: always energized
        let dark    = Theme.deEnergized.opacity(0.6)             // de-energized conductor
        let okC     = Theme.isFlat ? Theme.textHdr : Theme.statusNormal

        let cx   = fx(0.852)
        let xfR  = fx(0.013)
        let xfY  = gen.midY - 0.8 * xfR   // lower winding centred on the gen-lead line → clean side entry
        let gbY  = fy(0.148)
        let busY = fy(0.118)
        let busL = fx(0.828), busR = fx(0.966)
        let fA   = fx(0.873), fB = fx(0.938)
        let brkY = fy(0.085)
        let topY = fy(0.054)

        // Generator output lead → GSU low-voltage winding. One straight run into
        // the lower circle's edge — no bend, so no chamfer artifact at the joint.
        pipe(ctx, [(gen.maxX, gen.midY), (cx - xfR, gen.midY)], genCol, 2)
        // GSU transformer (two-winding).
        transformer(ctx, CGPoint(x: cx, y: xfY), r: xfR, live: !g52Open)
        txt(ctx, "GSU 21/400kV", CGPoint(x: cx + xfR + fx(0.005), y: xfY), .leading, Theme.textDim, 7)
        // GSU HV winding → generator breaker 52G → HV bus (single vertical run).
        pipe(ctx, [(cx, xfY - xfR), (cx, busY)], g52Open ? dark : bus, 2)
        breaker(ctx, CGPoint(x: cx, y: gbY), closed: !g52Open, color: g52Open ? Theme.deEnergized : bus)
        txt(ctx, "52G", CGPoint(x: cx + fx(0.009), y: gbY), .leading, Theme.textDim, 7)
        // HV busbar — the 400 kV grid stays energized even on a unit trip.
        ctx.stroke(Path { p in p.move(to: CGPoint(x: busL, y: busY)); p.addLine(to: CGPoint(x: busR, y: busY)) },
                   with: .color(bus), style: StrokeStyle(lineWidth: 3, lineCap: .round))
        txt(ctx, "400 kV 3\u{03C6}", CGPoint(x: busL, y: busY - fy(0.017)), .leading, Theme.textHdr, 8)
        txt(ctx, "50.00 Hz", CGPoint(x: busR, y: busY - fy(0.017)), .trailing, okC, 8)
        // Two 400 kV circuits — click a breaker to open/close it. One is redundant
        // (the other then carries full load); open BOTH → full load reject → trip.
        // The 3φ tick marks show each single conductor carries all three phases.
        let lineOpen = [sup.line1BreakerOpen, sup.line2BreakerOpen]
        let closedN  = lineOpen.filter { !$0 }.count
        let iTot     = g52Open ? 0 : snap.electricPowerW / (1.732 * 400e3 * 0.92) / 1000   // kA exported
        for (i, x) in [fA, fB].enumerated() {
            let open = lineOpen[i]
            let lc   = open ? dark : bus
            pipe(ctx, [(x, busY), (x, topY + 5)], lc, 2)
            phaseTicks(ctx, CGPoint(x: x, y: (busY + brkY) / 2), color: lc)
            breaker(ctx, CGPoint(x: x, y: brkY), closed: !open, color: open ? Theme.deEnergized : bus)
            let kA = open ? 0 : (closedN > 0 ? iTot / Double(closedN) : 0)
            txt(ctx, open ? "OPEN" : String(format: "%.2f kA", kA),
                CGPoint(x: x + fx(0.009), y: brkY), .leading,
                open ? Theme.caution : Theme.textDim, 7)
            tower(ctx, CGPoint(x: x, y: topY), color: lc)
            txt(ctx, "LINE \(i + 1)", CGPoint(x: x, y: topY - fy(0.017)), .center, Theme.textDim, 7)
        }
    }

    /// Switchyard breaker (IEEE square): filled = closed, hollow = open.
    private func breaker(_ ctx: GraphicsContext, _ c: CGPoint, closed: Bool, color: Color) {
        let s: CGFloat = 7
        let p = Path(roundedRect: CGRect(x: c.x - s/2, y: c.y - s/2, width: s, height: s), cornerRadius: 1.5)
        if closed { ctx.fill(p, with: .color(color)) }
        else { ctx.fill(p, with: .color(Theme.equipFill)); ctx.stroke(p, with: .color(color), lineWidth: 1.3) }
    }

    /// Transmission-line terminal: an arrowhead toward the grid + a short cap.
    private func tower(_ ctx: GraphicsContext, _ c: CGPoint, color: Color) {
        let a: CGFloat = 5
        ctx.stroke(Path { p in
            p.move(to: CGPoint(x: c.x - a, y: c.y + a)); p.addLine(to: CGPoint(x: c.x, y: c.y)); p.addLine(to: CGPoint(x: c.x + a, y: c.y + a))
        }, with: .color(color), style: StrokeStyle(lineWidth: 1.6, lineJoin: .round))
        ctx.stroke(Path { p in p.move(to: CGPoint(x: c.x - a, y: c.y)); p.addLine(to: CGPoint(x: c.x + a, y: c.y)) },
                   with: .color(color.opacity(0.55)), lineWidth: 1.2)
    }

    /// Three-phase tick marks: three short diagonal strokes across a (vertical)
    /// conductor — the one-line convention that a single line carries 3 phases.
    private func phaseTicks(_ ctx: GraphicsContext, _ c: CGPoint, color: Color) {
        let gap: CGFloat = 3.5, len: CGFloat = 3.6
        for k in 0..<3 {
            let yy = c.y - 3.5 + CGFloat(k) * gap
            ctx.stroke(Path { p in
                p.move(to: CGPoint(x: c.x - len, y: yy + len)); p.addLine(to: CGPoint(x: c.x + len, y: yy - len))
            }, with: .color(color), lineWidth: 1.1)
        }
    }

    /// Right dock under the electrical picture: turbine-generator mechanical +
    /// excitation supervision. The physics model doesn't simulate bearings / H₂
    /// cooling, so these are representative instrument values that track load and
    /// trip state — labeled honestly, same spirit as the other supervisory reads.
    private func drawTurbineGen(_ ctx: GraphicsContext, _ r: CGRect, snap: PlantSnapshot, sup: PlantSupervisor) {
        dockField(ctx, r, "TURBINE-GENERATOR")
        let trip = sup.turbineTrip
        let mwe  = snap.electricPowerW / 1e6
        let load = max(0, min(1, mwe / 990))
        let brgT = 60 + 32 * load
        let statT = 55 + 42 * load
        let statI = trip ? 0 : mwe / (1.732 * 21.0 * 0.92)     // kA at the 21 kV terminals
        let items: [(String, String, Color)] = [
            ("SPEED",    trip ? "0 rpm" : "3000 rpm", trip ? Theme.alarm : Theme.ink),
            ("BRG MTL",  String(format: "%.0f°C", brgT), brgT > 105 ? Theme.caution : Theme.ink),
            ("O/SPD",    "3300 rpm", Theme.textDim),
            ("LUBE",     "2.1 bar", Theme.ink),
            ("VIB",      String(format: "%.1f mil", 1.1 + 0.7 * load), Theme.ink),
            ("H2 P",     "4.0 bar", Theme.ink),
            ("ECC",      String(format: "%.2f mil", 0.30 + 0.12 * load), Theme.ink),
            ("H2 PUR",   "98.5%", Theme.ink),
            ("DIFF EXP", String(format: "%+.1f mm", -0.4 + 1.6 * load), Theme.ink),
            ("STAT T",   String(format: "%.0f°C", statT), statT > 110 ? Theme.caution : Theme.ink),
            ("THRUST",   String(format: "%.0f%%", 22 + 46 * load), Theme.ink),
            ("STAT I",   trip ? "0 kA" : String(format: "%.1f kA", statI), Theme.ink),
            ("FIELD V",  trip ? "0 V"  : String(format: "%.0f V", 250 + 150 * load), Theme.ink),
            ("SEAL OIL", "4.5 bar", Theme.ink),
            ("FIELD A",  trip ? "0 A"  : String(format: "%.0f A", 2400 + 1400 * load), Theme.ink),
            ("GAS T",    String(format: "%.0f°C", 40 + 6 * load), Theme.ink),
        ]
        let cols = 2, rpc = (items.count + cols - 1) / cols
        let gx: CGFloat = 8
        let colW = (r.width - gx * 2) / CGFloat(cols)
        let top = r.minY + 20
        let rowH = (r.maxY - top - 6) / CGFloat(rpc)
        for (i, it) in items.enumerated() {
            let col = i % cols, row = i / cols
            let cx = r.minX + gx + CGFloat(col) * colW
            let cy = top + CGFloat(row) * rowH + rowH / 2
            ctx.draw(Text(it.0).font(.system(size: 8, design: .monospaced)).foregroundColor(Theme.textDim),
                     at: CGPoint(x: cx, y: cy), anchor: .leading)
            ctx.draw(Text(it.1).font(.system(size: 10, weight: .semibold, design: .monospaced)).foregroundColor(it.2),
                     at: CGPoint(x: cx + colW - 14, y: cy), anchor: .trailing)
        }
    }

    private func condenser(_ ctx: GraphicsContext, _ r: CGRect, vacOK: Bool) {
        equipBox(ctx, r, "", nil)
        ctx.draw(Text("CONDENSER").font(.system(size: 10, design: .monospaced)).foregroundColor(Theme.textHdr),
                 at: CGPoint(x: r.midX, y: r.minY + 9), anchor: .top)
        // hotwell band
        let hw = CGRect(x: r.minX + 2, y: r.maxY - r.height * 0.22, width: r.width - 4, height: r.height * 0.22 - 2)
        ctx.fill(Path(hw), with: .color(Theme.water.opacity(0.35)))
        ctx.stroke(Path(CGRect(x: hw.minX, y: hw.minY, width: hw.width, height: 1)),
                   with: .color(Theme.water.opacity(0.6)), lineWidth: 1)
        // cooling-water tube hint
        for i in 0..<3 {
            let y = r.minY + r.height * (0.34 + 0.13 * CGFloat(i))
            ctx.stroke(Path { p in p.move(to: CGPoint(x: r.minX + 6, y: y)); p.addLine(to: CGPoint(x: r.maxX - 6, y: y)) },
                       with: .color(Theme.ink.opacity(0.10)), lineWidth: 1)
        }
        txt(ctx, vacOK ? "VAC ●" : "VAC ▲", CGPoint(x: r.minX + 7, y: r.minY + 7), .leading,
            vacOK ? (Theme.isFlat ? Theme.textHdr : Theme.statusNormal) : Theme.caution, 8)
    }

    /// Feedwater heater: box with the crossed-tube heat-exchanger hint.
    private func heaterBox(_ ctx: GraphicsContext, _ r: CGRect, _ label: String) {
        equipBox(ctx, r, "", nil)
        ctx.stroke(Path { p in
            p.move(to: CGPoint(x: r.minX + 4, y: r.midY)); p.addLine(to: CGPoint(x: r.maxX - 4, y: r.midY))
        }, with: .color(Theme.ink.opacity(0.22)), lineWidth: 1)
        ctx.draw(Text(label).font(.system(size: 8, design: .monospaced)).foregroundColor(Theme.textDim),
                 at: CGPoint(x: r.midX, y: r.maxY + 8), anchor: .top)
    }

    private func txt(_ ctx: GraphicsContext, _ s: String, _ pos: CGPoint, _ anchor: UnitPoint, _ color: Color, _ size: CGFloat) {
        ctx.draw(Text(s).font(.system(size: size, design: .monospaced)).foregroundColor(color), at: pos, anchor: anchor)
    }

    private func reading(_ ctx: GraphicsContext, _ pos: CGPoint, _ s: String, _ color: Color, _ size: CGFloat, _ anchor: UnitPoint) {
        ctx.draw(Text(s).font(.system(size: size, weight: .semibold, design: .monospaced)).foregroundColor(color), at: pos, anchor: anchor)
    }
}
