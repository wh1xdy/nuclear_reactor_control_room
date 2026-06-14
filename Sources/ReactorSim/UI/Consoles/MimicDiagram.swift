// MimicDiagram.swift — detailed full-screen PWR plant mimic.
// A purpose-built schematic (not the small panel P&ID) with proper component
// shapes, level indication, valves, flow animation, and live values embedded
// at every system — sized to fill the whole console.

import SwiftUI

struct MimicDiagram: View {
    let supervisor: PlantSupervisor

    var body: some View {
        TimelineView(.animation) { context in
            let t = context.date.timeIntervalSinceReferenceDate
            Canvas { ctx, size in
                draw(ctx, size, snap: supervisor.snapshot, sup: supervisor, t: t)
            }
        }
        .background(Theme.schematicBg)
    }

    private func draw(_ ctx: GraphicsContext, _ size: CGSize,
                      snap: PlantSnapshot, sup: PlantSupervisor, t: Double) {
        let W = size.width, H = size.height
        func fx(_ v: CGFloat) -> CGFloat { v * W }
        func fy(_ v: CGFloat) -> CGFloat { v * H }

        let scram = snap.scrammed
        let flow  = scram ? 0 : max(0, Double(sup.primaryFlow) * Double(sup.omegaRCP))
        let steam = sup.turbineValve

        // ── Component frames ────────────────────────────────────────────────
        // Both reactor nozzles are on the RIGHT face of the vessel (hot upper,
        // cold lower), so the bottom and left of the vessel stay clear for labels.
        let vessel = CGRect(x: fx(0.05),  y: fy(0.30), width: fx(0.085), height: fy(0.46))
        let pzr    = CGRect(x: fx(0.16),  y: fy(0.05), width: fx(0.045), height: fy(0.22))
        let sg     = CGRect(x: fx(0.30),  y: fy(0.12), width: fx(0.09),  height: fy(0.60))
        let hpT    = CGRect(x: fx(0.50),  y: fy(0.16), width: fx(0.075), height: fy(0.10))
        let lpT    = CGRect(x: fx(0.585), y: fy(0.13), width: fx(0.115), height: fy(0.16))
        let gen    = CGRect(x: fx(0.72),  y: fy(0.155), width: fx(0.075), height: fy(0.105))
        let cond   = CGRect(x: fx(0.62),  y: fy(0.56), width: fx(0.20),  height: fy(0.17))
        let fwHtr  = CGRect(x: fx(0.47),  y: fy(0.82), width: fx(0.07),  height: fy(0.07))

        // Key node points
        let hotY   = vessel.minY + fy(0.06)        // hot nozzle (upper-right)
        let coldY  = vessel.maxY - fy(0.12)        // cold nozzle (lower-right)
        let midX   = (vessel.maxX + sg.minX) / 2
        let rcp    = CGRect(x: midX - fx(0.03), y: coldY - fy(0.05), width: fx(0.06), height: fy(0.10))
        let sgHotIn   = CGPoint(x: sg.minX, y: sg.minY + fy(0.07))
        let sgColdOut = CGPoint(x: sg.minX, y: sg.maxY - fy(0.06))
        let sgSteam   = CGPoint(x: sg.maxX, y: sg.minY + fy(0.06))
        let sgFeed    = CGPoint(x: sg.maxX, y: sg.maxY - fy(0.10))

        // ── Pipes (primary) ─────────────────────────────────────────────────
        // Hot leg out the upper-right nozzle, up and over to the SG inlet.
        let hot = [(vessel.maxX, hotY), (midX, hotY), (midX, sgHotIn.y), (sg.minX, sgHotIn.y)]
        // Cold leg from SG bottom, straight across through the RCP into the
        // lower-right nozzle — one clean horizontal run, nothing under the vessel.
        let coldA = [(sgColdOut.x, sgColdOut.y), (sgColdOut.x, coldY), (rcp.maxX, coldY)]
        let coldB = [(rcp.minX, coldY), (vessel.maxX, coldY)]
        pipe(ctx, hot,   Theme.hotLeg, 4)
        pipe(ctx, coldA, Theme.water, 4)
        pipe(ctx, coldB, Theme.water, 4)
        // Surge line PZR → hot leg: a single vertical drop that tees into the run.
        pipe(ctx, [(pzr.midX, pzr.maxY), (pzr.midX, hotY)], Theme.ink.opacity(0.4), 2)

        // ── Pipes (secondary) ───────────────────────────────────────────────
        let msivX = sgSteam.x + fx(0.05)
        let steamHdr = [(sgSteam.x, sgSteam.y), (msivX, sgSteam.y), (hpT.minX, sgSteam.y), (hpT.minX, hpT.midY)]
        pipe(ctx, steamHdr, steam > 0.05 ? Theme.twophase : Theme.steam, 3)
        // turbine exhaust → condenser
        pipe(ctx, [(lpT.midX, lpT.maxY), (lpT.midX, cond.minY)], Theme.steam, 3)
        // condenser → condensate/feed → SG
        let fwY = fwHtr.midY
        pipe(ctx, [(cond.minX, cond.maxY - fy(0.03)), (cond.minX, fwY), (fwHtr.maxX, fwY)], Theme.water, 3)
        pipe(ctx, [(fwHtr.minX, fwY), (sgFeed.x + fx(0.02), fwY), (sgFeed.x + fx(0.02), sgFeed.y), (sgFeed.x, sgFeed.y)], Theme.water, 3)
        // steam dump (active post-trip)
        if sup.steamDumpValve > 0.001 {
            pipe(ctx, [(msivX, sgSteam.y), (msivX, cond.minY - fy(0.02)), (cond.midX, cond.minY - fy(0.02)), (cond.midX, cond.minY)],
                 Theme.steam.opacity(0.7), 2)
            txt(ctx, "DUMP", CGPoint(x: msivX + 4, y: (sgSteam.y + cond.minY)/2), .leading, Theme.caution, 8)
        }

        // ── Flow animation ──────────────────────────────────────────────────
        let ph = -(t * flow * 22)
        flowDash(ctx, hot, ph); flowDash(ctx, coldA, ph); flowDash(ctx, coldB, ph)
        if steam > 0.02 && !scram {
            flowDash(ctx, steamHdr, -(t * Double(steam) * 22))
        }

        // ── Valves ──────────────────────────────────────────────────────────
        valve(ctx, CGPoint(x: msivX, y: sgSteam.y), open: steam > 0.02 && !scram)

        // ── Reactor vessel + core + rods ────────────────────────────────────
        let fuelAlarm = snap.fuelTempK > 1200
        box(ctx, vessel, "REACTOR", fuelAlarm ? Theme.alarm : nil)
        let core = vessel.insetBy(dx: fx(0.012), dy: fy(0.06))
        let frac = max(0, min(1, snap.powerFraction))
        ctx.fill(Path(roundedRect: core, cornerRadius: 2, style: .continuous),
                 with: .color(Theme.hotLeg.opacity(0.06 + 0.16 * frac)))
        for i in 0..<4 {
            let rx = core.minX + core.width * CGFloat(i) / 4 + core.width * 0.06
            let guide = CGRect(x: rx, y: core.minY, width: core.width * 0.13, height: core.height)
            ctx.stroke(Path(guide), with: .color(Theme.ink.opacity(0.18)), lineWidth: 1)
            let depth = guide.height * CGFloat(max(0, min(1, snap.rodPosition)))
            if depth > 1 {
                ctx.fill(Path(CGRect(x: guide.minX + 1, y: guide.minY, width: guide.width - 2, height: depth)),
                         with: .color(Theme.ink.opacity(0.5)))
            }
        }
        // Reactor data block stacked BELOW the vessel — the bottom is now clear.
        reading(ctx, CGPoint(x: vessel.midX, y: vessel.maxY + fy(0.04)),
                String(format: "%.1f%%", snap.powerFraction * 100),
                snap.powerFraction > 1.1 ? Theme.alarm : Theme.ink, 14, .center)
        txt(ctx, String(format: "FUEL %.0f K", snap.fuelTempK),
            CGPoint(x: vessel.midX, y: vessel.maxY + fy(0.078)), .center, Theme.textDim, 9)
        txt(ctx, "\(Int((228 * (1 - snap.rodPosition)).rounded())) SWD",
            CGPoint(x: vessel.midX, y: vessel.maxY + fy(0.105)), .center, Theme.textDim, 9)

        // ── Pressurizer with level + PORV ───────────────────────────────────
        let pressAlarm = sup.pressureMPa > 17
        levelTank(ctx, pzr, level: max(0, min(1, (sup.pressureMPa - 13) / 4)),
                  liquid: Theme.water, label: "PZR", alarm: pressAlarm ? Theme.alarm : nil)
        reading(ctx, CGPoint(x: pzr.maxX + fx(0.01), y: pzr.minY + fy(0.05)),
                String(format: "%.2f MPa", sup.pressureMPa), pressAlarm ? Theme.alarm : Theme.ink, 11, .leading)
        if sup.porvOpen { txt(ctx, "● PORV", CGPoint(x: pzr.maxX + fx(0.01), y: pzr.minY + fy(0.10)), .leading, Theme.caution, 9) }

        // ── Steam generator with secondary level ────────────────────────────
        levelTank(ctx, sg, level: max(0.1, min(0.9, sup.steamInv * 0.5)),
                  liquid: Theme.water, label: "S/G", alarm: nil)
        reading(ctx, CGPoint(x: sg.midX, y: sg.maxY + fy(0.03)),
                String(format: "%.2f MPa", snap.steamPressureMPa), Theme.ink, 11, .center)
        txt(ctx, String(format: "%.0f K", snap.sgTempK), CGPoint(x: sg.midX, y: sg.maxY + fy(0.06)), .center, Theme.textDim, 9)

        // ── Hot / cold leg labels — placed on clear pipe runs, not on top ────
        reading(ctx, CGPoint(x: vessel.maxX + fx(0.008), y: hotY + fy(0.022)),
                String(format: "T-HOT %.0f K", snap.hotLegTempK), Theme.hotLeg, 10, .leading)
        reading(ctx, CGPoint(x: (vessel.maxX + rcp.minX) / 2, y: coldY - fy(0.02)),
                String(format: "T-COLD %.0f K", snap.coldLegTempK), Theme.water, 9, .center)

        // ── RCP ─────────────────────────────────────────────────────────────
        box(ctx, rcp, "", sup.pumpDegraded ? Theme.alarm : nil)
        pump(ctx, CGPoint(x: rcp.midX, y: rcp.midY), r: min(rcp.width, rcp.height) * 0.32, running: flow > 0.1)
        txt(ctx, "RCP", CGPoint(x: rcp.midX, y: rcp.minY + fy(0.018)), .center, Theme.textHdr, 9)
        txt(ctx, String(format: "%.0f%%", sup.omegaRCP * 100), CGPoint(x: rcp.midX, y: rcp.maxY + fy(0.025)), .center, Theme.textDim, 9)

        // ── Turbine (HP + LP) + generator ───────────────────────────────────
        trapezoid(ctx, hpT, widenRight: true, tripped: sup.turbineTrip)
        trapezoid(ctx, lpT, widenRight: true, tripped: sup.turbineTrip)
        txt(ctx, "HP", CGPoint(x: hpT.midX, y: hpT.minY - fy(0.022)), .center, Theme.textHdr, 9)
        txt(ctx, "LP TURBINE", CGPoint(x: lpT.midX, y: lpT.minY - fy(0.022)), .center, Theme.textHdr, 9)
        // shaft to generator
        pipe(ctx, [(lpT.maxX, lpT.midY), (gen.minX, gen.midY)], Theme.ink.opacity(0.5), 3)
        box(ctx, gen, "", sup.turbineTrip ? Theme.alarm : nil)
        let genCx = CGPoint(x: gen.midX, y: gen.midY - fy(0.005))
        ctx.stroke(Path(ellipseIn: CGRect(x: genCx.x - fx(0.018), y: genCx.y - fx(0.018), width: fx(0.036), height: fx(0.036))),
                   with: .color(Theme.ink.opacity(0.5)), lineWidth: 1)
        txt(ctx, "G", genCx, .center, Theme.ink, 12)
        txt(ctx, "GEN", CGPoint(x: gen.midX, y: gen.minY - fy(0.022)), .center, Theme.textHdr, 9)
        reading(ctx, CGPoint(x: gen.midX, y: gen.maxY + fy(0.035)),
                sup.turbineTrip ? "TRIPPED" : String(format: "%.0f MWe", snap.electricPowerW / 1e6),
                sup.turbineTrip ? Theme.alarm : Theme.ink, 14, .center)

        // ── Condenser + feed train ──────────────────────────────────────────
        box(ctx, cond, "CONDENSER", nil)
        // cooling-water hatch + hotwell level
        let hw = CGRect(x: cond.minX + 2, y: cond.maxY - fy(0.05), width: cond.width - 4, height: fy(0.045))
        ctx.fill(Path(hw), with: .color(Theme.water.opacity(0.35)))
        reading(ctx, CGPoint(x: cond.midX, y: cond.minY + fy(0.05)),
                String(format: "%.0f K", sup.condTempK), Theme.ink, 11, .center)
        txt(ctx, "VAC", CGPoint(x: cond.minX + fx(0.02), y: cond.minY + fy(0.02)), .leading, Theme.textDim, 8)
        // feed components
        pump(ctx, CGPoint(x: cond.minX, y: fwY), r: fy(0.022), running: !sup.feedwaterFault && sup.feedwaterValve > 0.05)
        box(ctx, fwHtr, "FWH", nil)
        reading(ctx, CGPoint(x: fwHtr.midX, y: fwHtr.maxY + fy(0.03)),
                String(format: "FW %.0f%%", sup.feedwaterValve * 100),
                sup.feedwaterInv < 0.1 ? Theme.alarm : Theme.textDim, 9, .center)
    }

    // MARK: — primitives -----------------------------------------------------

    private func pipe(_ ctx: GraphicsContext, _ pts: [(CGFloat, CGFloat)], _ color: Color, _ w: CGFloat) {
        guard pts.count >= 2 else { return }
        var p = Path(); p.move(to: CGPoint(x: pts[0].0, y: pts[0].1))
        for q in pts.dropFirst() { p.addLine(to: CGPoint(x: q.0, y: q.1)) }
        ctx.stroke(p, with: .color(color), lineWidth: w)
    }

    private func flowDash(_ ctx: GraphicsContext, _ pts: [(CGFloat, CGFloat)], _ phase: Double) {
        guard pts.count >= 2 else { return }
        var p = Path(); p.move(to: CGPoint(x: pts[0].0, y: pts[0].1))
        for q in pts.dropFirst() { p.addLine(to: CGPoint(x: q.0, y: q.1)) }
        ctx.stroke(p, with: .color(Theme.ink.opacity(0.5)),
                   style: StrokeStyle(lineWidth: 2, lineCap: .butt, dash: [3, 9], dashPhase: CGFloat(phase)))
    }

    private func box(_ ctx: GraphicsContext, _ r: CGRect, _ label: String, _ alarm: Color?) {
        let s = Path(roundedRect: r, cornerRadius: 5, style: .continuous)
        ctx.fill(s, with: .color(Theme.equipFill))
        ctx.stroke(s, with: .color(alarm ?? Theme.ink.opacity(0.30)), lineWidth: alarm == nil ? 1 : 2)
        if !label.isEmpty {
            ctx.draw(Text(label).font(.system(size: 10, design: .monospaced)).foregroundColor(alarm ?? Theme.textHdr),
                     at: CGPoint(x: r.midX, y: r.minY + 9), anchor: .top)
        }
    }

    // Tank with a liquid level fill (0–1) and a single-line label.
    private func levelTank(_ ctx: GraphicsContext, _ r: CGRect, level: Double, liquid: Color, label: String, alarm: Color?) {
        let s = Path(roundedRect: r, cornerRadius: 5, style: .continuous)
        ctx.fill(s, with: .color(Theme.equipFill))
        let lv = max(0, min(1, level))
        let liq = CGRect(x: r.minX + 1.5, y: r.maxY - r.height * CGFloat(lv) + 1, width: r.width - 3, height: r.height * CGFloat(lv) - 2)
        ctx.fill(Path(roundedRect: liq, cornerRadius: 3, style: .continuous), with: .color(liquid.opacity(0.30)))
        ctx.stroke(Path(CGRect(x: liq.minX, y: liq.minY, width: liq.width, height: 1)), with: .color(liquid.opacity(0.7)), lineWidth: 1)
        ctx.stroke(s, with: .color(alarm ?? Theme.ink.opacity(0.30)), lineWidth: alarm == nil ? 1 : 2)
        ctx.draw(Text(label).font(.system(size: 10, weight: .semibold, design: .monospaced))
                    .foregroundColor(alarm ?? Theme.textHdr),
                 at: CGPoint(x: r.midX, y: r.minY + 12), anchor: .top)
    }

    private func pump(_ ctx: GraphicsContext, _ c: CGPoint, r: CGFloat, running: Bool) {
        let rect = CGRect(x: c.x - r, y: c.y - r, width: r*2, height: r*2)
        ctx.fill(Path(ellipseIn: rect), with: .color(Theme.equipFill))
        ctx.stroke(Path(ellipseIn: rect), with: .color(Theme.ink.opacity(0.3)), lineWidth: 1)
        var tri = Path()
        tri.move(to: CGPoint(x: c.x - r*0.4, y: c.y - r*0.5))
        tri.addLine(to: CGPoint(x: c.x - r*0.4, y: c.y + r*0.5))
        tri.addLine(to: CGPoint(x: c.x + r*0.55, y: c.y))
        tri.closeSubpath()
        ctx.fill(tri, with: .color(running ? Theme.accent : Theme.ink.opacity(0.3)))
    }

    private func valve(_ ctx: GraphicsContext, _ c: CGPoint, open: Bool) {
        let s: CGFloat = 6
        var p = Path()
        p.move(to: CGPoint(x: c.x - s, y: c.y - s)); p.addLine(to: CGPoint(x: c.x + s, y: c.y + s))
        p.addLine(to: CGPoint(x: c.x + s, y: c.y - s)); p.addLine(to: CGPoint(x: c.x - s, y: c.y + s)); p.closeSubpath()
        ctx.fill(p, with: .color(open ? Theme.accent.opacity(0.8) : Theme.ink.opacity(0.35)))
        ctx.stroke(p, with: .color(Theme.ink.opacity(0.5)), lineWidth: 1)
    }

    private func trapezoid(_ ctx: GraphicsContext, _ r: CGRect, widenRight: Bool, tripped: Bool) {
        var p = Path()
        let inset = r.height * 0.28
        p.move(to: CGPoint(x: r.minX, y: r.minY + inset))
        p.addLine(to: CGPoint(x: r.maxX, y: r.minY))
        p.addLine(to: CGPoint(x: r.maxX, y: r.maxY))
        p.addLine(to: CGPoint(x: r.minX, y: r.maxY - inset))
        p.closeSubpath()
        ctx.fill(p, with: .color(Theme.equipFill))
        ctx.stroke(p, with: .color(tripped ? Theme.alarm : Theme.ink.opacity(0.3)), lineWidth: tripped ? 2 : 1)
    }

    private func txt(_ ctx: GraphicsContext, _ s: String, _ pos: CGPoint, _ anchor: UnitPoint, _ color: Color, _ size: CGFloat) {
        ctx.draw(Text(s).font(.system(size: size, design: .monospaced)).foregroundColor(color), at: pos, anchor: anchor)
    }

    private func reading(_ ctx: GraphicsContext, _ pos: CGPoint, _ s: String, _ color: Color, _ size: CGFloat, _ anchor: UnitPoint) {
        ctx.draw(Text(s).font(.system(size: size, weight: .semibold, design: .monospaced)).foregroundColor(color), at: pos, anchor: anchor)
    }
}
