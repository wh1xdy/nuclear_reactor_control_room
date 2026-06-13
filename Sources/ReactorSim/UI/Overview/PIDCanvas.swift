// PIDCanvas.swift — PWR P&ID schematic drawn in SwiftUI Canvas.
// Layout uses fractional coordinates so it fills the full panel at any size.
// All pipes are orthogonal. Animated flow dots travel at speed ∝ flow rate.

import SwiftUI

struct PIDCanvas: View {
    let supervisor: PlantSupervisor

    var body: some View {
        // snapshot is read HERE (leaf), not in the parent panel, so the glass
        // panel subtree doesn't re-evaluate at 60 Hz with every physics tick.
        TimelineView(.animation) { context in
            let t = context.date.timeIntervalSinceReferenceDate
            Canvas { ctx, size in
                drawPWR(ctx: ctx, size: size, snap: supervisor.snapshot, sup: supervisor, t: t)
            }
        }
        // No .drawingGroup(): offscreen Metal rasterization races the Liquid Glass
        // backdrop and causes flicker during window drag.
        .background(Color(r: 8, g: 11, b: 16))
        .clipShape(.rect(cornerRadius: Theme.controlRadius, style: .continuous))
    }

    // MARK: — PWR layout (fractional of W × H so it always fills the panel)
    private func drawPWR(ctx: GraphicsContext, size: CGSize,
                         snap: PlantSnapshot, sup: PlantSupervisor, t: Double) {
        let W = size.width; let H = size.height
        // Fractional helpers
        let fx: (CGFloat) -> CGFloat = { $0 * W }
        let fy: (CGFloat) -> CGFloat = { $0 * H }

        // ── Component rects (fraction of W × H) ──────────────────────────────
        let rv  = CGRect(x: fx(0.04),  y: fy(0.18), width: fx(0.14), height: fy(0.65))
        let pz  = CGRect(x: fx(0.055), y: fy(0.02), width: fx(0.10), height: fy(0.14))
        let sg  = CGRect(x: fx(0.38),  y: fy(0.08), width: fx(0.12), height: fy(0.70))
        let rcp = CGRect(x: fx(0.22),  y: fy(0.80), width: fx(0.12), height: fy(0.14))
        let tb  = CGRect(x: fx(0.65),  y: fy(0.08), width: fx(0.13), height: fy(0.22))
        let cd  = CGRect(x: fx(0.65),  y: fy(0.68), width: fx(0.13), height: fy(0.20))

        // ── Pipe routing ──────────────────────────────────────────────────────
        let frac     = max(0, min(1, snap.powerFraction))
        let flow     = snap.scrammed ? 0.0 : max(0, Double(sup.primaryFlow) * Double(sup.omegaRCP))
        let hotColor = Theme.hotLeg          // primary always hot at operating temp
        let cldColor = Theme.water
        let stmColor = sup.turbineValve > 0.05 ? Theme.twophase : Theme.steam

        // Midpoint between reactor and SG for L-routing
        let midX  = (rv.maxX + sg.minX) / 2
        let coldY = rcp.midY
        let fwY   = cd.midY
        let fwMid = sg.maxX + (cd.minX - sg.maxX) * 0.5

        // Hot leg: reactor right → midX → down → SG top-left inlet
        pipe(ctx, [(rv.maxX, rv.minY + fy(0.12)),
                   (midX,    rv.minY + fy(0.12)),
                   (midX,    sg.minY + fy(0.08)),
                   (sg.minX, sg.minY + fy(0.08))], hotColor, w: 4)

        // Cold leg: SG bottom outlet → midX → down to RCP level → RCP right
        pipe(ctx, [(sg.minX, sg.maxY - fy(0.08)),
                   (midX,    sg.maxY - fy(0.08)),
                   (midX,    coldY),
                   (rcp.maxX, coldY)], cldColor, w: 4)

        // Cold leg: RCP left → reactor bottom inlet
        pipe(ctx, [(rcp.minX, coldY),
                   (rv.midX,  coldY),
                   (rv.midX,  rv.maxY)], cldColor, w: 4)

        // Pressurizer surge line — plain line work, not a color of its own
        pipe(ctx, [(pz.midX, pz.maxY),
                   (rv.midX, rv.minY + fy(0.05))], Color.white.opacity(0.25), w: 1.5)

        // Main steam: SG right → horizontal run → TURB left
        pipe(ctx, [(sg.maxX, sg.minY + fy(0.10)),
                   (tb.minX, sg.minY + fy(0.10)),
                   (tb.minX, tb.midY)], stmColor, w: 3)

        // Turbine exhaust: TURB bottom → COND top
        pipe(ctx, [(tb.midX, tb.maxY),
                   (tb.midX, cd.minY)], Theme.steam, w: 3)

        // Feedwater: COND left → horizontal → SG bottom right inlet
        pipe(ctx, [(cd.minX,   fwY),
                   (fwMid,     fwY),
                   (sg.maxX,   fwY),
                   (sg.maxX,   sg.maxY - fy(0.10))], cldColor, w: 3)

        // ── Flow animation — SCADA-style scrolling dash overlay ───────────────
        // Not marching pellets. A faint, lighter dash pattern slides ALONG each
        // pipe in the direction of flow — the standard control-room "this line is
        // live and moving" cue. Speed ∝ flow; freezes when pumps coast down.
        let hotPath: [(CGFloat, CGFloat)] = [
            (rv.maxX, rv.minY+fy(0.12)), (midX, rv.minY+fy(0.12)),
            (midX, sg.minY+fy(0.08)),    (sg.minX, sg.minY+fy(0.08))]
        let coldPathA: [(CGFloat, CGFloat)] = [
            (sg.minX, sg.maxY-fy(0.08)), (midX, sg.maxY-fy(0.08)),
            (midX, coldY),               (rcp.maxX, coldY)]
        let coldPathB: [(CGFloat, CGFloat)] = [
            (rcp.minX, coldY), (rv.midX, coldY), (rv.midX, rv.maxY)]
        // dashPhase decreases over time → dashes travel in the path's direction.
        let phase = -(t * flow * 22.0)
        flowDash(ctx, hotPath,   phase: phase)
        flowDash(ctx, coldPathA, phase: phase)
        flowDash(ctx, coldPathB, phase: phase)
        // Secondary side: steam out
        if sup.turbineValve > 0.02 && !snap.scrammed {
            flowDash(ctx, [(sg.maxX, sg.minY+fy(0.10)), (tb.minX, sg.minY+fy(0.10)),
                           (tb.minX, tb.midY)], phase: -(t * Double(sup.turbineValve) * 22.0))
        }

        // Hot/cold leg temperatures — the loop ΔT from the two-node coolant model
        drawAnnotation(ctx, pos: CGPoint(x: (rv.maxX + midX) / 2 - 16, y: rv.minY + fy(0.12) - 10),
                       text: snap.hotLegTempK.isSafe ? String(format: "%.0f K", snap.hotLegTempK) : "---",
                       color: Theme.hotLeg, fontSize: 9)
        drawAnnotation(ctx, pos: CGPoint(x: (rcp.maxX + midX) / 2 - 16, y: coldY + 7),
                       text: snap.coldLegTempK.isSafe ? String(format: "%.0f K", snap.coldLegTempK) : "---",
                       color: Theme.water, fontSize: 9)

        // ── Component boxes — industrial line work ────────────────────────────
        // Equipment is NEUTRAL: one fill, one hairline. Color belongs to the
        // process (pipes) and to live alarm states only.
        let fuelAlarm  = snap.fuelTempK > 1200
        let pressAlarm = sup.pressureMPa > 17.0

        // Reactor vessel: rod-position indicators (drop from top with
        // insertion) over a faint power tint — an instrument, not a lava lamp.
        drawBox(ctx, rv, label: "", color: fuelAlarm ? Theme.alarm : nil)
        if frac > 0.01 {
            let core = CGRect(x: rv.minX + fy(0.025), y: rv.minY + fy(0.10),
                              width: rv.width - fy(0.05), height: rv.height - fy(0.16))
            ctx.fill(Path(roundedRect: core, cornerRadius: 3, style: .continuous),
                     with: .color(Theme.hotLeg.opacity(0.05 + 0.10 * frac)))
        }
        for ci in 0..<4 {
            let rx = rv.minX + rv.width * CGFloat(ci + 1) / 5.0 - fy(0.025)
            let guide = CGRect(x: rx, y: rv.minY + fy(0.10),
                               width: fy(0.05), height: rv.height - fy(0.16))
            // Channel guide
            ctx.stroke(Path(roundedRect: guide, cornerRadius: 2, style: .continuous),
                       with: .color(.white.opacity(0.12)), lineWidth: 1)
            // Rod inserted from the top, depth = actual rod position
            let depth = guide.height * CGFloat(max(0, min(1, snap.rodPosition)))
            if depth > 1 {
                ctx.fill(Path(roundedRect: CGRect(x: guide.minX + 1.5, y: guide.minY,
                                                  width: guide.width - 3, height: depth),
                              cornerRadius: 1.5, style: .continuous),
                         with: .color(.white.opacity(0.45)))
            }
        }
        drawLabel(ctx, snap.scrammed ? "REACTOR — TRIP" : "REACTOR", rect: rv, offsetY: fy(0.015),
                  color: snap.scrammed ? Theme.alarm : (fuelAlarm ? Theme.alarm : Theme.textHdr),
                  fontSize: 11)
        drawAnnotation(ctx, pos: CGPoint(x: rv.maxX + 8, y: rv.maxY - fy(0.10)),
                       text: snap.powerFraction.isSafe ? String(format: "%5.1f %%", snap.powerFraction*100) : "---",
                       color: snap.powerFraction > 1.1 ? Theme.alarm : .white, fontSize: 10)
        drawAnnotation(ctx, pos: CGPoint(x: rv.maxX + 8, y: rv.maxY - fy(0.04)),
                       text: snap.fuelTempK.isSafe ? String(format: "%5.0f K", snap.fuelTempK) : "---",
                       color: fuelAlarm ? Theme.alarm : Theme.textDim, fontSize: 9)

        // PZR
        drawBox(ctx, pz, label: "PZR", color: pressAlarm ? Theme.alarm : nil)
        drawAnnotation(ctx, pos: CGPoint(x: pz.maxX + 6, y: pz.midY - 6),
                       text: sup.pressureMPa.isSafe ? String(format: "%.2f MPa", sup.pressureMPa) : "---",
                       color: pressAlarm ? Theme.alarm : Theme.textDim, fontSize: 9)

        // S/G
        drawBox(ctx, sg, label: "S/G", color: nil)
        drawLabel(ctx, snap.sgTempK.isSafe ? String(format: "%.0f K", snap.sgTempK) : "---",
                  rect: sg, offsetY: fy(0.18), color: Theme.textDim, fontSize: 9)

        // RCP
        drawBox(ctx, rcp, label: "RCP", color: nil)
        drawPump(ctx, center: CGPoint(x: rcp.midX, y: rcp.midY),
                 r: min(rcp.width, rcp.height) * 0.28, running: flow > 0.1)

        // TURB — tripped state is a real alarm; valve % is just a number
        drawBox(ctx, tb, label: "TURB", color: sup.turbineTrip ? Theme.alarm : nil)
        drawLabel(ctx, sup.turbineTrip ? "TRIPPED"
                                       : String(format: "%.0f MWe", snap.electricPowerW / 1e6),
                  rect: tb, offsetY: fy(0.14),
                  color: sup.turbineTrip ? Theme.alarm : .white, fontSize: 9)

        // COND
        drawBox(ctx, cd, label: "COND", color: nil)
        drawLabel(ctx, sup.condTempK.isSafe ? String(format: "%.0f K", sup.condTempK) : "---",
                  rect: cd, offsetY: fy(0.12), color: Theme.textDim, fontSize: 9)

        // FW pump
        let fwPumpX = (cd.minX + sg.maxX) / 2
        drawPump(ctx, center: CGPoint(x: fwPumpX, y: fwY), r: H * 0.055,
                 running: sup.feedwaterValve > 0.1 && !sup.feedwaterFault)
    }

    // MARK: — Primitives

    private func pipe(_ ctx: GraphicsContext, _ pts: [(CGFloat, CGFloat)],
                      _ color: Color, w: CGFloat) {
        guard pts.count >= 2 else { return }
        var path = Path()
        path.move(to: CGPoint(x: pts[0].0, y: pts[0].1))
        for p in pts.dropFirst() { path.addLine(to: CGPoint(x: p.0, y: p.1)) }
        ctx.stroke(path, with: .color(color), lineWidth: w)
    }

    // Scrolling-dash flow overlay: a lighter dashed stroke laid over the pipe,
    // its dashPhase animated so the gaps slide along the line. Reads as moving
    // fluid without the toy "pellet" look.
    private func flowDash(_ ctx: GraphicsContext, _ pts: [(CGFloat, CGFloat)],
                          phase: Double) {
        guard pts.count >= 2 else { return }
        var path = Path()
        path.move(to: CGPoint(x: pts[0].0, y: pts[0].1))
        for p in pts.dropFirst() { path.addLine(to: CGPoint(x: p.0, y: p.1)) }
        let style = StrokeStyle(lineWidth: 2, lineCap: .butt,
                                dash: [3, 9], dashPhase: CGFloat(phase))
        ctx.stroke(path, with: .color(.white.opacity(0.55)), style: style)
    }

    // Neutral equipment body: one fill, one hairline. `color` overrides the
    // outline ONLY for a live alarm state (red); nil = standard line work.
    private func drawBox(_ ctx: GraphicsContext, _ rect: CGRect,
                         label: String, color: Color?) {
        let shape = Path(roundedRect: rect, cornerRadius: 6, style: .continuous)
        ctx.fill(shape, with: .color(Color(r: 15, g: 17, b: 21)))
        ctx.stroke(shape, with: .color(color ?? .white.opacity(0.22)),
                   lineWidth: color == nil ? 1 : 1.5)
        if !label.isEmpty {
            ctx.draw(Text(label).font(.system(size: 11, design: .monospaced))
                        .foregroundColor(color ?? Theme.textHdr),
                     at: CGPoint(x: rect.midX, y: rect.minY + 8), anchor: .top)
        }
    }

    private func drawLabel(_ ctx: GraphicsContext, _ text: String, rect: CGRect,
                           offsetY: CGFloat, color: Color, fontSize: CGFloat) {
        ctx.draw(Text(text).font(.system(size: fontSize, design: .monospaced))
                    .foregroundColor(color),
                 at: CGPoint(x: rect.midX, y: rect.minY + offsetY), anchor: .top)
    }

    private func drawAnnotation(_ ctx: GraphicsContext, pos: CGPoint,
                                text: String, color: Color, fontSize: CGFloat = 9) {
        ctx.draw(Text(text).font(.system(size: fontSize, design: .monospaced))
                    .foregroundColor(color),
                 at: pos, anchor: .leading)
    }

    // ISA-style pump symbol: circle + impeller triangle. Triangle fills accent
    // when running, stays hollow grey when stopped. No orbiting decoration.
    private func drawPump(_ ctx: GraphicsContext, center: CGPoint, r: CGFloat,
                          running: Bool) {
        let rect = CGRect(x: center.x-r, y: center.y-r, width: r*2, height: r*2)
        ctx.fill(Path(ellipseIn: rect), with: .color(Color(r: 15, g: 17, b: 21)))
        ctx.stroke(Path(ellipseIn: rect), with: .color(.white.opacity(0.25)), lineWidth: 1)
        var tri = Path()
        tri.move(to: CGPoint(x: center.x - r*0.40, y: center.y - r*0.50))
        tri.addLine(to: CGPoint(x: center.x - r*0.40, y: center.y + r*0.50))
        tri.addLine(to: CGPoint(x: center.x + r*0.55, y: center.y))
        tri.closeSubpath()
        if running {
            ctx.fill(tri, with: .color(Theme.accent.opacity(0.85)))
        } else {
            ctx.stroke(tri, with: .color(.white.opacity(0.30)), lineWidth: 1)
        }
    }
}
