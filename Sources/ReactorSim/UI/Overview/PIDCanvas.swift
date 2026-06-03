// PIDCanvas.swift — PWR P&ID schematic drawn in SwiftUI Canvas.
// Layout uses fractional coordinates so it fills the full panel at any size.
// All pipes are orthogonal. Animated flow dots travel at speed ∝ flow rate.

import SwiftUI

struct PIDCanvas: View {
    let snapshot: PlantSnapshot
    let supervisor: PlantSupervisor

    var body: some View {
        TimelineView(.animation) { context in
            let t = context.date.timeIntervalSinceReferenceDate
            Canvas { ctx, size in
                drawPWR(ctx: ctx, size: size, snap: snapshot, sup: supervisor, t: t)
            }
        }
        .drawingGroup()
        .background(Color(r: 8, g: 11, b: 16))
        .clipShape(RoundedRectangle(cornerRadius: 6))
        .overlay(CornerBrackets())
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

        // Pressurizer surge line
        pipe(ctx, [(pz.midX, pz.maxY),
                   (rv.midX, rv.minY + fy(0.05))], Color(r: 140, g: 90, b: 210), w: 2)

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

        // ── Animated flow dots (primary loop) ────────────────────────────────
        let segs: [(CGPoint, CGPoint, Color)] = [
            (CGPoint(x: rv.maxX,  y: rv.minY+fy(0.12)), CGPoint(x: midX,     y: rv.minY+fy(0.12)), hotColor),
            (CGPoint(x: midX,     y: rv.minY+fy(0.12)), CGPoint(x: midX,     y: sg.minY+fy(0.08)), hotColor),
            (CGPoint(x: midX,     y: sg.minY+fy(0.08)), CGPoint(x: sg.minX,  y: sg.minY+fy(0.08)), hotColor),
            (CGPoint(x: sg.minX,  y: sg.maxY-fy(0.08)), CGPoint(x: midX,     y: sg.maxY-fy(0.08)), cldColor),
            (CGPoint(x: midX,     y: sg.maxY-fy(0.08)), CGPoint(x: midX,     y: coldY),             cldColor),
            (CGPoint(x: midX,     y: coldY),             CGPoint(x: rcp.maxX, y: coldY),             cldColor),
            (CGPoint(x: rcp.minX, y: coldY),             CGPoint(x: rv.midX,  y: coldY),             cldColor),
            (CGPoint(x: rv.midX,  y: coldY),             CGPoint(x: rv.midX,  y: rv.maxY),           cldColor),
        ]
        for (si, (p0, p1, dotCol)) in segs.enumerated() {
            for di in 0..<3 {
                let phase = t * flow * 0.5 + Double(si) / Double(segs.count) + Double(di) / 3.0
                let frac2 = CGFloat(phase.truncatingRemainder(dividingBy: 1.0))
                let pt = CGPoint(x: p0.x + (p1.x - p0.x) * frac2,
                                 y: p0.y + (p1.y - p0.y) * frac2)
                var dp = Path()
                dp.addEllipse(in: CGRect(x: pt.x-3, y: pt.y-3, width: 6, height: 6))
                ctx.fill(dp, with: .color(dotCol.opacity(0.85)))
            }
        }

        // ── Component boxes ───────────────────────────────────────────────────
        let fuelAlarm  = snap.fuelTempK > 1200
        let pressAlarm = sup.pressureMPa > 17.0

        // Reactor vessel — draw fuel glow then rods, label last
        drawBox(ctx, rv, label: "", color: fuelAlarm ? Theme.alarm : Theme.border,
                bg: fuelAlarm ? Color(r: 50, g: 8, b: 8) : Color(r: 14, g: 20, b: 34))
        if frac > 0.01 {
            let gh = (rv.height - fy(0.10)) * CGFloat(frac)
            let gr = min(255, 60 + Int(180 * frac))
            ctx.fill(Path(roundedRect: CGRect(x: rv.minX + fy(0.02),
                                              y: rv.maxY - gh - fy(0.02),
                                              width: rv.width - fy(0.04), height: gh),
                          cornerRadius: 3),
                     with: .color(Color(r: gr, g: Int(60*(1-frac)), b: 0).opacity(0.65)))
        }
        for ci in 0..<4 {
            let rx = rv.minX + rv.width * CGFloat(ci + 1) / 5.0 - fy(0.03)
            let rodCol = frac > 0.05
                ? Color(r: min(255, 60 + Int(160*frac)), g: 25, b: 8)
                : Color(r: 35, g: 45, b: 65)
            ctx.fill(Path(roundedRect: CGRect(x: rx, y: rv.minY + fy(0.10),
                                              width: fy(0.05), height: rv.height - fy(0.16)),
                          cornerRadius: 2), with: .color(rodCol))
        }
        // Reactor label on top, value below
        drawLabel(ctx, "REACTOR", rect: rv, offsetY: fy(0.015),
                  color: fuelAlarm ? Theme.alarm : Theme.textHdr, fontSize: 11)
        drawLabel(ctx, snap.powerFraction.isSafe ? String(format: "%.1f %%", snap.powerFraction*100) : "---",
                  rect: rv, offsetY: fy(0.13),
                  color: snap.powerFraction > 1.1 ? Theme.alarm : Theme.normal, fontSize: 10)
        drawLabel(ctx, snap.fuelTempK.isSafe ? String(format: "%.0f K", snap.fuelTempK) : "---",
                  rect: rv, offsetY: fy(0.23),
                  color: snap.fuelTempK > 1200 ? Theme.alarm : Theme.textDim, fontSize: 9)

        // PZR
        drawBox(ctx, pz, label: "PZR",
                color: pressAlarm ? Theme.alarm : Color(r: 80, g: 60, b: 165),
                bg: Color(r: 14, g: 12, b: 30))
        drawAnnotation(ctx, pos: CGPoint(x: pz.maxX + 6, y: pz.midY - 6),
                       text: sup.pressureMPa.isSafe ? String(format: "%.2f MPa", sup.pressureMPa) : "---",
                       color: pressAlarm ? Theme.alarm : Theme.textDim, fontSize: 9)

        // S/G with value
        drawBox(ctx, sg, label: "S/G", color: Color(r: 55, g: 95, b: 185),
                bg: Color(r: 10, g: 16, b: 34))
        drawLabel(ctx, snap.sgTempK.isSafe ? String(format: "%.0f K", snap.sgTempK) : "---",
                  rect: sg, offsetY: fy(0.18),
                  color: Theme.twophase, fontSize: 9)

        // RCP
        drawBox(ctx, rcp, label: "RCP", color: Color(r: 65, g: 105, b: 185),
                bg: Color(r: 10, g: 16, b: 34))
        drawPump(ctx, center: CGPoint(x: rcp.midX, y: rcp.midY), r: min(rcp.width, rcp.height)*0.28,
                 running: flow > 0.1, color: Color(r: 70, g: 120, b: 210))

        // TURB
        drawBox(ctx, tb, label: "TURB", color: Color(r: 55, g: 135, b: 75),
                bg: Color(r: 10, g: 22, b: 14))
        drawLabel(ctx, String(format: "%.0f %%", sup.turbineValve * 100),
                  rect: tb, offsetY: fy(0.14), color: Theme.normal, fontSize: 9)

        // COND
        drawBox(ctx, cd, label: "COND", color: Color(r: 45, g: 85, b: 115),
                bg: Color(r: 8, g: 14, b: 22))
        drawLabel(ctx, sup.condTempK.isSafe ? String(format: "%.0f K", sup.condTempK) : "---",
                  rect: cd, offsetY: fy(0.12), color: Theme.textDim, fontSize: 9)

        // FW pump
        let fwPumpX = (cd.minX + sg.maxX) / 2
        drawPump(ctx, center: CGPoint(x: fwPumpX, y: fwY), r: H * 0.055,
                 running: sup.feedwaterValve > 0.1 && !sup.feedwaterFault,
                 color: Theme.water)
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

    private func drawBox(_ ctx: GraphicsContext, _ rect: CGRect,
                         label: String, color: Color, bg: Color) {
        ctx.fill(Path(roundedRect: rect, cornerRadius: 5), with: .color(bg))
        // Top highlight
        var hl = Path()
        hl.move(to: CGPoint(x: rect.minX+5, y: rect.minY+1))
        hl.addLine(to: CGPoint(x: rect.maxX-5, y: rect.minY+1))
        ctx.stroke(hl, with: .color(.white.opacity(0.07)), lineWidth: 1)
        ctx.stroke(Path(roundedRect: rect, cornerRadius: 5), with: .color(color), lineWidth: 1.5)
        if !label.isEmpty {
            ctx.draw(Text(label).font(.system(size: 11, design: .monospaced))
                        .foregroundColor(Theme.textHdr),
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

    private func drawPump(_ ctx: GraphicsContext, center: CGPoint, r: CGFloat,
                          running: Bool, color: Color) {
        let bg = running ? Color(r: 18, g: 28, b: 48) : Color(r: 12, g: 14, b: 18)
        let rect = CGRect(x: center.x-r, y: center.y-r, width: r*2, height: r*2)
        ctx.fill(Path(ellipseIn: rect), with: .color(bg))
        ctx.stroke(Path(ellipseIn: rect), with: .color(color), lineWidth: 1.5)
        var imp = Path()
        imp.move(to: center)
        imp.addLine(to: CGPoint(x: center.x + r*0.65, y: center.y - r*0.65))
        ctx.stroke(imp, with: .color(color), lineWidth: 1.5)
        if running {
            for a in stride(from: 0.0, to: 360.0, by: 90.0) {
                let rad = CGFloat(a) * .pi / 180
                let dx = (r * 0.75) * cos(rad); let dy = -(r * 0.75) * sin(rad)
                var dp = Path()
                dp.addEllipse(in: CGRect(x: center.x+dx-2, y: center.y+dy-2, width: 4, height: 4))
                ctx.fill(dp, with: .color(color.opacity(0.6)))
            }
        }
    }
}
