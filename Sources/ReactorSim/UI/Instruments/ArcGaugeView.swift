// ArcGaugeView.swift — Circular arc gauge (270° sweep, ISA color bands).

import SwiftUI

struct ArcGaugeView: View {
    let value: Double
    let lo: Double
    let hi: Double
    let label: String
    let unit: String
    var tripHi: Double? = nil

    private var fraction: Double { max(0, min(1, (value - lo) / max(1e-9, hi - lo))) }

    private var statusColor: Color {
        if value > (hi * 0.9) { return Theme.alarm }
        if value > (hi * 0.75) { return Theme.caution }
        return .white            // normal range: neutral, not ISA green
    }

    var body: some View {
        GeometryReader { geo in
            let r = min(geo.size.width, geo.size.height) / 2 - 18   // room for tick labels
            let cx = geo.size.width / 2
            let cy = geo.size.height / 2
            Canvas { ctx, _ in
                drawGauge(ctx: ctx, cx: cx, cy: cy, r: r)
            }
        }
        .aspectRatio(1, contentMode: .fit)
    }

    private func drawGauge(ctx: GraphicsContext, cx: CGFloat, cy: CGFloat, r: CGFloat) {
        let startDeg: Double = 225
        let sweepDeg: Double = 270
        let nSegs = 60
        let rOuter = r
        let rInner = r * 0.68

        // Background arc — neutral track (no decorative color zones)
        for i in 0..<nSegs {
            let a0 = Angle.degrees(startDeg - Double(i)   * sweepDeg / Double(nSegs))
            let a1 = Angle.degrees(startDeg - Double(i+1) * sweepDeg / Double(nSegs))
            ctx.fill(arcSegPath(cx: cx, cy: cy, rO: rOuter, rI: rInner, a0: a0, a1: a1),
                     with: .color(.white.opacity(0.07)))
        }

        // Filled arc up to current value — accent normally; amber/red only in the
        // caution/trip bands themselves (live abnormal indication).
        let nFill = max(1, Int(fraction * Double(nSegs)))
        for i in 0..<nFill {
            let a0 = Angle.degrees(startDeg - Double(i)   * sweepDeg / Double(nSegs))
            let a1 = Angle.degrees(startDeg - Double(i+1) * sweepDeg / Double(nSegs))
            let sf = Double(i) / Double(nSegs)
            let fillColor: Color = sf < 0.85 ? Theme.accent : Theme.alarm
            ctx.fill(arcSegPath(cx: cx, cy: cy, rO: rOuter, rI: rInner, a0: a0, a1: a1),
                     with: .color(fillColor.opacity(0.85)))
        }

        // Calibrated scale: 5 major numbered ticks + minors (a gauge without a
        // scale is a battery indicator, not an instrument)
        let tickFmt: String = {
            let step = (hi - lo) / 4
            return abs(step - step.rounded()) < 0.01 ? "%.0f" : "%.1f"
        }()
        for i in 0...4 {
            let f = Double(i) / 4
            let a = Angle.degrees(startDeg - f * sweepDeg)
            let x0 = cx + (rOuter + 2) * cos(a.radians)
            let y0 = cy - (rOuter + 2) * sin(a.radians)
            let x1 = cx + (rOuter + 7) * cos(a.radians)
            let y1 = cy - (rOuter + 7) * sin(a.radians)
            var p = Path(); p.move(to: .init(x: x0, y: y0)); p.addLine(to: .init(x: x1, y: y1))
            ctx.stroke(p, with: .color(.white.opacity(0.35)), lineWidth: 1)
            let lx = cx + (rOuter + 14) * cos(a.radians)
            let ly = cy - (rOuter + 14) * sin(a.radians)
            ctx.draw(
                Text(String(format: tickFmt, lo + (hi - lo) * f))
                    .font(.system(size: 8, design: .monospaced))
                    .foregroundColor(Theme.textDim),
                at: .init(x: lx, y: ly), anchor: .center)
        }
        for i in 0..<20 where i % 5 != 0 {
            let a = Angle.degrees(startDeg - Double(i) / 20 * sweepDeg)
            let x0 = cx + (rOuter + 2) * cos(a.radians)
            let y0 = cy - (rOuter + 2) * sin(a.radians)
            let x1 = cx + (rOuter + 5) * cos(a.radians)
            let y1 = cy - (rOuter + 5) * sin(a.radians)
            var p = Path(); p.move(to: .init(x: x0, y: y0)); p.addLine(to: .init(x: x1, y: y1))
            ctx.stroke(p, with: .color(.white.opacity(0.15)), lineWidth: 0.5)
        }

        // Trip line
        if let trip = tripHi {
            let tf = max(0, min(1, (trip - lo) / max(1e-9, hi - lo)))
            let ta = Angle.degrees(startDeg - tf * sweepDeg)
            let x0 = cx + rInner * cos(ta.radians)
            let y0 = cy - rInner * sin(ta.radians)
            let x1 = cx + rOuter * cos(ta.radians)
            let y1 = cy - rOuter * sin(ta.radians)
            var tp = Path(); tp.move(to: CGPoint(x: x0, y: y0)); tp.addLine(to: CGPoint(x: x1, y: y1))
            ctx.stroke(tp, with: .color(Theme.alarm), lineWidth: 2)
        }

        // Needle
        let needleA = Angle.degrees(startDeg - fraction * sweepDeg)
        let rTip = r * 0.62
        var needle = Path()
        needle.move(to: CGPoint(x: cx, y: cy))
        needle.addLine(to: CGPoint(x: cx + rTip * cos(needleA.radians),
                                   y: cy - rTip * sin(needleA.radians)))
        ctx.stroke(needle, with: .color(.white.opacity(0.9)), lineWidth: 2)
        ctx.fill(Path(ellipseIn: CGRect(x: cx-5, y: cy-5, width: 10, height: 10)),
                 with: .color(Color(r: 40, g: 50, b: 72)))
        ctx.fill(Path(ellipseIn: CGRect(x: cx-3, y: cy-3, width: 6, height: 6)),
                 with: .color(Color(r: 170, g: 185, b: 215)))

        // Digital readout
        let valStr = abs(value) < 100 ? String(format: "%.1f", value) : String(format: "%.0f", value)
        ctx.draw(
            Text(valStr).font(Theme.readoutSm).foregroundColor(statusColor),
            at: CGPoint(x: cx, y: cy + r * 0.22), anchor: .top
        )
        if !unit.isEmpty {
            ctx.draw(
                Text(unit).font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                at: CGPoint(x: cx, y: cy + r * 0.22 + 14), anchor: .top
            )
        }
        if !label.isEmpty {
            ctx.draw(
                Text(label).font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                at: CGPoint(x: cx, y: cy - r - 2), anchor: .bottom
            )
        }
    }

    private func arcSegPath(cx: CGFloat, cy: CGFloat,
                            rO: CGFloat, rI: CGFloat,
                            a0: Angle, a1: Angle) -> Path {
        var p = Path()
        p.move(to: CGPoint(x: cx + rO * cos(a0.radians), y: cy - rO * sin(a0.radians)))
        p.addLine(to: CGPoint(x: cx + rO * cos(a1.radians), y: cy - rO * sin(a1.radians)))
        p.addLine(to: CGPoint(x: cx + rI * cos(a1.radians), y: cy - rI * sin(a1.radians)))
        p.addLine(to: CGPoint(x: cx + rI * cos(a0.radians), y: cy - rI * sin(a0.radians)))
        p.closeSubpath()
        return p
    }
}
