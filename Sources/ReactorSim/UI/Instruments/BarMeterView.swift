// BarMeterView.swift — Vertical bar meter with a calibrated tick scale and
// setpoint markers. ISA colors appear only when the value is in an abnormal
// band; the normal-range fill is the accent.

import SwiftUI

struct BarMeterView: View {
    let value: Double
    let lo: Double
    let hi: Double
    let label: String
    let unit: String
    var tripHi: Double? = nil
    var tripLo: Double? = nil
    var warnHi: Double? = nil
    var warnLo: Double? = nil

    private var fraction: Double { max(0, min(1, (value - lo) / max(1e-9, hi - lo))) }

    // Color comes from REAL setpoints only — never from "% of scale", which
    // painted normal operating points amber (PZR at 15.6/18, RCP at 100/110).
    private var fillColor: Color {
        if let t = tripHi, value > t { return Theme.alarm }
        if let t = tripLo, value < t { return Theme.alarm }
        if let w = warnHi, value > w { return Theme.caution }
        if let w = warnLo, value < w { return Theme.caution }
        return Theme.accent
    }

    // Fixed decimals chosen from the SCALE, not the value — readouts never
    // jitter between widths.
    private var decimals: Int {
        let span = hi - lo
        if span >= 50 { return 1 }
        if span >= 5  { return 2 }
        return 3
    }

    private var tickFormat: String {
        let step = (hi - lo) / 4
        return abs(step - step.rounded()) < 0.01 ? "%.0f" : "%.1f"
    }

    var body: some View {
        VStack(spacing: 5) {
            Text(label)
                .font(.system(size: 9, weight: .medium, design: .monospaced))
                .foregroundStyle(Theme.textDim)
                .tracking(0.5)
                .lineLimit(1)
                .minimumScaleFactor(0.6)

            let flat = Theme.isFlat
            Canvas { ctx, size in
                let w = size.width
                let barW: CGFloat = min(16, max(10, w * 0.30))
                let barX = w - barW - 6              // room for setpoint flags
                let yTop: CGFloat = 4
                let yBot = size.height - 4
                let span = yBot - yTop
                let trackR: CGFloat = flat ? 0 : 5   // authentic: square sight-glass

                func yFor(_ v: Double) -> CGFloat {
                    let f = max(0, min(1, (v - lo) / max(1e-9, hi - lo)))
                    return yBot - span * CGFloat(f)
                }

                // Track
                ctx.fill(
                    Path(roundedRect: CGRect(x: barX, y: yTop, width: barW, height: span),
                         cornerRadius: trackR, style: .continuous),
                    with: .color(.white.opacity(flat ? 0.09 : 0.06)))
                if flat {
                    ctx.stroke(Path(CGRect(x: barX, y: yTop, width: barW, height: span)),
                               with: .color(.white.opacity(0.20)), lineWidth: 1)
                }

                // Calibrated scale: 5 major ticks with values, minor ticks between
                for i in 0...4 {
                    let v = lo + (hi - lo) * Double(i) / 4
                    let y = yFor(v)
                    var p = Path()
                    p.move(to: .init(x: barX - 8, y: y)); p.addLine(to: .init(x: barX - 2, y: y))
                    ctx.stroke(p, with: .color(.white.opacity(flat ? 0.5 : 0.35)), lineWidth: 1)
                    ctx.draw(
                        Text(String(format: tickFormat, v))
                            .font(.system(size: 8, design: .monospaced))
                            .foregroundColor(Theme.textDim),
                        at: .init(x: barX - 11, y: y), anchor: .trailing)
                }
                for i in 0..<20 where i % 5 != 0 {
                    let y = yFor(lo + (hi - lo) * Double(i) / 20)
                    var p = Path()
                    p.move(to: .init(x: barX - 5, y: y)); p.addLine(to: .init(x: barX - 2, y: y))
                    ctx.stroke(p, with: .color(.white.opacity(flat ? 0.25 : 0.15)), lineWidth: 0.5)
                }

                // Fill bar
                let fh = span * CGFloat(fraction)
                if fh > 0.5 {
                    ctx.fill(
                        Path(roundedRect: CGRect(x: barX + 2, y: yBot - fh, width: barW - 4, height: fh),
                             cornerRadius: flat ? 0 : 3, style: .continuous),
                        with: .color(fillColor.opacity(flat ? 0.95 : 0.85)))
                }

                // Setpoint markers: line across the bar + pointer flag at right
                func setpoint(_ v: Double, _ color: Color) {
                    guard v > lo, v < hi else { return }
                    let y = yFor(v)
                    var p = Path()
                    p.move(to: .init(x: barX, y: y)); p.addLine(to: .init(x: barX + barW, y: y))
                    ctx.stroke(p, with: .color(color), lineWidth: 1.2)
                    var tri = Path()
                    tri.move(to: .init(x: barX + barW + 6, y: y - 3.5))
                    tri.addLine(to: .init(x: barX + barW + 6, y: y + 3.5))
                    tri.addLine(to: .init(x: barX + barW + 1, y: y))
                    tri.closeSubpath()
                    ctx.fill(tri, with: .color(color))
                }
                if let w = warnHi { setpoint(w, Theme.caution) }
                if let w = warnLo { setpoint(w, Theme.caution) }
                if let t = tripHi { setpoint(t, Theme.alarm) }
                if let t = tripLo { setpoint(t, Theme.alarm) }
            }
            .frame(maxWidth: .infinity)

            // Numeric readout — fixed decimals from scale span
            Text(String(format: "%.\(decimals)f", value))
                .font(.system(size: 11, weight: .semibold, design: .monospaced))
                .foregroundStyle(fillColor == Theme.accent ? .white : fillColor)
            Text(unit.isEmpty ? "—" : unit)
                .font(.system(size: 8, design: .monospaced))
                .foregroundStyle(Theme.textDim)
        }
    }
}
