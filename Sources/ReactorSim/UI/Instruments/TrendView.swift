// TrendView.swift — industrial strip chart: framed plot area, calibrated
// grid with axis labels on the gridlines they describe, no decorative glow.

import SwiftUI

struct TrendView: View {
    let values: [Double]
    let yLo: Double
    let yHi: Double
    let color: Color
    let label: String

    // Axis format derived from the SCALE so widths never jitter
    private func axisStr(_ v: Double) -> String {
        let span = abs(yHi - yLo)
        if span >= 50  { return String(format: "%.0f", v) }
        if span >= 1   { return String(format: "%.1f", v) }
        return String(format: "%+.4f", v)
    }

    var body: some View {
        // Authentic draws a brighter, more present graticule like a real trend recorder.
        let flat   = Theme.isFlat
        let frameO = flat ? 0.32 : 0.13
        let gridO  = flat ? 0.11 : 0.05
        let vgridO = flat ? 0.075 : 0.035
        let tickO  = flat ? 0.42 : 0.25
        Canvas { ctx, size in
            let W = size.width; let H = size.height
            let compact = H < 70          // mini-trend: only label top & bottom
            let axisW: CGFloat = 38
            let plot = CGRect(x: 0.5, y: 0.5, width: W - axisW, height: H - 1)

            // Plot area + frame
            ctx.fill(Path(plot), with: .color(Theme.ink.opacity(flat ? 0.02 : 0.03)))
            ctx.stroke(Path(plot), with: .color(Theme.ink.opacity(frameO)), lineWidth: 1)

            // Horizontal grid: 4 divisions, each labeled with ITS OWN value
            for i in 0...4 {
                let fy = CGFloat(i) / 4
                let y = plot.minY + plot.height * fy
                if i > 0 && i < 4 {
                    var gp = Path()
                    gp.move(to: .init(x: plot.minX + 1, y: y))
                    gp.addLine(to: .init(x: plot.maxX - 1, y: y))
                    ctx.stroke(gp, with: .color(Theme.ink.opacity(gridO)), lineWidth: 0.5)
                }
                var tick = Path()
                tick.move(to: .init(x: plot.maxX, y: y))
                tick.addLine(to: .init(x: plot.maxX + 3, y: y))
                ctx.stroke(tick, with: .color(Theme.ink.opacity(tickO)), lineWidth: 1)

                if compact && i != 0 && i != 4 { continue }   // skip middle labels when tiny
                let v = yHi - (yHi - yLo) * Double(fy)
                let ly = min(max(y, 5), H - 5)
                ctx.draw(
                    Text(axisStr(v))
                        .font(.system(size: 8, design: .monospaced))
                        .foregroundColor(Theme.textDim),
                    at: .init(x: plot.maxX + 5, y: ly), anchor: .leading)
            }
            // Vertical grid: 6 divisions
            for i in 1..<6 {
                let x = plot.minX + plot.width * CGFloat(i) / 6
                var gp = Path()
                gp.move(to: .init(x: x, y: plot.minY + 1))
                gp.addLine(to: .init(x: x, y: plot.maxY - 1))
                ctx.stroke(gp, with: .color(Theme.ink.opacity(vgridO)), lineWidth: 0.5)
            }

            // Channel label (top-left) + current value (top-right, inside plot)
            ctx.draw(
                Text(label).font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                at: CGPoint(x: plot.minX + 5, y: plot.minY + 4), anchor: .topLeading)
            if let current = values.last {
                ctx.draw(
                    Text(axisStr(current))
                        .font(.system(size: 9, weight: .semibold, design: .monospaced))
                        .foregroundColor(Theme.ink),
                    at: CGPoint(x: plot.maxX - 5, y: plot.minY + 4), anchor: .topTrailing)
            }

            // Series — clipped to the plot, no end-dot decoration
            guard values.count >= 2 else { return }
            let n = values.count
            var lp = Path()
            for (i, v) in values.enumerated() {
                let frac = max(0, min(1, (v - yLo) / max(1e-9, yHi - yLo)))
                let pt = CGPoint(
                    x: plot.minX + 2 + CGFloat(i) * (plot.width - 4) / CGFloat(n - 1),
                    y: plot.maxY - 2 - CGFloat(frac) * (plot.height - 4))
                if i == 0 { lp.move(to: pt) } else { lp.addLine(to: pt) }
            }
            ctx.clip(to: Path(plot))
            ctx.stroke(lp, with: .color(color), lineWidth: 1.2)

            // Current-value pointer on the right frame edge
            if let current = values.last {
                let frac = max(0, min(1, (current - yLo) / max(1e-9, yHi - yLo)))
                let y = plot.maxY - 2 - CGFloat(frac) * (plot.height - 4)
                var marker = Path()
                marker.move(to: .init(x: plot.maxX - 6, y: y))
                marker.addLine(to: .init(x: plot.maxX - 1, y: y - 3))
                marker.addLine(to: .init(x: plot.maxX - 1, y: y + 3))
                marker.closeSubpath()
                ctx.fill(marker, with: .color(color))
            }
        }
    }
}
