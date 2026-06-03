// TrendView.swift — Strip-chart trend with grid lines.

import SwiftUI

struct TrendView: View {
    let values: [Double]
    let yLo: Double
    let yHi: Double
    let color: Color
    let label: String

    var body: some View {
        Canvas { ctx, size in
            let W = size.width; let H = size.height

            // Background
            ctx.fill(Path(roundedRect: CGRect(origin: .zero, size: size), cornerRadius: 3),
                     with: .color(Color(r: 8, g: 10, b: 12)))
            ctx.stroke(Path(roundedRect: CGRect(origin: .zero, size: size), cornerRadius: 3),
                       with: .color(Theme.border), lineWidth: 1)

            // Grid
            for i in 1..<4 {
                let gy = H * CGFloat(i) / 4
                var gp = Path(); gp.move(to: CGPoint(x: 2, y: gy)); gp.addLine(to: CGPoint(x: W-2, y: gy))
                ctx.stroke(gp, with: .color(Color(white: 1, opacity: 0.04)), lineWidth: 0.5)
            }
            for i in 1..<4 {
                let gx = W * CGFloat(i) / 4
                var gp = Path(); gp.move(to: CGPoint(x: gx, y: 2)); gp.addLine(to: CGPoint(x: gx, y: H-2))
                ctx.stroke(gp, with: .color(Color(white: 1, opacity: 0.03)), lineWidth: 0.5)
            }

            // Label
            ctx.draw(
                Text(label).font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                at: CGPoint(x: 5, y: 4), anchor: .topLeading
            )

            // Trend line
            guard values.count >= 2 else { return }
            let n = values.count
            var pts: [CGPoint] = []
            for (i, v) in values.enumerated() {
                let frac = max(0, min(1, (v - yLo) / max(1e-9, yHi - yLo)))
                pts.append(CGPoint(x: 2 + CGFloat(i) * (W-4) / CGFloat(n-1),
                                   y: H - 2 - CGFloat(frac) * (H-4)))
            }
            var lp = Path(); lp.move(to: pts[0])
            for pt in pts.dropFirst() { lp.addLine(to: pt) }
            ctx.stroke(lp, with: .color(color), lineWidth: 1.5)

            // Latest value dot
            if let last = pts.last {
                ctx.fill(Path(ellipseIn: CGRect(x: last.x-3, y: last.y-3, width: 6, height: 6)),
                         with: .color(color))
            }
        }
    }
}
