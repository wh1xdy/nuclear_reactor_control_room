// StartupPanel.swift — startup-physics popup (tap the NEUTRONICS dock).
// Three instruments for the approach to critical:
//  • ECP — estimated critical position from the live reactivity balance
//    (boron + xenon + moderator + Doppler inverted through the rod S-curve).
//  • 1/M plot — inverse multiplication vs rod withdrawal, with the straight-
//    line extrapolation to zero that predicts criticality.
//  • Rod worth — integral and differential curves with the current position.

import SwiftUI

struct StartupPanel: View {
    let supervisor: PlantSupervisor

    var body: some View {
        ZStack {
            Rectangle().fill(.black.opacity(0.55)).ignoresSafeArea()
                .contentShape(Rectangle())
                .onTapGesture { supervisor.startupPanelOpen = false }
            HStack(alignment: .top, spacing: 16) {
                ecpBlock.frame(width: 190)
                VStack(alignment: .leading, spacing: 6) {
                    Text("1/M · APPROACH TO CRITICAL")
                        .font(.system(size: 9, weight: .semibold, design: .monospaced))
                        .foregroundStyle(Theme.textHdr)
                    invMPlot.frame(width: 250, height: 210)
                    Text("plot clears on scram · points log as rods withdraw below 10% power")
                        .font(.system(size: 7, design: .monospaced)).foregroundStyle(Theme.textDim)
                }
                VStack(alignment: .leading, spacing: 6) {
                    HStack {
                        Text("ROD WORTH").font(.system(size: 9, weight: .semibold, design: .monospaced))
                            .foregroundStyle(Theme.textHdr)
                        Spacer()
                        Button { supervisor.startupPanelOpen = false } label: {
                            Text("CLOSE  esc").font(.system(size: 9, design: .monospaced))
                                .foregroundStyle(Theme.textDim)
                                .padding(.horizontal, 8).padding(.vertical, 4)
                                .contentShape(Rectangle())
                        }
                        .buttonStyle(.plain).controlSurface()
                    }
                    worthPlot.frame(width: 250, height: 210)
                    Text("integral (solid) · differential (dashed) · ● current position")
                        .font(.system(size: 7, design: .monospaced)).foregroundStyle(Theme.textDim)
                }
            }
            .padding(20)
            .background(RoundedRectangle(cornerRadius: 14).fill(Theme.panel))
            .overlay(RoundedRectangle(cornerRadius: 14).stroke(Theme.border, lineWidth: 1))
        }
    }

    // ── ECP ──────────────────────────────────────────────────────────────────

    private var ecpBlock: some View {
        let c = supervisor.rhoComponents
        let ecp = supervisor.ecpSWD
        return VStack(alignment: .leading, spacing: 7) {
            Text("ECP · REACTIVITY BALANCE")
                .font(.system(size: 9, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.textHdr)
            row("BORON",  c.boron)
            row("XENON",  c.xenon)
            row("MOD-T",  c.mod)
            row("DOPPLER", c.dop)
            Divider().background(Theme.sep)
            HStack {
                Text("ROD POS").font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                Spacer()
                Text("\(Int((228 * (1 - supervisor.snapshot.rodPosition)).rounded())) SWD")
                    .font(.system(size: 11, weight: .semibold, design: .monospaced)).foregroundStyle(Theme.ink)
            }
            HStack {
                Text("ECP").font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                Spacer()
                Text(ecp.map { "\($0) SWD" } ?? "PRECLUDED — XENON")
                    .font(.system(size: 12, weight: .bold, design: .monospaced))
                    .foregroundStyle(ecp == nil ? Theme.caution : (Theme.isFlat ? Theme.ink : Theme.statusNormal))
            }
            Text(ecp == nil ? "no critical position exists at the present balance — wait out the xenon or dilute"
                            : "withdraw toward the ECP on the 1/M plot, not past it")
                .font(.system(size: 7, design: .monospaced)).foregroundStyle(Theme.textDim)
        }
    }

    private func row(_ l: String, _ rho: Double) -> some View {
        HStack {
            Text(l).font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
            Spacer()
            Text(String(format: "%+.0f pcm", rho * 1e5))
                .font(.system(size: 10, weight: .semibold, design: .monospaced))
                .foregroundStyle(abs(rho) * 1e5 > 800 ? Theme.caution : Theme.ink)
        }
    }

    // ── 1/M ──────────────────────────────────────────────────────────────────

    private var invMPlot: some View {
        Canvas { ctx, size in
            let plot = CGRect(x: 26, y: 6, width: size.width - 32, height: size.height - 24)
            ctx.fill(Path(plot), with: .color(Theme.dockTint))
            ctx.stroke(Path(plot), with: .color(Theme.ink.opacity(0.18)), lineWidth: 1)
            func px(_ x: Double) -> CGFloat { plot.minX + plot.width * CGFloat(max(0, min(1, x))) }
            func py(_ m: Double) -> CGFloat { plot.maxY - plot.height * CGFloat(max(0, min(1.2, m)) / 1.2) }
            for v in [0.0, 0.5, 1.0] {
                ctx.stroke(Path { p in p.move(to: CGPoint(x: plot.minX, y: py(v))); p.addLine(to: CGPoint(x: plot.maxX, y: py(v))) },
                           with: .color(Theme.ink.opacity(v == 0 ? 0.3 : 0.08)), lineWidth: 0.75)
                ctx.draw(Text(String(format: "%.1f", v)).font(.system(size: 7, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: CGPoint(x: plot.minX - 3, y: py(v)), anchor: .trailing)
            }
            let pts = supervisor.invMPoints
            for p in pts {
                let c = CGPoint(x: px(p.x), y: py(p.invM))
                ctx.fill(Path(ellipseIn: CGRect(x: c.x - 2, y: c.y - 2, width: 4, height: 4)),
                         with: .color(Theme.accent))
            }
            // Straight-line extrapolation of the last few points to 1/M = 0.
            if pts.count >= 3 {
                let tail = pts.suffix(4)
                let n = Double(tail.count)
                let sx = tail.reduce(0) { $0 + $1.x }, sy = tail.reduce(0) { $0 + $1.invM }
                let sxx = tail.reduce(0) { $0 + $1.x * $1.x }, sxy = tail.reduce(0) { $0 + $1.x * $1.invM }
                let denom = n * sxx - sx * sx
                if abs(denom) > 1e-9 {
                    let slope = (n * sxy - sx * sy) / denom
                    let icept = (sy - slope * sx) / n
                    if slope < -1e-6 {
                        let xCrit = -icept / slope
                        let a = CGPoint(x: px(tail.first!.x), y: py(slope * tail.first!.x + icept))
                        let b = CGPoint(x: px(xCrit), y: py(0))
                        ctx.stroke(Path { p in p.move(to: a); p.addLine(to: b) },
                                   with: .color(Theme.caution.opacity(0.8)),
                                   style: StrokeStyle(lineWidth: 1, dash: [4, 3]))
                        ctx.draw(Text("crit ≈ \(Int((228 * (1 - min(1, max(0, xCrit)))).rounded())) SWD")
                                    .font(.system(size: 8, weight: .semibold, design: .monospaced))
                                    .foregroundColor(Theme.caution),
                                 at: CGPoint(x: min(b.x, plot.maxX - 4), y: py(0) - 8), anchor: .trailing)
                    }
                }
            }
            if pts.isEmpty {
                ctx.draw(Text("no approach in progress").font(.system(size: 8, design: .monospaced))
                            .foregroundColor(Theme.textDim),
                         at: CGPoint(x: plot.midX, y: plot.midY), anchor: .center)
            }
            ctx.draw(Text("withdrawn →").font(.system(size: 7, design: .monospaced)).foregroundColor(Theme.textDim),
                     at: CGPoint(x: plot.maxX, y: plot.maxY + 8), anchor: .trailing)
        }
    }

    // ── Rod worth ────────────────────────────────────────────────────────────

    private var worthPlot: some View {
        Canvas { ctx, size in
            let plot = CGRect(x: 26, y: 6, width: size.width - 32, height: size.height - 24)
            ctx.fill(Path(plot), with: .color(Theme.dockTint))
            ctx.stroke(Path(plot), with: .color(Theme.ink.opacity(0.18)), lineWidth: 1)
            let worth = 5000.0                                     // |rodWorth| in pcm
            func px(_ x: Double) -> CGFloat { plot.minX + plot.width * CGFloat(x) }
            func py(_ pcm: Double) -> CGFloat { plot.maxY - plot.height * CGFloat(pcm / worth) }
            // Integral worth W(x) = 3x²−2x³ (×|rodWorth|).
            var integral = Path()
            var differential = Path()
            for i in 0...60 {
                let x = Double(i) / 60
                let w = (3 * x * x - 2 * x * x * x) * worth
                let d = 6 * x * (1 - x) * worth / 4                // scaled to share the axis
                let pi = CGPoint(x: px(x), y: py(w))
                let pd = CGPoint(x: px(x), y: py(d))
                if i == 0 { integral.move(to: pi); differential.move(to: pd) }
                else { integral.addLine(to: pi); differential.addLine(to: pd) }
            }
            ctx.stroke(integral, with: .color(Theme.accent), lineWidth: 1.5)
            ctx.stroke(differential, with: .color(Theme.accent.opacity(0.6)),
                       style: StrokeStyle(lineWidth: 1, dash: [4, 3]))
            let xNow = supervisor.snapshot.rodPosition
            let wNow = (3 * xNow * xNow - 2 * xNow * xNow * xNow) * worth
            let dot = CGPoint(x: px(xNow), y: py(wNow))
            ctx.fill(Path(ellipseIn: CGRect(x: dot.x - 3, y: dot.y - 3, width: 6, height: 6)),
                     with: .color(Theme.isFlat ? Theme.ink : Theme.statusNormal))
            for v in [0.0, 2500.0, 5000.0] {
                ctx.draw(Text(String(format: "%.0f", v)).font(.system(size: 7, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: CGPoint(x: plot.minX - 3, y: py(v)), anchor: .trailing)
            }
            ctx.draw(Text("inserted →").font(.system(size: 7, design: .monospaced)).foregroundColor(Theme.textDim),
                     at: CGPoint(x: plot.maxX, y: plot.maxY + 8), anchor: .trailing)
        }
    }
}
