// ReactivityTab.swift — F4 Reactivity: budget bars, xenon trend, rod worth curve.

import SwiftUI

struct ReactivityTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 8) {
            VStack(spacing: 8) {
                ReactivityBudgetPanel(supervisor: supervisor)
                XenonTrendPanel(supervisor: supervisor)
            }
            .frame(maxWidth: .infinity)

            VStack(spacing: 8) {
                RodWorthCurvePanel(supervisor: supervisor)
                ReactivityReadoutsPanel(supervisor: supervisor)
            }
            .frame(width: 300)
        }
        .padding(8)
    }
}

// MARK: — Reactivity budget
private struct ReactivityBudgetPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        let p = PlantParams()

        // Approximate component breakdown
        let rodPos  = s.rodPosition
        let w       = 3*rodPos*rodPos - 2*rodPos*rodPos*rodPos
        let rhoRods = w * p.rodWorth + (s.scrammed ? p.scramExtraWorth : 0)
        let rhoFuel = p.fuelTempCoeff * (s.fuelTempK - p.nominalFuelTemp)
        let rhoCool = p.coolantTempCoeff * (s.coolantTempK - p.nominalCoolantTemp)
        let rhoXe   = -p.xenonReactivityCoeff * s.xenonInventory
        let rhoB    = -8e-5 * (supervisor.boronPPM - 800.0)
        let rhoTot  = s.reactivity

        let components: [(String, Double, Color)] = [
            ("Total",          rhoTot,  Theme.accent),
            ("Rod position",   rhoRods, Theme.twophase),
            ("Doppler (fuel)", rhoFuel, Theme.hotLeg),
            ("Moderator",      rhoCool, Theme.water),
            ("Xenon-135",      rhoXe,   Theme.caution),
            ("Boron",          rhoB,    Color(r: 130, g: 80, b: 200)),
        ]

        VStack(spacing: 0) {
            PanelHeader(title: "REACTIVITY BUDGET")
            Canvas { ctx, size in
                let W = size.width; let H = size.height
                let maxAbs = components.map { abs($0.1) }.max() ?? 0.05
                let scale  = (W / 2 - 160) / max(maxAbs, 0.001)
                let zeroX  = W / 2
                let rowH   = (H - 20) / CGFloat(components.count)

                // Zero line
                var zl = Path(); zl.move(to: .init(x: zeroX, y: 4)); zl.addLine(to: .init(x: zeroX, y: H - 4))
                ctx.stroke(zl, with: .color(Theme.textDim.opacity(0.4)), lineWidth: 1)

                for (i, (name, val, color)) in components.enumerated() {
                    let y = 10 + CGFloat(i) * rowH
                    // Label
                    ctx.draw(Text(name).font(Theme.readoutSm).foregroundColor(Theme.textDim),
                             at: .init(x: zeroX - 8, y: y + rowH/2), anchor: .trailing)
                    // Value label
                    ctx.draw(Text(String(format: "%+.5f", val)).font(Theme.readoutSm).foregroundColor(color),
                             at: .init(x: zeroX + 8, y: y + rowH/2), anchor: .leading)
                    // Bar
                    let bLen   = CGFloat(abs(val)) * scale
                    let labelW: CGFloat = 85
                    let barX = val >= 0 ? zeroX + labelW : zeroX - labelW - bLen
                    ctx.fill(
                        Path(roundedRect: CGRect(x: barX, y: y + 4, width: bLen, height: rowH - 8),
                             cornerRadius: 2),
                        with: .color(color.opacity(0.7))
                    )
                }
            }
            .padding(8)
            .frame(maxHeight: .infinity)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxWidth: .infinity)
    }
}

// MARK: — Xenon trend
private struct XenonTrendPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "XENON / IODINE INVENTORY TREND")
            HStack(spacing: 8) {
                TrendView(values: supervisor.orderedHistory(supervisor.histPower),
                          yLo: 0, yHi: 130, color: Theme.normal, label: "Core power (%)")
                TrendView(values: supervisor.orderedHistory(supervisor.histReact),
                          yLo: -0.05, yHi: 0.005, color: Theme.caution, label: "Reactivity (Δk/k)")
                TrendView(values: supervisor.orderedHistory(supervisor.histCoolT),
                          yLo: 400, yHi: 650, color: Theme.water, label: "Coolant temp (K)")
            }
            .padding(8)
            .frame(maxHeight: .infinity)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxWidth: .infinity)
    }
}

// MARK: — Rod worth S-curve
private struct RodWorthCurvePanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let rodPos = supervisor.snapshot.rodPosition
        VStack(spacing: 0) {
            PanelHeader(title: "ROD WORTH CURVE  (S-CURVE)")
            Canvas { ctx, size in
                let W = size.width; let H = size.height
                let pad: CGFloat = 20

                // Axes
                var axes = Path()
                axes.move(to: .init(x: pad, y: pad)); axes.addLine(to: .init(x: pad, y: H - pad))
                axes.move(to: .init(x: pad, y: H - pad)); axes.addLine(to: .init(x: W - pad, y: H - pad))
                ctx.stroke(axes, with: .color(Theme.border), lineWidth: 1)

                // S-curve: w(x) = 3x²-2x³, worth = w * rod_worth
                var curve = Path()
                for i in 0...100 {
                    let x = Double(i) / 100.0
                    let wx = 3*x*x - 2*x*x*x
                    let rho = wx * PlantParams().rodWorth
                    let px = pad + CGFloat(x)    * (W - 2*pad)
                    let py = H - pad - CGFloat((rho - PlantParams().rodWorth) / (-PlantParams().rodWorth)) * (H - 2*pad)
                    if i == 0 { curve.move(to: .init(x: px, y: py)) }
                    else { curve.addLine(to: .init(x: px, y: py)) }
                }
                ctx.stroke(curve, with: .color(Theme.accent), lineWidth: 2)

                // Current rod position marker
                let wx = 3*rodPos*rodPos - 2*rodPos*rodPos*rodPos
                let rho = wx * PlantParams().rodWorth
                let mpx = pad + CGFloat(rodPos) * (W - 2*pad)
                let mpy = H - pad - CGFloat((rho - PlantParams().rodWorth) / (-PlantParams().rodWorth)) * (H - 2*pad)
                ctx.fill(Path(ellipseIn: .init(x: mpx-5, y: mpy-5, width: 10, height: 10)),
                         with: .color(Theme.alarm))
                ctx.draw(
                    Text(String(format: "%.0f%%", rodPos * 100)).font(Theme.readoutSm).foregroundColor(Theme.alarm),
                    at: .init(x: mpx + 8, y: mpy), anchor: .leading
                )

                // Axis labels
                ctx.draw(Text("0%").font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: .init(x: pad, y: H - pad + 4), anchor: .top)
                ctx.draw(Text("100%").font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: .init(x: W - pad, y: H - pad + 4), anchor: .top)
            }
            .padding(4)
            .frame(maxHeight: .infinity)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
    }
}

// MARK: — Reactivity readouts
private struct ReactivityReadoutsPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        VStack(spacing: 0) {
            PanelHeader(title: "NEUTRON PHYSICS")
            VStack(alignment: .leading, spacing: 6) {
                readRow("Total reactivity",  String(format: "%+.5f Δk/k", s.reactivity), .reactivityStatus(s.reactivity))
                readRow("Neutron pop.",      String(format: "%.4f (n/n₀)", s.powerFraction), .powerStatus(s.powerFraction))
                readRow("Xenon inventory",   String(format: "%.4f (rel.)", s.xenonInventory), Theme.caution)
                readRow("Iodine inventory",  String(format: "%.4f (rel.)", s.iodineInventory), Theme.textDim)
                readRow("Decay heat",        String(format: "%.3f %%", s.decayHeatFraction * 100),
                        s.decayHeatFraction > 0.03 ? Theme.caution : Theme.text)
                readRow("Rod position",      String(format: "%.1f %% ins.", s.rodPosition * 100), Theme.text)
                readRow("Boron",             String(format: "%.1f ppm", supervisor.boronPPM),
                        supervisor.boronPPM < 100 ? Theme.caution : Theme.text)
                readRow("Sim time",          String(format: "%.1f s", s.time), Theme.textDim)
            }
            .padding(10)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxHeight: .infinity)
    }
    private func readRow(_ l: String, _ v: String, _ c: Color) -> some View {
        HStack {
            Text(l).font(Theme.readoutSm).foregroundStyle(Theme.textDim)
            Spacer()
            Text(v).font(Theme.readoutSm).foregroundStyle(c)
        }
    }
}
