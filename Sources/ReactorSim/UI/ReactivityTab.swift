// ReactivityTab.swift — F4 Reactivity: budget bars, xenon trend, rod worth curve.

import SwiftUI

struct ReactivityTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 12) {
            VStack(spacing: 12) {
                ReactivityBudgetPanel(supervisor: supervisor)
                XenonTrendPanel(supervisor: supervisor)
            }
            .frame(maxWidth: .infinity)

            VStack(spacing: 12) {
                RodWorthCurvePanel(supervisor: supervisor)
                ReactivityReadoutsPanel(supervisor: supervisor)
            }
            .frame(width: 300)
        }
        .padding(12)
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

        // Restrained palette: total in accent, components in neutral grey.
        let grey = Theme.ink.opacity(0.55)
        let components: [(String, Double, Color)] = [
            ("Total",          rhoTot,  Theme.accent),
            ("Rod position",   rhoRods, grey),
            ("Doppler (fuel)", rhoFuel, grey),
            ("Moderator",      rhoCool, grey),
            ("Xenon-135",      rhoXe,   grey),
            ("Boron",          rhoB,    grey),
        ]

        VStack(spacing: 0) {
            PanelHeader(title: "REACTIVITY BUDGET")
            // Fixed columns: label | bar field (zero-centered) | value.
            // Bars stay inside the field — they can never run under the text.
            Canvas { ctx, size in
                let W = size.width; let H = size.height
                let labelW: CGFloat = 118
                let valueW: CGFloat = 92
                let fieldX0 = labelW + 12
                let fieldX1 = W - valueW - 12
                let zeroX   = (fieldX0 + fieldX1) / 2
                let halfW   = (fieldX1 - fieldX0) / 2 - 4
                let maxAbs  = max(components.map { abs($0.1) }.max() ?? 0.001, 1e-4)
                let scale   = halfW / maxAbs
                let rowH    = (H - 16) / CGFloat(components.count)

                // Field frame + zero line
                ctx.stroke(Path(CGRect(x: fieldX0, y: 6, width: fieldX1 - fieldX0, height: H - 12)),
                           with: .color(Theme.ink.opacity(0.08)), lineWidth: 1)
                var zl = Path()
                zl.move(to: .init(x: zeroX, y: 6)); zl.addLine(to: .init(x: zeroX, y: H - 6))
                ctx.stroke(zl, with: .color(Theme.ink.opacity(0.25)), lineWidth: 1)

                for (i, (name, val, color)) in components.enumerated() {
                    let yMid = 8 + CGFloat(i) * rowH + rowH / 2
                    ctx.draw(Text(name).font(Theme.readoutSm).foregroundColor(Theme.textDim),
                             at: .init(x: labelW, y: yMid), anchor: .trailing)
                    ctx.draw(Text(String(format: "%+7.0f pcm", val * 1e5))
                                .font(Theme.readoutSm).foregroundColor(color == Theme.accent ? .white : color),
                             at: .init(x: W - 6, y: yMid), anchor: .trailing)

                    let bLen = min(halfW, CGFloat(abs(val)) * scale)
                    guard bLen > 0.5 else { continue }
                    let barX = val >= 0 ? zeroX : zeroX - bLen
                    ctx.fill(
                        Path(roundedRect: CGRect(x: barX, y: yMid - (rowH - 12) / 2,
                                                 width: bLen, height: rowH - 12),
                             cornerRadius: 2, style: .continuous),
                        with: .color(color.opacity(0.75)))
                }
            }
            .padding(.horizontal, Theme.panelPadding)
            .padding(.vertical, 8)
            .frame(maxHeight: .infinity)
        }
        .panelSurface()
        .frame(maxWidth: .infinity)
    }
}

// MARK: — Xenon trend
private struct XenonTrendPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "XENON / IODINE INVENTORY TREND")
            // Stacked WIDE strip charts — tall narrow towers read as broken
            VStack(spacing: 10) {
                TrendView(values: supervisor.orderedHistory(supervisor.histPower),
                          yLo: 0, yHi: 130, color: Theme.accent, label: "RX POWER %")
                TrendView(values: supervisor.orderedHistory(supervisor.histReact).map { $0 * 1e5 },
                          yLo: -18000, yHi: 1000, color: Theme.accent, label: "REACTIVITY pcm")
                TrendView(values: supervisor.orderedHistory(supervisor.histCoolT),
                          yLo: 400, yHi: 650, color: Theme.accent, label: "RCS T-AVG K")
            }
            .padding(Theme.panelPadding)
            .frame(maxHeight: .infinity)
        }
        .panelSurface()
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
                ctx.stroke(axes, with: .color(Theme.ink.opacity(0.15)), lineWidth: 1)

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
                // Position marker is a normal indication, not an alarm — accent, not red.
                // Label flips side near full insertion so it never clips the panel edge.
                ctx.fill(Path(ellipseIn: .init(x: mpx-4, y: mpy-4, width: 8, height: 8)),
                         with: .color(Theme.accent))
                let steps = Int((228 * (1 - rodPos)).rounded())
                ctx.draw(
                    Text("\(steps) SWD").font(Theme.readoutSm).foregroundColor(Theme.ink),
                    at: .init(x: rodPos > 0.72 ? mpx - 10 : mpx + 10, y: mpy),
                    anchor: rodPos > 0.72 ? .trailing : .leading
                )

                // Axis calibration: insertion on x, worth in pcm on y
                ctx.draw(Text("0%").font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: .init(x: pad, y: H - pad + 4), anchor: .top)
                ctx.draw(Text("INSERTION").font(.system(size: 8, design: .monospaced)).foregroundColor(Theme.textDim.opacity(0.7)),
                         at: .init(x: (W) / 2, y: H - pad + 4), anchor: .top)
                ctx.draw(Text("100%").font(.system(size: 9, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: .init(x: W - pad, y: H - pad + 4), anchor: .top)
                ctx.draw(Text("0 pcm").font(.system(size: 8, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: .init(x: pad + 5, y: pad + 2), anchor: .topLeading)
                ctx.draw(Text("−5000").font(.system(size: 8, design: .monospaced)).foregroundColor(Theme.textDim),
                         at: .init(x: pad + 5, y: H - pad - 4), anchor: .bottomLeading)
            }
            .padding(4)
            .frame(maxHeight: .infinity)
        }
        .panelSurface()
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
                readRow("REACTIVITY",    String(format: "%+7.0f pcm", s.reactivity * 1e5), .reactivityStatus(s.reactivity))
                readRow("RX POWER",      String(format: "%7.2f %%", s.powerFraction * 100), .powerStatus(s.powerFraction))
                readRow("XE-135 INV",    String(format: "%7.1f rel", s.xenonInventory), Theme.text)
                readRow("I-135 INV",     String(format: "%7.1f rel", s.iodineInventory), Theme.textDim)
                readRow("XENON WORTH",   String(format: "%+7.0f pcm", s.xenonInventory * -1.6e-5 * 1e5), Theme.text)
                readRow("DECAY HEAT",    String(format: "%7.3f %%", s.decayHeatFraction * 100),
                        s.decayHeatFraction > 0.03 ? Theme.caution : Theme.text)
                readRow("RODS D-BANK",   String(format: "%4d SWD", Int((228 * (1 - s.rodPosition)).rounded())), Theme.text)
                readRow("RCS BORON",     String(format: "%7.1f ppm", supervisor.boronPPM),
                        supervisor.boronPPM < 100 ? Theme.caution : Theme.text)
            }
            .padding(Theme.panelPadding)
        }
        .panelSurface()
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
