// PrimaryTab.swift — F2 Primary Circuit: pressurizer, RCP, coolant detail.

import SwiftUI

struct PrimaryTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 8) {
            // Left: large bar meters
            VStack(spacing: 0) {
                PanelHeader(title: "PRIMARY SYSTEM — INSTRUMENTATION")
                HStack(alignment: .top, spacing: 10) {
                    ForEach(barDefs, id: \.label) { def in
                        BarMeterView(value: def.value, lo: def.lo, hi: def.hi,
                                     label: def.label, unit: def.unit,
                                     tripHi: def.tripHi, tripLo: def.tripLo)
                            .frame(maxWidth: .infinity)
                    }
                }
                .padding(16)
                .frame(maxHeight: .infinity)
            }
            .background(Theme.panel)
            .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
            .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
            .overlay(CornerBrackets())
            .frame(maxWidth: .infinity)

            // Right: DCS readout grid + status
            VStack(spacing: 8) {
                DCSGrid(supervisor: supervisor)
                RCPPanel(supervisor: supervisor)
                PressPanel(supervisor: supervisor)
            }
            .frame(width: 320)
        }
        .padding(8)
    }

    private var barDefs: [(label: String, unit: String, value: Double,
                          lo: Double, hi: Double, tripHi: Double?, tripLo: Double?)] {
        let s = supervisor.snapshot
        let nomP: Double = 15.5
        return [
            ("PRESSURE",  "MPa", supervisor.pressureMPa, 0, 18, nomP * 1.08, nil),
            ("FUEL T",    "K",   s.fuelTempK,   300, 1800, 1400, nil),
            ("COOLANT T", "K",   s.coolantTempK, 300, 700, 616, nil),
            ("RCP SPEED", "%",   supervisor.omegaRCP * 100, 0, 110, nil, 87),
            ("POWER",     "%",   s.powerFraction * 100, 0, 130, 120, nil),
            ("REACTIVITY","pcm", s.reactivity * 1e5, -500, 500, nil, nil),
        ]
    }
}

private struct DCSGrid: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        let nomP = 15.5
        VStack(spacing: 0) {
            PanelHeader(title: "PRIMARY PARAMETERS")
            LazyVGrid(columns: [GridItem(.flexible()), GridItem(.flexible())], spacing: 6) {
                ForEach(readouts(s: s, nomP: nomP), id: \.0) { (lbl, val, col) in
                    DCSReadout(label: lbl, value: val, color: col)
                }
            }
            .padding(8)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
    }

    private func readouts(s: PlantSnapshot, nomP: Double) -> [(String, String, Color)] {[
        ("Core power",    String(format: "%.2f %%",     s.powerFraction * 100),   .powerStatus(s.powerFraction)),
        ("Thermal power", String(format: "%.1f MWt",   s.thermalPowerW / 1e6),   Theme.text),
        ("Pressure",      String(format: "%.3f MPa",   supervisor.pressureMPa),
            supervisor.pressureMPa > nomP * 1.08 ? Theme.alarm : supervisor.pressureMPa > nomP * 1.04 ? Theme.caution : Theme.normal),
        ("Fuel temp",     String(format: "%.1f K",     s.fuelTempK),
            s.fuelTempK > 1400 ? Theme.alarm : s.fuelTempK > 1200 ? Theme.caution : Theme.normal),
        ("Coolant temp",  String(format: "%.1f K",     s.coolantTempK),
            s.coolantTempK > 616 ? Theme.alarm : Theme.normal),
        ("Reactivity",    String(format: "%+.5f",      s.reactivity),             .reactivityStatus(s.reactivity)),
        ("Xenon worth",   String(format: "%+.0f pcm",  s.reactivity * 1e5),       Theme.textDim),
        ("Boron",         String(format: "%.1f ppm",   supervisor.boronPPM),
            supervisor.boronPPM < 100 ? Theme.caution : Theme.text),
    ]}
}

private struct RCPPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        HStack(spacing: 12) {
            VStack(spacing: 4) {
                LEDIndicator(label: "RCP RUN", on: supervisor.omegaRCP > 0.9, color: Theme.normal)
                LEDIndicator(label: "PORV",    on: supervisor.porvOpen,        color: Theme.caution)
                LEDIndicator(label: "ECCS",    on: supervisor.eccsActuated,    color: Theme.alarm)
            }
            VStack(alignment: .leading, spacing: 4) {
                readoutRow("RCP speed", String(format: "%.1f %%", supervisor.omegaRCP * 100),
                           supervisor.omegaRCP < 0.87 ? Theme.caution : Theme.normal)
                readoutRow("Flow", String(format: "%.1f %%", supervisor.primaryFlow * supervisor.omegaRCP * 100), Theme.text)
                readoutRow("FW inv.", String(format: "%.3f", supervisor.feedwaterInv),
                           supervisor.feedwaterInv < 0.1 ? Theme.alarm : Theme.text)
            }
        }
        .padding(10)
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxWidth: .infinity)
    }
    private func readoutRow(_ l: String, _ v: String, _ c: Color) -> some View {
        HStack {
            Text(l).font(Theme.readoutSm).foregroundStyle(Theme.textDim)
            Spacer()
            Text(v).font(Theme.readout).foregroundStyle(c)
        }
    }
}

private struct PressPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "PRESSURIZER")
            TrendView(values: supervisor.orderedHistory(supervisor.histPress),
                      yLo: 13, yHi: 17, color: Theme.accent, label: "Pressure (MPa)")
                .padding(8)
                .frame(height: 120)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
    }
}

// MARK: — Shared small components

struct DCSReadout: View {
    let label: String
    let value: String
    let color: Color
    var body: some View {
        VStack(alignment: .leading, spacing: 2) {
            Text(label).font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
            Text(value).font(Theme.readoutSm).foregroundStyle(color).lineLimit(1).minimumScaleFactor(0.7)
        }
        .padding(6)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(Color(r: 8, g: 10, b: 14), in: RoundedRectangle(cornerRadius: 4))
        .overlay(RoundedRectangle(cornerRadius: 4).stroke(Theme.border, lineWidth: 0.5))
    }
}

struct LEDIndicator: View {
    let label: String
    let on: Bool
    let color: Color
    var body: some View {
        HStack(spacing: 6) {
            Circle()
                .fill(on ? color : Theme.border)
                .frame(width: 8, height: 8)
                .shadow(color: on ? color.opacity(0.8) : .clear, radius: 4)
            Text(label)
                .font(Theme.readoutSm)
                .foregroundStyle(on ? color : Theme.textDim)
        }
    }
}
