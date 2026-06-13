// PrimaryTab.swift — F2 Primary Circuit: pressurizer, RCP, coolant detail.

import SwiftUI

struct PrimaryTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 12) {
            // Left: large bar meters
            VStack(spacing: 0) {
                PanelHeader(title: "PRIMARY SYSTEM — INSTRUMENTATION")
                HStack(alignment: .top, spacing: 14) {
                    ForEach(barDefs, id: \.label) { def in
                        BarMeterView(value: def.value, lo: def.lo, hi: def.hi,
                                     label: def.label, unit: def.unit,
                                     tripHi: def.tripHi, tripLo: def.tripLo,
                                     warnHi: def.warnHi, warnLo: def.warnLo)
                            .frame(maxWidth: .infinity)
                    }
                }
                .padding(20)
                .frame(maxHeight: .infinity)
            }
            .panelSurface()
            .frame(maxWidth: .infinity)

            // Right: DCS readout grid + status
            VStack(spacing: 12) {
                DCSGrid(supervisor: supervisor)
                RCPPanel(supervisor: supervisor)
                PressPanel(supervisor: supervisor)
            }
            .frame(width: 320)
        }
        .padding(12)
    }

    private var barDefs: [(label: String, unit: String, value: Double,
                          lo: Double, hi: Double, tripHi: Double?, tripLo: Double?,
                          warnHi: Double?, warnLo: Double?)] {
        let s = supervisor.snapshot
        // Warn/trip markers match the protection setpoints in updateAlarms()
        return [
            ("PZR PRESS", "MPa", supervisor.pressureMPa, 0, 18, 17.0, nil, 16.3, nil),
            ("T-FUEL",    "K",   s.fuelTempK,   300, 1800, 1500, nil, 1400, nil),
            ("RCS T-AVG", "K",   s.coolantTempK, 300, 700, 620, nil, 610, nil),
            ("RCP SPEED", "%",   supervisor.omegaRCP * 100, 0, 110, nil, nil, nil, 87),
            ("RX POWER",  "%",   s.powerFraction * 100, 0, 130, 120, nil, 115, nil),
            ("REACTIVITY","pcm", s.reactivity * 1e5, -500, 500, nil, nil, nil, nil),
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
        .panelSurface()
    }

    private func readouts(s: PlantSnapshot, nomP: Double) -> [(String, String, Color)] {[
        ("RX POWER",    String(format: "%6.2f %%",    s.powerFraction * 100),   .powerStatus(s.powerFraction)),
        ("THERMAL PWR", String(format: "%6.1f MWt",   s.thermalPowerW / 1e6),   Theme.text),
        ("PZR PRESS",   String(format: "%6.3f MPa",   supervisor.pressureMPa),
            supervisor.pressureMPa > 17.0 ? Theme.alarm : supervisor.pressureMPa > 16.3 ? Theme.caution : Theme.text),
        ("T-FUEL AVG",  String(format: "%6.1f K",     s.fuelTempK),
            s.fuelTempK > 1500 ? Theme.alarm : s.fuelTempK > 1200 ? Theme.caution : Theme.text),
        ("RCS T-AVG",   String(format: "%6.1f K",     s.coolantTempK),
            s.coolantTempK > 620 ? Theme.alarm : Theme.text),
        ("REACTIVITY",  String(format: "%+6.0f pcm",  s.reactivity * 1e5),       .reactivityStatus(s.reactivity)),
        // Xenon worth = −coeff·X (the old row mistakenly showed TOTAL reactivity)
        ("XENON WORTH", String(format: "%+6.0f pcm",  s.xenonInventory * -1.6e-5 * 1e5), Theme.textDim),
        ("RCS BORON",   String(format: "%6.1f ppm",   supervisor.boronPPM),
            supervisor.boronPPM < 100 ? Theme.caution : Theme.text),
    ]}
}

private struct RCPPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        HStack(spacing: 12) {
            VStack(spacing: 6) {
                LEDIndicator(label: "RCP RUN", on: supervisor.omegaRCP > 0.9, color: Theme.accent)
                LEDIndicator(label: "PORV",    on: supervisor.porvOpen,        color: Theme.caution)
                LEDIndicator(label: "ECCS",    on: supervisor.eccsActuated,    color: Theme.alarm)
            }
            VStack(alignment: .leading, spacing: 6) {
                readoutRow("RCP SPEED", String(format: "%5.1f %%", supervisor.omegaRCP * 100),
                           supervisor.omegaRCP < 0.87 ? Theme.caution : Theme.text)
                readoutRow("RCS FLOW", String(format: "%5.1f %%", supervisor.primaryFlow * supervisor.omegaRCP * 100), Theme.text)
                readoutRow("FW INV", String(format: "%5.3f", supervisor.feedwaterInv),
                           supervisor.feedwaterInv < 0.1 ? Theme.alarm : Theme.text)
            }
        }
        .padding(10)
        .panelSurface()
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
                      yLo: 13, yHi: 17, color: Theme.accent, label: "PZR PRESS MPa")
                .padding(8)
                .frame(height: 120)
        }
        .panelSurface()
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
        .padding(8)
        .frame(maxWidth: .infinity, alignment: .leading)
        .readoutSurface()
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
