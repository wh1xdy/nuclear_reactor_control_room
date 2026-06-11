// OverviewTab.swift — F1 Overview: P&ID schematic + instruments + status readouts.
// Glass panels (macOS 26) over Canvas drawing layer.

import SwiftUI

struct OverviewTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 12) {
            // Center: P&ID + wide strip charts (2×2 — proper chart proportions)
            VStack(spacing: 12) {
                PIDPanel(supervisor: supervisor)
                TrendsPanel(supervisor: supervisor)
            }
            .frame(maxWidth: .infinity)

            // Right: instruments + dense status column
            VStack(spacing: 12) {
                InstrumentPanel(supervisor: supervisor)
                StatusReadoutsPanel(supervisor: supervisor)
            }
            .frame(width: 330)
        }
        .padding(12)
    }
}

// MARK: — P&ID Panel
private struct PIDPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "PLANT OVERVIEW — PWR")
            PIDCanvas(supervisor: supervisor)
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .padding(8)
        }
        .glassEffect(.regular, in: .rect(cornerRadius: Theme.panelRadius, style: .continuous))
        .frame(maxWidth: .infinity)
        .frame(height: 320)
    }
}

// MARK: — Status readouts panel
private struct StatusReadoutsPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let snap = supervisor.snapshot
        // Rod position in steps withdrawn (0–228), real CRDM convention
        let rodSWD = Int((Double(228) * (1.0 - snap.rodPosition)).rounded())
        VStack(spacing: 0) {
            PanelHeader(title: "PLANT STATUS")
            Grid(alignment: .leading, horizontalSpacing: 14, verticalSpacing: 4) {
                statusRow("RX POWER",     (snap.powerFraction*100).fmt("%6.2f"), "%",
                          .powerStatus(snap.powerFraction))
                statusRow("THERMAL PWR",  (snap.thermalPowerW/1e6).fmt("%6.1f"), "MWt", Theme.text)
                statusRow("GROSS ELEC",   (snap.electricPowerW/1e6).fmt("%6.1f"), "MWe", Theme.text)
                Divider().background(Theme.sep)
                statusRow("REACTIVITY",   (snap.reactivity*1e5).fmt("%+6.0f"), "pcm",
                          .reactivityStatus(snap.reactivity))
                statusRow("XENON WORTH",  (snap.xenonInventory * -1.6e-5 * 1e5).fmt("%+6.0f"), "pcm", Theme.textDim)
                statusRow("RODS D-BANK",  "\(rodSWD)", "SWD", Theme.text)
                Divider().background(Theme.sep)
                statusRow("T-FUEL AVG",   snap.fuelTempK.fmt("%6.1f"), "K",
                          snap.fuelTempK > 1400 ? Theme.alarm : snap.fuelTempK > 1200 ? Theme.caution : Theme.text)
                statusRow("RCS T-AVG",    snap.coolantTempK.fmt("%6.1f"), "K",
                          snap.coolantTempK > 616 ? Theme.alarm : snap.coolantTempK > 580 ? Theme.caution : Theme.text)
                statusRow("SG TEMP",      snap.sgTempK.fmt("%6.1f"), "K", Theme.text)
                Divider().background(Theme.sep)
                statusRow("PZR PRESS",    supervisor.pressureMPa.fmt("%6.3f"), "MPa",
                          supervisor.pressureMPa > 17.0 ? Theme.alarm : supervisor.pressureMPa > 16.3 ? Theme.caution : Theme.text)
                statusRow("RCS BORON",    supervisor.boronPPM.fmt("%6.1f"), "ppm",
                          supervisor.boronPPM < 100 ? Theme.caution : Theme.text)
                statusRow("DECAY HEAT",   (snap.decayHeatFraction*100).fmt("%6.2f"), "%",
                          snap.decayHeatFraction > 0.03 ? Theme.caution : Theme.text)
                Divider().background(Theme.sep)
                statusRow("RCP SPEED",    (supervisor.omegaRCP*100).fmt("%6.1f"), "%",
                          supervisor.omegaRCP < 0.87 ? Theme.caution : Theme.text)
                statusRow("RCS FLOW",     (supervisor.primaryFlow*supervisor.omegaRCP*100).fmt("%6.1f"), "%", Theme.text)
                statusRow("FW INV",       supervisor.feedwaterInv.fmt("%6.3f"), "rel",
                          supervisor.feedwaterInv < 0.1 ? Theme.alarm : Theme.text)
                statusRow("STEAM INV",    supervisor.steamInv.fmt("%6.3f"), "rel",
                          supervisor.steamInv < 0.3 ? Theme.alarm : Theme.text)
                statusRow("COND TEMP",    supervisor.condTempK.fmt("%6.1f"), "K",
                          supervisor.condTempK > 330 ? Theme.caution : Theme.text)
                statusRow("STM DUMP",     (supervisor.steamDumpValve*100).fmt("%6.1f"), "%",
                          supervisor.steamDumpValve > 0 ? Theme.caution : Theme.text)
                statusRow("TBN BREAKER",  supervisor.turbineTrip ? "OPEN" : "CLOSED", "",
                          supervisor.turbineTrip ? Theme.alarm : Theme.text)
            }
            .padding(.horizontal, Theme.panelPadding)
            .padding(.vertical, 10)
            Spacer(minLength: 0)
        }
        .glassEffect(.regular, in: .rect(cornerRadius: Theme.panelRadius, style: .continuous))
        .frame(maxWidth: .infinity, maxHeight: .infinity)
    }

    // Label / value / unit columns — values right-aligned, fixed decimals
    private func statusRow(_ label: String, _ value: String, _ unit: String,
                           _ color: Color) -> some View {
        GridRow {
            Text(label)
                .font(.system(size: 11, design: .monospaced))
                .foregroundStyle(Theme.textDim)
                .tracking(0.5)
            Text(value)
                .font(.system(size: 13, weight: .medium, design: .monospaced))
                .foregroundStyle(color)
                .frame(maxWidth: .infinity, alignment: .trailing)
            Text(unit)
                .font(.system(size: 10, design: .monospaced))
                .foregroundStyle(Theme.textDim.opacity(0.8))
                .frame(width: 34, alignment: .leading)
        }
    }
}

// MARK: — Instruments panel (arc gauge + bar meters)
private struct InstrumentPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let snap = supervisor.snapshot
        VStack(spacing: 0) {
            PanelHeader(title: "INSTRUMENTS")
            HStack(alignment: .top, spacing: 8) {
                // Power arc gauge
                VStack(spacing: 4) {
                    Text("RX POWER")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.textDim)
                    ArcGaugeView(value: snap.powerFraction * 100, lo: 0, hi: 130,
                                 label: "", unit: "%", tripHi: 120)
                }
                .frame(maxWidth: .infinity)

                // Bar meters
                HStack(alignment: .top, spacing: 8) {
                    BarMeterView(value: supervisor.pressureMPa, lo: 0, hi: 18,
                                 label: "PZR PRESS", unit: "MPa", tripHi: 17.0, warnHi: 16.3)
                    BarMeterView(value: snap.fuelTempK, lo: 300, hi: 1800,
                                 label: "T-FUEL", unit: "K",   tripHi: 1500, warnHi: 1400)
                    BarMeterView(value: snap.coolantTempK, lo: 300, hi: 700,
                                 label: "RCS T-AVG", unit: "K",   tripHi: 620, warnHi: 610)
                }
                .frame(maxWidth: .infinity)
            }
            .padding(12)
        }
        .glassEffect(.regular, in: .rect(cornerRadius: Theme.panelRadius, style: .continuous))
        .frame(height: 180)
    }
}

// MARK: — Trends panel — 2×2 grid of WIDE strip charts (≈2.5:1, not towers)
private struct TrendsPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "TRENDS")
            VStack(spacing: 10) {
                HStack(spacing: 10) {
                    TrendView(values: supervisor.orderedHistory(supervisor.histPower),
                              yLo: 0, yHi: 130, color: Theme.accent, label: "RX POWER %")
                    TrendView(values: supervisor.orderedHistory(supervisor.histReact).map { $0 * 1e5 },
                              yLo: -1000, yHi: 1000, color: Theme.accent, label: "REACTIVITY pcm")
                }
                HStack(spacing: 10) {
                    TrendView(values: supervisor.orderedHistory(supervisor.histFuelT),
                              yLo: 400, yHi: 1800, color: Theme.accent, label: "T-FUEL K")
                    TrendView(values: supervisor.orderedHistory(supervisor.histDecay),
                              yLo: 0, yHi: 10, color: Theme.accent, label: "DECAY HEAT %")
                }
            }
            .padding(Theme.panelPadding)
        }
        .glassEffect(.regular, in: .rect(cornerRadius: Theme.panelRadius, style: .continuous))
        .frame(maxHeight: .infinity)
    }
}
