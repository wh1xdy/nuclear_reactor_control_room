// OverviewTab.swift — F1 Overview: P&ID schematic + instruments + status readouts.
// Glass panels (macOS 26) over Canvas drawing layer.

import SwiftUI

struct OverviewTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 8) {
            // Center: P&ID + status readouts
            VStack(spacing: 8) {
                PIDPanel(supervisor: supervisor)
                StatusReadoutsPanel(supervisor: supervisor)
            }
            .frame(maxWidth: .infinity)

            // Right: instruments + trends
            VStack(spacing: 8) {
                InstrumentPanel(supervisor: supervisor)
                TrendsPanel(supervisor: supervisor)
            }
            .frame(width: 280)
        }
        .padding(8)
    }
}

// MARK: — P&ID Panel
private struct PIDPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "PLANT OVERVIEW — PWR")
            PIDCanvas(snapshot: supervisor.snapshot, supervisor: supervisor)
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .padding(8)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxWidth: .infinity)
        .frame(height: 320)
    }
}

// MARK: — Status readouts panel
private struct StatusReadoutsPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let snap = supervisor.snapshot
        VStack(spacing: 0) {
            PanelHeader(title: "PLANT STATUS")
            Grid(alignment: .leading, horizontalSpacing: 20, verticalSpacing: 6) {
                statusRow("Core power",    (snap.powerFraction*100).fmt("%.2f %%"),
                          .powerStatus(snap.powerFraction))
                statusRow("Thermal power", (snap.thermalPowerW/1e6).fmt("%.1f MWt"), Theme.text)
                statusRow("Electric power",(snap.electricPowerW/1e6).fmt("%.1f MWe"), Theme.text)
                Divider().background(Theme.sep)
                statusRow("Reactivity",    snap.reactivity.fmt("%+.5f Δk/k"),
                          .reactivityStatus(snap.reactivity))
                statusRow("Xenon worth",   (snap.xenonInventory * -1.6e-5).fmt("%+.5f Δk/k"), Theme.textDim)
                Divider().background(Theme.sep)
                statusRow("Fuel temp",     snap.fuelTempK.fmt("%.1f K"),
                          snap.fuelTempK > 1400 ? Theme.alarm : snap.fuelTempK > 1200 ? Theme.caution : Theme.normal)
                statusRow("Coolant temp",  snap.coolantTempK.fmt("%.1f K"),
                          snap.coolantTempK > 616 ? Theme.alarm : snap.coolantTempK > 580 ? Theme.caution : Theme.normal)
                statusRow("Steam temp",    snap.sgTempK.fmt("%.1f K"), Theme.text)
                Divider().background(Theme.sep)
                statusRow("Pressure",      supervisor.pressureMPa.fmt("%.3f MPa"),
                          supervisor.pressureMPa > 17.0 ? Theme.alarm : supervisor.pressureMPa > 16.3 ? Theme.caution : Theme.normal)
                statusRow("Boron",         supervisor.boronPPM.fmt("%.1f ppm"),
                          supervisor.boronPPM < 100 ? Theme.caution : Theme.text)
                statusRow("Decay heat",    (snap.decayHeatFraction*100).fmt("%.2f %%"),
                          snap.decayHeatFraction > 0.03 ? Theme.caution : Theme.text)
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 8)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxWidth: .infinity)
    }

    private func statusRow(_ label: String, _ value: String, _ color: Color) -> some View {
        GridRow {
            Text(label)
                .font(Theme.readout)
                .foregroundStyle(Theme.textDim)
            Text(value)
                .font(Theme.readout)
                .foregroundStyle(color)
                .frame(maxWidth: .infinity, alignment: .trailing)
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
                    Text("POWER")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.textDim)
                    ArcGaugeView(value: snap.powerFraction * 100, lo: 0, hi: 130,
                                 label: "", unit: "%", tripHi: 120)
                }
                .frame(maxWidth: .infinity)

                // Bar meters
                HStack(alignment: .top, spacing: 6) {
                    BarMeterView(value: supervisor.pressureMPa, lo: 0, hi: 18,
                                 label: "PRESS", unit: "MPa", tripHi: 16.7)
                    BarMeterView(value: snap.fuelTempK, lo: 300, hi: 1800,
                                 label: "FUEL T", unit: "K",   tripHi: 1400)
                    BarMeterView(value: snap.coolantTempK, lo: 300, hi: 700,
                                 label: "COOL T", unit: "K",   tripHi: 616)
                }
                .frame(maxWidth: .infinity)
            }
            .padding(10)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(height: 180)
    }
}

// MARK: — Trends panel
private struct TrendsPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "TRENDS")
            VStack(spacing: 6) {
                TrendView(values: supervisor.orderedHistory(supervisor.histPower),
                          yLo: 0, yHi: 130, color: Theme.normal, label: "Core power (%)")
                TrendView(values: supervisor.orderedHistory(supervisor.histReact),
                          yLo: -0.01, yHi: 0.01, color: Theme.accent, label: "Reactivity (Δk/k)")
                TrendView(values: supervisor.orderedHistory(supervisor.histFuelT),
                          yLo: 400, yHi: 1800, color: Theme.warning, label: "Fuel temp (K)")
                TrendView(values: supervisor.orderedHistory(supervisor.histDecay),
                          yLo: 0, yHi: 10, color: Theme.caution, label: "Decay heat (%)")
            }
            .padding(8)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxHeight: .infinity)
    }
}
