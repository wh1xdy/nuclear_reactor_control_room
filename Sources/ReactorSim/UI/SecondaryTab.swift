// SecondaryTab.swift — F3 Secondary: steam generator, turbine, condenser, feedwater.

import SwiftUI

struct SecondaryTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 12) {
            // Bar meters
            VStack(spacing: 0) {
                PanelHeader(title: "SECONDARY SYSTEM — INSTRUMENTATION")
                HStack(alignment: .top, spacing: 14) {
                    ForEach(barDefs, id: \.label) { def in
                        BarMeterView(value: def.value, lo: def.lo, hi: def.hi,
                                     label: def.label, unit: def.unit,
                                     tripHi: def.tripHi, warnHi: def.warnHi)
                            .frame(maxWidth: .infinity)
                    }
                }
                .padding(20)
                .frame(maxHeight: .infinity)
            }
            .panelSurface()
            .frame(maxWidth: .infinity)

            // Right: DCS boxes + electric output
            VStack(spacing: 12) {
                SecondaryDCSGrid(supervisor: supervisor)
                ElectricOutputPanel(supervisor: supervisor)
                SecondaryTrends(supervisor: supervisor)
            }
            .frame(width: 320)
        }
        .padding(12)
    }

    private var barDefs: [(label: String, unit: String, value: Double,
                          lo: Double, hi: Double, tripHi: Double?, warnHi: Double?)] {
        let s = supervisor.snapshot
        return [
            // BWR is direct-cycle with no steam generator — this reads dome temperature.
            (supervisor.hasSteamGenerator ? "SG TEMP" : "DOME TEMP", "K", s.sgTempK, 300, 700, 620, 600),
            ("COND TEMP",  "K",   supervisor.condTempK, 280, 380, 340, 330),
            ("STEAM INV",  "rel", supervisor.steamInv,  0,   2.0, nil, nil),
            ("FW INV",     "rel", supervisor.feedwaterInv, 0, 1.2, nil, nil),
            ("GROSS ELEC", "MWe", s.electricPowerW / 1e6, 0, supervisor.nominalMWe * 1.15, nil, nil),
            ("TBN GOV",    "%",   supervisor.turbineValve * 100, 0, 105, nil, nil),
        ]
    }
}

private struct SecondaryDCSGrid: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        VStack(spacing: 0) {
            PanelHeader(title: "SECONDARY PARAMETERS")
            LazyVGrid(columns: [GridItem(.flexible()), GridItem(.flexible())], spacing: 6) {
                DCSReadout(label: "THERMAL PWR",  value: String(format: "%6.1f MWt", s.thermalPowerW / 1e6), color: Theme.text)
                DCSReadout(label: "GROSS ELEC",   value: String(format: "%6.1f MWe", s.electricPowerW / 1e6), color: Theme.text)
                DCSReadout(label: supervisor.hasSteamGenerator ? "SG TEMP" : "DOME TEMP",
                           value: String(format: "%6.1f K", s.sgTempK),
                           color: s.sgTempK > 600 ? Theme.caution : Theme.text)
                DCSReadout(label: "STEAM PRESS",  value: String(format: "%6.3f MPa", s.steamPressureMPa), color: Theme.text)
                DCSReadout(label: "COND TEMP",    value: String(format: "%6.1f K", supervisor.condTempK),
                           color: supervisor.condTempK > 330 ? Theme.caution : Theme.text)
                DCSReadout(label: "STEAM INV",    value: String(format: "%6.3f", supervisor.steamInv),
                           color: supervisor.steamInv < 0.3 ? Theme.alarm : Theme.text)
                DCSReadout(label: "FW INV",       value: String(format: "%6.3f", supervisor.feedwaterInv),
                           color: supervisor.feedwaterInv < 0.1 ? Theme.alarm : Theme.text)
                DCSReadout(label: "TBN GOV VLV",  value: String(format: "%6.1f %%", supervisor.turbineValve * 100), color: Theme.text)
                DCSReadout(label: "STEAM DUMP",   value: String(format: "%6.1f %%", supervisor.steamDumpValve * 100),
                           color: supervisor.steamDumpValve > 0 ? Theme.caution : Theme.text)
                DCSReadout(label: "DECAY HEAT",   value: String(format: "%6.2f %%", s.decayHeatFraction * 100),
                           color: s.decayHeatFraction > 0.03 ? Theme.caution : Theme.text)
                DCSReadout(label: "TBN STATUS",   value: supervisor.turbineTrip ? "TRIPPED" : "ON LINE",
                           color: supervisor.turbineTrip ? Theme.alarm : Theme.text)
            }
            .padding(10)
        }
        .panelSurface()
    }
}

private struct ElectricOutputPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        let efficiency = s.thermalPowerW > 0 ? s.electricPowerW / s.thermalPowerW * 100 : 0
        VStack(spacing: 0) {
            PanelHeader(title: "TURBINE-GENERATOR")
            HStack(spacing: 16) {
                // Big electric power readout
                VStack(spacing: 2) {
                    Text(String(format: "%.1f", s.electricPowerW / 1e6))
                        .font(Theme.readoutXl)
                        .foregroundStyle(Theme.ink)
                    Text("MWe")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.textDim)
                }
                VStack(alignment: .leading, spacing: 6) {
                    readoutRow("THERMAL",    String(format: "%6.1f MWt", s.thermalPowerW / 1e6), Theme.text)
                    readoutRow("EFFICIENCY", String(format: "%6.1f %%",  efficiency), Theme.text)
                    readoutRow("STM DUMP",   String(format: "%6.1f %%",  supervisor.steamDumpValve * 100),
                               supervisor.steamDumpValve > 0 ? Theme.caution : Theme.text)
                    readoutRow("BREAKER",    supervisor.turbineTrip ? "OPEN" : "CLOSED",
                               supervisor.turbineTrip ? Theme.alarm : Theme.text)
                }
            }
            .padding(Theme.panelPadding)
        }
        .panelSurface()
    }
    private func readoutRow(_ l: String, _ v: String, _ c: Color) -> some View {
        HStack { Text(l).font(Theme.readoutSm).foregroundStyle(Theme.textDim); Spacer()
            Text(v).font(Theme.readoutSm).foregroundStyle(c) }
    }
}

private struct SecondaryTrends: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "SECONDARY TRENDS")
            VStack(spacing: 8) {
                TrendView(values: supervisor.orderedHistory(supervisor.histElec),
                          yLo: 0, yHi: supervisor.nominalMWe * 1.1, color: Theme.accent, label: "GROSS ELEC MWe")
                TrendView(values: supervisor.orderedHistory(supervisor.histSteamT),
                          yLo: 400, yHi: 650, color: Theme.accent,
                          label: supervisor.hasSteamGenerator ? "SG TEMP K" : "DOME TEMP K")
                TrendView(values: supervisor.orderedHistory(supervisor.histDecay),
                          yLo: 0, yHi: 10, color: Theme.accent, label: "DECAY HEAT %")
            }
            .padding(12)
        }
        .panelSurface()
        .frame(maxHeight: .infinity)
    }
}
