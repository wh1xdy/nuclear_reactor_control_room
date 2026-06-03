// SecondaryTab.swift — F3 Secondary: steam generator, turbine, condenser, feedwater.

import SwiftUI

struct SecondaryTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 8) {
            // Bar meters
            VStack(spacing: 0) {
                PanelHeader(title: "SECONDARY SYSTEM — INSTRUMENTATION")
                HStack(alignment: .top, spacing: 10) {
                    ForEach(barDefs, id: \.label) { def in
                        BarMeterView(value: def.value, lo: def.lo, hi: def.hi,
                                     label: def.label, unit: def.unit,
                                     tripHi: def.tripHi)
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

            // Right: DCS boxes + electric output
            VStack(spacing: 8) {
                SecondaryDCSGrid(supervisor: supervisor)
                ElectricOutputPanel(supervisor: supervisor)
                SecondaryTrends(supervisor: supervisor)
            }
            .frame(width: 320)
        }
        .padding(8)
    }

    private var barDefs: [(label: String, unit: String, value: Double,
                          lo: Double, hi: Double, tripHi: Double?)] {
        let s = supervisor.snapshot
        return [
            ("SG TEMP",    "K",   s.sgTempK,           300, 700, 620),
            ("COND TEMP",  "K",   supervisor.condTempK, 280, 380, 340),
            ("STEAM INV",  "",    supervisor.steamInv,  0,   2.0, nil),
            ("FW INV",     "",    supervisor.feedwaterInv, 0, 1.2, nil),
            ("ELEC PWR",   "MWe", s.electricPowerW / 1e6, 0, 1200, nil),
            ("TURB VALVE", "%",   supervisor.turbineValve * 100, 0, 105, nil),
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
                DCSReadout(label: "Thermal power",  value: String(format: "%.1f MWt", s.thermalPowerW / 1e6), color: Theme.text)
                DCSReadout(label: "Electric power", value: String(format: "%.1f MWe", s.electricPowerW / 1e6), color: Theme.normal)
                DCSReadout(label: "SG temp",        value: String(format: "%.1f K", s.sgTempK),
                           color: s.sgTempK > 600 ? Theme.caution : Theme.text)
                DCSReadout(label: "Condenser",      value: String(format: "%.1f K", supervisor.condTempK),
                           color: supervisor.condTempK > 330 ? Theme.caution : Theme.text)
                DCSReadout(label: "Steam inv.",     value: String(format: "%.3f", supervisor.steamInv),
                           color: supervisor.steamInv < 0.3 ? Theme.alarm : Theme.text)
                DCSReadout(label: "FW inv.",        value: String(format: "%.3f", supervisor.feedwaterInv),
                           color: supervisor.feedwaterInv < 0.1 ? Theme.alarm : Theme.text)
                DCSReadout(label: "Turbine valve",  value: String(format: "%.0f %%", supervisor.turbineValve * 100), color: Theme.text)
                DCSReadout(label: "Decay heat",     value: String(format: "%.2f %%", s.decayHeatFraction * 100),
                           color: s.decayHeatFraction > 0.03 ? Theme.caution : Theme.text)
            }
            .padding(8)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
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
                        .foregroundStyle(Theme.normal)
                    Text("MWe")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.textDim)
                }
                VStack(alignment: .leading, spacing: 4) {
                    readoutRow("Thermal",    String(format: "%.1f MWt", s.thermalPowerW / 1e6), Theme.text)
                    readoutRow("Efficiency", String(format: "%.1f %%",  efficiency), Theme.text)
                    readoutRow("Trip",       supervisor.turbineTrip ? "TRIPPED" : "NORMAL",
                               supervisor.turbineTrip ? Theme.alarm : Theme.normal)
                }
            }
            .padding(12)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
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
            VStack(spacing: 6) {
                TrendView(values: supervisor.orderedHistory(supervisor.histElec),
                          yLo: 0, yHi: 1100, color: Theme.normal, label: "Electric power (MWe)")
                TrendView(values: supervisor.orderedHistory(supervisor.histSteamT),
                          yLo: 400, yHi: 650, color: Theme.twophase, label: "SG temperature (K)")
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
