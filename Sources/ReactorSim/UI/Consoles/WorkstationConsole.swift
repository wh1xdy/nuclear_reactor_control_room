// WorkstationConsole.swift — tiled DCS operator workstation.
// Four fixed panes, all visible at once (no tabs): alarm summary, plant mimic,
// a multi-pen trend group, and control faceplates with SP/PV/OUT. Reuses the
// same components as the other consoles over the shared supervisor.

import SwiftUI

struct WorkstationConsole: View {
    let supervisor: PlantSupervisor

    var body: some View {
        VStack(spacing: 10) {
            HStack(spacing: 10) {
                AlarmSummaryPane(supervisor: supervisor).frame(width: 420)
                MimicPane(supervisor: supervisor).frame(maxWidth: .infinity)
            }
            .frame(maxHeight: .infinity)
            HStack(spacing: 10) {
                TrendGroupPane(supervisor: supervisor).frame(maxWidth: .infinity)
                FaceplatePane(supervisor: supervisor).frame(width: 460)
            }
            .frame(height: 280)
        }
        .padding(10)
    }
}

// MARK: — Alarm summary pane
private struct AlarmSummaryPane: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            HStack {
                PanelHeader(title: "ALARM SUMMARY — CHRONOLOGICAL")
                Spacer()
                Button { supervisor.acknowledgeAllAlarms() } label: {
                    Text("ACK [C]").font(.system(size: 9, weight: .semibold, design: .monospaced))
                        .foregroundStyle(Theme.accent).padding(.horizontal, 10).padding(.vertical, 5)
                        .contentShape(Rectangle())
                }.buttonStyle(.plain).controlSurface()
                .padding(.trailing, 8)
            }
            if supervisor.alarms.isEmpty {
                Spacer()
                Text("NO ACTIVE ALARMS")
                    .font(.system(size: 11, design: .monospaced)).foregroundStyle(Theme.textDim.opacity(0.7))
                Spacer()
            } else {
                ScrollView {
                    VStack(spacing: 0) {
                        ForEach(supervisor.alarms) { a in
                            HStack(spacing: 8) {
                                Text(String(format: "T+%05.0f", a.simTime))
                                    .frame(width: 70, alignment: .leading).foregroundStyle(Theme.textDim)
                                Text("P\(a.priority)").frame(width: 26).foregroundStyle(a.isTrip ? Theme.alarm : Theme.caution)
                                Text(a.state == "unack" ? "●" : "○").frame(width: 16)
                                    .foregroundStyle(a.state == "unack" ? Theme.caution : Theme.textDim)
                                Text(a.message).frame(maxWidth: .infinity, alignment: .leading)
                                    .foregroundStyle(a.isTrip ? Theme.alarm : Theme.text).lineLimit(1)
                            }
                            .font(.system(size: 10, design: .monospaced))
                            .padding(.horizontal, 12).padding(.vertical, 6)
                            Divider().background(Theme.sep.opacity(0.5))
                        }
                    }
                }
            }
        }
        .panelSurface()
    }
}

// MARK: — Mimic pane
private struct MimicPane: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "PLANT MIMIC")
            MimicDiagram(supervisor: supervisor)
                .clipShape(.rect(cornerRadius: Theme.controlRadius, style: .continuous))
                .padding(8)
        }
        .panelSurface()
    }
}

// MARK: — Trend group pane (4 pens)
private struct TrendGroupPane: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "TREND GROUP — 4 PEN")
            VStack(spacing: 6) {
                HStack(spacing: 6) {
                    TrendView(values: supervisor.orderedHistory(supervisor.histPower),
                              yLo: 0, yHi: 130, color: Theme.accent, label: "RX PWR %")
                    TrendView(values: supervisor.orderedHistory(supervisor.histCoolT),
                              yLo: 400, yHi: 650, color: Theme.accent, label: "T-AVG K")
                }
                HStack(spacing: 6) {
                    // Pressure trend: label + axis follow the kind (BWR dome ≈ 7 MPa vs PWR/SMR 15.5).
                    TrendView(values: supervisor.orderedHistory(supervisor.histPress),
                              yLo: supervisor.nominalPressureMPa - 2.5, yHi: supervisor.nominalPressureMPa + 1.5,
                              color: Theme.accent,
                              label: supervisor.hasPressurizer ? "PZR P MPa" : "DOME P MPa")
                    // Electric axis scales from the kind's nominal so an SMR isn't a flatline.
                    TrendView(values: supervisor.orderedHistory(supervisor.histElec),
                              yLo: 0, yHi: supervisor.nominalMWe * 1.1, color: Theme.accent, label: "GROSS MWe")
                }
            }
            .padding(10)
        }
        .panelSurface()
    }
}

// MARK: — Control faceplates (SP / PV / OUT)
private struct FaceplatePane: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        VStack(spacing: 0) {
            PanelHeader(title: "CONTROL FACEPLATES")
            VStack(spacing: 8) {
                faceplate("ROD CONTROL", auto: supervisor.rodAutoEnabled,
                          sp: "550.0 K", pv: String(format: "%.1f K", s.coolantTempK),
                          out: "\(Int((228 * (1 - supervisor.rodPosition)).rounded())) SWD") {
                    supervisor.rodAutoEnabled.toggle()
                }
                // PWR/SMR have a controllable pressurizer; a BWR holds dome pressure with
                // no PZR controller, so show a read-only STEAM DOME faceplate instead.
                if supervisor.hasPressurizer {
                    faceplate("PRESSURIZER", auto: supervisor.pzrAutoEnabled,
                              sp: String(format: "%.2f MPa", supervisor.nominalPressureMPa),
                              pv: String(format: "%.2f MPa", supervisor.pressureMPa),
                              out: supervisor.pressureMPa < supervisor.nominalPressureMPa ? "HEATERS" : "SPRAY") {
                        supervisor.pzrAutoEnabled.toggle()
                    }
                } else {
                    faceplate("STEAM DOME", auto: false,
                              sp: String(format: "%.2f MPa", supervisor.nominalPressureMPa),
                              pv: String(format: "%.2f MPa", supervisor.pressureMPa),
                              out: supervisor.porvOpen ? "SRV OPEN" : "TURB RELIEF",
                              showAuto: false) { }
                }
                faceplate("FEEDWATER", auto: supervisor.fwAutoEnabled,
                          sp: supervisor.hasSteamGenerator ? "SG LEVEL" : "RPV LVL",
                          pv: String(format: "%.3f", supervisor.feedwaterInv),
                          out: String(format: "%.0f%% VLV", supervisor.feedwaterValve * 100)) {
                    supervisor.fwAutoEnabled.toggle()
                }
            }
            .padding(10)
            Spacer(minLength: 0)
        }
        .panelSurface()
    }

    private func faceplate(_ name: String, auto: Bool, sp: String, pv: String, out: String,
                           showAuto: Bool = true,
                           toggle: @escaping () -> Void) -> some View {
        HStack(spacing: 10) {
            VStack(alignment: .leading, spacing: 4) {
                Text(name).font(.system(size: 11, weight: .semibold, design: .monospaced)).foregroundStyle(Theme.ink)
                HStack(spacing: 12) {
                    col("SP", sp); col("PV", pv); col("OUT", out)
                }
            }
            Spacer()
            // No AUTO/MAN toggle for read-only faceplates (e.g. BWR steam dome — no PZR controller).
            if showAuto {
                Button(action: toggle) {
                    Text(auto ? "AUTO" : "MAN")
                        .font(.system(size: 10, weight: .bold, design: .monospaced))
                        .foregroundStyle(auto ? Theme.ink : Theme.textDim)
                        .frame(width: 56).padding(.vertical, 8)
                        .contentShape(Rectangle())
                }
                .buttonStyle(.plain)
                .controlSurface(tint: auto ? Theme.accent : nil)
            }
        }
        .padding(10)
        .background(Theme.ink.opacity(0.04))
        .clipShape(.rect(cornerRadius: Theme.controlRadius, style: .continuous))
        .overlay(RoundedRectangle(cornerRadius: Theme.controlRadius, style: .continuous)
            .strokeBorder(Theme.border.opacity(0.6), lineWidth: 1))
    }

    private func col(_ label: String, _ value: String) -> some View {
        VStack(alignment: .leading, spacing: 1) {
            Text(label).font(.system(size: 8, design: .monospaced)).foregroundStyle(Theme.textDim)
            Text(value).font(.system(size: 11, weight: .medium, design: .monospaced)).foregroundStyle(Theme.text)
        }
    }
}
