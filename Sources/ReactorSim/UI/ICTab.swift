// ICTab.swift — F6 I&C: 2-of-3 protection channels + setpoints + auto controller status.

import SwiftUI

struct ICTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        VStack(spacing: 12) {
            ProtectionChannelsPanel(supervisor: supervisor)
            AutoControllersPanel(supervisor: supervisor)
        }
        .padding(12)
    }
}

// MARK: — 2-of-3 voting display
private struct ProtectionChannelsPanel: View {
    let supervisor: PlantSupervisor

    // Small gaussian noise to differentiate channels (seeded per-channel)
    private func channelVal(_ base: Double, seed: Int) -> Double {
        let noise = sin(Double(seed) * 12.9898 + supervisor.snapshot.time * 0.1) * 0.003
        return base * (1.0 + noise)
    }

    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "PROTECTION CHANNELS  —  2-OF-3 VOTING LOGIC")
            let s = supervisor.snapshot

            // Warn/trip setpoints match the protection system in updateAlarms()
            // Pressure channel + SG channel labels/setpoints follow the reactor kind.
            let nomP = supervisor.nominalPressureMPa
            let params: [(String, Double, String, Double, Double?)] = [
                ("NEUTRON FLUX",  s.powerFraction * 100, "%",    115.0, 120.0),
                ("T-FUEL AVG",    s.fuelTempK,           "K",   1400.0, 1500.0),
                (supervisor.hasPressurizer ? "PZR PRESSURE" : "DOME PRESSURE",
                                  supervisor.pressureMPa,"MPa",   nomP * 1.05, nomP * 1.097),
                ("RCS T-AVG",     s.coolantTempK,        "K",    610.0, 620.0),
                (supervisor.hasSteamGenerator ? "SG TEMP" : "DOME TEMP",
                                  s.sgTempK,             "K",    600.0, 620.0),
            ]

            VStack(spacing: 0) {
                // Column headers
                HStack(spacing: 0) {
                    Text("PARAMETER").font(Theme.readoutSm).foregroundStyle(Theme.textDim).frame(width: 200, alignment: .leading)
                    Text("VALUE").font(Theme.readoutSm).foregroundStyle(Theme.textDim).frame(width: 140, alignment: .leading)
                    Text("VOTES").font(Theme.readoutSm).foregroundStyle(Theme.textDim).frame(width: 80, alignment: .center)
                    Text("CH-A").font(Theme.readoutSm).foregroundStyle(Theme.textDim).frame(maxWidth: .infinity)
                    Text("CH-B").font(Theme.readoutSm).foregroundStyle(Theme.textDim).frame(maxWidth: .infinity)
                    Text("CH-C").font(Theme.readoutSm).foregroundStyle(Theme.textDim).frame(maxWidth: .infinity)
                    Text("SETPOINT").font(Theme.readoutSm).foregroundStyle(Theme.textDim).frame(width: 160, alignment: .trailing)
                }
                .padding(.horizontal, 12)
                .padding(.vertical, 8)

                Divider().background(Theme.sep)

                ForEach(Array(params.enumerated()), id: \.offset) { (i, param) in
                    let (name, value, unit, warnSP, tripSP) = param
                    let chVals = [channelVal(value, seed: i*3+1),
                                  channelVal(value, seed: i*3+2),
                                  channelVal(value, seed: i*3+3)]
                    let votes  = chVals.filter { tripSP != nil && $0 > tripSP! }.count
                    let warnOK = value < warnSP
                    let valCol = warnOK ? Theme.text : (tripSP == nil || value < tripSP! ? Theme.caution : Theme.alarm)
                    let voteCol = votes >= 2 ? Theme.alarm : votes == 1 ? Theme.caution : Theme.textDim

                    HStack(spacing: 0) {
                        Text(name).font(Theme.readout).foregroundStyle(Theme.text).frame(width: 200, alignment: .leading)
                        Text(String(format: "%.3f %@", value, unit))
                            .font(Theme.readout).foregroundStyle(valCol).frame(width: 140, alignment: .leading)
                        Text("\(votes)/3")
                            .font(Theme.readoutMd).foregroundStyle(voteCol).frame(width: 80, alignment: .center)
                        ForEach(Array(chVals.enumerated()), id: \.offset) { (ci, cv) in
                            let chCol: Color = {
                                if let t = tripSP, cv > t { return Theme.alarm }
                                if cv > warnSP             { return Theme.caution }
                                return Theme.text
                            }()
                            HStack(spacing: 4) {
                                Circle().fill(chCol).frame(width: 6, height: 6)
                                    .shadow(color: chCol.opacity(0.8), radius: 3)
                                Text(String(format: "%.2f", cv)).font(Theme.readoutSm).foregroundStyle(chCol)
                            }
                            .frame(maxWidth: .infinity)
                        }
                        // Static setpoint labels are reference info, not live alarms — grey.
                        VStack(alignment: .trailing, spacing: 2) {
                            Text(String(format: "Warn: %.2f %@", warnSP, unit))
                                .font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                            if let t = tripSP {
                                Text(String(format: "Trip: %.2f %@", t, unit))
                                    .font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                            }
                        }
                        .frame(width: 160, alignment: .trailing)
                    }
                    .padding(.horizontal, 12)
                    .padding(.vertical, 8)
                    .background((i % 2 == 0) ? Color.clear : Color(white: 1, opacity: 0.015))

                    Divider().background(Theme.sep)
                }
            }
            .padding(.bottom, 4)
            Spacer(minLength: 0)
            // Footer: coincidence logic status line — fills the panel honestly
            HStack {
                Text("2-OF-3 COINCIDENCE LOGIC  ·  CHANNELS IN TEST: NONE  ·  BYPASS: NONE")
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundStyle(Theme.textDim.opacity(0.7))
                Spacer()
                Text("RPS NORMAL")
                    .font(.system(size: 9, weight: .semibold, design: .monospaced))
                    .foregroundStyle(Theme.textDim)
            }
            .padding(.horizontal, Theme.panelPadding)
            .padding(.bottom, 10)
        }
        .frame(maxHeight: .infinity)
        .panelSurface()
    }
}

// MARK: — Auto controllers (live, clickable)
private struct AutoControllersPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        VStack(spacing: 0) {
            PanelHeader(title: "AUTO CONTROLLERS  —  CLICK TO TOGGLE")
            HStack(spacing: 12) {
                autoCard(name: "ROD AUTO-CONTROL",
                         enabled: supervisor.rodAutoEnabled,
                         setpoint: "550.0 K",
                         measurement: String(format: "%6.1f K", s.coolantTempK),
                         output: String(format: "OUT %4d SWD", Int((228 * (1 - supervisor.rodPosition)).rounded())),
                         description: "T-avg program — drives rod demand at CRDM rate, 0.5 K deadband") {
                    supervisor.rodAutoEnabled.toggle()
                }
                // A BWR has no pressurizer — dome pressure is held by the turbine EHC / bypass valves.
                let nomP = supervisor.nominalPressureMPa
                autoCard(name: supervisor.hasPressurizer ? "PRESSURIZER AUTO" : "PRESS CTRL (EHC)",
                         enabled: supervisor.pzrAutoEnabled,
                         setpoint: String(format: "%.3f MPa", nomP),
                         measurement: String(format: "%6.3f MPa", supervisor.pressureMPa),
                         output: supervisor.hasPressurizer
                            ? (supervisor.pressureMPa < nomP ? "OUT HEATERS" : "OUT SPRAY")
                            : "OUT BYPASS VLV",
                         description: supervisor.hasPressurizer
                            ? "Heaters + spray hold primary pressure at \(String(format: "%.1f", nomP)) MPa"
                            : "Turbine inlet/bypass holds dome pressure at \(String(format: "%.1f", nomP)) MPa") {
                    supervisor.pzrAutoEnabled.toggle()
                }
                autoCard(name: "FEEDWATER AUTO",
                         enabled: supervisor.fwAutoEnabled,
                         setpoint: "1.000",
                         measurement: String(format: "%6.3f", supervisor.feedwaterInv),
                         output: String(format: "OUT %5.1f %%", supervisor.feedwaterValve * 100),
                         description: "Power feedforward + inventory trim on the FW reg valve") {
                    supervisor.fwAutoEnabled.toggle()
                }
                autoCard(name: "STARTUP SEQUENCER",
                         enabled: supervisor.autoStartup,
                         setpoint: "100.0 %",
                         measurement: String(format: "%6.2f %%", s.powerFraction * 100),
                         output: supervisor.startupPhase,
                         description: "Semi-auto ascension: rod withdrawal → T-avg ramp → at power [U]") {
                    supervisor.autoStartup.toggle()
                }
            }
            .padding(Theme.panelPadding)
        }
        .panelSurface()
        .frame(maxHeight: .infinity)
    }

    private func autoCard(name: String, enabled: Bool, setpoint: String,
                          measurement: String, output: String,
                          description: String, action: @escaping () -> Void) -> some View {
        Button(action: action) {
            VStack(alignment: .leading, spacing: 6) {
                HStack {
                    Circle().fill(enabled ? Theme.accent : Theme.border)
                        .frame(width: 8, height: 8)
                        .shadow(color: enabled ? Theme.accent.opacity(0.8) : .clear, radius: 4)
                    Text(name).font(Theme.readoutSm).foregroundStyle(enabled ? .white : Theme.textDim)
                    Spacer()
                    Text(enabled ? "AUTO" : "MANUAL").font(Theme.readoutSm)
                        .foregroundStyle(enabled ? Theme.accent : Theme.textDim)
                }
                Divider().background(Theme.sep)
                HStack {
                    VStack(alignment: .leading, spacing: 2) {
                        Text("SP").font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                        Text(setpoint).font(Theme.readoutSm).foregroundStyle(Theme.text)
                    }
                    Spacer()
                    VStack(alignment: .trailing, spacing: 2) {
                        Text("PV").font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                        Text(measurement).font(Theme.readoutSm).foregroundStyle(Theme.text)
                    }
                }
                Text(output)
                    .font(.system(size: 10, weight: .semibold, design: .monospaced))
                    .foregroundStyle(enabled ? Theme.accent : Theme.textDim.opacity(0.7))
                Text(description)
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundStyle(Theme.textDim)
                    .fixedSize(horizontal: false, vertical: true)
            }
            .padding(12)
            .frame(maxWidth: .infinity)
            .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .controlSurface(tint: enabled ? Theme.accent : nil)
    }
}
