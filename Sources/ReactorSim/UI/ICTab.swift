// ICTab.swift — F6 I&C: 2-of-3 protection channels + setpoints + auto controller status.

import SwiftUI

struct ICTab: View {
    let supervisor: PlantSupervisor

    var body: some View {
        VStack(spacing: 8) {
            ProtectionChannelsPanel(supervisor: supervisor)
            AutoControllersPanel(supervisor: supervisor)
        }
        .padding(8)
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
            let nomP = 15.5

            let params: [(String, Double, String, Double, Double?)] = [
                ("Neutron flux",    s.powerFraction * 100, "%",    120.0, 130.0),
                ("Fuel temperature",s.fuelTempK,           "K",   1400.0, nil),
                ("Reactor pressure",supervisor.pressureMPa,"MPa", nomP * 1.08, nomP * 1.1),
                ("Coolant temp",    s.coolantTempK,        "K",    610.0, 616.0),
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
                .padding(.vertical, 6)
                .background(Color(r: 14, g: 16, b: 20))

                Divider().background(Theme.sep)

                ForEach(Array(params.enumerated()), id: \.offset) { (i, param) in
                    let (name, value, unit, warnSP, tripSP) = param
                    let chVals = [channelVal(value, seed: i*3+1),
                                  channelVal(value, seed: i*3+2),
                                  channelVal(value, seed: i*3+3)]
                    let votes  = chVals.filter { tripSP != nil && $0 > tripSP! }.count
                    let warnOK = value < warnSP
                    let valCol = warnOK ? Theme.normal : (tripSP == nil || value < tripSP! ? Theme.caution : Theme.alarm)
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
                                return Theme.normal
                            }()
                            HStack(spacing: 4) {
                                Circle().fill(chCol).frame(width: 6, height: 6)
                                    .shadow(color: chCol.opacity(0.8), radius: 3)
                                Text(String(format: "%.2f", cv)).font(Theme.readoutSm).foregroundStyle(chCol)
                            }
                            .frame(maxWidth: .infinity)
                        }
                        VStack(alignment: .trailing, spacing: 2) {
                            Text(String(format: "Warn: %.2f %@", warnSP, unit))
                                .font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.caution)
                            if let t = tripSP {
                                Text(String(format: "Trip: %.2f %@", t, unit))
                                    .font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.alarm)
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
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
    }
}

// MARK: — Auto controller status
private struct AutoControllersPanel: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "AUTO CONTROLLERS")
            HStack(spacing: 12) {
                autoCard(name: "ROD AUTO-CONTROL", enabled: false,
                         setpoint: "550.0 K", measurement: String(format: "%.1f K", supervisor.snapshot.coolantTempK),
                         description: "Maintains avg coolant temp at 550 K by adjusting rod position")
                autoCard(name: "PRESSURIZER AUTO", enabled: false,
                         setpoint: "15.500 MPa", measurement: String(format: "%.3f MPa", supervisor.pressureMPa),
                         description: "Maintains primary pressure at 15.5 MPa via heater + spray")
                autoCard(name: "FEEDWATER AUTO", enabled: false,
                         setpoint: "1.000", measurement: String(format: "%.3f", supervisor.feedwaterInv),
                         description: "Maintains feedwater inventory via feedwater valve position")
            }
            .padding(10)
        }
        .background(Theme.panel)
        .clipShape(RoundedRectangle(cornerRadius: Theme.cornerRadius))
        .overlay(RoundedRectangle(cornerRadius: Theme.cornerRadius).stroke(Theme.border, lineWidth: 1))
        .overlay(CornerBrackets())
        .frame(maxHeight: .infinity)
    }

    private func autoCard(name: String, enabled: Bool, setpoint: String,
                          measurement: String, description: String) -> some View {
        VStack(alignment: .leading, spacing: 6) {
            HStack {
                Circle().fill(enabled ? Theme.normal : Theme.border)
                    .frame(width: 8, height: 8)
                    .shadow(color: enabled ? Theme.normal.opacity(0.8) : .clear, radius: 4)
                Text(name).font(Theme.readoutSm).foregroundStyle(enabled ? Theme.normal : Theme.textDim)
                Spacer()
                Text(enabled ? "AUTO" : "MANUAL").font(Theme.readoutSm)
                    .foregroundStyle(enabled ? Theme.accent : Theme.textDim)
            }
            Divider().background(Theme.sep)
            HStack {
                VStack(alignment: .leading, spacing: 2) {
                    Text("SP").font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                    Text(setpoint).font(Theme.readoutSm).foregroundStyle(Theme.accent)
                }
                Spacer()
                VStack(alignment: .trailing, spacing: 2) {
                    Text("PV").font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
                    Text(measurement).font(Theme.readoutSm).foregroundStyle(Theme.text)
                }
            }
            Text(description)
                .font(.system(size: 9, design: .monospaced))
                .foregroundStyle(Theme.textDim)
                .fixedSize(horizontal: false, vertical: true)
        }
        .padding(10)
        .frame(maxWidth: .infinity)
        .background(Color(r: 10, g: 12, b: 16), in: RoundedRectangle(cornerRadius: 6))
        .overlay(RoundedRectangle(cornerRadius: 6).stroke(Theme.border, lineWidth: 1))
    }
}
