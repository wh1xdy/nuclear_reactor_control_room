// AlarmsTab.swift — F5 Alarms: fixed annunciator matrix + chronological summary.
// The matrix shows EVERY defined alarm window at all times — dark when clear,
// lit when active, blinking until acknowledged. An empty alarm screen is an
// amateur tell; a dark-but-present matrix is how real control rooms look.

import SwiftUI

private struct AnnWindow: Identifiable {
    let id: String
    let label: String
    let trip: Bool          // trip-class (red) vs warning-class (amber)
    let active: Bool
    let unack: Bool
}

struct AlarmsTab: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.5, on: .main, in: .common).autoconnect()

    // Fixed 5×3 window layout — annunciators don't rearrange themselves.
    private var windows: [AnnWindow] {
        func win(_ id: String, _ label: String, trip: Bool) -> AnnWindow {
            let a = supervisor.alarms.first { $0.id == id }
            return AnnWindow(id: id, label: label, trip: trip,
                             active: a != nil, unack: a?.state == "unack")
        }
        func state(_ id: String, _ label: String, _ active: Bool, trip: Bool) -> AnnWindow {
            AnnWindow(id: id, label: label, trip: trip, active: active, unack: false)
        }
        return [
            state("RX_TRIP", "REACTOR\nTRIP", supervisor.scrammed, trip: true),
            win("HIGH_FLUX",   "HI NEUTRON\nFLUX",   trip: true),
            win("HIGH_FUEL_T", "HI FUEL\nTEMP",      trip: true),
            win("HIGH_PRESS",  "HI RCS\nPRESS",      trip: true),
            win("HIGH_COOL_T", "HI RCS\nT-AVG",      trip: true),

            win("ECCS_ACT",    "ECCS\nACTUATION",    trip: true),
            win("PORV_OPEN",   "PZR PORV\nOPEN",     trip: false),
            win("WARN_PRESS",  "PZR PRESS\nHI WARN", trip: false),
            win("LOW_FEED",    "LO FW\nINVENTORY",   trip: false),
            win("XE_TRANSIENT","XENON\nTRANSIENT",   trip: false),

            state("TBN_TRIP",  "TURBINE\nTRIP",      supervisor.turbineTrip, trip: false),
            state("STM_DUMP",  "STEAM DUMP\nOPEN",   supervisor.steamDumpValve > 0.001, trip: false),
            state("RCP_DEG",   "RCP\nDEGRADED",      supervisor.pumpDegraded, trip: false),
            state("FW_FAULT",  "FEEDWATER\nFAULT",   supervisor.feedwaterFault, trip: false),
            state("SPARE_1",   "SPARE",              false, trip: false),
        ]
    }

    var body: some View {
        VStack(spacing: 12) {
            // ── Annunciator matrix ────────────────────────────────────────────
            VStack(spacing: 0) {
                PanelHeader(title: "ANNUNCIATOR")
                let cols = [GridItem](repeating: GridItem(.flexible(), spacing: 8), count: 5)
                LazyVGrid(columns: cols, spacing: 8) {
                    ForEach(windows) { w in
                        AnnunciatorTile(window: w, blink: blink)
                    }
                }
                .padding(Theme.panelPadding)
            }
            .glassEffect(.regular, in: .rect(cornerRadius: Theme.panelRadius, style: .continuous))

            // ── Chronological alarm summary ───────────────────────────────────
            VStack(spacing: 0) {
                PanelHeader(title: "ALARM SUMMARY  —  CHRONOLOGICAL")
                // Column headers
                summaryRow("TIME", "PRI", "ST", "FO", "MESSAGE",
                           color: Theme.textDim, header: true)
                Divider().background(Theme.sep)
                if supervisor.alarms.isEmpty {
                    summaryRow("—", "—", "—", "—", "NO ACTIVE ALARMS",
                               color: Theme.textDim.opacity(0.6), header: false)
                } else {
                    ScrollView(.vertical, showsIndicators: false) {
                        VStack(spacing: 0) {
                            ForEach(supervisor.alarms) { alarm in
                                summaryRow(String(format: "T+%05.0fs", alarm.simTime),
                                           "P\(alarm.priority)",
                                           alarm.state == "unack" ? "UNACK" : "ACK",
                                           alarm.isFirstOut ? "FO" : "",
                                           alarm.message,
                                           color: alarm.isTrip ? Theme.alarm :
                                                  alarm.state == "unack" ? Theme.caution : Theme.text,
                                           header: false)
                                Divider().background(Theme.sep.opacity(0.5))
                            }
                        }
                    }
                }
                Spacer(minLength: 0)
                HStack {
                    Spacer()
                    Button {
                        supervisor.acknowledgeAllAlarms()
                    } label: {
                        Text("ACKNOWLEDGE ALL [C]")
                            .font(Theme.readoutSm)
                            .foregroundStyle(Theme.accent)
                            .padding(10)
                            .contentShape(Rectangle())
                    }
                    .buttonStyle(.plain)
                }
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity)
            .glassEffect(.regular, in: .rect(cornerRadius: Theme.panelRadius, style: .continuous))
        }
        .padding(12)
        .onReceive(timer) { _ in blink.toggle() }
    }

    private func summaryRow(_ time: String, _ pri: String, _ st: String,
                            _ fo: String, _ msg: String,
                            color: Color, header: Bool) -> some View {
        HStack(spacing: 0) {
            Text(time).frame(width: 90, alignment: .leading)
            Text(pri).frame(width: 44, alignment: .leading)
            Text(st).frame(width: 64, alignment: .leading)
            Text(fo).frame(width: 36, alignment: .leading)
                .foregroundStyle(fo == "FO" ? Theme.caution : color)
            Text(msg).frame(maxWidth: .infinity, alignment: .leading)
        }
        .font(.system(size: header ? 9 : 11, design: .monospaced))
        .foregroundStyle(color)
        .padding(.horizontal, Theme.panelPadding)
        .padding(.vertical, 5)
    }
}

// MARK: — Annunciator tile
private struct AnnunciatorTile: View {
    let window: AnnWindow
    let blink: Bool

    private var fill: Color {
        guard window.active else { return Color.white.opacity(0.035) }
        let base = window.trip ? Theme.alarm : Theme.caution
        if window.unack { return base.opacity(blink ? 0.95 : 0.45) }
        return base.opacity(0.75)
    }

    private var textColor: Color {
        if window.active { return .white }
        return Color.white.opacity(window.label == "SPARE" ? 0.15 : 0.30)
    }

    var body: some View {
        Text(window.label)
            .font(.system(size: 9, weight: window.active ? .bold : .medium,
                          design: .monospaced))
            .multilineTextAlignment(.center)
            .lineLimit(2)
            .minimumScaleFactor(0.7)
            .foregroundStyle(textColor)
            .frame(maxWidth: .infinity)
            .frame(height: 46)
            .background(fill, in: .rect(cornerRadius: 8, style: .continuous))
    }
}
