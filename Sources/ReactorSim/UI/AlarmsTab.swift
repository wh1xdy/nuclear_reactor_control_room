// AlarmsTab.swift — F5 Alarms: annunciator tile grid.

import SwiftUI

struct AlarmsTab: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.5, on: .main, in: .common).autoconnect()

    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "ALARM ANNUNCIATOR")
            if supervisor.alarms.isEmpty {
                VStack {
                    Spacer()
                    HStack(spacing: 12) {
                        Circle().fill(Theme.normal).frame(width: 10, height: 10)
                            .shadow(color: Theme.normal.opacity(0.8), radius: 6)
                        Text("All systems nominal")
                            .font(Theme.readoutMd)
                            .foregroundStyle(Theme.normal)
                    }
                    Spacer()
                }
            } else {
                LazyVGrid(columns: [GridItem(.adaptive(minimum: 200, maximum: 300))], spacing: 6) {
                    ForEach(supervisor.alarms) { alarm in
                        AlarmTile(alarm: alarm, blink: blink)
                    }
                }
                .padding(10)

                Divider().background(Theme.sep)

                Button("ACKNOWLEDGE ALL [C]") { supervisor.acknowledgeAllAlarms() }
                    .font(Theme.readoutSm)
                    .foregroundStyle(Theme.textDim)
                    .buttonStyle(.plain)
                    .padding(8)
            }
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        .background(Theme.panel)
        .overlay(CornerBrackets())
        .padding(8)
        .onReceive(timer) { _ in blink.toggle() }
    }
}

private struct AlarmTile: View {
    let alarm: ReactorAlarm
    let blink: Bool

    private var bg: Color {
        if alarm.isTrip  { return blink ? Color(r: 100, g: 10, b: 10) : Color(r: 60, g: 8, b: 8) }
        if alarm.state == "unack" { return blink ? Color(r: 80, g: 50, b: 0) : Color(r: 50, g: 30, b: 0) }
        return Color(r: 30, g: 40, b: 20)
    }

    private var borderColor: Color {
        if alarm.isTrip  { return Theme.alarm }
        if alarm.state == "unack" { return Theme.caution }
        return Theme.normal
    }

    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack {
                Text(alarm.isTrip ? "TRIP" : alarm.state == "unack" ? "WARN" : "ACK")
                    .font(Theme.readoutSm)
                    .foregroundStyle(borderColor)
                if alarm.isFirstOut {
                    Text("FO")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.caution)
                        .padding(.horizontal, 4)
                        .background(Theme.caution.opacity(0.2), in: RoundedRectangle(cornerRadius: 2))
                }
                Spacer()
                Text(String(format: "T+%.0fs", alarm.simTime))
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundStyle(Theme.textDim)
            }
            Text(alarm.message)
                .font(Theme.readoutSm)
                .foregroundStyle(Theme.text)
                .lineLimit(2)
                .minimumScaleFactor(0.8)
        }
        .padding(8)
        .background(bg, in: RoundedRectangle(cornerRadius: 5))
        .overlay(RoundedRectangle(cornerRadius: 5).stroke(borderColor, lineWidth: 1))
    }
}
