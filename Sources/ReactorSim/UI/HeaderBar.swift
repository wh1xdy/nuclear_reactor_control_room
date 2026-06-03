// HeaderBar.swift — Top bar: facility ID, sim clock, speed, alarm indicator.

import SwiftUI

struct HeaderBar: View {
    let supervisor: PlantSupervisor
    let timeSpeed: Int
    let onSpeedCycle: () -> Void

    var body: some View {
        ZStack {
            Theme.panelHdr

            HStack(spacing: 16) {
                // Facility name
                VStack(alignment: .leading, spacing: 2) {
                    Text("UNIT 1")
                        .font(Theme.readoutLg)
                        .foregroundStyle(Theme.text)
                    Text("PRESSURISED WATER REACTOR PLANT")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.textDim)
                }
                .padding(.leading, 20)

                Spacer()

                // Sim clock
                VStack(alignment: .center, spacing: 1) {
                    Text("SIM TIME")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.textDim)
                    Text(clockString)
                        .font(Theme.readoutMd)
                        .foregroundStyle(Theme.accent)
                }

                // Speed control
                Button(action: onSpeedCycle) {
                    Text("×\(timeSpeed)")
                        .font(Theme.readoutMd)
                        .foregroundStyle(timeSpeed > 1 ? Theme.caution : Theme.textDim)
                        .frame(width: 60, height: 30)
                        .background(Theme.panel, in: RoundedRectangle(cornerRadius: 6))
                        .overlay(RoundedRectangle(cornerRadius: 6).stroke(Theme.border, lineWidth: 1))
                }
                .buttonStyle(.plain)

                Spacer()

                // Alarm/trip indicator
                AlarmIndicator(supervisor: supervisor)
                    .frame(width: 220, height: 40)
                    .padding(.trailing, 20)
            }
        }
        .overlay(alignment: .bottom) {
            Theme.border.frame(height: 1)
        }
    }

    private var clockString: String {
        let t = Int(supervisor.snapshot.time)
        let h = t / 3600; let m = (t % 3600) / 60; let s = t % 60
        return String(format: "%02d:%02d:%02d", h, m, s)
    }
}

private struct AlarmIndicator: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.5, on: .main, in: .common).autoconnect()

    var body: some View {
        let hasTrips  = !supervisor.trips.isEmpty
        let hasAlarms = !supervisor.alarms.isEmpty

        Group {
            if hasTrips {
                RoundedRectangle(cornerRadius: 6)
                    .fill(blink ? Theme.alarm : Theme.alarm.opacity(0.5))
                    .overlay(
                        Text("▶  REACTOR TRIP")
                            .font(Theme.readoutMd)
                            .foregroundStyle(.white)
                    )
            } else if hasAlarms {
                RoundedRectangle(cornerRadius: 6)
                    .fill(blink ? Theme.caution : Theme.caution.opacity(0.4))
                    .overlay(
                        Text("▶  ALARM")
                            .font(Theme.readoutMd)
                            .foregroundStyle(.white)
                    )
            } else {
                RoundedRectangle(cornerRadius: 6)
                    .fill(Color(r: 14, g: 36, b: 14))
                    .overlay(
                        RoundedRectangle(cornerRadius: 6)
                            .stroke(Color(r: 30, g: 80, b: 30), lineWidth: 1)
                    )
                    .overlay(
                        Text("ALL SYSTEMS NORMAL")
                            .font(Theme.readoutSm)
                            .foregroundStyle(Theme.normal)
                    )
            }
        }
        .onReceive(timer) { _ in blink.toggle() }
    }
}
