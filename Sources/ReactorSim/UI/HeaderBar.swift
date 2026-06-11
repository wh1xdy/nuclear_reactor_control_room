// HeaderBar.swift — Liquid Glass header strip.

import SwiftUI

struct HeaderBar: View {
    let supervisor: PlantSupervisor
    let timeSpeed: Int
    let onSpeedCycle: () -> Void

    var body: some View {
        ZStack {
            // Transparent — inherits window background so no seam at top
            Color.clear

            HStack(spacing: 0) {
                // Facility nameplate
                VStack(alignment: .leading, spacing: 1) {
                    Text("UNIT 1")
                        .font(.system(size: 22, weight: .bold, design: .monospaced))
                        .foregroundStyle(.white)
                    Text("PRESSURISED WATER REACTOR PLANT")
                        .font(Theme.readoutSm)
                        .foregroundStyle(Theme.textDim)
                }
                .padding(.leading, 22)
                .padding(.vertical, 10)

                Spacer()

                // Sim clock + speed — glass pill
                HStack(spacing: 16) {
                    // Leaf view: isolates the 60 Hz snapshot.time re-evaluation
                    // so the glass pill subtree isn't rebuilt every physics tick.
                    SimClock(supervisor: supervisor)

                    Button(action: onSpeedCycle) {
                        Text("×\(timeSpeed)")
                            .font(.system(size: 13, weight: .semibold, design: .monospaced))
                            // Accelerated time is an active selection, not a warning.
                            .foregroundStyle(timeSpeed > 1 ? Theme.accent : Theme.textDim)
                            .padding(.horizontal, 14)
                            .padding(.vertical, 6)
                    }
                    .buttonStyle(.plain)
                    .glassEffect(.regular.interactive(), in: Capsule())
                }
                .padding(.horizontal, 24)
                .padding(.vertical, 8)
                .glassEffect(.regular, in: .rect(cornerRadius: Theme.controlRadius, style: .continuous))
                .padding(.horizontal, 20)

                Spacer()

                // Alarm/trip indicator
                AlarmIndicator(supervisor: supervisor)
                    .frame(width: 230, height: 36)
                    .padding(.trailing, 20)
            }
        }
    }
}

private struct SimClock: View {
    let supervisor: PlantSupervisor

    var body: some View {
        VStack(spacing: 1) {
            Text("SIM TIME")
                .font(.system(size: 9, weight: .medium, design: .monospaced))
                .foregroundStyle(Theme.textDim)
            Text(clockString)
                .font(.system(size: 20, weight: .semibold, design: .monospaced))
                .foregroundStyle(.white)
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
                HStack(spacing: 8) {
                    Circle().fill(blink ? Theme.alarm : Theme.alarm.opacity(0.5))
                        .frame(width: 8, height: 8)
                        .shadow(color: Theme.alarm, radius: blink ? 6 : 2)
                    Text("REACTOR TRIP")
                        .font(.system(size: 13, weight: .bold, design: .monospaced))
                        .foregroundStyle(.white)
                }
                .padding(.horizontal, 16)
                .glassEffect(.regular.tint(Theme.alarm.opacity(0.3)).interactive(),
                             in: .rect(cornerRadius: Theme.controlRadius, style: .continuous))
            } else if hasAlarms {
                HStack(spacing: 8) {
                    Circle().fill(blink ? Theme.caution : Theme.caution.opacity(0.5))
                        .frame(width: 8, height: 8)
                        .shadow(color: Theme.caution, radius: 4)
                    Text("ALARM ACTIVE")
                        .font(.system(size: 13, weight: .semibold, design: .monospaced))
                        .foregroundStyle(.white)
                }
                .padding(.horizontal, 16)
                .glassEffect(.regular.tint(Theme.caution.opacity(0.2)),
                             in: .rect(cornerRadius: Theme.controlRadius, style: .continuous))
            } else {
                // Normal state stays quiet — no ISA green on a non-alarm condition.
                HStack(spacing: 8) {
                    Circle().fill(Color.white.opacity(0.35))
                        .frame(width: 7, height: 7)
                    Text("ALL SYSTEMS NORMAL")
                        .font(.system(size: 11, weight: .medium, design: .monospaced))
                        .foregroundStyle(Theme.textDim)
                }
                .padding(.horizontal, 16)
                .glassEffect(.regular,
                             in: .rect(cornerRadius: Theme.controlRadius, style: .continuous))
            }
        }
        .frame(maxWidth: .infinity)
        .onReceive(timer) { _ in blink.toggle() }
    }
}
