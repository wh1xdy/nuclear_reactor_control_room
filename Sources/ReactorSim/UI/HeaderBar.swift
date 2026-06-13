// HeaderBar.swift — Liquid Glass header strip.

import SwiftUI

struct HeaderBar: View {
    let supervisor: PlantSupervisor
    let timeSpeed: Int
    let skin: Skin
    let onSpeedCycle: () -> Void
    let onSkinToggle: () -> Void

    var body: some View {
        ZStack {
            // Transparent — inherits window background so no seam at top
            Color.clear

            HStack(spacing: 0) {
                // Facility nameplate — utility bar scale, not a billboard
                HStack(spacing: 10) {
                    Text("UNIT 1")
                        .font(.system(size: 14, weight: .bold, design: .monospaced))
                        .foregroundStyle(.white)
                    Text("PWR  3000 MWt")
                        .font(.system(size: 10, design: .monospaced))
                        .foregroundStyle(Theme.textDim)
                }
                .padding(.leading, 18)

                // Skin switch — GUIDED (Liquid Glass) ⇄ AUTHENTIC (flat DCS)
                Button(action: onSkinToggle) {
                    HStack(spacing: 6) {
                        Image(systemName: skin == .guided ? "sparkles" : "square.grid.3x3.fill")
                            .font(.system(size: 9))
                            .foregroundStyle(Theme.accent)
                        Text(skin.label)
                            .font(.system(size: 10, weight: .semibold, design: .monospaced))
                            .foregroundStyle(.white)
                        Text("[M]")
                            .font(.system(size: 9, design: .monospaced))
                            .foregroundStyle(Theme.textDim)
                    }
                    .padding(.horizontal, 10)
                    .padding(.vertical, 4)
                    .contentShape(Rectangle())
                }
                .buttonStyle(.plain)
                .controlSurface(shape: .capsule)
                .padding(.leading, 14)

                Spacer()

                // Sim clock + speed — glass pill
                HStack(spacing: 12) {
                    // Leaf view: isolates the 60 Hz snapshot.time re-evaluation
                    // so the glass pill subtree isn't rebuilt every physics tick.
                    SimClock(supervisor: supervisor)

                    Button(action: onSpeedCycle) {
                        Text("×\(timeSpeed)")
                            .font(.system(size: 11, weight: .semibold, design: .monospaced))
                            // Accelerated time is an active selection, not a warning.
                            .foregroundStyle(timeSpeed > 1 ? Theme.accent : Theme.textDim)
                            .padding(.horizontal, 10)
                            .padding(.vertical, 4)
                            .contentShape(Rectangle())
                    }
                    .buttonStyle(.plain)
                    .controlSurface(shape: .capsule)
                }
                .padding(.horizontal, 14)
                .padding(.vertical, 5)
                .controlSurface()
                .padding(.horizontal, 16)

                Spacer()

                // Alarm/trip indicator
                AlarmIndicator(supervisor: supervisor)
                    .frame(width: 210, height: 28)
                    .padding(.trailing, 16)
            }
        }
    }
}

private struct SimClock: View {
    let supervisor: PlantSupervisor

    var body: some View {
        HStack(spacing: 8) {
            Text("SIM TIME")
                .font(.system(size: 8, weight: .medium, design: .monospaced))
                .foregroundStyle(Theme.textDim)
            Text(clockString)
                .font(.system(size: 15, weight: .semibold, design: .monospaced))
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
                .controlSurface(tint: Theme.alarm)
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
                .controlSurface(tint: Theme.caution)
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
                .controlSurface()
            }
        }
        .frame(maxWidth: .infinity)
        .onReceive(timer) { _ in blink.toggle() }
    }
}
