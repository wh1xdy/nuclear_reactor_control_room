// SystemBar.swift — slim global top strip shared by every console.
// Carries the nameplate, sim clock, plant alarm state, a quick console
// switcher, and the settings gear. Consoles render below it and bring their
// own domain chrome (mimic banner, panels, etc.).

import SwiftUI
import AppKit

struct SystemBar: View {
    let supervisor: PlantSupervisor
    let timeSpeed: Int
    @Binding var console: Console
    let onSpeedCycle: () -> Void
    let onOpenSettings: () -> Void

    var body: some View {
        HStack(spacing: 12) {
            // Nameplate
            HStack(spacing: 8) {
                Text("UNIT 1")
                    .font(.system(size: 13, weight: .bold, design: .monospaced))
                    .foregroundStyle(Theme.ink)
                Text("\(supervisor.reactorKind.rawValue.uppercased()) \(Int(supervisor.nominalMWt)) MWt")
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundStyle(Theme.textDim)
            }
            .padding(.leading, 14)

            // Plant alarm state
            SystemAlarmPill(supervisor: supervisor)

                // Global acknowledge — works in every console, not just the mimic.
            if !supervisor.alarms.isEmpty {
                Button { supervisor.acknowledgeAllAlarms() } label: {
                    Text("ACK ALL [C]")
                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                        .foregroundStyle(Theme.accent)
                        .padding(.horizontal, 10).padding(.vertical, 5)
                        .contentShape(Rectangle())
                }
                .buttonStyle(.plain)
                .controlSurface()
            }

            Spacer(minLength: 8)

            // Sim clock + speed
            SystemClock(supervisor: supervisor)
            Button(action: onSpeedCycle) {
                Text("×\(timeSpeed)")
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundStyle(timeSpeed > 1 ? Theme.accent : Theme.textDim)
                    .padding(.horizontal, 9).padding(.vertical, 3)
                    .contentShape(Rectangle())
            }
            .buttonStyle(.plain)
            .controlSurface(shape: .capsule)

            // Console quick switcher
            ConsoleSwitcher(console: $console)

            // Settings
            Button(action: onOpenSettings) {
                Image(systemName: "gearshape.fill")
                    .font(.system(size: 13))
                    .foregroundStyle(Theme.textHdr)
                    .padding(.horizontal, 8).padding(.vertical, 5)
                    .contentShape(Rectangle())
            }
            .buttonStyle(.plain)
            .controlSurface()
            .padding(.trailing, 14)
        }
        .frame(maxWidth: .infinity)
        .frame(height: 36)
        .background(Theme.bg)
        // The window uses a hidden title bar, so the OS title-bar double-click
        // has no target. This strip stands in for it: double-click empty bar
        // space to zoom the window (single taps still hit the buttons).
        .contentShape(Rectangle())
        .onTapGesture(count: 2) {
            (NSApp.keyWindow ?? NSApp.mainWindow ?? NSApp.windows.first)?.zoom(nil)
        }
    }
}

private struct ConsoleSwitcher: View {
    @Binding var console: Console
    var body: some View {
        HStack(spacing: 2) {
            ForEach(Console.allCases) { c in
                Button { console = c } label: {
                    Text(c.short)
                        .font(.system(size: 9, weight: .semibold, design: .monospaced))
                        .foregroundStyle(console == c ? Theme.ink : Theme.textDim)
                        .padding(.horizontal, 9).padding(.vertical, 5)
                        .contentShape(Rectangle())
                }
                .buttonStyle(.plain)
                .controlSurface(tint: console == c ? Theme.accent : nil)
            }
        }
    }
}

private struct SystemClock: View {
    let supervisor: PlantSupervisor
    var body: some View {
        HStack(spacing: 6) {
            Text("SIM")
                .font(.system(size: 8, weight: .medium, design: .monospaced))
                .foregroundStyle(Theme.textDim)
            Text(clockString)
                .font(.system(size: 13, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.ink)
        }
    }
    private var clockString: String {
        let t = Int(supervisor.snapshot.time)
        return String(format: "%02d:%02d:%02d", t / 3600, (t % 3600) / 60, t % 60)
    }
}

private struct SystemAlarmPill: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.5, on: .main, in: .common).autoconnect()
    var body: some View {
        let hasTrips  = !supervisor.trips.isEmpty
        let hasAlarms = !supervisor.alarms.isEmpty
        let (dot, text, tint): (Color, String, Color?) =
            hasTrips  ? (Theme.alarm,   "REACTOR TRIP",       Theme.alarm) :
            hasAlarms ? (Theme.caution, "ALARM ACTIVE",       Theme.caution) :
                        (Theme.ink.opacity(0.4), "ALL SYSTEMS NORMAL", nil)
        return HStack(spacing: 7) {
            Circle().fill(hasTrips && !blink ? dot.opacity(0.4) : dot)
                .frame(width: 8, height: 8)
            Text(text)
                .font(.system(size: 11, weight: hasTrips ? .bold : .medium, design: .monospaced))
                .foregroundStyle(tint == nil ? Theme.textDim : Theme.ink)
        }
        .padding(.horizontal, 12).padding(.vertical, 5)
        .controlSurface(tint: tint)
        .onReceive(timer) { _ in blink.toggle() }
    }
}
