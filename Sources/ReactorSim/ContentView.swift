// ContentView.swift
// Physics runs on a 60 Hz Timer (separate from rendering).
// TimelineView(.animation) only READS state — no mutations inside view body.

import SwiftUI
import AppKit

struct ContentView: View {
    @State private var supervisor   = PlantSupervisor()
    @State private var activeTab    = 0
    @State private var timeSpeed    = 1
    @State private var physicsTimer: Timer?
    @State private var keyMonitor:   Any?

    // Visual skin — persisted across launches. Mirrored into Theme.skin so the
    // static design tokens resolve correctly; .id(skin) rebuilds the tree on switch.
    @AppStorage("reactorSkin") private var skinRaw = Skin.guided.rawValue
    private var skin: Skin { Skin(rawValue: skinRaw) ?? .guided }
    private var skinBinding: Binding<Skin> {
        Binding(get: { skin }, set: { skinRaw = $0.rawValue })
    }

    // Selected operator console + reactor type (persisted).
    @AppStorage("console") private var consoleRaw = Console.mimic.rawValue
    private var console: Binding<Console> {
        Binding(get: { Console(rawValue: consoleRaw) ?? .mimic }, set: { consoleRaw = $0.rawValue })
    }
    @AppStorage("reactorType") private var reactorRaw = ReactorType.pwr.rawValue
    private var reactor: Binding<ReactorType> {
        Binding(get: { ReactorType(rawValue: reactorRaw) ?? .pwr }, set: { reactorRaw = $0.rawValue })
    }
    @State private var showSettings = false

    private let speeds = [1, 10, 60, 600]

    var body: some View {
        // Resolve the active skin for this build pass before any surface reads it.
        let _ = (Theme.skin = skin)
        // No root TimelineView — @Observable supervisor drives redraws via the Timer.
        // Scoped TimelineView lives inside PIDCanvas only (flow animation).
        layout
            .id(skin)                       // rebuild only the content subtree on skin switch
            .background(Theme.bg)
            .frame(minWidth: 1100, minHeight: 700)
            .onAppear {                     // stable identity — physics/keys start once
                // Honor a persisted reactor selection from a previous launch.
                if supervisor.reactorKind != reactor.wrappedValue.kind {
                    supervisor = PlantSupervisor(kind: reactor.wrappedValue.kind)
                }
                startPhysics()
                startKeyMonitor()
            }
            .onChange(of: reactorRaw) { _, newVal in
                // Switching reactor type rebuilds the plant model from scratch.
                supervisor = PlantSupervisor(kind: (ReactorType(rawValue: newVal) ?? .pwr).kind)
                startPhysics()
            }
            .onDisappear {
                physicsTimer?.invalidate()
                if let m = keyMonitor { NSEvent.removeMonitor(m) }
            }
    }

    private func toggleSkin() { skinRaw = skin.next.rawValue }

    // MARK: — Layout (read-only from supervisor)

    private var layout: some View {
        VStack(spacing: 0) {
            SystemBar(supervisor: supervisor, timeSpeed: timeSpeed, console: console,
                      onSpeedCycle: cycleSpeed, onOpenSettings: { showSettings = true })
                .frame(height: 36)
            Divider().background(Theme.sep)

            consoleView
                .frame(maxWidth: .infinity, maxHeight: .infinity)
        }
        .sheet(isPresented: $showSettings) {
            SettingsView(console: console, skin: skinBinding, reactor: reactor,
                         supervisor: supervisor, onClose: { showSettings = false })
        }
        .overlay {
            if supervisor.simPaused {
                PauseMenu(
                    speed: timeSpeed,
                    onResume:   { supervisor.simPaused = false },
                    onSettings: { supervisor.simPaused = false; showSettings = true },
                    onReset:    { supervisor = PlantSupervisor(kind: reactor.wrappedValue.kind)
                                  startPhysics() })
            }
        }
    }

    @ViewBuilder
    private var consoleView: some View {
        switch console.wrappedValue {
        case .mimic:       MimicConsole(supervisor: supervisor)
        case .workstation: WorkstationConsole(supervisor: supervisor)
        case .benchboard:  BenchboardConsole(supervisor: supervisor)
        case .dashboard:   DashboardConsole(supervisor: supervisor, activeTab: $activeTab)
        }
    }

    // MARK: — Physics timer (60 Hz, main thread — physics is < 0.5 ms/step)

    private func startPhysics() {
        physicsTimer?.invalidate()
        // .common mode keeps the timer firing during scroll/drag (default mode pauses).
        // Step inline — the timer already fires on the main thread. A Task{@MainActor}
        // hop would queue on the main dispatch queue, which stalls during event
        // tracking and would freeze physics during drags.
        let t = Timer(timeInterval: 1.0 / 60.0, repeats: true) { [self] _ in
            MainActor.assumeIsolated {
                guard !self.supervisor.simPaused else { return }   // pause freezes the physics clock (read live off the supervisor)
                self.supervisor.step(dt: 1.0 / 60.0 * Double(self.timeSpeed))
            }
        }
        RunLoop.main.add(t, forMode: .common)
        physicsTimer = t
    }

    private func cycleSpeed() {
        let idx = ((speeds.firstIndex(of: timeSpeed) ?? 0) + 1) % speeds.count
        timeSpeed = speeds[idx]
        // Restart timer so new speed takes effect cleanly
        startPhysics()
    }

    // MARK: — Keyboard (single monitor, cleaned up on disappear)

    private func startKeyMonitor() {
        if keyMonitor != nil { return }          // already registered
        keyMonitor = NSEvent.addLocalMonitorForEvents(matching: .keyDown) { [self] event in
            // Esc: close the core map first if open; otherwise toggle the pause
            // menu (space is SCRAM). Skip while the settings sheet is up.
            if event.keyCode == 53 && !showSettings {
                if supervisor.coreMapOpen { supervisor.coreMapOpen = false }
                else { supervisor.simPaused.toggle() }
                return nil
            }
            handleKey(event)
            return event
        }
    }

    private func handleKey(_ event: NSEvent) {
        let fast = event.modifierFlags.contains(.shift)
        let step = fast ? 0.05 : 0.01

        switch event.characters?.lowercased() {
        case "w": supervisor.rodPosition    = max(0, supervisor.rodPosition    - step)
        case "s": supervisor.rodPosition    = min(1, supervisor.rodPosition    + step)
        case "a": supervisor.primaryFlow    = max(0, supervisor.primaryFlow    - step)
        case "d": supervisor.primaryFlow    = min(1, supervisor.primaryFlow    + step)
        case "q": supervisor.turbineValve   = max(0, supervisor.turbineValve   - step)
        case "e": supervisor.turbineValve   = min(1, supervisor.turbineValve   + step)
        case "f": supervisor.feedwaterValve = min(1, supervisor.feedwaterValve + step)
        case "v": supervisor.feedwaterValve = max(0, supervisor.feedwaterValve - step)
        case "b": supervisor.borationRate   = min(1, supervisor.borationRate   + step)
        case "g": supervisor.borationRate   = max(0, supervisor.borationRate   - step)
        case "p": supervisor.startupPermit  = !supervisor.startupPermit
        case "t": supervisor.turbineTrip    = !supervisor.turbineTrip
        case "u": supervisor.autoStartup    = !supervisor.autoStartup
        case "o": supervisor.rodAutoEnabled = !supervisor.rodAutoEnabled
        case "z": supervisor.pumpDegraded   = !supervisor.pumpDegraded
        case "x": supervisor.feedwaterFault = !supervisor.feedwaterFault
        case "m": toggleSkin()
        case ",": showSettings.toggle()
        case "c": supervisor.acknowledgeAllAlarms()
        case "l": supervisor.resetScram()
        case "r": supervisor = PlantSupervisor(kind: reactor.wrappedValue.kind); startPhysics()
        case "+", "=": cycleSpeed()
        case "-":
            let idx = max(0, (speeds.firstIndex(of: timeSpeed) ?? 0) - 1)
            timeSpeed = speeds[idx]; startPhysics()
        default:
            switch event.keyCode {
            case 122: activeTab = 0   // F1
            case 120: activeTab = 1   // F2
            case 99:  activeTab = 2   // F3
            case 118: activeTab = 3   // F4
            case 96:  activeTab = 4   // F5
            case 97:  activeTab = 5   // F6
            case 49:                  // Space — SCRAM
                if supervisor.scrammed { supervisor.resetScram() }
                else { supervisor.triggerScram() }
            default: break
            }
        }
    }
}

// MARK: — Pause overlay (Esc)

private struct PauseMenu: View {
    let speed: Int
    let onResume:   () -> Void
    let onSettings: () -> Void
    let onReset:    () -> Void

    var body: some View {
        ZStack {
            Rectangle().fill(.black.opacity(0.55)).ignoresSafeArea()
                .contentShape(Rectangle())
                .onTapGesture(perform: onResume)          // click backdrop to resume
            VStack(spacing: 4) {
                Text("PAUSED").font(.system(size: 24, weight: .bold, design: .monospaced))
                    .foregroundStyle(Theme.textHdr).tracking(3)
                Text("simulation clock frozen · was ×\(speed)")
                    .font(.system(size: 10, design: .monospaced)).foregroundStyle(Theme.textDim)
                VStack(spacing: 8) {
                    row("RESUME", "esc", onResume)
                    row("SETTINGS", ",", onSettings)
                    row("RESET PLANT", "R", onReset)
                }
                .padding(.top, 16)
            }
            .padding(30)
            .frame(width: 320)
            .background(RoundedRectangle(cornerRadius: 16).fill(Theme.panel))
            .overlay(RoundedRectangle(cornerRadius: 16).stroke(Theme.border, lineWidth: 1))
        }
    }

    private func row(_ title: String, _ key: String, _ action: @escaping () -> Void) -> some View {
        Button(action: action) {
            HStack {
                Text(title).font(.system(size: 12, weight: .semibold, design: .monospaced))
                    .foregroundStyle(Theme.ink)
                Spacer()
                Text(key).font(.system(size: 10, design: .monospaced)).foregroundStyle(Theme.textDim)
            }
            .padding(.horizontal, 14).padding(.vertical, 11)
            .frame(maxWidth: .infinity)
            .contentShape(Rectangle())
            .background(RoundedRectangle(cornerRadius: 9).fill(Theme.dockTint))
            .overlay(RoundedRectangle(cornerRadius: 9).stroke(Theme.border.opacity(0.6), lineWidth: 1))
        }
        .buttonStyle(.plain)
    }
}
