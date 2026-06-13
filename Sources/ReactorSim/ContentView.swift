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
                startPhysics()
                startKeyMonitor()
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
            HeaderBar(supervisor: supervisor, timeSpeed: timeSpeed,
                      skin: skin, onSpeedCycle: cycleSpeed, onSkinToggle: toggleSkin)
                .frame(height: Theme.headerHeight)

            TabBar(activeTab: $activeTab)
                .frame(height: Theme.tabHeight)

            HStack(spacing: 0) {
                ControlsPanel(supervisor: supervisor)
                    .frame(width: Theme.controlsWidth)

                tabContent
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
            }
            .frame(maxWidth: .infinity, maxHeight: .infinity)
        }
    }

    @ViewBuilder
    private var tabContent: some View {
        switch activeTab {
        case 0: OverviewTab(supervisor: supervisor)
        case 1: PrimaryTab(supervisor: supervisor)
        case 2: SecondaryTab(supervisor: supervisor)
        case 3: ReactivityTab(supervisor: supervisor)
        case 4: AlarmsTab(supervisor: supervisor)
        case 5: ICTab(supervisor: supervisor)
        default: EmptyView()
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
        case "c": supervisor.acknowledgeAllAlarms()
        case "l": supervisor.resetScram()
        case "r": supervisor = PlantSupervisor(); startPhysics()
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
