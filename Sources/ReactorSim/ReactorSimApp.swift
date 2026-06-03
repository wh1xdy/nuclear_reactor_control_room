// ReactorSimApp.swift — App entry point with proper macOS activation.

import SwiftUI
import AppKit

@main
struct ReactorSimApp: App {
    @NSApplicationDelegateAdaptor(AppDelegate.self) var appDelegate

    var body: some Scene {
        WindowGroup {
            ContentView()
                .preferredColorScheme(.dark)
        }
        .windowStyle(.hiddenTitleBar)
        .defaultSize(width: 1440, height: 900)
        .windowResizability(.contentSize)
        .commands {
            // Remove default menu items that don't apply
            CommandGroup(replacing: .newItem) {}
        }
    }
}

// MARK: — AppDelegate: forces window to front and sets up dock presence
final class AppDelegate: NSObject, NSApplicationDelegate {
    func applicationDidFinishLaunching(_ notification: Notification) {
        NSApp.setActivationPolicy(.regular)
        NSApp.activate(ignoringOtherApps: true)
        // Ensure window comes to front
        DispatchQueue.main.async {
            NSApp.windows.first?.makeKeyAndOrderFront(nil)
        }
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool {
        true
    }
}
