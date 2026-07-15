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
        // contentMinSize (not contentSize) so the window can be freely resized,
        // zoomed with the green button, and double-clicked to fill the screen —
        // contentSize locked it to the layout size, defeating maximize.
        .windowResizability(.contentMinSize)
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
        // Ensure the window comes to front and opens filling the screen (a
        // control room wants the whole display). setFrame to the visible frame
        // rather than zoom(nil) so it's deterministic regardless of the restored
        // autosave frame; the user can still resize/zoom afterwards.
        DispatchQueue.main.async {
            guard let win = NSApp.windows.first else { return }
            win.makeKeyAndOrderFront(nil)
            win.collectionBehavior.insert(.fullScreenPrimary)   // green button → full screen too
            if let vis = win.screen?.visibleFrame { win.setFrame(vis, display: true, animate: false) }
        }
    }

    func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool {
        true
    }
}
