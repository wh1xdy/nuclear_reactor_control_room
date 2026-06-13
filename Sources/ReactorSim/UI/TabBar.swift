// TabBar.swift — Liquid Glass tab bar with morphing-merge behaviour.

import SwiftUI

private let tabs: [(String, String)] = [
    ("F1", "OVERVIEW"),
    ("F2", "PRIMARY"),
    ("F3", "SECONDARY"),
    ("F4", "REACTIVITY"),
    ("F5", "ALARMS"),
    ("F6", "I&C"),
]

struct TabBar: View {
    @Binding var activeTab: Int

    var body: some View {
        // Tabs inside a GlassEffectContainer so they merge when adjacent —
        // pure glass, no background strip, no separator rule.
        GlassEffectContainer(spacing: 2) {
            HStack(spacing: 4) {
                ForEach(Array(tabs.enumerated()), id: \.offset) { i, tab in
                    TabButton(fkey: tab.0, name: tab.1, active: activeTab == i) {
                        withAnimation(.spring(duration: 0.3)) { activeTab = i }
                    }
                }
            }
            .padding(.horizontal, 12)
        }
        .frame(height: Theme.tabHeight)
    }
}

private struct TabButton: View {
    let fkey: String
    let name: String
    let active: Bool
    let action: () -> Void

    var body: some View {
        Button(action: action) {
            HStack(spacing: 5) {
                Text(fkey)
                    .font(Theme.readoutSm)
                    .foregroundStyle(active ? Theme.accent : Theme.textDim)
                Text(name)
                    .font(Theme.readoutSm)
                    .foregroundStyle(active ? .white : Theme.textDim)
            }
            .padding(.horizontal, 16)
            .padding(.vertical, 7)
            .frame(maxWidth: .infinity)
            .contentShape(Rectangle())   // full rect is hittable, not just the text glyphs
        }
        .buttonStyle(.plain)
        .controlSurface(tint: active ? Theme.accent : nil)
    }
}
