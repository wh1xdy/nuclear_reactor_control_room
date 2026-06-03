// TabBar.swift — Six-tab bar with electric-blue active underline.

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
        HStack(spacing: 0) {
            ForEach(Array(tabs.enumerated()), id: \.offset) { i, tab in
                TabButton(fkey: tab.0, name: tab.1, active: activeTab == i) {
                    activeTab = i
                }
            }
        }
        .background(Color(r: 12, g: 14, b: 17))
        .overlay(alignment: .bottom) { Theme.border.frame(height: 1) }
    }
}

private struct TabButton: View {
    let fkey: String
    let name: String
    let active: Bool
    let action: () -> Void

    var body: some View {
        Button(action: action) {
            HStack(spacing: 4) {
                Text(fkey)
                    .font(Theme.readout)
                    .foregroundStyle(active ? Theme.accent : Theme.textDim)
                Text(name)
                    .font(Theme.readout)
                    .foregroundStyle(active ? Theme.text : Theme.textDim)
            }
            .frame(maxWidth: .infinity)
            .frame(height: Theme.tabHeight)
            .background(active ? Color(r: 22, g: 28, b: 36) : Color.clear)
            .overlay(alignment: .bottom) {
                if active {
                    Theme.accent.frame(height: 3)
                }
            }
            .overlay {
                if active {
                    RoundedRectangle(cornerRadius: 0)
                        .stroke(Theme.border, lineWidth: 0.5)
                }
            }
        }
        .buttonStyle(.plain)
    }
}
