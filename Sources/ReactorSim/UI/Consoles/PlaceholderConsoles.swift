// PlaceholderConsoles.swift — stubs for consoles still being built.
// Each is a real selectable layout that currently shows a "coming soon" card so
// the console switcher + settings work end-to-end while the real layout lands.

import SwiftUI

struct BenchboardConsole: View {
    let supervisor: PlantSupervisor
    var body: some View { ComingSoonConsole(console: .benchboard) }
}

private struct ComingSoonConsole: View {
    let console: Console
    var body: some View {
        VStack(spacing: 14) {
            Image(systemName: console.sfSymbol)
                .font(.system(size: 44))
                .foregroundStyle(Theme.textDim)
            Text(console.label)
                .font(.system(size: 16, weight: .bold, design: .monospaced))
                .foregroundStyle(Theme.ink).tracking(1.5)
            Text("UNDER CONSTRUCTION")
                .font(.system(size: 11, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.caution).tracking(2)
            Text(console.blurb)
                .font(.system(size: 11, design: .monospaced))
                .foregroundStyle(Theme.textDim)
                .multilineTextAlignment(.center)
                .frame(maxWidth: 420)
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        .background(Theme.bg)
    }
}
