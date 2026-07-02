// SettingsView.swift — picks the console layout, theme skin, and reactor type.
// Presented as a sheet from the system bar's gear.

import SwiftUI

struct SettingsView: View {
    @Binding var console: Console
    @Binding var skin: Skin
    @Binding var reactor: ReactorType
    var supervisor: PlantSupervisor? = nil   // scenario toggles + data export
    let onClose: () -> Void
    @State private var exportMsg: String? = nil

    var body: some View {
        VStack(alignment: .leading, spacing: 0) {
            // Title bar
            HStack {
                Text("CONTROL ROOM SETTINGS")
                    .font(.system(size: 13, weight: .bold, design: .monospaced))
                    .foregroundStyle(Theme.ink)
                    .tracking(1.5)
                Spacer()
                Button(action: onClose) {
                    Text("DONE")
                        .font(.system(size: 11, weight: .semibold, design: .monospaced))
                        .foregroundStyle(Theme.accent)
                        .padding(.horizontal, 14).padding(.vertical, 7)
                        .contentShape(Rectangle())
                }
                .buttonStyle(.plain)
                .controlSurface(tint: Theme.accent)
            }
            .padding(16)
            Divider().background(Theme.sep)

            ScrollView {
                VStack(alignment: .leading, spacing: 22) {
                    // ── Console layout ──
                    section("CONSOLE LAYOUT", "Choose how the control room is laid out. All layouts drive the same plant.")
                    VStack(spacing: 8) {
                        ForEach(Console.allCases) { c in
                            choiceRow(symbol: c.sfSymbol, title: c.label, blurb: c.blurb,
                                      selected: console == c, enabled: true) { console = c }
                        }
                    }

                    // ── Theme ──
                    section("THEME", "Visual styling applied across the chosen layout.")
                    HStack(spacing: 8) {
                        ForEach(Skin.allCases, id: \.self) { s in
                            skinChip(s)
                        }
                    }

                    // ── Reactor ──
                    section("REACTOR TYPE", "Switching rebuilds the plant model. Each shares the same six-group kinetics, decay heat and xenon, with kind-specific thermal-hydraulics and balance-of-plant.")
                    VStack(spacing: 8) {
                        ForEach(ReactorType.allCases) { r in
                            choiceRow(symbol: "atom", title: r.label,
                                      blurb: r.modelBlurb,
                                      selected: reactor == r, enabled: true) { reactor = r }
                        }
                    }

                    // ── Scenario & data (needs the live supervisor) ──
                    if let sup = supervisor {
                        section("SCENARIO & DATA", "Core-cycle scenario and run-data export.")
                        choiceRow(symbol: "waveform.path.ecg", title: "End-of-cycle core",
                                  blurb: "Strengthens the axial-xenon feedback toward the divergent EOL oscillation — hold ΔI in the flyspeck band with rods, at speed, or it runs away.",
                                  selected: sup.eolCore, enabled: true) { sup.eolCore.toggle() }
                        choiceRow(symbol: "square.and.arrow.down", title: "Export run history (CSV)",
                                  blurb: exportMsg ?? "Writes the slow historian (5 s cadence, up to ~6 sim-hours: power, temps, pressures, ΔI, DNBR, xenon, boron) to ~/Downloads.",
                                  selected: false, enabled: true) {
                            exportMsg = sup.exportCSV().map { "Saved: \($0)" } ?? "EXPORT FAILED"
                        }
                        choiceRow(symbol: "speaker.wave.2", title: "Control-room audio",
                                  blurb: "Procedural turbine hum (tracks shaft speed through a coastdown), annunciator chime on new alarms, breaker clunk.",
                                  selected: sup.soundEnabled, enabled: true) { sup.soundEnabled.toggle() }
                    }
                }
                .padding(16)
            }
        }
        .frame(width: 560, height: 620)
        .background(Theme.panel)
    }

    // MARK: — pieces

    private func section(_ title: String, _ sub: String) -> some View {
        VStack(alignment: .leading, spacing: 3) {
            Text(title)
                .font(.system(size: 10, weight: .bold, design: .monospaced))
                .foregroundStyle(Theme.textHdr).tracking(1.5)
            Text(sub)
                .font(.system(size: 10, design: .monospaced))
                .foregroundStyle(Theme.textDim)
        }
    }

    private func choiceRow(symbol: String, title: String, blurb: String,
                           selected: Bool, enabled: Bool, action: @escaping () -> Void) -> some View {
        Button(action: action) {
            HStack(spacing: 12) {
                Image(systemName: symbol)
                    .font(.system(size: 16))
                    .foregroundStyle(selected ? Theme.accent : Theme.textDim)
                    .frame(width: 26)
                VStack(alignment: .leading, spacing: 2) {
                    Text(title)
                        .font(.system(size: 12, weight: .semibold, design: .monospaced))
                        .foregroundStyle(enabled ? Theme.ink : Theme.textDim)
                    Text(blurb)
                        .font(.system(size: 10, design: .monospaced))
                        .foregroundStyle(Theme.textDim)
                        .fixedSize(horizontal: false, vertical: true)
                }
                Spacer()
                ZStack {
                    Circle().strokeBorder(selected ? Theme.accent : Theme.border, lineWidth: 1.5)
                        .frame(width: 16, height: 16)
                    if selected { Circle().fill(Theme.accent).frame(width: 8, height: 8) }
                }
            }
            .padding(12)
            .frame(maxWidth: .infinity, alignment: .leading)
            .contentShape(Rectangle())
            .opacity(enabled ? 1 : 0.5)
        }
        .buttonStyle(.plain)
        .controlSurface(tint: selected ? Theme.accent : nil)
    }

    private func skinChip(_ s: Skin) -> some View {
        let sub: String
        switch s {
        case .guided:        sub = "Liquid Glass · dark"
        case .authentic:     sub = "ISA-101 · light gray"
        case .authenticDark: sub = "ISA-101 · dark"
        }
        return Button { skin = s } label: {
            VStack(spacing: 3) {
                Text(s.label)
                    .font(.system(size: 10, weight: .semibold, design: .monospaced))
                    .foregroundStyle(skin == s ? Theme.ink : Theme.textDim)
                Text(sub)
                    .font(.system(size: 9, design: .monospaced))
                    .foregroundStyle(Theme.textDim)
            }
            .frame(maxWidth: .infinity)
            .padding(.vertical, 10)
            .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .controlSurface(tint: skin == s ? Theme.accent : nil)
    }
}
