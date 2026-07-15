// StartMenu.swift — launch screen. Pick the reactor and how the plant starts:
// AT POWER (drop into steady full-power operation) or COLD STARTUP (boot in
// hot standby, subcritical — approach criticality, roll the turbine, sync the
// generator, and load the unit yourself).

import SwiftUI

struct StartMenu: View {
    @Binding var reactor: ReactorType
    @State private var mode: PlantSupervisor.StartMode = .atPower
    let onStart: (ReactorType, PlantSupervisor.StartMode) -> Void

    var body: some View {
        ZStack {
            Theme.bg.ignoresSafeArea()
            VStack(spacing: 26) {
                VStack(spacing: 6) {
                    Text("REACTORSIM")
                        .font(.system(size: 34, weight: .bold, design: .monospaced))
                        .foregroundStyle(Theme.ink).tracking(6)
                    Text("nuclear plant control-room simulator")
                        .font(.system(size: 11, design: .monospaced))
                        .foregroundStyle(Theme.textDim).tracking(1)
                }
                .padding(.bottom, 4)

                section("REACTOR") {
                    HStack(spacing: 10) {
                        ForEach(ReactorType.allCases) { t in
                            pill(t.short, sub: subLabel(t), on: reactor == t) { reactor = t }
                        }
                    }
                }

                section("START CONDITION") {
                    HStack(spacing: 10) {
                        modeCard("AT POWER", "Steady 100 % — jump straight into operating the plant.",
                                 on: mode == .atPower) { mode = .atPower }
                        modeCard("COLD STARTUP", "Hot standby, subcritical. Pull rods to critical on the 1/M plot, roll the turbine, sync the generator, load up.",
                                 on: mode == .coldStartup) { mode = .coldStartup }
                    }
                }

                Button { onStart(reactor, mode) } label: {
                    Text("START  ▸")
                        .font(.system(size: 15, weight: .bold, design: .monospaced))
                        .foregroundStyle(Theme.isFlat ? Theme.ink : .white)
                        .tracking(2)
                        .padding(.horizontal, 40).padding(.vertical, 12)
                        .background(RoundedRectangle(cornerRadius: 10).fill(Theme.accent))
                        .contentShape(Rectangle())
                }
                .buttonStyle(.plain)
                .padding(.top, 6)
            }
            .frame(maxWidth: 620)
            .padding(40)
        }
    }

    private func subLabel(_ t: ReactorType) -> String {
        switch t { case .pwr: return "3000 MWt"; case .bwr: return "3000 MWt"; case .smr: return "200 MWt" }
    }

    private func section<Content: View>(_ title: String, @ViewBuilder _ content: () -> Content) -> some View {
        VStack(alignment: .leading, spacing: 10) {
            Text(title).font(.system(size: 10, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.textDim).tracking(2)
            content()
        }
    }

    private func pill(_ label: String, sub: String, on: Bool, _ act: @escaping () -> Void) -> some View {
        Button(action: act) {
            VStack(spacing: 3) {
                Text(label).font(.system(size: 15, weight: .semibold, design: .monospaced))
                    .foregroundStyle(on ? (Theme.isFlat ? Theme.ink : .white) : Theme.ink)
                Text(sub).font(.system(size: 8, design: .monospaced))
                    .foregroundStyle(on ? (Theme.isFlat ? Theme.ink : .white).opacity(0.8) : Theme.textDim)
            }
            .frame(maxWidth: .infinity).padding(.vertical, 14)
            .background(RoundedRectangle(cornerRadius: 10)
                .fill(on ? Theme.accent : Theme.panel))
            .overlay(RoundedRectangle(cornerRadius: 10)
                .stroke(Theme.ink.opacity(on ? 0 : 0.12), lineWidth: 1))
            .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
    }

    private func modeCard(_ title: String, _ desc: String, on: Bool, _ act: @escaping () -> Void) -> some View {
        Button(action: act) {
            VStack(alignment: .leading, spacing: 6) {
                Text(title).font(.system(size: 13, weight: .bold, design: .monospaced))
                    .foregroundStyle(on ? (Theme.isFlat ? Theme.ink : .white) : Theme.ink)
                Text(desc).font(.system(size: 9, design: .monospaced))
                    .foregroundStyle(on ? (Theme.isFlat ? Theme.ink : .white).opacity(0.85) : Theme.textDim)
                    .fixedSize(horizontal: false, vertical: true)
            }
            .frame(maxWidth: .infinity, minHeight: 92, alignment: .topLeading)
            .padding(14)
            .background(RoundedRectangle(cornerRadius: 10).fill(on ? Theme.accent : Theme.panel))
            .overlay(RoundedRectangle(cornerRadius: 10).stroke(Theme.ink.opacity(on ? 0 : 0.12), lineWidth: 1))
            .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
    }
}
