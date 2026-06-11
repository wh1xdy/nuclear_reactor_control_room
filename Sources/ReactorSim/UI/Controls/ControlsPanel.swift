// ControlsPanel.swift — Left panel with Liquid Glass throughout.

import SwiftUI

struct ControlsPanel: View {
    let supervisor: PlantSupervisor

    var body: some View {
        ScrollView(.vertical, showsIndicators: false) {
            VStack(alignment: .leading, spacing: 0) {
                // Section header
                sectionLabel("CONTROLS")

                // Restrained palette: all sliders use the single electric-blue accent.
                // Rod demand reads in steps withdrawn (0–228), real CRDM convention.
                VStack(spacing: 18) {
                    DCSSlider(label: "ROD DEMAND",
                              value: Binding(get: { supervisor.rodPosition },
                                            set: { supervisor.rodPosition = $0 }),
                              displayStr: { "\(Int((228 * (1 - $0)).rounded())) SWD" })
                    DCSSlider(label: "RCS FLOW",
                              value: Binding(get: { supervisor.primaryFlow },
                                            set: { supervisor.primaryFlow = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" })
                    DCSSlider(label: "TBN GOV VALVE",
                              value: Binding(get: { supervisor.turbineValve },
                                            set: { supervisor.turbineValve = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" })
                    DCSSlider(label: "FW REG VALVE",
                              value: Binding(get: { supervisor.feedwaterValve },
                                            set: { supervisor.feedwaterValve = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" })
                    DCSSlider(label: "BORATION RATE",
                              value: Binding(get: { supervisor.borationRate },
                                            set: { supervisor.borationRate = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" })
                }
                .padding(.horizontal, 16)
                .padding(.top, 12)

                sectionLabel("STATUS")
                    .padding(.top, 6)

                GlassEffectContainer(spacing: 4) {
                    VStack(spacing: 4) {
                        ToggleButton(label: "STARTUP PERMIT",
                                     state: supervisor.startupPermit,
                                     statusColor: Theme.accent, keyHint: "P") {
                            supervisor.startupPermit = !supervisor.startupPermit
                        }
                        ToggleButton(label: "AUTO STARTUP SEQ",
                                     state: supervisor.autoStartup,
                                     statusColor: Theme.accent, keyHint: "U") {
                            supervisor.autoStartup = !supervisor.autoStartup
                        }
                        ToggleButton(label: "ROD AUTO (T-AVG)",
                                     state: supervisor.rodAutoEnabled,
                                     statusColor: Theme.accent, keyHint: "O") {
                            supervisor.rodAutoEnabled = !supervisor.rodAutoEnabled
                        }
                        ToggleButton(label: "TURBINE TRIP",
                                     state: supervisor.turbineTrip,
                                     statusColor: Theme.caution, keyHint: "T") {
                            supervisor.turbineTrip = !supervisor.turbineTrip
                        }
                        ToggleButton(label: "PUMP DEGRADED",
                                     state: supervisor.pumpDegraded,
                                     statusColor: Theme.alarm, keyHint: "Z") {
                            supervisor.pumpDegraded = !supervisor.pumpDegraded
                        }
                        ToggleButton(label: "FEEDWATER FAULT",
                                     state: supervisor.feedwaterFault,
                                     statusColor: Theme.alarm, keyHint: "X") {
                            supervisor.feedwaterFault = !supervisor.feedwaterFault
                        }
                    }
                }
                .padding(.horizontal, 12)

                ScramButton(supervisor: supervisor)
                    .padding(.horizontal, 12)
                    .padding(.top, 12)
                    .padding(.bottom, 20)
            }
        }
        .background(.ultraThinMaterial)
    }

    private func sectionLabel(_ text: String) -> some View {
        Text(text)
            .font(.system(size: 9, weight: .semibold, design: .monospaced))
            .foregroundStyle(Theme.textDim)
            .tracking(2)
            .padding(.leading, 16)
            .padding(.top, 16)
            .padding(.bottom, 4)
    }
}

// MARK: — Slim PanelHeader used across all tabs
// Pure glass design: plain tracked label, no material strip, no accent rule.
struct PanelHeader: View {
    let title: String
    var body: some View {
        HStack {
            Text(title)
                .font(.system(size: 9, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.textHdr)
                .tracking(1.6)
            Spacer()
        }
        .padding(.horizontal, Theme.panelPadding)
        .padding(.top, 10)
        .padding(.bottom, 3)
    }
}

// MARK: — DCS Slider — industrial fader: thin groove, rectangular handle with
// center notch, calibrated tick row. No glow, no glass ball.
struct DCSSlider: View {
    let label: String
    @Binding var value: Double
    let displayStr: (Double) -> String
    var color: Color = Theme.accent   // single accent — restrained palette

    private let thumbW: CGFloat = 9
    private let thumbH: CGFloat = 20

    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack {
                Text(label)
                    .font(.system(size: 10, design: .monospaced))
                    .foregroundStyle(Theme.textDim)
                    .tracking(0.5)
                Spacer()
                Text(displayStr(value))
                    .font(.system(size: 11, weight: .semibold, design: .monospaced))
                    .foregroundStyle(.white)
            }
            GeometryReader { geo in
                let w = geo.size.width
                VStack(spacing: 3) {
                    ZStack(alignment: .leading) {
                        // Groove
                        Capsule()
                            .fill(Color.black.opacity(0.45))
                            .frame(height: 3)
                        // Fill up to handle
                        Capsule()
                            .fill(color.opacity(0.55))
                            .frame(width: max(0, (w - thumbW) * value + thumbW / 2), height: 3)
                        // Fader handle: rectangle with index notch
                        RoundedRectangle(cornerRadius: 1.5, style: .continuous)
                            .fill(Color(white: 0.82))
                            .frame(width: thumbW, height: thumbH)
                            .overlay(
                                Rectangle()
                                    .fill(Color.black.opacity(0.55))
                                    .frame(height: 1.2))
                            .offset(x: (w - thumbW) * value)
                    }
                    .frame(height: thumbH)
                    // Tick row: majors at 0/25/50/75/100%
                    Canvas { ctx, size in
                        for i in 0...20 {
                            let major = i % 5 == 0
                            let x = thumbW / 2 + (size.width - thumbW) * CGFloat(i) / 20
                            var p = Path()
                            p.move(to: .init(x: x, y: 0))
                            p.addLine(to: .init(x: x, y: major ? size.height : size.height * 0.5))
                            ctx.stroke(p, with: .color(.white.opacity(major ? 0.30 : 0.14)),
                                       lineWidth: 1)
                        }
                    }
                    .frame(height: 5)
                }
                .contentShape(Rectangle())
                .gesture(DragGesture(minimumDistance: 0)
                    .onChanged { g in
                        let x = g.location.x - thumbW / 2
                        value = max(0, min(1, x / max(1, w - thumbW)))
                    }
                )
            }
            .frame(height: 28)
        }
    }
}

// MARK: — Toggle Button (glass)
struct ToggleButton: View {
    let label: String
    let state: Bool
    let statusColor: Color
    let keyHint: String
    let action: () -> Void

    var body: some View {
        Button(action: action) {
            HStack(spacing: 8) {
                Circle()
                    .fill(state ? statusColor : Color.white.opacity(0.2))
                    .frame(width: 6, height: 6)
                    .shadow(color: state ? statusColor : .clear, radius: 4)
                Text(label)
                    .font(.system(size: 11, design: .monospaced))
                    .foregroundStyle(state ? .white : Theme.textDim)
                Spacer()
                Text("[\(keyHint)]")
                    .font(.system(size: 10, design: .monospaced))
                    .foregroundStyle(Theme.textDim.opacity(0.6))
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 8)
        }
        .buttonStyle(.plain)
        .glassEffect(
            state
                ? .regular.tint(statusColor.opacity(0.2)).interactive()
                : .regular.interactive(),
            in: .rect(cornerRadius: Theme.controlRadius, style: .continuous)
        )
    }
}

// MARK: — SCRAM Button (glass with red tint)
struct ScramButton: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.4, on: .main, in: .common).autoconnect()

    var body: some View {
        Button {
            if supervisor.scrammed { supervisor.resetScram() }
            else { supervisor.triggerScram() }
        } label: {
            VStack(spacing: 3) {
                if supervisor.scrammed {
                    Text("⚠  SCRAM ACTIVE")
                        .font(.system(size: 14, weight: .bold, design: .monospaced))
                        .foregroundStyle(.white)
                    Text("tap to reset")
                        .font(.system(size: 10, design: .monospaced))
                        .foregroundStyle(.white.opacity(0.6))
                } else {
                    Text("SCRAM")
                        .font(.system(size: 22, weight: .bold, design: .monospaced))
                        .foregroundStyle(.white)
                }
            }
            .frame(maxWidth: .infinity)
            .frame(height: 64)
        }
        .buttonStyle(.plain)
        .glassEffect(
            supervisor.scrammed
                ? .regular.tint(Theme.alarm.opacity(blink ? 0.5 : 0.2)).interactive()
                : .regular.tint(Theme.alarm.opacity(0.15)).interactive(),
            in: .rect(cornerRadius: Theme.controlRadius, style: .continuous)
        )
        .onReceive(timer) { _ in blink.toggle() }
    }
}
