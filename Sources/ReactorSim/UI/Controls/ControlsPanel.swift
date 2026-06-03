// ControlsPanel.swift — Left panel: sliders, toggles, SCRAM button.
// Uses Liquid Glass for the panel surface on macOS 26+.

import SwiftUI

struct ControlsPanel: View {
    let supervisor: PlantSupervisor

    var body: some View {
        ScrollView(.vertical, showsIndicators: false) {
            VStack(alignment: .leading, spacing: 0) {
                PanelHeader(title: "CONTROLS")

                VStack(spacing: 12) {
                    DCSSlider(label: "ROD POSITION",
                              value: Binding(get: { supervisor.rodPosition },
                                            set: { supervisor.rodPosition = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" },
                              color: Theme.accent)

                    DCSSlider(label: "PRIMARY FLOW",
                              value: Binding(get: { supervisor.primaryFlow },
                                            set: { supervisor.primaryFlow = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" },
                              color: Theme.water)

                    DCSSlider(label: "TURBINE VALVE",
                              value: Binding(get: { supervisor.turbineValve },
                                            set: { supervisor.turbineValve = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" },
                              color: Theme.steam)

                    DCSSlider(label: "FEEDWATER VALVE",
                              value: Binding(get: { supervisor.feedwaterValve },
                                            set: { supervisor.feedwaterValve = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" },
                              color: Theme.twophase)

                    DCSSlider(label: "BORATION RATE",
                              value: Binding(get: { supervisor.borationRate },
                                            set: { supervisor.borationRate = $0 }),
                              displayStr: { "\(Int($0 * 100)) %" },
                              color: Theme.caution)
                }
                .padding(.horizontal, 12)
                .padding(.top, 8)

                Divider().background(Theme.sep).padding(.vertical, 10).padding(.horizontal, 12)

                VStack(spacing: 6) {
                    ToggleButton(label: "STARTUP PERMIT",
                                 state: supervisor.startupPermit,
                                 statusColor: Theme.normal, keyHint: "P") {
                        supervisor.startupPermit = !supervisor.startupPermit
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
                .padding(.horizontal, 12)

                Divider().background(Theme.sep).padding(.vertical, 10).padding(.horizontal, 12)

                ScramButton(supervisor: supervisor)
                    .padding(.horizontal, 12)
                    .padding(.bottom, 16)
            }
        }
        .background {
            // Liquid Glass: enable in Xcode 17 + macOS 26 SDK once confirmed:
            // Rectangle().glassEffect(in: Rectangle())
            Theme.panel
        }
        .overlay(alignment: .trailing) { Theme.border.frame(width: 1) }
    }
}

// MARK: — Panel header with accent underline + corner brackets
struct PanelHeader: View {
    let title: String

    var body: some View {
        ZStack(alignment: .leading) {
            Theme.panelHdr
            Text(title)
                .font(Theme.readoutSm)
                .foregroundStyle(Theme.textHdr)
                .padding(.leading, 12)
        }
        .frame(height: 32)
        .overlay(alignment: .bottom) { Theme.accent.frame(height: 1) }
        .overlay(CornerBrackets())
    }
}

// MARK: — DCS Slider
struct DCSSlider: View {
    let label: String
    @Binding var value: Double
    let displayStr: (Double) -> String
    let color: Color

    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack {
                Text(label)
                    .font(Theme.readoutSm)
                    .foregroundStyle(Theme.textDim)
                Spacer()
                Text(displayStr(value))
                    .font(Theme.readoutSm)
                    .foregroundStyle(color)
            }
            GeometryReader { geo in
                ZStack(alignment: .leading) {
                    // Track
                    RoundedRectangle(cornerRadius: 4)
                        .fill(Theme.sliderBg)
                        .frame(height: 14)

                    // Fill
                    RoundedRectangle(cornerRadius: 4)
                        .fill(color)
                        .frame(width: max(0, geo.size.width * value), height: 14)

                    // Thumb
                    Circle()
                        .fill(Theme.text)
                        .frame(width: 18, height: 18)
                        .shadow(color: color.opacity(0.6), radius: 4)
                        .offset(x: max(0, geo.size.width * value - 9))
                }
                .contentShape(Rectangle())
                .gesture(DragGesture(minimumDistance: 0)
                    .onChanged { g in
                        value = max(0, min(1, g.location.x / geo.size.width))
                    }
                )
            }
            .frame(height: 18)
        }
    }
}

// MARK: — Toggle Button (bordered, glass-style)
struct ToggleButton: View {
    let label: String
    let state: Bool
    let statusColor: Color
    let keyHint: String
    let action: () -> Void

    var body: some View {
        Button(action: action) {
            HStack {
                // LED
                Circle()
                    .fill(state ? statusColor : Theme.border)
                    .frame(width: 8, height: 8)
                    .shadow(color: state ? statusColor.opacity(0.8) : .clear, radius: 4)

                Text(label)
                    .font(Theme.readoutSm)
                    .foregroundStyle(state ? statusColor : Theme.textDim)

                Spacer()

                Text("[\(keyHint)]")
                    .font(Theme.readoutSm)
                    .foregroundStyle(Theme.textDim)
            }
            .padding(.horizontal, 10)
            .frame(height: 30)
            .background(
                state
                    ? statusColor.opacity(0.12)
                    : Theme.panel,
                in: RoundedRectangle(cornerRadius: 5)
            )
            .overlay(
                RoundedRectangle(cornerRadius: 5)
                    .stroke(state ? statusColor : Theme.border, lineWidth: 1)
            )
        }
        .buttonStyle(.plain)
    }
}

// MARK: — SCRAM Button
struct ScramButton: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.4, on: .main, in: .common).autoconnect()

    var body: some View {
        Button {
            if supervisor.scrammed { supervisor.resetScram() }
            else { supervisor.triggerScram() }
        } label: {
            VStack(spacing: 4) {
                if supervisor.scrammed {
                    Text("⚠  SCRAM ACTIVE")
                        .font(Theme.readoutMd)
                        .foregroundStyle(.white)
                    Text("tap to reset")
                        .font(Theme.readoutSm)
                        .foregroundStyle(.white.opacity(0.7))
                } else {
                    Text("SCRAM")
                        .font(Theme.readoutLg)
                        .foregroundStyle(.white)
                }
            }
            .frame(maxWidth: .infinity)
            .frame(height: 60)
            .background(
                supervisor.scrammed
                    ? (blink ? Theme.alarm : Theme.alarm.opacity(0.6))
                    : Color(r: 80, g: 18, b: 18),
                in: RoundedRectangle(cornerRadius: 8)
            )
            .overlay(
                RoundedRectangle(cornerRadius: 8)
                    .stroke(Theme.alarm, lineWidth: 2)
            )
        }
        .buttonStyle(.plain)
        .onReceive(timer) { _ in blink.toggle() }
    }
}

// MARK: — Corner brackets decoration
struct CornerBrackets: View {
    var color: Color = Theme.accent
    var arm: CGFloat = 8

    var body: some View {
        GeometryReader { geo in
            let w = geo.size.width; let h = geo.size.height
            Canvas { ctx, _ in
                let c = GraphicsContext.Shading.color(color)
                func bracket(x: CGFloat, y: CGFloat, dx: CGFloat, dy: CGFloat) {
                    var h1 = Path(); h1.move(to: CGPoint(x: x, y: y)); h1.addLine(to: CGPoint(x: x + dx*arm, y: y))
                    var v1 = Path(); v1.move(to: CGPoint(x: x, y: y)); v1.addLine(to: CGPoint(x: x, y: y + dy*arm))
                    ctx.stroke(h1, with: c, lineWidth: 1)
                    ctx.stroke(v1, with: c, lineWidth: 1)
                }
                bracket(x: 0,   y: 0,   dx:  1, dy:  1)
                bracket(x: w,   y: 0,   dx: -1, dy:  1)
                bracket(x: 0,   y: h,   dx:  1, dy: -1)
                bracket(x: w,   y: h,   dx: -1, dy: -1)
            }
        }
        .allowsHitTesting(false)
    }
}
