// MimicConsole.swift — full-screen plant mimic.
// Alarm banner → key-parameter row → detailed plant mimic → control strip.
// The mimic IS the screen; no tabs, no side panels.

import SwiftUI

struct MimicConsole: View {
    let supervisor: PlantSupervisor

    var body: some View {
        VStack(spacing: 0) {
            MimicAlarmBanner(supervisor: supervisor)
            MimicKeyParams(supervisor: supervisor)
            MimicDiagram(supervisor: supervisor)
                .frame(maxWidth: .infinity, maxHeight: .infinity)
            MimicControlStrip(supervisor: supervisor)
        }
    }
}

// MARK: — Key parameter row (no overlap; lives under the banner)
private struct MimicKeyParams: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        HStack(spacing: 0) {
            cell("REACTOR POWER", String(format: "%.1f", s.powerFraction * 100), "% RTP",
                 .powerStatus(s.powerFraction))
            sep
            cell("THERMAL", String(format: "%.0f", s.thermalPowerW / 1e6), "MWt", Theme.ink)
            sep
            cell("GROSS ELECTRIC", String(format: "%.0f", s.electricPowerW / 1e6), "MWe", Theme.ink)
            sep
            cell("RCS T-AVG", String(format: "%.1f", s.coolantTempK), "K",
                 s.coolantTempK > 616 ? Theme.alarm : Theme.ink)
            sep
            cell(supervisor.hasPressurizer ? "PZR PRESS" : "DOME PRESS",
                 String(format: "%.2f", supervisor.pressureMPa), "MPa",
                 supervisor.pressureMPa > supervisor.nominalPressureMPa * 1.097 ? Theme.alarm : Theme.ink)
            sep
            cell("REACTIVITY", pcmString(s.reactivity), "pcm",
                 .reactivityStatus(s.reactivity))
        }
        .frame(height: 62)
        .background(Theme.panel)
        .overlay(alignment: .bottom) { Rectangle().fill(Theme.border).frame(height: 1) }
    }
    // Avoids the ugly "−0" that "%+.0f" prints for tiny negative reactivity.
    private func pcmString(_ rho: Double) -> String {
        let pcm = rho * 1e5
        return abs(pcm) < 0.5 ? "0" : String(format: "%+.0f", pcm)
    }
    private var sep: some View { Rectangle().fill(Theme.sep).frame(width: 1).padding(.vertical, 12) }
    private func cell(_ label: String, _ value: String, _ unit: String, _ color: Color) -> some View {
        VStack(spacing: 2) {
            Text(label).font(.system(size: 9, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.textDim).tracking(0.5)
            HStack(alignment: .lastTextBaseline, spacing: 4) {
                Text(value).font(.system(size: 22, weight: .bold, design: .monospaced)).foregroundStyle(color)
                Text(unit).font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
            }
        }
        .frame(maxWidth: .infinity)
    }
}

// MARK: — Alarm banner
private struct MimicAlarmBanner: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.5, on: .main, in: .common).autoconnect()

    var body: some View {
        let trips  = supervisor.trips
        let alarms = supervisor.alarms
        let (tint, lead): (Color?, String) =
            !trips.isEmpty  ? (Theme.alarm,   "⚠ REACTOR TRIP") :
            !alarms.isEmpty ? (Theme.caution, "● ALARM") :
                              (nil,           "● PLANT NORMAL")
        let detail = trips.first ?? alarms.first?.message ?? "all parameters within limits"

        return HStack(spacing: 12) {
            Text(lead)
                .font(.system(size: 13, weight: .bold, design: .monospaced))
                .foregroundStyle(tint == nil ? Theme.textDim : Theme.ink)
                .opacity(!trips.isEmpty && blink ? 0.5 : 1)
            Rectangle().fill(Theme.sep).frame(width: 1, height: 16)
            Text(detail.uppercased())
                .font(.system(size: 11, design: .monospaced))
                .foregroundStyle(Theme.text).lineLimit(1)
            Spacer()
            Text("\(alarms.count) ACTIVE")
                .font(.system(size: 10, weight: .semibold, design: .monospaced))
                .foregroundStyle(Theme.textDim)
            Button { supervisor.acknowledgeAllAlarms() } label: {
                Text("ACK [C]")
                    .font(.system(size: 10, weight: .semibold, design: .monospaced))
                    .foregroundStyle(Theme.accent)
                    .padding(.horizontal, 12).padding(.vertical, 6)
                    .contentShape(Rectangle())
            }
            .buttonStyle(.plain)
            .controlSurface()
        }
        .padding(.horizontal, 16).padding(.vertical, 8)
        .frame(maxWidth: .infinity)
        .background((tint ?? Theme.panelHdr).opacity(tint == nil ? 1 : 0.22))
        .overlay(alignment: .bottom) { Rectangle().fill(tint ?? Theme.sep).frame(height: tint == nil ? 1 : 2) }
        .onReceive(timer) { _ in blink.toggle() }
    }
}

// MARK: — Bottom control strip
// Faders with per-loop AUTO/MANUAL stations (real A/M stations) + mode toggles.
private struct MimicControlStrip: View {
    let supervisor: PlantSupervisor
    var body: some View {
        HStack(alignment: .top, spacing: 12) {
            // ── Faders. ROD and FW have a real auto controller → A/M station. ──
            faderAM("ROD DEMAND",
                    Binding(get: { supervisor.rodPosition },
                            set: { supervisor.rodPosition = $0; supervisor.rodAutoEnabled = false }),
                    { "\(Int((228 * (1 - $0)).rounded())) SWD" },
                    auto: Binding(get: { supervisor.rodAutoEnabled }, set: { supervisor.rodAutoEnabled = $0 }))
            // BWR: recirc flow IS the power lever; SMR: natural circulation means
            // there is no flow lever at all (the fader is hidden, not disabled).
            if !supervisor.isNaturalCirc {
                fader(supervisor.reactorKind == .bwr ? "RECIRC" : "RCS FLOW",
                      Binding(get: { supervisor.primaryFlow }, set: { supervisor.primaryFlow = $0 }),
                      { "\(Int($0 * 100)) %" })
            }
            fader("TBN GOV",
                  Binding(get: { supervisor.turbineValve }, set: { supervisor.turbineValve = $0 }),
                  { "\(Int($0 * 100)) %" })
            faderAM("FW REG",
                    Binding(get: { supervisor.feedwaterValve },
                            set: { supervisor.feedwaterValve = $0; supervisor.fwAutoEnabled = false }),
                    { "\(Int($0 * 100)) %" },
                    auto: Binding(get: { supervisor.fwAutoEnabled }, set: { supervisor.fwAutoEnabled = $0 }))
            // No chemical shim on a BWR — the fader would be a dead control.
            if supervisor.hasBoron {
                fader("BORATION",
                      Binding(get: { supervisor.borationRate }, set: { supervisor.borationRate = $0 }),
                      { "\(Int($0 * 100)) %" })
            }

            // ── Mode / permissive station ──
            VStack(spacing: 5) {
                HStack(spacing: 5) {
                    chip("MASTER AUTO", on: allAuto, tint: Theme.accent) { setAllAuto(!allAuto) }
                    chip("PZR AUTO", on: supervisor.pzrAutoEnabled, tint: Theme.accent) { supervisor.pzrAutoEnabled.toggle() }
                }
                HStack(spacing: 5) {
                    chip("PERMIT", on: supervisor.startupPermit, tint: Theme.accent) { supervisor.startupPermit.toggle() }
                    chip("STARTUP", on: supervisor.autoStartup, tint: Theme.accent) { supervisor.autoStartup.toggle() }
                }
                HStack(spacing: 5) {
                    chip("TBN TRIP", on: supervisor.turbineTrip, tint: Theme.caution) { supervisor.turbineTrip.toggle() }
                    chip("RESET [L]", on: false, tint: Theme.accent) { supervisor.resetScram() }
                }
                HStack(spacing: 5) {
                    chip("RCP FLT", on: supervisor.pumpDegraded, tint: Theme.alarm) { supervisor.pumpDegraded.toggle() }
                    chip("FW FLT", on: supervisor.feedwaterFault, tint: Theme.alarm) { supervisor.feedwaterFault.toggle() }
                }
            }
            .frame(width: 210)

            ScramButton(supervisor: supervisor).frame(width: 140)

            VStack(spacing: 6) {
                TrendView(values: supervisor.orderedHistory(supervisor.histPower),
                          yLo: 0, yHi: 130, color: Theme.accent, label: "RX PWR %").frame(height: 40)
                TrendView(values: supervisor.orderedHistory(supervisor.histCoolT),
                          yLo: 400, yHi: 650, color: Theme.accent, label: "T-AVG K").frame(height: 40)
            }
            .frame(minWidth: 180)
        }
        .padding(.horizontal, 16).padding(.vertical, 10)
        .frame(maxWidth: .infinity).frame(height: 164)
        .background(Theme.panel)
        .overlay(alignment: .top) { Rectangle().fill(Theme.border).frame(height: 1) }
    }

    private var allAuto: Bool { supervisor.rodAutoEnabled && supervisor.fwAutoEnabled && supervisor.pzrAutoEnabled }
    private func setAllAuto(_ on: Bool) {
        supervisor.rodAutoEnabled = on; supervisor.fwAutoEnabled = on; supervisor.pzrAutoEnabled = on
    }

    private func fader(_ label: String, _ value: Binding<Double>, _ display: @escaping (Double) -> String) -> some View {
        VStack(spacing: 0) {
            DCSSlider(label: label, value: value, displayStr: display)
            Spacer(minLength: 0)
        }.frame(width: 128)
    }

    // Fader with an AUTO/MANUAL station beneath it.
    private func faderAM(_ label: String, _ value: Binding<Double>,
                         _ display: @escaping (Double) -> String, auto: Binding<Bool>) -> some View {
        VStack(spacing: 6) {
            DCSSlider(label: label, value: value, displayStr: display)
            HStack(spacing: 4) {
                amButton("MAN", active: !auto.wrappedValue) { auto.wrappedValue = false }
                amButton("AUTO", active: auto.wrappedValue) { auto.wrappedValue = true }
            }
        }.frame(width: 128)
    }

    private func amButton(_ label: String, active: Bool, action: @escaping () -> Void) -> some View {
        Button(action: action) {
            Text(label)
                .font(.system(size: 9, weight: .semibold, design: .monospaced))
                .foregroundStyle(active ? Theme.ink : Theme.textDim)
                .frame(maxWidth: .infinity).padding(.vertical, 4)
                .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .controlSurface(tint: active ? Theme.accent : nil)
    }

    private func chip(_ label: String, on: Bool, tint: Color, action: @escaping () -> Void) -> some View {
        Button(action: action) {
            HStack(spacing: 5) {
                Circle().fill(on ? tint : Theme.ink.opacity(0.25)).frame(width: 6, height: 6)
                Text(label).font(.system(size: 9, weight: .semibold, design: .monospaced))
                    .foregroundStyle(on ? Theme.ink : Theme.textDim)
                Spacer(minLength: 0)
            }
            .padding(.horizontal, 8).padding(.vertical, 7)
            .frame(maxWidth: .infinity)
            .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .controlSurface(tint: on ? tint : nil)
    }
}
