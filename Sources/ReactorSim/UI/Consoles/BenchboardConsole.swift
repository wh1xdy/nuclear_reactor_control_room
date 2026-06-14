// BenchboardConsole.swift — hardware control-panel layout.
// An annunciator lamp wall across the top, a row of bezeled round gauges and
// strip recorders, and a switch deck (hand switches + faders + SCRAM) at the
// bottom — like standing at a physical benchboard rather than a software app.

import SwiftUI

struct BenchboardConsole: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 10) {
            AnnunciatorWall(supervisor: supervisor)
            HStack(spacing: 10) {
                GaugeBay(supervisor: supervisor).frame(maxWidth: .infinity)
                RecorderBay(supervisor: supervisor).frame(width: 320)
            }
            .frame(maxHeight: .infinity)
            SwitchDeck(supervisor: supervisor)
        }
        .padding(10)
    }
}

// MARK: — Annunciator lamp wall
private struct AnnunciatorWall: View {
    let supervisor: PlantSupervisor
    @State private var blink = false
    let timer = Timer.publish(every: 0.5, on: .main, in: .common).autoconnect()

    private struct Win: Identifiable { let id: String; let label: String; let active: Bool; let trip: Bool; let unack: Bool }

    private var windows: [Win] {
        func a(_ id: String) -> ReactorAlarm? { supervisor.alarms.first { $0.id == id } }
        func w(_ id: String, _ l: String, trip: Bool) -> Win {
            let al = a(id); return Win(id: id, label: l, active: al != nil, trip: trip, unack: al?.state == "unack")
        }
        func st(_ id: String, _ l: String, _ on: Bool, trip: Bool) -> Win { Win(id: id, label: l, active: on, trip: trip, unack: false) }
        return [
            st("RX_TRIP", "REACTOR\nTRIP", supervisor.scrammed, trip: true),
            w("HIGH_FLUX", "HI NEUTRON\nFLUX", trip: true),
            w("HIGH_FUEL_T", "HI FUEL\nTEMP", trip: true),
            w("HIGH_PRESS", "HI RCS\nPRESS", trip: true),
            w("HIGH_COOL_T", "HI RCS\nT-AVG", trip: true),
            w("ECCS_ACT", "ECCS\nACTUATION", trip: true),
            w("PORV_OPEN", "PZR PORV\nOPEN", trip: false),
            w("WARN_PRESS", "PZR PRESS\nHI", trip: false),
            w("LOW_FEED", "LO FW\nINVENTORY", trip: false),
            w("XE_TRANSIENT", "XENON\nTRANSIENT", trip: false),
            st("TBN_TRIP", "TURBINE\nTRIP", supervisor.turbineTrip, trip: false),
            st("STM_DUMP", "STEAM DUMP\nOPEN", supervisor.steamDumpValve > 0.001, trip: false),
            st("RCP_DEG", "RCP\nDEGRADED", supervisor.pumpDegraded, trip: false),
            st("FW_FAULT", "FEEDWATER\nFAULT", supervisor.feedwaterFault, trip: false),
            st("SPARE", "SPARE", false, trip: false),
        ]
    }

    var body: some View {
        let cols = [GridItem](repeating: GridItem(.flexible(), spacing: 6), count: 5)
        LazyVGrid(columns: cols, spacing: 6) {
            ForEach(windows) { win in
                let fill: Color = win.active
                    ? (win.trip ? Theme.alarm : Theme.caution).opacity(win.unack ? (blink ? 0.95 : 0.4) : 0.7)
                    : Theme.ink.opacity(0.05)
                Text(win.label)
                    .font(.system(size: 9, weight: win.active ? .bold : .medium, design: .monospaced))
                    .multilineTextAlignment(.center).lineLimit(2).minimumScaleFactor(0.7)
                    .foregroundStyle(win.active ? .white : Theme.textDim)
                    .frame(maxWidth: .infinity).frame(height: 40)
                    .background(fill)
                    .overlay(RoundedRectangle(cornerRadius: 3).strokeBorder(Theme.border, lineWidth: 1))
                    .clipShape(.rect(cornerRadius: 3))
            }
        }
        .padding(10)
        .panelSurface()
        .frame(height: 120)
        .onReceive(timer) { _ in blink.toggle() }
    }
}

// MARK: — Round gauge bay
private struct GaugeBay: View {
    let supervisor: PlantSupervisor
    var body: some View {
        let s = supervisor.snapshot
        VStack(spacing: 0) {
            PanelHeader(title: "PROCESS GAUGES")
            HStack(spacing: 6) {
                gauge("RX POWER", ArcGaugeView(value: s.powerFraction * 100, lo: 0, hi: 130, label: "", unit: "%", tripHi: 120))
                gauge("PZR PRESS", ArcGaugeView(value: supervisor.pressureMPa, lo: 0, hi: 18, label: "", unit: "MPa", tripHi: 17))
                gauge("RCS T-AVG", ArcGaugeView(value: s.coolantTempK, lo: 500, hi: 600, label: "", unit: "K", tripHi: 593))
                gauge("GROSS", ArcGaugeView(value: s.electricPowerW / 1e6, lo: 0, hi: 1100, label: "", unit: "MWe", tripHi: nil))
            }
            .padding(12)
            .frame(maxHeight: .infinity)
        }
        .panelSurface()
    }
    private func gauge(_ label: String, _ g: ArcGaugeView) -> some View {
        VStack(spacing: 4) {
            Text(label).font(.system(size: 9, weight: .semibold, design: .monospaced)).foregroundStyle(Theme.textDim)
            g
        }
        .frame(maxWidth: .infinity)
        .padding(10)
        .background(Theme.ink.opacity(0.04))
        .overlay(RoundedRectangle(cornerRadius: Theme.controlRadius).strokeBorder(Theme.border.opacity(0.7), lineWidth: 1))
        .clipShape(.rect(cornerRadius: Theme.controlRadius))
    }
}

// MARK: — Strip recorders
private struct RecorderBay: View {
    let supervisor: PlantSupervisor
    var body: some View {
        VStack(spacing: 0) {
            PanelHeader(title: "STRIP RECORDERS")
            VStack(spacing: 8) {
                TrendView(values: supervisor.orderedHistory(supervisor.histPower),
                          yLo: 0, yHi: 130, color: Theme.accent, label: "RX PWR %").frame(maxHeight: .infinity)
                TrendView(values: supervisor.orderedHistory(supervisor.histCoolT),
                          yLo: 400, yHi: 650, color: Theme.accent, label: "T-AVG K").frame(maxHeight: .infinity)
                TrendView(values: supervisor.orderedHistory(supervisor.histPress),
                          yLo: 13, yHi: 17, color: Theme.accent, label: "PZR P MPa").frame(maxHeight: .infinity)
            }
            .padding(10)
        }
        .panelSurface()
    }
}

// MARK: — Switch deck (hand switches + faders + SCRAM)
private struct SwitchDeck: View {
    let supervisor: PlantSupervisor
    var body: some View {
        HStack(alignment: .top, spacing: 14) {
            fader("ROD", Binding(get: { supervisor.rodPosition },
                                 set: { supervisor.rodPosition = $0; supervisor.rodAutoEnabled = false }),
                  { "\(Int((228 * (1 - $0)).rounded())) SWD" })
            fader("FLOW", Binding(get: { supervisor.primaryFlow }, set: { supervisor.primaryFlow = $0 }), { "\(Int($0*100))%" })
            fader("TBN", Binding(get: { supervisor.turbineValve }, set: { supervisor.turbineValve = $0 }), { "\(Int($0*100))%" })
            fader("FW", Binding(get: { supervisor.feedwaterValve },
                                set: { supervisor.feedwaterValve = $0; supervisor.fwAutoEnabled = false }),
                  { "\(Int($0*100))%" })

            // Hand switches
            VStack(spacing: 6) {
                HStack(spacing: 6) {
                    handSwitch("ROD AUTO", on: supervisor.rodAutoEnabled) { supervisor.rodAutoEnabled.toggle() }
                    handSwitch("FW AUTO", on: supervisor.fwAutoEnabled) { supervisor.fwAutoEnabled.toggle() }
                }
                HStack(spacing: 6) {
                    handSwitch("PZR AUTO", on: supervisor.pzrAutoEnabled) { supervisor.pzrAutoEnabled.toggle() }
                    handSwitch("STARTUP", on: supervisor.autoStartup) { supervisor.autoStartup.toggle() }
                }
            }
            .frame(width: 220)

            ScramButton(supervisor: supervisor).frame(width: 150)
        }
        .padding(.horizontal, 14).padding(.vertical, 12)
        .frame(maxWidth: .infinity).frame(height: 130)
        .panelSurface()
    }

    private func fader(_ label: String, _ value: Binding<Double>, _ display: @escaping (Double) -> String) -> some View {
        DCSSlider(label: label, value: value, displayStr: display).frame(width: 130)
    }

    // Two-position panel switch (vertical rocker look).
    private func handSwitch(_ label: String, on: Bool, action: @escaping () -> Void) -> some View {
        Button(action: action) {
            VStack(spacing: 5) {
                Text(label).font(.system(size: 8, weight: .semibold, design: .monospaced)).foregroundStyle(Theme.textDim)
                ZStack(alignment: on ? .top : .bottom) {
                    RoundedRectangle(cornerRadius: 4).fill(Theme.ink.opacity(0.10)).frame(width: 26, height: 40)
                    RoundedRectangle(cornerRadius: 3)
                        .fill(on ? Theme.accent : Theme.ink.opacity(0.4))
                        .frame(width: 22, height: 18).padding(2)
                }
                Text(on ? "AUTO" : "MAN").font(.system(size: 8, weight: .bold, design: .monospaced))
                    .foregroundStyle(on ? Theme.accent : Theme.textDim)
            }
            .frame(maxWidth: .infinity).padding(.vertical, 6)
            .contentShape(Rectangle())
        }
        .buttonStyle(.plain)
        .background(Theme.ink.opacity(0.03))
        .overlay(RoundedRectangle(cornerRadius: Theme.controlRadius).strokeBorder(Theme.border.opacity(0.6), lineWidth: 1))
        .clipShape(.rect(cornerRadius: Theme.controlRadius))
    }
}
