// MimicConsole.swift — full-screen plant mimic.
// Alarm banner → key-parameter row → detailed plant mimic → control strip.
// The mimic IS the screen; no tabs, no side panels.

import SwiftUI

struct MimicConsole: View {
    let supervisor: PlantSupervisor

    var body: some View {
        ZStack {
            VStack(spacing: 0) {
                MimicAlarmBanner(supervisor: supervisor)
                MimicKeyParams(supervisor: supervisor)
                MimicDiagram(supervisor: supervisor)
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
                MimicControlStrip(supervisor: supervisor)
            }
            if supervisor.coreMapOpen {
                CoreMapView(supervisor: supervisor)
            }
            if supervisor.malfMenuOpen {
                MalfunctionMenu(supervisor: supervisor)
            }
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
                // Individual fault chips moved into the instructor panel.
                chip("MALFUNCTIONS [I]", on: supervisor.anyMalfunctionActive, tint: Theme.alarm) {
                    supervisor.malfMenuOpen.toggle()
                }
            }
            .frame(width: 210)

            ScramButton(supervisor: supervisor).frame(width: 140)

            // ΔI–power flyspeck: the CAOC operating-band chart. Keep the point
            // in the corridor with rods as xenon walks the axial shape around.
            FlyspeckChart(supervisor: supervisor).frame(width: 148)

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

// MARK: — ΔI–power flyspeck (CAOC operating band)

/// The classic PWR axial-offset control chart: ΔI on x, power on y, with the
/// allowed operating corridor around the target axial offset. The live point
/// must be held inside the band with rods while xenon walks the shape around
/// (the axial-xenon oscillation makes this a real exercise at speed). The band
/// is centred on the kind's natural full-power AO (a BWR runs bottom-peaked).
private struct FlyspeckChart: View {
    let supervisor: PlantSupervisor

    var body: some View {
        Canvas { ctx, size in
            let refAO: Double = supervisor.reactorKind == .bwr ? -20 : -5
            let plot = CGRect(x: 26, y: 14, width: size.width - 32, height: size.height - 26)
            // Both axes CLAMP to the plot — a post-scram ΔI swing must pin the
            // point at the frame edge, never let it wander out of the window.
            func px(_ ao: Double) -> CGFloat {
                let f = max(0, min(1, (ao - (refAO - 20)) / 40))   // ref ± 20 %
                return plot.minX + plot.width * CGFloat(f)
            }
            func py(_ p: Double) -> CGFloat {    // power % → y (0…110)
                plot.maxY - plot.height * CGFloat(max(0, min(110, p)) / 110)
            }
            func halfW(_ p: Double) -> Double {  // corridor half-width vs power
                p >= 90 ? 6 : min(18, 6 + (90 - p) * 0.25)
            }

            ctx.fill(Path(plot), with: .color(Theme.dockTint))
            // Operating band (trapezoid, wider at low power).
            var band = Path()
            let ps = stride(from: 20.0, through: 110.0, by: 10.0).map { $0 }
            band.move(to: CGPoint(x: px(refAO - halfW(ps[0])), y: py(ps[0])))
            for p in ps { band.addLine(to: CGPoint(x: px(refAO - halfW(p)), y: py(p))) }
            for p in ps.reversed() { band.addLine(to: CGPoint(x: px(refAO + halfW(p)), y: py(p))) }
            band.closeSubpath()
            let bandTint = Theme.isFlat ? Theme.ink.opacity(0.07) : Theme.statusNormal.opacity(0.10)
            ctx.fill(band, with: .color(bandTint))
            ctx.stroke(band, with: .color(Theme.ink.opacity(0.20)), lineWidth: 0.75)
            // Target AO reference line.
            ctx.stroke(Path { p in p.move(to: CGPoint(x: px(refAO), y: plot.minY)); p.addLine(to: CGPoint(x: px(refAO), y: plot.maxY)) },
                       with: .color(Theme.ink.opacity(0.18)), style: StrokeStyle(lineWidth: 0.75, dash: [3, 3]))

            // Trail: recent ΔI/power history, fading toward the past.
            let aos = supervisor.orderedHistory(supervisor.histAO)
            let pws = supervisor.orderedHistory(supervisor.histPower)
            let n = min(aos.count, pws.count)
            if n > 8 {
                let step = max(1, n / 90)
                var prev: CGPoint? = nil
                var i = 0
                while i < n {
                    let pt = CGPoint(x: px(aos[i]), y: py(pws[i]))
                    if let a = prev {
                        let age = Double(i) / Double(n)          // 0 old → 1 new
                        ctx.stroke(Path { p in p.move(to: a); p.addLine(to: pt) },
                                   with: .color(Theme.accent.opacity(0.08 + 0.42 * age)), lineWidth: 1)
                    }
                    prev = pt
                    i += step
                }
            }

            // Live point — calm in the band, caution outside, HOLLOW when the
            // true value is pinned at the frame (off-scale).
            let ao = supervisor.snapshot.axialOffsetPct
            let pw = supervisor.snapshot.powerFraction * 100
            let inBand = abs(ao - refAO) <= halfW(pw)
            let offScale = abs(ao - refAO) > 20 || pw > 110
            let dotC: Color = inBand ? (Theme.isFlat ? Theme.ink : Theme.statusNormal) : Theme.caution
            let dot = CGPoint(x: px(ao), y: py(pw))
            if offScale {
                ctx.stroke(Path(ellipseIn: CGRect(x: dot.x - 3, y: dot.y - 3, width: 6, height: 6)),
                           with: .color(dotC), lineWidth: 1.5)
            } else {
                ctx.fill(Path(ellipseIn: CGRect(x: dot.x - 3, y: dot.y - 3, width: 6, height: 6)), with: .color(dotC))
            }
            ctx.stroke(Path(ellipseIn: CGRect(x: dot.x - 5.5, y: dot.y - 5.5, width: 11, height: 11)),
                       with: .color(dotC.opacity(0.4)), lineWidth: 1)

            // Frame + labels.
            ctx.stroke(Path(plot), with: .color(Theme.ink.opacity(0.18)), lineWidth: 1)
            func lbl(_ s: String, _ at: CGPoint, _ anchor: UnitPoint = .center) {
                ctx.draw(Text(s).font(.system(size: 7, design: .monospaced)).foregroundColor(Theme.textDim), at: at, anchor: anchor)
            }
            lbl("ΔI-PWR", CGPoint(x: plot.minX, y: 5), .topLeading)
            lbl(String(format: "%+.0f", refAO), CGPoint(x: px(refAO), y: plot.maxY + 6))
            lbl(String(format: "%+.0f", refAO - 15), CGPoint(x: px(refAO - 15), y: plot.maxY + 6))
            lbl(String(format: "%+.0f", refAO + 15), CGPoint(x: px(refAO + 15), y: plot.maxY + 6))
            lbl("100", CGPoint(x: plot.minX - 3, y: py(100)), .trailing)
            lbl("50", CGPoint(x: plot.minX - 3, y: py(50)), .trailing)
        }
    }
}
