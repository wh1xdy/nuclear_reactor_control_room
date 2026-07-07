// MalfunctionMenu.swift — instructor panel.
// Opened with the MALFUNCTIONS chip (or I); Esc / backdrop closes. Each fault
// is a latched toggle with real physics behind it — clearing a fault restores
// the SYSTEM, the transient it caused is still the operator's problem.

import SwiftUI

struct MalfunctionMenu: View {
    let supervisor: PlantSupervisor

    var body: some View {
        ZStack {
            Rectangle().fill(.black.opacity(0.55)).ignoresSafeArea()
                .contentShape(Rectangle())
                .onTapGesture { supervisor.malfMenuOpen = false }
            VStack(alignment: .leading, spacing: 10) {
                HStack {
                    Text("INSTRUCTOR · MALFUNCTIONS")
                        .font(.system(size: 12, weight: .bold, design: .monospaced))
                        .foregroundStyle(Theme.textHdr).tracking(1.2)
                    Spacer()
                    Button { supervisor.malfMenuOpen = false } label: {
                        Text("CLOSE  esc").font(.system(size: 9, design: .monospaced))
                            .foregroundStyle(Theme.textDim)
                            .padding(.horizontal, 8).padding(.vertical, 4)
                            .contentShape(Rectangle())
                    }
                    .buttonStyle(.plain).controlSurface()
                }

                row("SMALL LOCA", "50 kg/s primary break → depressurization, SI, containment pressurizes",
                    on: supervisor.primaryLeakKgs > 0) {
                    supervisor.primaryLeakKgs = supervisor.primaryLeakKgs > 0 ? 0 : 50
                }
                row("SG TUBE RUPTURE", "primary → secondary leak: SG level rises uninvited, air-ejector radiation",
                    on: supervisor.sgtrLeakKgs > 0) {
                    supervisor.sgtrLeakKgs = supervisor.sgtrLeakKgs > 0 ? 0 : 15
                }
                row("STUCK ROD", "one RCCA frozen: flux tilt (watch QPTR), degraded scram worth",
                    on: supervisor.stuckRod) { supervisor.stuckRod.toggle() }
                row("DROPPED ROD", "one RCCA on the bottom: −350 pcm + tilt — recover with rods/boron",
                    on: supervisor.droppedRod) { supervisor.droppedRod.toggle() }
                row("ATWS", "reactor protection fails to scram — emergency boration is the drill",
                    on: supervisor.atwsFault) { supervisor.atwsFault.toggle() }
                row("MSIV CLOSURE", "spurious main-steam isolation: turbine trips, SRVs are the only heat path",
                    on: supervisor.msivClosed) { supervisor.msivClosed.toggle() }
                row("LOSS OF OFFSITE POWER", "opens both 400 kV lines — house load or the diesel ladder",
                    on: supervisor.line1BreakerOpen && supervisor.line2BreakerOpen) {
                    if !(supervisor.line1BreakerOpen && supervisor.line2BreakerOpen) {
                        if !supervisor.line1BreakerOpen { supervisor.toggleLineBreaker(0) }
                        if !supervisor.line2BreakerOpen { supervisor.toggleLineBreaker(1) }
                    } else {
                        if supervisor.line1BreakerOpen { supervisor.toggleLineBreaker(0) }
                        if supervisor.line2BreakerOpen { supervisor.toggleLineBreaker(1) }
                    }
                }
                row("DIESEL FAILURE", "EDGs refuse to start — a grid loss becomes a station blackout",
                    on: supervisor.dieselFault) { supervisor.dieselFault.toggle() }
                row("RCP DEGRADATION", "reactor coolant pumps coast down", on: supervisor.pumpDegraded) {
                    supervisor.pumpDegraded.toggle()
                }
                row("FEEDWATER FAULT", "main feedwater lost", on: supervisor.feedwaterFault) {
                    supervisor.feedwaterFault.toggle()
                }

                Button {
                    supervisor.primaryLeakKgs = 0; supervisor.sgtrLeakKgs = 0
                    supervisor.stuckRod = false; supervisor.droppedRod = false
                    supervisor.atwsFault = false; supervisor.msivClosed = false
                    supervisor.dieselFault = false; supervisor.pumpDegraded = false
                    supervisor.feedwaterFault = false
                } label: {
                    Text("CLEAR ALL MALFUNCTIONS")
                        .font(.system(size: 10, weight: .semibold, design: .monospaced))
                        .foregroundStyle(Theme.ink)
                        .frame(maxWidth: .infinity).padding(.vertical, 9)
                        .contentShape(Rectangle())
                }
                .buttonStyle(.plain).controlSurface(tint: Theme.accent)
                .padding(.top, 4)
            }
            .padding(20)
            .frame(width: 470)
            .background(RoundedRectangle(cornerRadius: 14).fill(Theme.panel))
            .overlay(RoundedRectangle(cornerRadius: 14).stroke(Theme.border, lineWidth: 1))
        }
    }

    private func row(_ title: String, _ blurb: String, on: Bool, action: @escaping () -> Void) -> some View {
        Button(action: action) {
            HStack(spacing: 8) {
                Circle().fill(on ? Theme.alarm : Theme.ink.opacity(0.22)).frame(width: 7, height: 7)
                VStack(alignment: .leading, spacing: 1) {
                    Text(title).font(.system(size: 10, weight: .semibold, design: .monospaced))
                        .foregroundStyle(on ? Theme.ink : Theme.text)
                    Text(blurb).font(.system(size: 8, design: .monospaced))
                        .foregroundStyle(Theme.textDim).lineLimit(1)
                }
                Spacer()
                Text(on ? "ACTIVE" : "OFF")
                    .font(.system(size: 8, weight: .bold, design: .monospaced))
                    .foregroundStyle(on ? Theme.alarm : Theme.textDim)
            }
            .padding(.horizontal, 10).padding(.vertical, 6)
            .frame(maxWidth: .infinity)
            .contentShape(Rectangle())
            .background(RoundedRectangle(cornerRadius: 8).fill(on ? Theme.alarm.opacity(0.10) : Theme.dockTint))
            .overlay(RoundedRectangle(cornerRadius: 8).stroke(Theme.border.opacity(0.6), lineWidth: 1))
        }
        .buttonStyle(.plain)
    }
}
