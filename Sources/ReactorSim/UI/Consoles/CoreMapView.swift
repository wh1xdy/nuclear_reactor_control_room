// CoreMapView.swift — BEACON/GARDEL-style core-map popup.
// Opened by tapping the reactor vessel on the mimic; Esc / backdrop closes.
//
// The assembly grid combines the REAL model outputs — radial rings (3-node
// solver with ring-wise xenon), the azimuthal first-harmonic tilt (QPTR), and
// the 1-D axial shape — with a fixed ±3% fuel-batch checkerboard (real cores
// are loaded fresh/once/twice-burnt; the pattern is representative, the rest
// is live physics). Values are relative assembly power ×100 (100 = core avg).

import SwiftUI

struct CoreMapView: View {
    let supervisor: PlantSupervisor

    private let lattice = 15          // 15×15 positions, circle-masked

    var body: some View {
        let snap = supervisor.snapshot
        ZStack {
            Rectangle().fill(.black.opacity(0.55)).ignoresSafeArea()
                .contentShape(Rectangle())
                .onTapGesture { supervisor.coreMapOpen = false }
            HStack(alignment: .top, spacing: 18) {
                VStack(alignment: .leading, spacing: 8) {
                    Text("CORE MAP · RELATIVE ASSEMBLY POWER ×100")
                        .font(.system(size: 11, weight: .semibold, design: .monospaced))
                        .foregroundStyle(Theme.textHdr).tracking(0.5)
                    grid(snap: snap)
                        .frame(width: 430, height: 430)
                    Text("radial rings × axial × azimuthal tilt · fuel-batch checkerboard ±3% representative")
                        .font(.system(size: 7, design: .monospaced)).foregroundStyle(Theme.textDim)
                }
                VStack(alignment: .leading, spacing: 0) {
                    HStack {
                        Text("LIMITS").font(.system(size: 10, weight: .semibold, design: .monospaced))
                            .foregroundStyle(Theme.textHdr)
                        Spacer()
                        Button { supervisor.coreMapOpen = false } label: {
                            Text("CLOSE  esc").font(.system(size: 9, design: .monospaced))
                                .foregroundStyle(Theme.textDim)
                                .padding(.horizontal, 8).padding(.vertical, 4)
                                .contentShape(Rectangle())
                        }
                        .buttonStyle(.plain).controlSurface()
                    }
                    .padding(.bottom, 10)
                    readouts(snap: snap)
                    Text("AXIAL SHAPE").font(.system(size: 9, weight: .semibold, design: .monospaced))
                        .foregroundStyle(Theme.textHdr).padding(.top, 14).padding(.bottom, 6)
                    axialBars(snap: snap).frame(width: 190, height: 150)
                }
                .frame(width: 200)
            }
            .padding(22)
            .background(RoundedRectangle(cornerRadius: 16).fill(Theme.panel))
            .overlay(RoundedRectangle(cornerRadius: 16).stroke(Theme.border, lineWidth: 1))
        }
    }

    // ── Assembly grid ────────────────────────────────────────────────────────

    /// Assembly factor at lattice position (i, j): rings interpolated at the
    /// assembly radius × azimuthal tilt × batch checkerboard.
    private func factor(_ i: Int, _ j: Int, snap: PlantSnapshot) -> Double? {
        let c = Double(lattice - 1) / 2
        let dx = (Double(i) - c) / c, dy = (Double(j) - c) / c
        let r = (dx * dx + dy * dy).squareRoot()
        guard r <= 1.02 else { return nil }                    // circle mask
        let rings = snap.radialRings.count == 3 ? snap.radialRings : [1.14, 1.05, 0.81]
        // Equal-area ring centres at r ≈ 0.41 / 0.71 / 0.91.
        let radial: Double
        if r < 0.41      { radial = rings[0] }
        else if r < 0.71 { radial = rings[0] + (rings[1] - rings[0]) * (r - 0.41) / 0.30 }
        else if r < 0.91 { radial = rings[1] + (rings[2] - rings[1]) * (r - 0.71) / 0.20 }
        else             { radial = rings[2] * (1.0 - (r - 0.91) * 1.2) }
        let tilt = 1.0 + (snap.tiltX * dx + snap.tiltY * dy)
        let batch = (i + j) % 2 == 0 ? 1.03 : 0.97
        return max(0, radial * tilt * batch)
    }

    private func grid(snap: PlantSnapshot) -> some View {
        Canvas { ctx, size in
            let cell = size.width / CGFloat(lattice)
            let pf = max(0, snap.powerFraction)
            for i in 0..<lattice {
                for j in 0..<lattice {
                    guard let f = factor(i, j, snap: snap) else { continue }
                    let rect = CGRect(x: CGFloat(i) * cell + 1, y: CGFloat(j) * cell + 1,
                                      width: cell - 2, height: cell - 2)
                    let g = max(0, min(1, pf * f / 1.10))
                    let glow = Theme.fluxGlow(g)
                    ctx.fill(Path(roundedRect: rect, cornerRadius: 2), with: .color(glow.center))
                    ctx.stroke(Path(roundedRect: rect, cornerRadius: 2),
                               with: .color(Theme.ink.opacity(0.10)), lineWidth: 0.5)
                    let txtC: Color = Theme.isFlat ? Theme.ink.opacity(0.75)
                                    : (g > 0.45 ? .black.opacity(0.72) : Theme.ink.opacity(0.55))
                    ctx.draw(Text("\(Int((f * 100).rounded()))")
                                .font(.system(size: 8, weight: .medium, design: .monospaced))
                                .foregroundColor(txtC),
                             at: CGPoint(x: rect.midX, y: rect.midY), anchor: .center)
                }
            }
        }
    }

    // ── Right column ─────────────────────────────────────────────────────────

    private func readouts(snap: PlantSnapshot) -> some View {
        VStack(spacing: 7) {
            row("QPTR", String(format: "%.3f", snap.qptr), snap.qptr > 1.02 ? Theme.caution : Theme.ink)
            row("Fz",   String(format: "%.2f", snap.fz), Theme.ink)
            row("Fq",   String(format: "%.2f", snap.fq), snap.fq > 2.3 ? Theme.caution : Theme.ink)
            row("FΔH",  String(format: "%.2f", snap.fdh), Theme.ink)
            row("ΔI",   String(format: "%+.1f %%", snap.axialOffsetPct), Theme.ink)
            row(supervisor.reactorKind == .bwr ? "MCPR" : "DNBR",
                String(format: "%.2f", snap.minDNBR), Theme.ink)
            row("PK CLAD", String(format: "%.0f K", snap.peakCladTempK), Theme.ink)
            row("RINGS", snap.radialRings.count == 3
                ? String(format: "%.2f/%.2f/%.2f", snap.radialRings[0], snap.radialRings[1], snap.radialRings[2])
                : "—", Theme.textDim)
        }
    }

    private func row(_ l: String, _ v: String, _ c: Color) -> some View {
        HStack {
            Text(l).font(.system(size: 9, design: .monospaced)).foregroundStyle(Theme.textDim)
            Spacer()
            Text(v).font(.system(size: 11, weight: .semibold, design: .monospaced)).foregroundStyle(c)
        }
    }

    /// Horizontal axial-shape bars, node 0 (bottom of core) at the bottom.
    private func axialBars(snap: PlantSnapshot) -> some View {
        Canvas { ctx, size in
            let prof = snap.axialProfile
            guard !prof.isEmpty else { return }
            let n = prof.count
            let bh = size.height / CGFloat(n)
            let maxV = max(1.4, prof.max() ?? 1.4)
            for i in 0..<n {
                let y = size.height - CGFloat(i + 1) * bh
                let wFrac = CGFloat(prof[i] / maxV)
                let bar = CGRect(x: 0, y: y + 1, width: size.width * wFrac, height: bh - 2)
                let g = max(0, min(1, snap.powerFraction * prof[i] / 1.10))
                ctx.fill(Path(roundedRect: bar, cornerRadius: 1.5),
                         with: .color(Theme.isFlat ? Theme.ink.opacity(0.35) : Theme.fluxGlow(g).center))
            }
            // Mean line at 1.0.
            let xm = size.width * CGFloat(1.0 / maxV)
            ctx.stroke(Path { p in p.move(to: CGPoint(x: xm, y: 0)); p.addLine(to: CGPoint(x: xm, y: size.height)) },
                       with: .color(Theme.ink.opacity(0.25)), style: StrokeStyle(lineWidth: 0.75, dash: [3, 3]))
        }
    }
}
