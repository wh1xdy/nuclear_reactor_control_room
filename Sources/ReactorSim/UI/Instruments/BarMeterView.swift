// BarMeterView.swift — Vertical bar meter with ISA color bands + trip lines.

import SwiftUI

struct BarMeterView: View {
    let value: Double
    let lo: Double
    let hi: Double
    let label: String
    let unit: String
    var tripHi: Double? = nil
    var tripLo: Double? = nil

    private var fraction: Double { max(0, min(1, (value - lo) / max(1e-9, hi - lo))) }

    private var fillColor: Color {
        if let t = tripHi, value > t { return Theme.alarm }
        if value > hi * 0.85         { return Theme.caution }
        return Theme.normal
    }

    var body: some View {
        VStack(spacing: 4) {
            Text(label)
                .font(.system(size: 9, design: .monospaced))
                .foregroundStyle(Theme.textDim)
                .lineLimit(1)
                .minimumScaleFactor(0.6)

            GeometryReader { geo in
                ZStack(alignment: .bottom) {
                    // Bezel
                    RoundedRectangle(cornerRadius: 4)
                        .fill(Color(r: 8, g: 11, b: 18))
                        .overlay(RoundedRectangle(cornerRadius: 4).stroke(Theme.border, lineWidth: 1))

                    // Fill bar
                    let barH = max(0, (geo.size.height - 4) * CGFloat(fraction))
                    RoundedRectangle(cornerRadius: 3)
                        .fill(fillColor)
                        .frame(width: geo.size.width - 4, height: barH)
                        .padding(.bottom, 2)
                        .overlay(alignment: .top) {
                            // Bright top edge
                            fillColor.opacity(0.6)
                                .frame(height: 1)
                                .offset(y: -barH + 1)
                        }

                    // Trip lines
                    if let t = tripHi {
                        let ty = geo.size.height - (geo.size.height - 4) * CGFloat((t - lo) / max(1e-9, hi - lo)) - 2
                        Rectangle()
                            .fill(Theme.alarm)
                            .frame(height: 1.5)
                            .offset(y: ty - geo.size.height / 2)
                    }
                }
            }
            .frame(maxWidth: .infinity)

            // Numeric readout
            Text(formattedValue)
                .font(.system(size: 10, design: .monospaced))
                .foregroundStyle(fillColor)
            if !unit.isEmpty {
                Text(unit)
                    .font(.system(size: 8, design: .monospaced))
                    .foregroundStyle(Theme.textDim)
            }
        }
    }

    private var formattedValue: String {
        abs(value) < 1000 ? String(format: "%.1f", value) : String(format: "%.0f", value)
    }
}
