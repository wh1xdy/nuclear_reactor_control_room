// Theme.swift — ISA DCS color palette + design constants.
// All colors match the approved palette from the design brief.

import SwiftUI

enum Theme {
    // MARK: — Background
    static let bg         = Color(r: 10,  g: 12,  b: 14)
    static let panel      = Color(r: 16,  g: 19,  b: 22)
    static let panelHdr   = Color(r: 20,  g: 24,  b: 28)
    static let border     = Color(r: 44,  g: 50,  b: 56)
    static let sep        = Color(r: 30,  g: 34,  b: 38)

    // MARK: — Text
    static let text       = Color(r: 215, g: 218, b: 215)
    static let textDim    = Color(r: 85,  g: 94,  b: 88)
    static let textHdr    = Color(r: 135, g: 145, b: 140)

    // MARK: — ISA status colors
    static let normal     = Color(r: 55,  g: 185, b: 75)    // green
    static let caution    = Color(r: 205, g: 160, b: 18)    // amber
    static let warning    = Color(r: 195, g: 95,  b: 18)    // orange
    static let alarm      = Color(r: 205, g: 38,  b: 38)    // red

    // MARK: — Accent (electric blue)
    static let accent     = Color(r: 48,  g: 144, b: 200)

    // MARK: — P&ID fluid colors
    static let water      = Color(r: 38,  g: 95,  b: 185)
    static let twophase   = Color(r: 28,  g: 145, b: 185)
    static let steam      = Color(r: 105, g: 118, b: 128)
    static let hotLeg     = Color(r: 185, g: 65,  b: 28)

    // MARK: — Slider
    static let sliderBg   = Color(r: 28,  g: 32,  b: 36)
    static let sliderFg   = Color(r: 48,  g: 144, b: 200)

    // MARK: — Sizing
    // Locked-in design: squircle (.continuous) corners — panels 16, controls 12.
    static let panelRadius: CGFloat     = 16
    static let controlRadius: CGFloat   = 12
    static let panelPadding: CGFloat    = 12
    static let headerHeight: CGFloat    = 44   // thin utility bar, not a billboard
    static let tabHeight: CGFloat       = 30
    static let controlsWidth: CGFloat   = 300

    // MARK: — Fonts (monospace DCS readouts)
    static func monoFont(_ size: CGFloat, weight: Font.Weight = .regular) -> Font {
        .system(size: size, weight: weight, design: .monospaced)
    }
    static let readoutSm:  Font = monoFont(12)
    static let readout:    Font = monoFont(14)
    static let readoutMd:  Font = monoFont(16, weight: .semibold)
    static let readoutLg:  Font = monoFont(22, weight: .bold)
    static let readoutXl:  Font = monoFont(32, weight: .bold)
}

// MARK: — Helpers

// MARK: — NaN-safe formatting

func safeFormat(_ value: Double, _ fmt: String, fallback: String = "---") -> String {
    guard value.isFinite else { return fallback }
    return String(format: fmt, value)
}

extension Double {
    var isSafe: Bool { isFinite }
    func fmt(_ format: String, fallback: String = "---") -> String {
        safeFormat(self, format, fallback: fallback)
    }
}

extension Color {
    init(r: Int, g: Int, b: Int) {
        self.init(red: Double(r)/255, green: Double(g)/255, blue: Double(b)/255)
    }

    /// Status color from power fraction.
    /// Restrained palette: white in normal range — ISA colors only when abnormal.
    static func powerStatus(_ pf: Double) -> Color {
        if pf > 1.1   { return Theme.alarm }
        if pf > 1.05  { return Theme.caution }   // 100% is normal — amber only above
        return Theme.text
    }

    /// Status color from reactivity — white unless abnormal.
    static func reactivityStatus(_ rho: Double) -> Color {
        if abs(rho) > 0.005 { return Theme.alarm }
        if abs(rho) > 0.001 { return Theme.caution }
        return Theme.text
    }
}
