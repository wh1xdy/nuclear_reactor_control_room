// Theme.swift — ISA DCS color palette + design constants.
// All colors match the approved palette from the design brief.

import SwiftUI

// MARK: — Skin
// Three visual identities the operator can switch between at runtime:
//  • .guided        — Liquid Glass, soft squircles, dark. The "premium" look.
//  • .authentic     — flat ISA-101 High-Performance HMI, LIGHT desaturated gray.
//  • .authenticDark — the same flat ISA-101 layout on a DARK desk.
// `isFlat` (surfaces/corners) is true for both authentic variants; `isLight`
// (palette) is true only for the light one — so guided and authenticDark share
// the dark palette and differ only in glass-vs-flat surface treatment.
enum Skin: String, CaseIterable {
    case guided
    case authentic
    case authenticDark

    var label: String {
        switch self {
        case .guided:        return "GUIDED"
        case .authentic:     return "AUTHENTIC LIGHT"
        case .authenticDark: return "AUTHENTIC DARK"
        }
    }
    var next: Skin {
        switch self {
        case .guided:        return .authentic
        case .authentic:     return .authenticDark
        case .authenticDark: return .guided
        }
    }
}

enum Theme {
    // Active skin — set once by ContentView before the tree builds.
    static var skin: Skin = .guided
    static var isFlat:  Bool { skin != .guided }       // flat surfaces + square corners
    static var isLight: Bool { skin == .authentic }    // light ISA-101 palette

    // MARK: — Background
    // GUIDED is the dark Liquid-Glass theme. AUTHENTIC is ISA-101 High-Performance
    // HMI: a desaturated medium-gray desk with light-gray panels, near-black
    // process text, and color RESERVED for alarms — exactly how a modern nuclear
    // control room is mandated to look.
    static var bg         : Color { isLight ? Color(r: 126, g: 132, b: 138) : Color(r: 10, g: 12, b: 14) }
    static var panel      : Color { isLight ? Color(r: 178, g: 183, b: 188) : Color(r: 16, g: 19, b: 22) }
    static var panelHdr   : Color { isLight ? Color(r: 150, g: 156, b: 162) : Color(r: 20, g: 24, b: 28) }
    static var border     : Color { isLight ? Color(r: 92, g: 99, b: 107) : Color(r: 44, g: 50, b: 56) }
    static var sep        : Color { isLight ? Color(r: 120, g: 127, b: 134) : Color(r: 30, g: 34, b: 38) }

    // MARK: — Text
    static var text       : Color { isLight ? Color(r: 26, g: 30, b: 34) : Color(r: 215, g: 218, b: 215) }
    static var textDim    : Color { isLight ? Color(r: 74, g: 80, b: 86) : Color(r: 85, g: 94, b: 88) }
    static var textHdr    : Color { isLight ? Color(r: 40, g: 46, b: 52) : Color(r: 135, g: 145, b: 140) }

    /// Primary "ink" — bright foreground on dark / dark ink on the light skin.
    /// Replaces every hardcoded white so canvases invert correctly.
    static var ink        : Color { isLight ? Color(r: 26, g: 30, b: 34) : .white }
    /// Equipment body fill and schematic field for the P&ID, per skin.
    static var equipFill  : Color { isLight ? Color(r: 196, g: 201, b: 206) : Color(r: 15, g: 17, b: 21) }
    static var schematicBg: Color { isLight ? Color(r: 162, g: 168, b: 174) : Color(r: 8, g: 11, b: 16) }

    // MARK: — ISA status colors (saturated so they POP against the desk)
    static let normal     = Color(r: 40,  g: 150, b: 60)    // green
    static var caution    : Color { isLight ? Color(r: 200, g: 130, b: 0) : Color(r: 205, g: 160, b: 18) }   // amber
    static let warning    = Color(r: 195, g: 95,  b: 18)    // orange
    static var alarm      : Color { isLight ? Color(r: 196, g: 22, b: 22) : Color(r: 205, g: 38, b: 38) }    // red

    // MARK: — Accent
    // Light ISA-101 uses a muted dark slate (grayscale philosophy). The dark
    // skins (guided + authenticDark) keep the electric blue.
    static var accent     : Color { isLight ? Color(r: 60, g: 76, b: 94) : Color(r: 48, g: 144, b: 200) }

    // MARK: — P&ID fluid colors
    static let water      = Color(r: 38,  g: 95,  b: 185)
    static let twophase   = Color(r: 28,  g: 145, b: 185)
    static let steam      = Color(r: 105, g: 118, b: 128)
    static let hotLeg     = Color(r: 185, g: 65,  b: 28)

    // MARK: — Slider
    static let sliderBg   = Color(r: 28,  g: 32,  b: 36)
    static let sliderFg   = Color(r: 48,  g: 144, b: 200)

    // MARK: — Sizing
    // Guided uses squircle (.continuous) corners; authentic uses near-square
    // corners like a real DCS mimic. Both read the same token names.
    static var panelRadius: CGFloat     { isFlat ? 0 : 16 }
    static var controlRadius: CGFloat   { isFlat ? 0 : 12 }
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

// MARK: — Skin-aware surfaces
// One choke point for the glass-vs-flat decision. Every panel and control routes
// through these so a skin switch is a single branch, not a per-view rewrite.

enum CtrlShape { case rounded, capsule, panel }

extension View {
    /// Panel container surface. Guided → Liquid Glass; Authentic → flat fill + hard 1px border.
    @ViewBuilder
    func panelSurface() -> some View {
        switch Theme.skin {
        case .guided:
            self.glassEffect(.regular,
                             in: .rect(cornerRadius: Theme.panelRadius, style: .continuous))
        case .authentic, .authenticDark:
            self.background(Theme.panel)
                .clipShape(.rect(cornerRadius: Theme.panelRadius, style: .continuous))
                .overlay(
                    RoundedRectangle(cornerRadius: Theme.panelRadius, style: .continuous)
                        .strokeBorder(Theme.border, lineWidth: 1))
        }
    }

    /// Lightweight static readout chip (non-interactive value tiles).
    @ViewBuilder
    func readoutSurface() -> some View {
        switch Theme.skin {
        case .guided:
            self.glassEffect(.clear, in: .rect(cornerRadius: Theme.controlRadius, style: .continuous))
        case .authentic, .authenticDark:
            self.background(Theme.ink.opacity(0.06))
                .clipShape(.rect(cornerRadius: Theme.controlRadius, style: .continuous))
                .overlay(
                    RoundedRectangle(cornerRadius: Theme.controlRadius, style: .continuous)
                        .strokeBorder(Theme.border.opacity(0.6), lineWidth: 1))
        }
    }

    /// Interactive control surface (buttons, pills, tiles).
    /// `tint` is the PURE status color (nil = neutral). Each skin decides how to express it.
    @ViewBuilder
    func controlSurface(tint: Color? = nil, shape: CtrlShape = .rounded) -> some View {
        let r: CGFloat = {
            switch shape {
            case .capsule: return 999
            case .panel:   return Theme.panelRadius
            case .rounded: return Theme.controlRadius
            }
        }()
        switch Theme.skin {
        case .guided:
            if let tint {
                self.glassEffect(.regular.tint(tint.opacity(0.20)).interactive(),
                                 in: .rect(cornerRadius: r, style: .continuous))
            } else {
                self.glassEffect(.regular.interactive(),
                                 in: .rect(cornerRadius: r, style: .continuous))
            }
        case .authentic, .authenticDark:
            self.background(tint == nil ? AnyShapeStyle(Theme.panelHdr)
                                        : AnyShapeStyle(tint!.opacity(0.22)))
                .clipShape(.rect(cornerRadius: r, style: .continuous))
                .overlay(
                    RoundedRectangle(cornerRadius: r, style: .continuous)
                        .strokeBorder(tint ?? Theme.border, lineWidth: tint == nil ? 1 : 1.5))
        }
    }
}
