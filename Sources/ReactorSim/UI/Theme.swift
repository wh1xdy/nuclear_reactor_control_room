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
    // Dark-skin dim lifted ~25% so small scanned readouts (kg/s, K) survive glass blur.
    static var textDim    : Color { isLight ? Color(r: 74, g: 80, b: 86) : Color(r: 108, g: 117, b: 111) }
    static var textHdr    : Color { isLight ? Color(r: 40, g: 46, b: 52) : Color(r: 135, g: 145, b: 140) }

    /// Primary "ink" — bright foreground on dark / dark ink on the light skin.
    /// Replaces every hardcoded white so canvases invert correctly.
    static var ink        : Color { isLight ? Color(r: 26, g: 30, b: 34) : .white }
    /// Equipment body fill and schematic field for the P&ID, per skin.
    static var equipFill  : Color { isLight ? Color(r: 196, g: 201, b: 206) : Color(r: 15, g: 17, b: 21) }
    static var schematicBg: Color { isLight ? Color(r: 162, g: 168, b: 174) : Color(r: 8, g: 11, b: 16) }

    // MARK: — ISA status colors (saturated so they POP against the desk)
    // Three-way: light ISA / authentic-dark ISA (both kept disciplined) / GUIDED
    // brightened so reserved colors survive the Liquid-Glass blur on near-black.
    static let normal     = Color(r: 40,  g: 150, b: 60)    // green (legacy app-wide constant)
    static var caution    : Color { isLight ? Color(r: 200, g: 130, b: 0) : (isFlat ? Color(r: 205, g: 160, b: 18) : Color(r: 240, g: 180, b: 41)) }   // amber
    static var warning    : Color { isFlat ? Color(r: 195, g: 95, b: 18) : Color(r: 242, g: 114, b: 43) }    // orange
    static var alarm      : Color { isLight ? Color(r: 196, g: 22, b: 22) : (isFlat ? Color(r: 205, g: 38, b: 38) : Color(r: 255, g: 69, b: 58)) }    // red

    // MARK: — Accent
    // Light ISA-101 uses a muted dark slate (grayscale philosophy). The dark
    // skins (guided + authenticDark) keep the electric blue.
    // Branch on isFlat so BOTH authentic skins desaturate (light→slate, dark→graphite);
    // only GUIDED keeps the electric blue. Fixes pump triangles reading vivid blue on authentic-dark.
    static var accent     : Color { isFlat ? (isLight ? Color(r: 60, g: 76, b: 94) : Color(r: 120, g: 130, b: 142)) : Color(r: 48, g: 144, b: 200) }

    // MARK: — P&ID fluid colors
    // Branched like the newer fluid tokens so liquid fills / hot-leg nozzles read
    // as translucent ink (monochrome) in the AUTHENTIC skins instead of leaking
    // saturated blue/red onto a reserved-color HMI.
    static var water      : Color { isFlat ? ink.opacity(0.18) : Color(r: 38,  g: 95,  b: 185) }
    static let twophase   = Color(r: 28,  g: 145, b: 185)
    static let steam      = Color(r: 105, g: 118, b: 128)
    static var hotLeg     : Color { isFlat ? ink.opacity(0.5)  : Color(r: 185, g: 65,  b: 28) }

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

// MARK: — Purposeful color system (Plant Mimic)
// Two laws: (1) HUE encodes physical medium + reactor state via continuous
// functions of live values; (2) every saturated pixel maps to a REAL setpoint,
// so any warm/orange/red pixel is a genuine off-normal measurement. The
// anti-cliché inversion — operating temperature is a CALM teal, not warm — is
// what gives the warning colors meaning. Every hue collapses to graphite/ink in
// the AUTHENTIC skins (ISA-101 discipline); color erupts there only off-normal.
extension Theme {
    /// RGB interpolation across (knot,r,g,b) stops, knots ascending in [0,1].
    private static func lerpRGB(_ stops: [(k: Double, r: Double, g: Double, b: Double)], _ t: Double) -> Color {
        let x = max(0, min(1, t))
        var a = stops[0], b = stops[stops.count - 1]
        for i in 0..<(stops.count - 1) where x >= stops[i].k && x <= stops[i + 1].k {
            a = stops[i]; b = stops[i + 1]; break
        }
        let f = max(0, min(1, (x - a.k) / max(1e-9, b.k - a.k)))
        return Color(r: Int((a.r + (b.r - a.r) * f).rounded()),
                     g: Int((a.g + (b.g - a.g) * f).rounded()),
                     b: Int((a.b + (b.b - a.b) * f).rounded()))
    }

    private static let thermalStops: [(k: Double, r: Double, g: Double, b: Double)] = [
        (0.00,  46, 107, 214),   // #2E6BD6 cold blue
        (0.05,  70, 199, 192),   // #46C7C0 operating teal — the calm signature
        (0.28, 224, 169,  46),   // #E0A92E amber
        (0.64, 232,  99,  28),   // #E8631C orange
        (1.00, 255,  59,  48),   // #FF3B30 incandescent red
    ]

    /// Temperature → color. GUIDED: 5-stop ramp over 525–1500 K. AUTHENTIC:
    /// graphite ink unless K exceeds the caller's real trip (then caution/alarm).
    static func colorFor(_ K: Double, trip: Double = .infinity) -> Color {
        if isFlat {
            if K >= trip { return alarm }
            if K >= trip * 0.97 { return caution }
            return ink
        }
        return lerpRGB(thermalStops, (K - 525) / 975)
    }

    private static let fuelStops: [(k: Double, r: Double, g: Double, b: Double)] = [
        (0.00,  70, 199, 192),   // #46C7C0 teal at the ~900 K normal centerline (calm)
        (0.50, 224, 169,  46),   // #E0A92E amber midway to trip
        (1.00, 255,  59,  48),   // #FF3B30 at the 1500 K fuel trip
    ]
    /// Fuel temperature → color, anchored to FUEL's own normal band (≈900 K) so a
    /// warm fuel readout means abnormal — not a normal-state amber like colorFor.
    static func colorForFuel(_ K: Double, trip: Double = 1500) -> Color {
        if isFlat { return K >= trip ? alarm : (K >= trip * 0.93 ? caution : ink) }
        return lerpRGB(fuelStops, (K - 900) / 600)
    }

    private static let fluxStops: [(k: Double, r: Double, g: Double, b: Double)] = [
        (0.00,  21,  36,  58),   // #15243A subcritical dark
        (0.33, 122,  74, 107),   // #7A4A6B low-power ember (narrow transition only)
        (0.66, 232, 144,  42),   // #E8902A power amber
        (1.00, 255, 210, 122),   // #FFD27A full-power white-hot
    ]
    /// Core incandescence by power fraction g∈[0,1] → (center, edge) for a radial
    /// fill. AUTHENTIC stays flat (equipFill) — flux is shown numerically there.
    static func fluxGlow(_ g: Double) -> (center: Color, edge: Color) {
        if isFlat { return (equipFill, equipFill) }
        return (lerpRGB(fluxStops, g), Color(r: 232, g: 144, b: 42))
    }

    /// Setpoint deviation → color. Within ±band reads calm; drift warms.
    /// GUIDED: teal → amber (dev 1–2) → red. AUTHENTIC: ink → caution → alarm.
    static func setpointDev(_ v: Double, _ sp: Double, _ band: Double) -> Color {
        let dev = abs(v - sp) / max(1e-9, band)
        if isFlat { return dev <= 1 ? ink : (dev <= 2 ? caution : alarm) }
        if dev <= 1 { return Color(r: 70, g: 199, b: 192) }                 // #46C7C0
        if dev <= 2 { return lerpRGB([(0, 70, 199, 192), (1, 240, 180, 41)], dev - 1) }
        return Color(r: 255, g: 69, b: 58)                                   // #FF453A
    }

    // ── Fluid / medium tokens (GUIDED hue; AUTHENTIC graphite process line) ──
    static var fluidSubcooled: Color { isFlat ? text : Color(r: 47,  g: 111, b: 224) }  // #2F6FE0
    static var fluidTwoPhase:  Color { isFlat ? text : Color(r: 52,  g: 198, b: 200) }  // #34C6C8
    static var fluidExhaust:   Color { isFlat ? text : Color(r: 138, g: 147, b: 176) }  // #8A93B0
    static var fluidFeedwater: Color { isFlat ? text : Color(r: 63,  g: 169, b: 201) }  // #3FA9C9
    static var fluidHotLeg:    Color { isFlat ? text : Color(r: 70,  g: 199, b: 192) }  // #46C7C0

    // ── Status ladder: alive / off / electrical (off ≠ failed) ──
    static var statusNormal: Color { isFlat ? normal : Color(r: 39, g: 194, b: 129) }   // #27C281
    static var deEnergized:  Color { isFlat ? textDim : Color(r: 88, g: 97, b: 115) }   // #586173
    // Cooler brass-gold, a clear hue step off caution amber (#F0B429) so a synced
    // generator is never mistaken for a warning.
    static var elecGold:     Color { isFlat ? text : Color(r: 206, g: 178, b: 96) }     // brass
    /// Faint structural tint behind instrument docks — the only background tint, toward blue.
    static var dockTint:     Color { isFlat ? ink.opacity(0.04) : Color(r: 26, g: 42, b: 56).opacity(0.5) }
}
