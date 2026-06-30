// Console.swift — selectable operator-console layouts + reactor type.
// Every console is just a different "shell" over the SAME PlantSupervisor, so
// the physics is shared and new layouts (or reactor types) plug in here.

import SwiftUI

/// The on-screen layout / operating environment. Picked in Settings, persisted.
enum Console: String, CaseIterable, Identifiable {
    case mimic        // full-screen plant mimic with embedded values
    case workstation  // tiled multi-pane DCS operator workstation
    case benchboard   // skeuomorphic hardware panel
    case dashboard    // the original tabbed dashboard

    var id: String { rawValue }

    var label: String {
        switch self {
        case .mimic:       return "PLANT MIMIC"
        case .workstation: return "DCS WORKSTATION"
        case .benchboard:  return "HARDWARE BENCHBOARD"
        case .dashboard:   return "ENGINEERING DASHBOARD"
        }
    }

    var short: String {
        switch self {
        case .mimic:       return "MIMIC"
        case .workstation: return "DCS"
        case .benchboard:  return "PANEL"
        case .dashboard:   return "DASH"
        }
    }

    var blurb: String {
        switch self {
        case .mimic:       return "One large P&ID with live values on the schematic. How operators actually work."
        case .workstation: return "Alarm list, mimic, trends and a faceplate — all tiled and visible at once."
        case .benchboard:  return "A physical control panel: lamp tiles, bezeled gauges, strip recorders, hand switches."
        case .dashboard:   return "The original tabbed engineering dashboard (F1–F6)."
        }
    }

    var sfSymbol: String {
        switch self {
        case .mimic:       return "point.topleft.down.to.point.bottomright.curvepath"
        case .workstation: return "rectangle.split.2x2"
        case .benchboard:  return "switch.2"
        case .dashboard:   return "square.grid.2x2"
        }
    }
}

/// Reactor type. Only PWR is modelled today; the rest are placeholders so the
/// settings UI (and the architecture) is ready for them.
enum ReactorType: String, CaseIterable, Identifiable {
    case pwr
    case bwr
    case smr

    var id: String { rawValue }
    var available: Bool { true }

    /// Maps the UI selection onto the physics-layer reactor profile.
    var kind: ReactorKind {
        switch self {
        case .pwr: return .pwr
        case .bwr: return .bwr
        case .smr: return .smr
        }
    }

    var label: String {
        switch self {
        case .pwr: return "PWR — Pressurized Water Reactor"
        case .bwr: return "BWR — Boiling Water Reactor"
        case .smr: return "SMR — Small Modular Reactor"
        }
    }

    /// One-line model description shown in Settings.
    var modelBlurb: String {
        switch self {
        case .pwr: return "3000 MWt · pressurizer + steam generators · soluble-boron shim."
        case .bwr: return "Direct cycle · 7 MPa steam dome · recirc-flow power control via void feedback."
        case .smr: return "Integral 200 MWt · natural circulation · passive — no reactor coolant pumps."
        }
    }
    var short: String { rawValue.uppercased() }
}
