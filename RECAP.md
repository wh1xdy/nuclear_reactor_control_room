# ReactorSim — Session Recap

**Project:** Native macOS SwiftUI nuclear reactor (PWR) training simulator.
**Dir:** `/Users/alexanderwessbladhardh/kärnreaktor/ReactorSim/`
**Repo:** github.com/wh1xdy/nuclear_reactor_control_room (branch `macos-swift`).
**Build:** `swift build` (debug) / `swift build -c release` (release for 600× speed).
**Test:** `swift test` — 16 headless tests, all passing (2026-06-15).
**Run:** `swift run ReactorSim`. macOS 26 (Tahoe) required for Liquid Glass.

---

## Architecture

One physics model (`PlantSupervisor` → `PWRPlant`) drives several swappable
**operator consoles**, chosen in **Settings** (gear in the system bar, or `,`)
and persisted via `@AppStorage`. `ContentView` is a thin router: a global
`SystemBar` (nameplate, sim clock, alarm state, **ACK ALL**, console switcher,
gear) over the selected console.

### Consoles (`Sources/ReactorSim/UI/Consoles/`)
- **Plant Mimic** (`MimicConsole` + `MimicDiagram`) — full-screen PWR schematic:
  vessel + core/rods, pressurizer (level + PORV), SG (secondary level), hot/cold
  legs, RCP, MSIV, HP/LP turbine, generator, condenser, feed train; live values
  embedded; alarm banner; control strip with **per-loop A/M stations**.
- **DCS Workstation** (`WorkstationConsole`) — tiled 4-pane: alarm summary,
  mimic, 4-pen trend group, control faceplates (SP/PV/OUT).
- **Hardware Benchboard** (`BenchboardConsole`) — annunciator lamp wall, bezeled
  round gauges, strip recorders, hand switches, SCRAM.
- **Engineering Dashboard** (`DashboardConsole`) — the original F1–F6 tabbed
  layout (Overview/Primary/Secondary/Reactivity/Alarms/I&C) + left control panel.

### Skins (`Theme.skin`, cycle with `M`)
- **GUIDED** — Liquid Glass, dark, soft squircles.
- **AUTHENTIC LIGHT** — flat ISA-101 High-Performance HMI, light desaturated gray.
- **AUTHENTIC DARK** — same flat ISA-101 layout on a dark desk.

`Theme.isFlat` (flat surfaces + square corners) is true for both authentic
variants; `Theme.isLight` (palette) only for the light one — so guided and
authentic-dark share the dark palette, differing only in glass-vs-flat surfaces.
Every panel routes through `panelSurface()` / `controlSurface()` /
`readoutSurface()`; every former hardcoded white goes through `Theme.ink` so
canvases invert correctly between light and dark.

## Controls / automation
- Per-loop **AUTO/MANUAL stations** on ROD and FW (real A/M stations); dragging a
  fader drops it to MANUAL. **MASTER AUTO** engages rod + feedwater + pressurizer.
- Turbine-following: with rods on AUTO the reactor follows turbine load and holds
  T-avg = 550 K hands-off (verified by `AutomationTests`). Manual = you babysit it.
- Global **ACK ALL** in the system bar (works in every console; `C` key too).

## Physics (`Sources/ReactorSim/Physics/`)
- **PointKinetics** — six-group U-235 delayed neutrons (β≈650 pcm, Keepin λ's,
  Λ=20 µs), RK2 with adaptive substepping + prompt-jump fast path.
- **DecayHeat** — full ANS/ANSI-5.1 23-group, exact exponential integrator,
  ~6.5% equilibrium.
- **XenonIodine** — standard I-135/Xe-135 two-ODE poison model.
- **ThermalHydraulics** — TWO-node primary loop (hot leg / cold leg) with a
  ~30 K core ΔT and a transport-delay lag (~4 s/leg, ∝1/flow); Dittus-Boelter
  (convective HTC ∝ flow^0.8, advection linear in flow); saturation steam
  pressure (Antoine); Carnot-scaled gross efficiency. Steady state: fuel 900 K,
  T-avg 550 K, ~990 MWe.
- Fidelity: neutronics + decay heat + xenon are genuine industry-standard 0-D
  models; thermal side is a calibrated reduced-order model (good intuition,
  approximate off-nominal numbers).

## Tests (`Tests/ReactorSimTests/PhysicsTests.swift`) — 16, all green
DecayHeat (3), Plant (3: steady-state, scram/decay, rod rate limit),
ThermalHydraulics (4: leg split, Dittus-Boelter, transport lag, saturation),
Automation (2: MASTER AUTO holds steady, rods follow turbine load),
Supervisor (4: BOP at 600×, alarm latching, startup sequencer, scram-reset sync).

## Known next steps / open items
- Visual polish of the new consoles (Workstation/Benchboard built but not yet
  eyeballed); remaining mimic label/pipe spacing.
- Optional higher fidelity: real pressurizer (two-phase + heaters/spray on
  saturation), boron worth vs temperature/burnup.
- Architecture is ready for more reactor types (`ReactorType` enum: PWR active,
  BWR/SMR stubbed).
