# ReactorSim — Session Recap

**Project:** Native macOS SwiftUI nuclear reactor (PWR) training simulator.
**Dir:** `/Users/alexanderwessbladhardh/kärnreaktor/ReactorSim/`
**Repo:** github.com/wh1xdy/nuclear_reactor_control_room (branch `macos-swift`).
**Build:** `swift build` (debug) / `swift build -c release` (release for 600× speed).
**Test:** `swift test` — 22 headless tests, all passing (2026-07-02).
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
  vessel + core (live axial-flux bands) + rods, pressurizer (level + PORV), SG
  (secondary level), hot/cold legs, RCP, MSIV, HP/MSR/LP turbine train,
  generator, **400 kV switchyard one-line with click-to-operate breakers**
  (52G + two 3φ circuits; load rejection trips the unit), condenser, feed
  train; instrument docks (neutronics, electrical + synchroscope, RCS primary,
  steam cycle, turbine-generator, CORE·THERMAL data page, margin-to-trip
  strip); 45° PCB-style pipe routing; alarm banner; control strip with
  **per-loop A/M stations**. **Esc** pauses the sim clock (pause menu).
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
  pressure (Antoine, with the exact inverse `satTempK`); Carnot-scaled gross
  efficiency (SG 553 K / 6.4 MPa, 12 K pinch off the hot leg). Steady state:
  fuel 900 K, T-avg 550 K, ~990 MWe.
- **AxialCore** — 1-D 15-node axial flux + per-node xenon layered on the
  point-kinetics (total power stays 0-D; the module solves the SHAPE). Gives
  REAL axial offset, Fz/Fq, peak clad (flow-dependent film ΔT) and min-DNBR
  (W-3-style CHF at the LIVE primary pressure). Rod kick at equilibrium xenon
  → damped ΔI oscillation, ~32 h period (verified). PWR/SMR rods top-entry;
  BWR blades bottom-entry + a void-profile tilt.
- **TurbineGenerator** — synchronous machine on an infinite bus (Xs = 1.8 pu)
  with a first-order AVR: real MVAr / power factor / load angle / field V-A /
  stator kA; thermal lags (bearing, stator, lube, H₂, τ = minutes); rpm
  coastdown on trip with a vibration bump through the rotor criticals.
- **Reactor kinds** (`PlantParams.pwr()/.bwr()/.smr()`) — BWR: void feedback,
  direct cycle, 7 MPa, no boron/PZR; SMR: integral PWR, ~200 MWt, natural
  circulation. All three selectable in Settings.
- Fidelity: neutronics + decay heat + xenon are genuine industry-standard 0-D
  models; axial shape / DNBR / clad are calibrated reduced-order correlations
  (right trends, approximate magnitudes); thermal side is a calibrated
  reduced-order model.

## Tests (`Tests/ReactorSimTests/`) — 22, all green
PhysicsTests: DecayHeat (3), Plant (3), ThermalHydraulics (4), Automation (2),
Supervisor (4), ReactorKind (4: BWR void feedback/direct cycle, SMR natcirc,
PWR unchanged). AxialCoreTests (2: nominal shape sane, rod insertion shifts
flux down).

## Known next steps / open items
- Visual polish of the new consoles (Workstation/Benchboard built but not yet
  eyeballed).
- Kind-aware mimic: BWR/SMR still draw the PWR plant picture (pressurizer/SG/
  boron shown for reactors that lack them) — needs gating on `sup.params`.
- Optional higher fidelity: real pressurizer (two-phase + heaters/spray on
  saturation), boron worth vs temperature/burnup, axial-xenon EOL divergence
  (raise `xeMult` toward ~2.5).
