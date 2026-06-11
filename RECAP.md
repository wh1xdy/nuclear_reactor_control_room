# ReactorSim — Session Recap

**Project:** Native macOS SwiftUI nuclear reactor (PWR) training simulator.
**Dir:** `/Users/alexanderwessbladhardh/kärnreaktor/ReactorSim/`
**Build:** `swift build` (debug) / `swift build -c release` (use release for 600× speed).
**Test:** `swift test` — 9 headless physics tests, all passing (2026-06-11).

---

## Design (locked in)

Gorgeous, minimalistic, native Apple polish: Liquid Glass, 120 Hz, squircles.
1. **Pure clean glass** — no scanlines, no grid, no corner brackets. ✅ DONE
2. **Spacious** — fewer elements, larger instruments, breathing room. ✅ DONE
3. **Restrained color** — white/grey text; electric blue accent ONLY for
   active/selected; ISA green/amber/red ONLY on live alarm/trip. ✅ DONE

### Redesign checklist — COMPLETE (2026-06-11)
- [x] Deleted `ScanlineOverlay`, `GridBackground`, `CornerBrackets`.
- [x] Removed every `.drawingGroup()` (flicker root cause — offscreen Metal
      texture races the Liquid Glass backdrop).
- [x] Squircles everywhere: `.rect(cornerRadius:style:.continuous)`.
- [x] Consistent radii: `Theme.panelRadius=16`, `Theme.controlRadius=12`.
- [x] No hard 1px strokes — glass edges define panels.
      (Exception kept: alarm tiles use ISA color FILLS, which is semantic.)
- [x] Restrained color: `Theme.powerStatus`/`reactivityStatus` return white in
      normal range; instruments fill accent-blue, amber/red only in abnormal
      bands; trend lines all accent; "ALL SYSTEMS NORMAL" is grey, not green.
- [x] Render-scope isolation: `SimClock` leaf in HeaderBar; `PIDCanvas` reads
      `supervisor.snapshot` itself — parent glass panels don't re-eval at 60 Hz.
- [x] Physics timer steps inline via `MainActor.assumeIsolated` — the old
      `Task { @MainActor }` hop queued on the main dispatch queue, which stalls
      during event tracking (would freeze physics during drags).

---

## Physics — verified by `swift test` (Tests/ReactorSimTests/PhysicsTests.swift)

- **Point kinetics:** six-group RK2, adaptive substepping with an explicit
  stability bound (`|ρ−β|/Λ·subDt ≤ 0.5`), plus a **prompt-jump fast path**
  for deep shutdown (n slaved to precursor source — O(1) instead of ~900
  substeps when scrammed).
- **Decay heat:** full **ANS-5.1 23-group** model, exact exponential update
  (stable at any dt). Builds toward ~6.5% during operation, decays along the
  ANS curve after shutdown. Rated power = fission share (93.5%) + decay share;
  `powerFraction` reports total thermal (reads 100% at steady full power).
- **Rod motion:** CRDM rate-limited at 0.0053/s (~3 min full travel) — rods
  walk toward demand, never teleport. Scram reset syncs demand to actual
  position (rods stay in until deliberately withdrawn).
- **Scram interlocks:** turbine auto-trips on scram; condenser **steam dump**
  (proportional valve) holds no-load T_avg ≈ 550 K post-trip; tripped turbine
  generates 0 MWe.
- **BOP:** steam/feedwater inventories balance at steady state (steam balance
  self-regulating); condenser + pressurizer ODEs solved exactly (600×-stable);
  pressurizer coupling asymmetric (0.05 MPa/K overheat → 17.0 MPa HIGH_PRESS
  trip reachable; 0.01 MPa/K cooldown — heaters hold pressure, no spurious
  ECCS post-scram).
- **Alarms:** annunciator latching — cleared conditions stay on the board
  until acknowledged; trips block scram reset until acked.
- `internalDt = 0.05 s` (coolant node τ ≈ 0.7 s ⇒ 14× margin); use a release
  build for smooth 600× time compression.

## Xenon startup
Fresh core (X=0, I=0), critical at rod 0 / boron 800 ppm. Xenon builds over
6–8 sim-hours; equilibrium worth ≈ −2500 pcm — operator compensates by boron
dilution (rods alone can't add positive reactivity). Decay heat also builds
from zero (fresh fuel), so the core self-adjusts slightly via temperature
feedback in the first hours — intentional realism, total thermal stays ~100%.

## Liquid Glass API notes (post-knowledge-cutoff, verified)
- `.glassEffect(.regular, in: <shape>)` ✅ — `Glass.regular/.clear/.identity`.
- `.glassEffect(.regular.tint(color).interactive(), in: <shape>)` for controls.
- `GlassEffectContainer(spacing:) { }` merges adjacent glass.
- `.glassEffectID(id, in: ns)` needs a STABLE `@Namespace`.
- `GlassEffect()` does NOT exist — compile error.

## Possible next steps
- Run the app and eyeball the new look (spacing/typography may want tuning).
- Auto controllers (ICTab cards are placeholders, all MANUAL).
- SG inventory ↔ core coupling (feedwater loss currently doesn't dry the SG).
- Multi-nuclide xenon display, boron dilution UI affordance.
