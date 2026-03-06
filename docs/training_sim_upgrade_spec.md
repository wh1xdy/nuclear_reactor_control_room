# Operator Training Simulator Upgrade Specification

This document analyses the current codebase and lists every change required
to make the simulator credible as an operator training tool. Each section
describes the current state, the gap, and what must be implemented. Items
are ordered by impact on training value.

---

## 1. Physics Engines

### 1.1 Boron / Chemical Shim (PWR) — **missing entirely**

Real PWR operators spend most of their shift adjusting boron concentration
(300–1200 ppm) to control long-term reactivity, compensate for xenon, and
perform load-following. Without it the PWR model is not trainable.

**What to add:**
- State variable `boron_ppm` (float, ~0–2000 ppm).
- Reactivity worth: `rho_boron = -boron_worth_per_ppm * (boron_ppm - nominal_boron_ppm)`
  where `boron_worth_per_ppm ≈ -8 pcm/ppm` (i.e. `−8e-5 Δk/k per ppm`).
- Two control inputs: `boration_rate` (kg/s of boric acid injection) and
  `dilution_rate` (kg/s of demineralized water). Both change boron_ppm
  through a first-order mixing model with primary coolant inventory.
- Display: show boron ppm as a numeric readout; add a trend.
- UI keys: e.g. `B`/`N` or mouse sliders for boration/dilution rate.
- Interlocks: startup without meeting boron requirement should block startup
  permit.

---

### 1.2 Pressurizer — **severely simplified**

Current model: `pressure_mpa` driven by a single linear ODE with heater
coefficient and turbine draw. No surge line, no spray, no PORV/SRV.

**What to add:**
- Pressurizer liquid inventory `V_surge` (m³, 0–40 m³).
- Surge line flow: `Q_surge = K_surge * sign(T_hot - T_sat(P)) * |delta_T|^0.5`
  representing thermal expansion/contraction of primary coolant.
- Spray: opening spray valve (`F_spray` 0–1) injects cold coolant into the
  steam space, condensing steam and reducing pressure.
- Heater banks (proportional+backup): electric power input heats liquid,
  vaporises water, raises pressure.
- PORV (Power-Operated Relief Valve): opens automatically at setpoint
  (e.g. 16.55 MPa for a 4-loop Westinghouse), relieves steam to the
  pressurizer relief tank (PRT). PRT rupture disk if PRT overpressure.
- Safety valves (SRV): open at 17.2 MPa; reclosing setpoint lower.
- Pressure-volume relationship: treat pressurizer steam space as saturated
  steam, use IAPWS specific volume to relate level to pressure.
- Control: proportional heater (raise P), proportional spray (lower P).
  This gives realistic pressure oscillations during transients.

---

### 1.3 Steam Generator Level / Feedwater Control (PWR) — **missing**

The current SG node is a single temperature. Real PWR operators control
SG water level (narrow range: 25–75%) via feedwater control valves
(one per SG, typically three-element control). Low SG level trips the
reactor; high SG level trips the turbine.

**What to add:**
- SG liquid inventory `m_sg` [kg] (or dimensionless level 0–1).
- Mass balance: `dm_sg/dt = m_dot_fw - m_dot_steam`
  - `m_dot_fw` = feedwater flow (function of feedwater valve + pump head)
  - `m_dot_steam` = steam demand (function of turbine valve and pressure)
- Level indication: narrow-range (25–75%) and wide-range (0–100%).
- Trips: low-low level (e.g. 17%) trips reactor; high-high level (e.g. 85%)
  trips turbine.
- Emergency feedwater (EFW): auto-start on low-low level from motor-driven
  and turbine-driven pumps.

---

### 1.4 BWR Steam Drum Level — **missing**

In a BWR the reactor water level (above top of core) is the most important
controlled variable. It is controlled by feedwater flow. Low level exposes
the core. A SCRAM on low level is the reactor's primary self-protection.

**What to add:**
- State variable `level_m` [m above top of active fuel, range −0.5 to +1.5 m].
- Mass balance:
  `dm/dt = m_dot_fw + m_dot_recirc_makeup - m_dot_steam - m_dot_boil`
- Level trips: Low Level 1 alarm (−0.3 m), Low Level 2 SCRAM (−0.6 m),
  High Level 8 turbine trip (+1.2 m).
- Display: vertical level bar on HMI (very prominent).
- Auto feedwater control: proportional+integral on level setpoint.

---

### 1.5 Reactor Coolant Pump (RCP) Model — **absent**

Current model: `flow` is a scalar 0–1 that scales heat transfer coefficients.
No pump head/flow curve, no inertia, no coast-down after trip.

**What to add:**
- RCP state: `omega_rcp` (fractional speed, 0–1), inertia time constant
  ~10–15 s for coast-down.
- On pump trip: `domega/dt = -omega / tau_coast` (10–15 s time constant).
- Flow as function of speed: approximately linear for large RCPs.
- Natural circulation: residual flow ~3–5% after RCP coast-down, sustained
  by buoyancy (temperature-driven); implement as a floor on effective flow.
- Trip conditions: loss of voltage, manual trip, bearing temperature high.
- Multiple pump model: e.g. four RCPs in a PWR (each contributing ~25% flow),
  so 2-pump / 1-pump operation is trainable.

---

### 1.6 Non-linear Rod Worth Curve — **linear only**

Current model: `rho_rods = rod_position * rod_worth` (linear).
Real reactivity worth is strongly S-shaped: highest differential worth at
mid-core, near-zero at full in / full out.

**What to add:**
- Replace linear scaling with a cubic or tabulated S-curve:
  `w(x) = 3x² − 2x³` (integrated cosine profile), then
  `rho_rods = rod_worth * w(rod_position)`.
- For PWR: separate worth for gray (absorber) rods and black (control) rods.
- For RBMK: additional term for the graphite displacer tip — positive
  reactivity spike for the first ~1 m of insertion (this was a key factor
  in the Chernobyl accident). Model with a bump `+0.8 beta` when rod moves
  from 0 → 0.2 insertion before turning negative.

---

### 1.7 Moderator Temperature Coefficient (MTC) — **missing from PWR**

The PWR model has `fuel_temp_coeff` (Doppler) and `coolant_temp_coeff` but
does not separate moderator density feedback from Doppler. At End of Cycle
(EOC) the MTC becomes more negative. The model conflates these.

**What to add:**
- Explicit `moderator_temp_coeff` (typically −10 to −40 pcm/K at EOC,
  more negative than Doppler alone, `−3e-4` to `−1e-3` Δk/k per K).
- Tie this to a `burnup_fraction` parameter (0 = BOC, 1 = EOC) that scales
  the MTC more negative as burnup increases.
- Makes power defect and temperature defect trainable.

---

### 1.8 Samarium (Sm-149) Poisoning — **not modeled**

Samarium-149 builds in after shutdown and burns off slowly on restart.
Less dramatic than xenon but visible and trainable (affects reactor startup
after long shutdown).

**What to add:**
- Two ODEs for Pm-149 (promethium precursor) and Sm-149:
  - `dPm/dt = gamma_Pm * P_fis * n - lambda_Pm * Pm`
  - `dSm/dt = lambda_Pm * Pm - sigma_a_Sm * phi * Sm`
  - where `gamma_Pm ≈ 0.011`, `lambda_Pm ≈ 3.6e-6 /s`, neutron absorption
    removes Sm at `sigma_a_Sm * phi`.
- Reactivity: `rho_Sm = -0.007 * (Sm / Sm_eq - 1)` relative to equilibrium.
- Show Sm inventory in HMI (small secondary numeric readout).

---

### 1.9 Integration Stability — **Euler in PWR kinetics**

The PWR `PointKinetics.step()` uses first-order Euler integration. The BWR
and RBMK use RK2 (midpoint). Euler is sufficient for the substep sizes used
(~10 ms) but can accumulate errors.

**Fix:**
- Upgrade PWR `PointKinetics.step()` to the same RK2 (midpoint) method used
  in BWR/RBMK for consistency.
- Add analytical precursor update (matrix exponential for stiff step) as an
  optional path for substeps > 50 ms, which would allow a time-acceleration
  mode without changing substep count.

---

### 1.10 Primary Pressure — **no hydraulic coupling**

Current model: pressure is driven by steam inventory / heater ODE in
`_update_bop`, completely decoupled from primary coolant temperature.
Real primary pressure rises when primary temperature rises (thermal expansion)
and drops when temperature drops.

**Fix:**
- Link primary pressure change to coolant temperature via coolant expansion
  coefficient: `dP/dt += K_exp * dT_coolant/dt`
  where `K_exp ≈ 0.08–0.12 MPa/K` (representative of a pressurized water
  system with a pressurizer).
- Retain the pressurizer heater/spray model from section 1.2 as the dominant
  long-term pressure control.

---

### 1.11 Loss of Coolant Accident (LOCA) — **not modelable**

The current plant_supervisor only models pump degradation (flow reduction).
A LOCA is the design-basis accident for all power reactors.

**What to add:**
- `fault_loca_break_area` parameter (m²): sets a flow leak out of the primary.
- Leak flow: `m_dot_leak = Cd * A * sqrt(2 * rho * (P_primary - P_atm))`
- Primary inventory `m_primary` state variable; as it drops, flow drops and
  eventually core uncovers.
- ECCS actuation: high-pressure injection (HPIS, >~9 MPa), accumulator
  injection (~4 MPa), low-pressure injection (<~1 MPa). Simplified as step
  inputs of makeup flow when thresholds are crossed.
- This enables Loss of Coolant, Loss of Flow, and Anticipated Transient
  Without SCRAM (ATWS) scenarios.

---

## 2. Protection and I&C System

### 2.1 Trip Logic — **only 2 out of ~15 needed**

Current trips: high fuel temperature (>1400 K), high pressure. No flux trip,
no flow trip, no level trip, no overpower ΔT trip.

**What to add (minimum credible set):**

| Trip | Setpoint (typical) | Reactor |
|---|---|---|
| High neutron flux | >120% rated (prompt trip at >130%) | All |
| Overpower ΔT (OPDT) | f(Tin, pressure, flow) | PWR |
| High fuel temperature | >1600 K (cladding damage) | All |
| High reactor pressure | 16.8 MPa / 8.3 MPa | PWR / BWR |
| Low reactor coolant flow | <87% (2/4 pumps) | PWR |
| High coolant temperature | >616 K Thot | PWR |
| Low SG level (2/3 SGs) | <17% | PWR |
| BWR high-high level | Level 8 | BWR |
| BWR low-low level | Level 2 | BWR |
| RBMK: low pump flow | <60% nominal | RBMK |
| Manual SCRAM | operator SPACE | All |
| Turbine trip (loss of load) | auto-coupled | All |

Each trip should:
- Be latching (once tripped, requires manual reset).
- Have a separate `bypass` flag per channel (bypassable only with admin key).
- Be annunciated separately (trip name, time stamp).

---

### 2.2 2-out-of-3 (2oo3) Voting — **documented but not implemented**

The README mentions 2oo3 voting; the code does not implement it.
`plant_supervisor._protection()` uses single-channel comparisons.

**What to implement:**
- For each protection parameter, maintain three independent
  pseudo-channel readings with additive Gaussian noise
  (`sigma ≈ 0.5–1% of setpoint`).
- Trip logic: trip if 2 of 3 channels exceed setpoint.
- Single-channel failure: alarm "Channel A disagreement", no trip.
- 1-of-2 trip logic for manual actions.
- Bypass one channel: now running on 1oo2 — annunciate "Protective system
  degraded".
- Display channels A/B/C on dedicated I&C screen.

---

### 2.3 Alarm Prioritisation — **flat list**

Current alarms are a plain list appended in `_protection()`. No priorities,
no acknowledgement state, no timestamps, no shelving.

**What to implement:**
- `Alarm` dataclass with fields: `id`, `message`, `priority` (1=emergency/
  2=warning/3=advisory), `time`, `state` (unacknowledged/acknowledged/clear).
- `AlarmManager` class:
  - `raise_alarm(id, priority, message)` — sets alarm if not already active.
  - `clear_alarm(id)` — clears it (annunciator returns to normal).
  - `acknowledge(id)` — operator acknowledge; still displayed until clear.
  - `shelve(id, duration_s)` — suppresses alarm for duration (with logging).
- Priority 1 alarms should flash in HMI; priority 2 steady on; priority 3
  advisory bar.
- Keep a circular log of last 500 alarm events with timestamps.

---

### 2.4 SCRAM Reset Logic — **not implemented**

Pressing `L` is listed in the README but the key handler for `L` is not in
`run_control_room.py`. The `PlantSupervisor` also has no reset-permissive
logic.

**What to implement:**
- `PlantSupervisor.reset_scram()` method:
  - Requires all trip conditions to be cleared.
  - Requires manual rod insertion > 95% (rods are in).
  - Clears `plant.scrammed` flag, allowing rod withdrawal.
- `C` key: acknowledge active alarms.
- `L` key: attempt SCRAM reset (fails with message if conditions not met).
- On reset, HMI shows "SCRAM RESET APPROVED" or "RESET BLOCKED: <reason>".

---

### 2.5 Auto Control Loops — **not implemented**

Real reactors have several automatic control loops that the operator monitors
and can place in manual. Their absence makes the simulator unphysically
unstable to operate.

**Minimum auto loops to add:**

| Loop | Controlled variable | Actuator | Notes |
|---|---|---|---|
| Pressurizer pressure control | Primary pressure | Heater + spray | Proportional + integral |
| Pressurizer level control | Pressurizer level | Charging/letdown flow | P+I |
| SG level control (per SG) | SG level | Feedwater control valve | 3-element |
| Turbine speed/load control | Generator output | Turbine governor valve | Governor characteristic |
| BWR feedwater control | Drum level | Feedwater valve | 3-element |
| Rod auto-control (PWR) | Average coolant temp | Control rod group | Deadband ±1°C |

Each loop should be selectable auto/manual from the HMI. When in manual the
operator directly drives the actuator.

---

## 3. HMI / GUI

### 3.1 Window Size and Layout — **1280×760, single screen**

Current screen is too small and uses only plain text. A training simulator
needs at minimum a two-panel layout with instruments and trends, ideally
with multiple screens or tabs.

**Minimum changes:**
- Increase resolution to 1920×1080 (full HD) or allow resizable window.
- Divide screen into three vertical columns:
  - Left (280 px): control rod bank display, manual control inputs.
  - Center (860 px): main plant diagram (simplified P&ID), digital readouts,
    analog gauges.
  - Right (780 px): trend plots (4–6 trends, each 180 s window).
- Add a tab bar or function keys to switch between sub-screens:
  - F1: Main plant overview
  - F2: Primary system (pressurizer, RCP, pressure detail)
  - F3: Secondary system (SG levels, feedwater, turbine)
  - F4: Alarm summary / event log
  - F5: I&C / protection channel status
  - F6: Rod worth / reactivity breakdown

---

### 3.2 Trend Recorder — **~7 seconds, 4 parameters only**

`TREND_LEN = 400` at 60 Hz = 6.7 seconds. This is useless for training.
Operators need to see minutes of history.

**What to implement:**
- Separate slow historian (1 Hz) storing 60 minutes of data per parameter.
- Fast historian (10 Hz) storing 2 minutes per parameter.
- Display: selectable time window (30 s / 5 min / 30 min / 60 min).
- At least 8 trend parameters, user-selectable from a list including: power,
  thermal power, fuel temperature, coolant temperature, pressure, xenon
  inventory, boron (PWR), drum level (BWR), void fraction.
- Trend pen colors match ISA/ANSI standard (green = normal, yellow = alert,
  red = trip).

---

### 3.3 Analog / Digital Instrument Displays — **text only**

Current HMI shows all values as formatted text strings. Real control rooms
use panel meters, bar graphs, and circular gauges.

**What to implement (using Pygame drawing):**

1. **Vertical bar meters** (drawn with `pygame.draw.rect`):
   - Core power bar (0–120%) with colored bands: green 0–100, yellow
     100–110, red 110–120.
   - Reactor coolant pressure bar with operating range / trip lines.
   - SG level bars (one per SG).
   - BWR drum level bar with level marks L2/L8.

2. **Circular gauge** (drawn with `pygame.draw.arc` + pointer line):
   - Turbine speed/load (0–105%).
   - Rod position indicator arc (0 = out, 100% = in).

3. **Digital readouts** in DCS-style boxes:
   - Fixed-width monospace digits, green on black.
   - Color changes: green (normal) → yellow (alarm) → flashing red (trip).

4. **Reactor rod bank display**:
   - Horizontal bar per rod bank (A, B, C, D, SA).
   - Shows step position (0–228 steps for Westinghouse).

---

### 3.4 Simplified P&ID Diagram — **not present**

No flow diagram exists. Operators must be able to see where water/steam is
going in the plant.

**What to implement:**
- A schematic view (~700×400 px) with drawn lines representing:
  - PWR: reactor vessel → hot leg → SG primary side → cold leg → RCP → back
    to vessel. Secondary: SG secondary → steam header → turbine → condenser
    → feedwater pump → feedwater header → SG.
  - BWR: reactor vessel → steam line → turbine → condenser → feedwater →
    vessel. Recirculation loop.
  - RBMK: simplified channel-tube diagram.
- Live color/thickness on lines: blue (subcooled water), cyan (two-phase),
  white/gray (steam).
- Component boxes that flash when in alarm.
- Key valve positions shown as small squares (green=open, red=closed).

---

### 3.5 Alarm Annunciator Panel — **single flat list**

Current code: `snapshot.alarms + snapshot.trips` printed as a list.

**What to implement:**
- Grid of alarm tiles (e.g. 8×6 = 48 tiles), each 120×40 px.
- Each tile: fixed label, background color: black=normal, yellow=acknowledged,
  flashing amber=unacknowledged, flashing red=trip/emergency.
- Clicking a tile acknowledges it.
- "First-out" indicator: the first alarm to activate after plant is normal
  gets a "FO" marker (helps identify initiating event).
- Sound: `pygame.mixer.Sound` plays a tone on new unacknowledged alarm.
- `C` key = acknowledge all visible alarms.

---

### 3.6 Event Sequence Recorder (ESR) — **not present**

Critical for post-exercise review.

**What to implement:**
- A `EventLog` list of `(sim_time, wall_time, category, message)` tuples.
- Automatically log:
  - All alarms raised/cleared/acknowledged.
  - All operator key presses with affected parameter and value.
  - All auto-SCRAM events with trip cause.
  - Protection channel disagreements.
- Display on F4 screen: scrollable table, most recent at top.
- Export to CSV on `ESC` / quit.

---

### 3.7 Xenon / Samarium Trend Screen — **xenon not visible enough**

Xenon oscillations occur over hours. The simulator should support viewing
the full xenon buildup after shutdown.

**What to implement:**
- On the reactivity/poisons sub-screen, show:
  - Xenon inventory vs. time (last 12 hours of slow historian).
  - Iodine inventory vs. time.
  - Samarium inventory vs. time (after Sm model is added).
  - Reactivity budget bar chart: rod worth, Doppler, MTC, xenon, samarium,
    boron (PWR) — stacked bars showing each component's Δk/k.
- Time-acceleration button (10×, 100×, 1000×) that runs the physics faster
  for slow processes like xenon recovery. Flag this clearly as "fast time".

---

### 3.8 SCRAM / Safety System Indication — **unclear**

Current: `SCRAM` word appears in alarms list. No clear top-of-screen
safety status indicator.

**What to implement:**
- Persistent safety function status bar at top of every screen:
  - `[REACTOR: CRITICAL]` / `[REACTOR: SHUTDOWN]` / `[REACTOR: SCRAMMED]`
    in large bold text with appropriate background color.
  - `[SAFETY FUNCTION: OK]` / `[SAFETY FUNCTION: CHALLENGED]`.
  - `[SCRAM LATCHED — rods in — press L to evaluate reset]`.
- Rod position summary: "All rods IN" / "Rods withdrawn X%" with visual bar.
- Flashing red banner for unresolved trip conditions.

---

## 4. Training Infrastructure

### 4.1 Initial Condition (IC) Save and Load — **not present**

Simulator training always uses saved initial conditions to reproduce scenarios.

**What to implement:**
- `save_ic(filepath)`: pickle (or JSON-serialize) the entire
  `PlantSupervisor` state including all physics engine internal variables,
  time, alarms, and historian.
- `load_ic(filepath)`: restore the state.
- Default ICs to ship:
  - `ic_100pct_pwr.pkl` — full-power steady state PWR
  - `ic_50pct_pwr.pkl` — part-power PWR (xenon built in)
  - `ic_100pct_bwr.pkl` — full-power BWR
  - `ic_xenon_peak.pkl` — PWR at xenon peak (post-shutdown 9h)
- Key bindings: `Ctrl+S` to save current state to `ic_quicksave.pkl`,
  `Ctrl+L` to reload last quicksave.

---

### 4.2 Time Acceleration — **not present**

Xenon transients occur over 8–24 hours of simulated time. Without time
acceleration a trainee cannot practice xenon management.

**What to implement:**
- `time_speed_factor` (default 1.0, options: 1, 10, 60, 600, 3600).
- When > 1: multiply each frame's `dt` by the factor before calling
  `supervisor.step(dt)`. Increase internal substep count proportionally
  to maintain physics stability.
- Display current time factor prominently: "Sim time: 1× real" or "10×".
- Disable time acceleration automatically if any alarm is active (prevents
  skipping over developing transients).

---

### 4.3 Instructor Override Panel — **not present**

A training simulator without instructor controls requires the student to
cause every scenario manually.

**What to implement (in-process, no separate window needed initially):**
- `I` key opens an overlay "Instructor panel" with:
  - Fault injection: pump trip, loss of feedwater, stuck valve, instrument
    failure (fixed reading), break size (LOCA 0–100%).
  - Parameter override: force a parameter to a specific value for one step
    (e.g. artificially spike fuel temperature to test protection).
  - Freeze/resume: pause all physics while the instructor explains.
  - Reset to IC: reload last saved initial condition.
  - Insert reactivity: inject a step Δρ for SRT (shutdown reactivity test).
- Log all instructor actions in the event recorder.

---

### 4.4 Scenario Scripting — **not present**

Allows predefined transients to be run repeatably for competency testing.

**What to implement:**
- Simple text-based scenario file format:
  ```
  t=0       load_ic ic_100pct_pwr.pkl
  t=30      set fault_pump_degraded True
  t=120     expect trip "AUTO SCRAM: high fuel temperature"
  t=300     end
  ```
- Scenario runner that executes the script in headless or GUI mode.
- Pass/fail evaluation based on `expect` conditions and operator responses.
- Integrate with `validation.py` framework.

---

### 4.5 Performance Assessment — **not present**

**What to implement:**
- Track operator actions: number of unnecessary control movements, time to
  respond to alarms, correct procedure adherence.
- Score rubric: configurable per scenario.
- Post-exercise report printed to console or saved as text file.

---

## 5. Reactor-Specific Physics Gaps

### 5.1 RBMK — Positive Scram Effect (AZ-5 / EPS)

The Chernobyl accident was exacerbated by a positive reactivity insertion
when the emergency protection rods (AZ-5) began inserting. The graphite
displacer tips entered the core before the absorber portion, adding ~2–3 β
of reactivity before reversing.

**What to add:**
- In `RBMKPlant.compute_reactivity`: detect when `rod_pos_effective` is
  increasing from < 0.2 (rod entering core tip first).
- Add a positive spike term:
  `rho_tip = +tip_worth * max(0, 0.2 - rod_pos_effective) / 0.2`
  where `tip_worth ≈ +2 * beta_total ≈ +0.014` Δk/k.
- This makes a realistic RBMK low-power SCRAM scenario possible.

---

### 5.2 BWR — Stability (Density Wave Oscillation)

At low flow / low power, BWRs can exhibit channel density wave oscillations
(Type I: core-wide, Type II: regional). Current model has no mechanism for
this.

**What to add:**
- Transport delay in void fraction response: `alpha(t) = f(power(t - tau))`
  where `tau ≈ channel_length / steam_velocity ≈ 1–3 s`.
- Implement as a ring buffer of power history to compute delayed void.
- This produces oscillatory void fraction response at low flow, visually
  demonstrable on the trend display.

---

### 5.3 PWR — Reactivity Feedback at Low Power

At low power (< ~5%) point kinetics initialised at zero xenon/samarium
(fresh startup from cold) show very different behaviour than warm restarts.
The current `XenonIodine.__init__` sets equilibrium at `n=1` even for a
cold startup.

**Fix:**
- Add `cold_startup` flag to `PWRPlant.__init__`. If True, initialise
  xenon/iodine to 0.
- Initialise decay heat to `initial_power=0.0` for a truly cold startup.

---

## 6. Parameter Tuning Required

Several parameters in the current code produce unrealistic dynamics on
inspection:

| Parameter | Current value | Realistic value | Note |
|---|---|---|---|
| `PWRPlant.rod_worth` | −0.05 Δk/k | −0.06 to −0.08 total | Too low; 5000 pcm for ~15 banks |
| `BWRParams.internal_dt` | 0.001 s | 0.001–0.005 s OK | Fine but slow; RK2 allows 0.005 |
| `RBMKParams.void_coeff` | +0.08 Δk/k | +0.03 to +0.05 total | +8 pcm/%void is correct order |
| `PlantParams.xenon_burn_coeff` | 0.08 | ~0.08 at full power | Correct; verify saturation Xe at EOL |
| `PlantParams.xenon_reactivity_coeff` | 0.02 | ~0.028 for PWR full core | Slightly low; −2800 pcm max Xe worth |
| `PlantParams.coolant_heat_capacity` | 4.0e7 J/K | ~1.0e8 J/K | Too small; slows thermal time constant |
| BWR `rod_drop_tau` | 2.5 s | 3–5 s for hydraulic blades | Close but toward lower end |
| RBMK `rod_drop_tau` | 3.0 s | 12–18 s for RBMK (pre-upgrade) | Significantly too fast |
| `h_fuel_to_coolant` (PWR) | 1.0e8 W/K | Calibrate to 20 s fuel-to-coolant τ | Check: τ = C_fuel / h |

RBMK rod insertion time of 3 s is particularly wrong — pre-Chernobyl RBMK
rods took 18–20 s to fully insert. This is trainably significant and should
be corrected to `rod_drop_tau = 18.0` for RBMK.

---

## 7. Implementation Priority Order

| Priority | Item | Files affected |
|---|---|---|
| 1 | RBMK rod insertion time fix (3 s → 18 s) | `physics_engine_rbm.py` |
| 2 | Non-linear rod worth S-curve | All three physics engines |
| 3 | RBMK graphite tip positive reactivity | `physics_engine_rbm.py` |
| 4 | SCRAM reset key `L` — add to `run_control_room.py` | `run_control_room.py` |
| 5 | Longer trend window (historian at 1 Hz, 60 min) | `run_control_room.py` |
| 6 | More trips (flux, flow, level) | `plant_supervisor.py` |
| 7 | Alarm annunciator panel (tiles) | `run_control_room.py` |
| 8 | BWR drum level state variable | `physics_engine_bwr.py`, `plant_supervisor.py` |
| 9 | Boron concentration (PWR) | `physics_engine.py`, `plant_supervisor.py` |
| 10 | Pressurizer dynamic model | `plant_supervisor.py` |
| 11 | RCP coast-down model | `plant_supervisor.py` |
| 12 | SG level model (PWR) | `plant_supervisor.py` |
| 13 | 2oo3 voting with noisy channels | `plant_supervisor.py` |
| 14 | IC save/load | `plant_supervisor.py`, `run_control_room.py` |
| 15 | Time acceleration | `run_control_room.py` |
| 16 | Auto control loops (P control) | `plant_supervisor.py` |
| 17 | Simplified P&ID diagram | `run_control_room.py` |
| 18 | Multi-screen HMI (F1–F6 tabs) | `run_control_room.py` |
| 19 | Event sequence recorder | `plant_supervisor.py`, `run_control_room.py` |
| 20 | Samarium model | All three physics engines |
| 21 | LOCA fault | `plant_supervisor.py` |
| 22 | BWR stability / density wave delay | `physics_engine_bwr.py` |
| 23 | Instructor override panel | `run_control_room.py` |

---

## 8. What This Simulator Is Already Doing Well

For context, the following are correctly implemented and do not need rework:

- Six-group point kinetics with correct U-235 β and λ values.
- Four-group ANS-5.1 decay heat (7% at shutdown, correct time constants).
- Xenon/iodine two-ODE model with burn-out.
- SCRAM rod-drop dynamics (except RBMK time constant).
- Positive void coefficient in RBMK, negative in BWR.
- Doppler fuel temperature feedback.
- Mechanistic void fraction ODE in PlantSupervisor (pressure-dependent
  bubble collapse and boil source terms).
- Decay heat included in thermal load to steam generator.
- Startup permit blocking rod withdrawal.
- Turbine trip auto-coupling to turbine valve.
- Fault injection for pump degradation and feedwater loss.
