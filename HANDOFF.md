# Nuclear Reactor Control Room — UI Rewrite Handoff

## What this project is

A Pygame-CE nuclear reactor training simulator running on macOS M5 Pro (arm64).
Three reactor types: PWR, BWR, RBMK-1000. Fully working physics engine.
The UI needs a **complete visual redesign** — not patches, a ground-up rewrite of all drawing code.

**Run it:**
```bash
source .venv/bin/activate && python run_control_room.py
```

---

## What the user asked for (approved design brief)

1. **Tabbed layout** — six tabs: Overview / Primary / Secondary / Reactivity / Alarms / I&C
2. **Dark industrial DCS aesthetic** — Foxboro/ABB style. Near-black background, muted borders, ISA status colors (green/amber/red). NOT green-tinted, NOT amateurish.
3. **Electric blue (#3090C0) as the accent** for selected/active elements
4. **Monospace engineering readouts** — Menlo font, looks like a real DCS terminal
5. **P&ID flow diagrams** — interactive, orthogonal pipes only (NO diagonals), actually pretty, animated flow dots showing fluid direction
6. **Draggable sliders** for rod position, flow, valves etc — left panel, always visible
7. **Both P&ID and bar meters** — integrated properly, not crammed together
8. **ALL of these decorations:**
   - CRT scanline texture (subtle, every other row)
   - Corner bracket accents on every panel (L-shaped, in accent color)
   - Animated flow dots on P&ID pipes (pulsing dots traveling along pipes, speed ∝ flow rate)
   - Glowing status LEDs (filled circle + larger alpha halo)
   - Grid-line background on the main body area (subtle, ~48px spacing)

---

## Current state of the code

The code has been partially modified but does NOT yet achieve the design brief. What exists now:
- Working physics and supervisor (DO NOT touch these)
- Partial decorations: scanlines, corner brackets, grid bg, animated flow dots on PWR primary loop, glowing LEDs, improved tab bar
- The overall layout and P&ID are still ugly — the user said "this is not near what I asked for"
- The P&ID still has some diagonal pipes (BWR, RBMK), labels that clip, boxes too small
- The left panel controls work (draggable sliders, button-style toggles) but look basic
- The right panel trends and bar meters are functional but visually weak

**The user wants a real redesign — not incremental patches. Start from scratch on the drawing layer.**

---

## File structure

| File | Purpose | Touch? |
|------|---------|--------|
| `run_control_room.py` | All UI drawing + main loop (~1804 lines) | YES — full rewrite of drawing code |
| `plant_supervisor.py` | Physics wrapper, BOP, protection system | NO |
| `physics_engine.py` | PWR point kinetics | NO |
| `physics_engine_bwr.py` | BWR physics | NO |
| `physics_engine_rbm.py` | RBMK physics | NO |
| `CLAUDE.md` | Project guidance | Read it |

---

## Key snapshot fields (from `UnifiedSnapshot`)

```python
snap.power_fraction        # 0.0–1.0+ (>1 = over power)
snap.rod_position          # 0.0 (withdrawn) – 1.0 (inserted)
snap.reactivity            # dk/k, ~0 at steady state
snap.fuel_temp_k           # Kelvin, trip >1400 K
snap.coolant_temp_k        # Kelvin, trip >616 K (PWR)
snap.outlet_temp_k
snap.primary_pressure_pa   # Pascals; snap.pressure_mpa also available
snap.steam_pressure_pa
snap.steam_temp_k
snap.condenser_temp_k
snap.steam_flow_kgs
snap.feedwater_inventory   # 0–1
snap.steam_inventory       # 0–1
snap.flow                  # primary flow fraction 0–1
snap.void_fraction         # BWR/RBMK, 0–1
snap.xenon_worth           # negative dk/k
snap.iodine_conc
snap.xenon_conc
snap.reactor_type          # "PWR", "BWR", "RBMK"
snap.boron_ppm             # PWR only
snap.thermal_mw
snap.electric_mw
snap.decay_heat_fraction
snap.time                  # simulation time in seconds
snap.alarms                # list of strings
snap.trips                 # list of strings
snap.alarm_objects         # list of Alarm objects with .message, .state, .is_trip, .first_out
snap.scrammed              # bool
snap.porv_open             # bool
snap.eccs_actuated         # bool
snap.turbine_valve         # 0–1
snap.channel_noise         # list of 3 floats for I&C channels
snap.orm                   # RBMK only — operational reactivity margin (rods)
snap.skala_orm             # RBMK only
snap.skala_age_s           # RBMK only
```

---

## Current color palette (keep this)

```python
C_BG        = (10,  12,  14)
C_PANEL     = (16,  19,  22)
C_PANEL_HDR = (20,  24,  28)
C_BORDER    = (44,  50,  56)
C_TEXT      = (215, 218, 215)
C_TEXT_DIM  = (85,  94,  88)
C_TEXT_HDR  = (135, 145, 140)
C_NORMAL    = (55,  185,  75)   # ISA green
C_CAUTION   = (205, 160,  18)   # ISA amber
C_WARNING   = (195,  95,  18)   # ISA orange
C_ALARM     = (205,  38,  38)   # ISA red
C_ACCENT    = (48,  125, 195)   # electric blue
C_WATER     = (38,   95, 185)
C_TWOPHASE  = (28,  145, 185)
C_STEAM     = (105, 118, 128)
C_HOT_LEG   = (185,  65,  28)
C_SLIDER_BG = (28,  32,  36)
C_SLIDER_FG = (48, 125, 195)
```

---

## Existing helpers to keep / improve

```python
_s(n)                   # scale design pixels to native — used everywhere, keep it
SCALE                   # set from screen.get_size()[0] / 1920.0
draw_arc_gauge(...)     # circular arc gauge — keep, it's good
draw_bar_meter(...)     # vertical bar — keep, improve visually
draw_dcs_box(...)       # DCS readout box — keep
draw_trend(...)         # strip chart trend — keep
draw_valve(...)         # butterfly valve symbol — keep
draw_pump_symbol(...)   # IEC pump circle — keep
draw_led(...)           # glowing LED — keep (just improved)
draw_corner_brackets(...)  # panel corner L-brackets — keep
draw_scanlines(...)     # CRT overlay — keep
draw_grid_bg(...)       # body grid lines — keep
```

---

## P&ID design rules (CRITICAL)

1. **All pipes must be horizontal or vertical only** — zero diagonal lines
2. **L-shaped routing** — two segments per connection (horizontal then vertical, or vice versa)
3. **Component boxes** — minimum 90px logical width, label drawn LAST (on top of fills)
4. **Pipe colors:**
   - Hot primary coolant: `C_HOT_LEG` (orange-red)
   - Cold/subcooled: `C_WATER` (blue)
   - Two-phase/boiling: `C_TWOPHASE` (cyan)
   - Steam/vapor: `C_STEAM` (grey)
5. **Animated flow dots** — dots travel along each pipe segment, speed ∝ flow rate
6. **Interactive** — clicking a component box highlights it and shows a detail overlay
7. **Labels** — always uppercase, drawn after fills so they're never occluded

### PWR recommended layout (logical px, relative to draw area)
```
REACTOR:  x=30,  y=60,  w=110, h=180
PZR:      x=30,  y=8,   w=80,  h=48
S/G:      x=240, y=40,  w=90,  h=180
RCP:      x=155, y=270, w=80,  h=44
TURB:     x=430, y=40,  w=90,  h=55
COND:     x=430, y=260, w=90,  h=50
```

---

## Slider & controls spec

- Left panel always visible, ~300px wide
- Sliders: track 16px tall, hit rect 32px centered on track, blue fill, white thumb when dragging
- Value in engineering units: `"75 %"` not `"0.750"`
- Toggle buttons: bordered rect per button, ON = faint colored bg + bright border, OFF = C_PANEL + C_BORDER
- SCRAM button: big red at bottom of left panel

---

## What "interactive P&ID" means

- Each component box is a clickable `pygame.Rect`
- Store rects in a list each frame: `pid_component_rects = [(rect, component_name), ...]`
- On click: set `selected_component = component_name`
- Draw a detail overlay (small panel) near the clicked component showing key values
- Example: click REACTOR → show fuel temp, power, rod pos, reactivity in overlay

---

## Architecture rule

**The physics layer is not to be touched.** The supervisor `.step(dt)` call produces a snapshot every frame. All UI work is read-only from the snapshot. The only writes from the UI to physics are through `supervisor.controls` fields (which the existing keyboard/slider handlers already update correctly).

---

## How to approach the rewrite

1. Read `run_control_room.py` in full — understand what exists
2. Keep: `main()` skeleton, event loop, keyboard handlers, supervisor wiring, Slider class, historian deques
3. Rewrite: all `draw_*` functions, panel layout, P&ID schematics, tab content functions
4. Test each tab renders without crash before moving to the next
5. Syntax check: `python -c "import ast; ast.parse(open('run_control_room.py').read()); print('OK')"`
6. Quick launch test: run for 4 seconds, check no crash

---

## The user's exact words on what they want

> "dark industrial dcs with cool looking decorations, and scientific stats"
> "do flowcharts, interactive and actually pretty"
> "both [P&ID and meters], but integrated properly and actually pretty"
> "engineering readouts [monospace]"
> Decorations: "yes all of them" (scanlines, corner brackets, animated flow dots, glowing LEDs, grid-line background)
