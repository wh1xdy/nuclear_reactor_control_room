# Nuclear Reactor Control Room Simulator

A real-time nuclear reactor control room simulator written in Python and Pygame. Supports three reactor types — PWR, BWR, and RBMK — each with physics-informed models of reactor kinetics, thermal feedback, xenon poisoning, and plant-specific behaviour.

## Features

**Three reactor models:**
- **PWR** – two-node thermal-hydraulics, steam generator, optional IAPWS water properties
- **BWR** – void fraction feedback, axial two-phase channel model, direct steam cycle, reactor pressure
- **RBMK** – graphite moderator node, positive void coefficient, three thermal nodes (fuel / graphite / coolant)

**All three engines share:**
- Six-group point kinetics (U-235 β and λ values)
- Xenon/iodine poison model (hours-long reactivity swings after power changes)
- Doppler fuel temperature feedback (negative, stabilising)

**Supervisor / protection layer (on top of physics engines):**
- Simplified balance-of-plant: pressure, steam inventory, condenser, feedwater
- Instrumentation & Control channels with 2-out-of-3 voting, filtering, and fault detection
- Alarm management and latched SCRAM with startup interlocks
- Injectable faults: pump degradation, feedwater loss
- Axial two-phase void and pressure-drop model for BWR/RBMK channels

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Optional – more accurate steam/water properties:

```bash
pip install iapws
```

## Running

```bash
nrcr
```

or directly:

```bash
python run_control_room.py
```

## Keyboard Controls

| Key | Action |
|---|---|
| `1` / `2` / `3` | Switch to PWR / BWR / RBMK |
| `W` / `S` | Control rods up / down |
| `A` / `D` | Primary flow / recirculation pump speed |
| `Q` / `E` | Turbine valve open / close |
| `F` / `V` | Feedwater valve open / close |
| `H` / `N` | Pressurizer heater on / off |
| `P` | Startup permit interlock toggle |
| `T` | Turbine trip toggle |
| `Z` | Inject / clear pump fault |
| `X` | Inject / clear feedwater fault |
| `C` | Acknowledge active alarms |
| `L` | Reset SCRAM latch (permissive matrix) |
| `B` / `G` | Bypass high-pressure / high-fuel-temperature trip |
| `I` | Inhibit auto-SCRAM |
| `SHIFT` | Fast adjustment (5× step size) |
| `SPACE` | SCRAM (emergency shutdown) |
| `R` | Reset current reactor |
| `ESC` | Quit |

## Architecture

The project is split into two independent layers:

### Physics engines (no Pygame dependency)

Each engine is a self-contained module that exposes:
- `*Params` dataclass – tunable physical constants
- `*ControlInputs` dataclass – operator inputs (rod position 0–1, flow 0–1, turbine valve 0–1, scram bool)
- `*Plant` class – integrates the ODE each call to `step(dt, controls) -> snapshot`
- Snapshot dataclass – observable state returned after every step

| Module | Reactor | Unique features |
|---|---|---|
| `physics_engine.py` | PWR | Two-node thermal-hydraulics, steam generator, iapws optional |
| `physics_engine_bwr.py` | BWR | Void fraction feedback, direct steam cycle, reactor pressure |
| `physics_engine_rbm.py` | RBMK | Graphite moderator node, positive void coefficient, three thermal nodes |

### UI layer (`run_control_room.py`)

A single-file Pygame application running at 60 Hz. It instantiates one `*Plant` object at a time, calls `plant.step(dt, controls)` on each frame, and drives the display. Switching reactor type (`1`/`2`/`3`) replaces the plant object and resets sliders.

### Supervisor layer (`plant_supervisor.py`)

Wraps the physics engines with a protection and control layer: balance-of-plant coupling, I&C channel logic, alarm/trip management, and fault injection.

## Validation

```bash
python validation.py
```

Verifies qualitatively that:

1. SCRAM reduces power across all three reactor types.
2. Protection systems and trips behave deterministically.
3. BWR/RBMK transients run correctly with the axial two-phase void/pressure-drop model.
4. Reference transients are compared against benchmark envelopes with uncertainty bands (`data/benchmark_envelopes.json`).

See `docs/model_limits_and_vv.md` for a discussion of model scope, known limitations, and recommended next V&V steps.

## Realism and Limitations

This is an educational simulator, not a certified training or licensing tool. It captures qualitative transient behaviour and relative trends for control actions around normal and upset conditions. It does not model full 3D neutronics, detailed CFD, complete plant hydraulic networks, or regulatory acceptance criteria.

## macOS Distribution

The repository includes a GitHub Actions workflow (`.github/workflows/build-macos.yml`) that triggers on `v*` tags and builds two `.app` bundles (x86_64 + arm64), merges them into a universal `.app`, and packages a DMG artefact.

To build locally:

```bash
pip install pyinstaller dmgbuild
pyinstaller --noconfirm --windowed --name "NuclearReactorControlRoom" --collect-all pygame run_control_room.py
```

Remaining steps for a full end-user distribution chain:

1. **Code signing** of `.app` and `.dmg` (`codesign`, Developer ID)
2. **Notarization** + stapling (`notarytool` + `stapler`)
3. **Release publishing** – attach DMG to a GitHub Release

## Adding a New Reactor Type

1. Create `physics_engine_<name>.py` following the `step(dt, controls) -> snapshot` API.
2. Register the module in `pyproject.toml` under `[tool.setuptools] py-modules`.
3. Add a `select()` branch and display section in `run_control_room.py`.

## Future Work

- Additional reactor models (CANDU, AP1000, etc.)
- Real-time trend charts and improved panel layouts
- More detailed thermal-hydraulic and pressure modelling
- Windows and Linux packaging workflows
- Web-based dashboard UI
