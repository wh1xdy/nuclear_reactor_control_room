# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Target hardware

**Primary development machine: Apple MacBook Pro 16" M5 Pro (2025)**

All performance decisions should be made with this machine as the baseline:
- Apple Silicon arm64 — no x86 assumptions; use `math` / `numpy` over `ctypes`/SIMD hacks
- Performance cores: 12 (6P + 6E); the main loop runs on one P-core at 60 Hz
- Unified memory: 48 GB — large historian deques and ring buffers are fine; don't micro-optimise allocation
- GPU: 20-core Apple GPU shared with CPU via unified memory — Pygame uses Metal via SDL2; avoid excessive `pygame.Surface` copies per frame
- Python: CPython arm64 (Homebrew or python.org universal2 build); numpy and scipy link against Accelerate framework automatically

**Performance targets on M5 Pro:**
- Physics step (`supervisor.step(dt)`) must complete in < 1 ms at 1× real-time (budget: ~16 ms/frame at 60 Hz, physics should use < 6%)
- At 600× time acceleration the physics loop runs ~600 substeps per frame; keep each `step()` call well under 10 µs
- Render pass (all `draw_*` calls) must stay under 8 ms per frame
- No busy-wait loops; use `clock.tick(60)` to yield to the OS scheduler between frames

**What to avoid:**
- `time.sleep()` in the main loop (use `clock.tick`)
- Per-frame Python `list` allocations inside the physics engines — mutate in-place or use `deque`
- Unnecessary numpy array creation inside `step()` — the engines are pure-Python scalar math by design; keep them that way unless profiling shows a bottleneck
- Spawning threads or subprocesses from the physics engines

## Running the simulator

```bash
# Install dependencies (in a venv)
python -m venv .venv && source .venv/bin/activate
pip install -e .                  # installs pygame, numpy (and registers the `nrcr` entry point)
pip install iapws                 # optional: more accurate steam/water properties

# Run
nrcr
# or
python run_control_room.py
```

## Building a macOS DMG (local)

```bash
pip install pyinstaller dmgbuild
pyinstaller --noconfirm --windowed --name "NuclearReactorControlRoom" --collect-all pygame run_control_room.py
```

The CI workflow (`.github/workflows/build-macos.yml`) triggers on `v*` tags and produces a universal (arm64 + x86_64) DMG via the `build-universal-dmg` job.

## Architecture

The project is split into two completely independent layers:

### Physics engines (no Pygame dependency)

Each engine is a self-contained module that exposes this pattern:
- `*Params` dataclass – tunable physical constants
- `*ControlInputs` dataclass – operator inputs (rod position 0–1, flow/pump speed 0–1, turbine valve 0–1, scram bool)
- `*Plant` class – integrates the ODE each call to `step(dt, controls) -> snapshot`
- Snapshot dataclass – observable state returned after every step

| Module | Reactor | Unique features |
|---|---|---|
| `physics_engine.py` | PWR | Two-node thermal-hydraulics, steam generator, iapws optional |
| `physics_engine_bwr.py` | BWR | Void fraction feedback, direct steam cycle, reactor pressure |
| `physics_engine_rbm.py` | RBMK | Graphite moderator node, positive void coefficient, three thermal nodes (fuel/graphite/coolant) |

All three engines share the same underlying dynamics:
- **Six-group point kinetics** (U-235 β_i and λ_i values)
- **Xenon/iodine poison model** (hours-long negative reactivity swings after power changes)
- **Doppler fuel temperature feedback** (negative, stabilising)

### UI layer (`run_control_room.py`)

A single-file Pygame application running at 1920×1080, 60 Hz. It instantiates one of the three `*Plant` objects via `PlantSupervisor` and calls `supervisor.step(dt)` once per frame. Reactor switching (`1/2/3`) replaces the plant object and resets sliders. The six F-key tabs (F1–F6) render different views of the same snapshot; only the active tab is drawn each frame.

## Keyboard controls (runtime)

| Key | Action |
|---|---|
| `1` / `2` / `3` | Switch to PWR / BWR / RBMK |
| `W` / `S` | Control rods up / down |
| `A` / `D` | Primary flow / recirc pump / pump speed |
| `Q` / `E` | Turbine valve open / close |
| `F` / `V` | Feedwater valve open / close |
| `H` / `N` | Pressurizer heater up / down |
| `B` / `G` | Boration rate up / down (PWR) |
| `U` / `I key` | Dilution rate up / down (PWR) |
| `SHIFT` | Fast adjustment (5× step) |
| `SPACE` | SCRAM (emergency shutdown) |
| `L` | Reset SCRAM (requires rods in, all trips clear) |
| `C` | Acknowledge all alarms |
| `P` | Toggle startup permit |
| `T` | Toggle turbine trip |
| `Z` | Toggle pump-degraded fault |
| `X` | Toggle feedwater-loss fault |
| `I` | Cycle LOCA break area (0 → 10 cm² → 100 cm²) |
| `O` | Toggle auto rod control — PWR (maintains 550 K avg coolant) |
| `M` | Toggle auto pressurizer control — PWR (maintains 15.5 MPa) |
| `[` / `]` | Cycle trend time window (10 s / 5 min / 30 min / 60 min) |
| `+` / `-` | Increase / decrease time acceleration (1× / 10× / 60× / 600×) |
| `F1`–`F6` | Switch display tab (Overview / Primary / Secondary / Alarms / I&C / Reactivity) |
| `Ctrl+S` | Save initial condition to `ic_quicksave.pkl` |
| `Ctrl+L` | Load initial condition from `ic_quicksave.pkl` |
| `R` | Reset current reactor to initial state |
| `ESC` | Quit (exports `event_log.csv`) |

## Adding a new reactor type

1. Create `physics_engine_<name>.py` following the `step(dt, controls) -> snapshot` API.
2. Register the module in `pyproject.toml` under `[tool.setuptools] py-modules`.
3. Add a `select()` branch and display section in `run_control_room.py`.
