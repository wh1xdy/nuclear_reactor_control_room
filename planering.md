# Nuclear Reactor Control Room Simulator

This repository provides a real-time nuclear power plant simulator with an interactive control room user interface and multiple realistic physics engines.

## Overview

This project simulates the dynamics of a nuclear reactor and lets users operate a virtual control room. It is designed to be educational and engaging without requiring a supercomputer. The core components are:

- **Physics engines** modelling different reactor types:
  - **Pressurized Water Reactor (PWR)** – includes six-group point kinetics, two-node thermal-hydraulics, xenon poisoning, and a steam generator with turbine/generator loop.
  - **Boiling Water Reactor (BWR)** – models two-phase coolant in the core, void fraction feedback, xenon poisoning, and a direct steam cycle with reactor pressure dynamics.
  - **RBMK-type Reactor** – graphite moderated with positive void coefficient, including separate fuel and graphite nodes, two-phase coolant, and xenon/iodine transients.

- **UI** built in Pygame that acts as a simplified control room:
  - Choose reactor type (PWR/BWR/RBMK).
  - Control rod position, pump/recirc flow, turbine valve opening.
  - Trigger a SCRAM (rapid shutdown).
  - View live readouts (power, temperature, pressure, reactivity, xenon).
  - Basic time trend graphs placeholder for future enhancements.

- **Packaging**: A `pyproject.toml` describes dependencies and entry points. A GitHub Actions workflow builds standalone macOS DMG installers for both Intel and Apple Silicon via PyInstaller and dmgbuild.

## Goals

This project aims to:

1. **Provide an educational nuclear simulator** that strikes a balance between realism and playability. We use point kinetics and lumped-parameter thermal models to capture the qualitative behaviour of real plants.

2. **Support multiple reactor types**. The current engines simulate PWR, BWR, and RBMK designs. Additional models can be added by implementing the common `step(dt, controls) -> snapshot` API.

3. **Allow headless physics testing**. The physics modules can be executed without a UI, making it easy to unit test and validate the models or run long transients quickly.

4. **Build cross-platform desktop applications**. With the provided GitHub Actions workflow, macOS binaries (.app) and DMG installers are generated automatically for both arm64 and x86_64 architectures. Support for Windows and Linux builds can be added later.

5. **Lay the foundation for a modern web UI**. While the Pygame control room is usable, we plan to integrate the physics engines with a React/TypeScript frontend and FastAPI backend to create a browser-based dashboard with charts and alarms.

## How to run

After cloning, install the dependencies and run the control room:

```bash
python -m pip install -e .
nrcr
```

The application launches a window with keyboard controls:
- `1` `2` `3` – switch between PWR, BWR and RBMK.
- `W/S` – move control rods up/down.
- `A/D` – change primary flow / recirc pump / pump speed.
- `Q/E` – adjust turbine valve.
- `SPACE` – SCRAM (emergency shutdown).
- `R` – reset the current reactor.

## Building installers

GitHub Actions is configured to build macOS DMG installers on tag pushes. To build locally:

```bash
pip install pyinstaller dmgbuild
pyinstaller --windowed --name NuclearReactorControlRoom run_control_room.py
python -c "import dmgbuild; dmgbuild.build_dmg('NRCR-intel.dmg', 'NuclearReactorControlRoom', settings)"
```

See `.github/workflows/build-macos.yml` for details.

## Contributing

This is a work-in-progress project. Contributions are welcome! Possible areas:
- Enhancing the user interface and adding graphs.
- Adding more realistic thermal-hydraulic models or new reactor types.
- Porting the UI to a web technology stack.
- Improving documentation.

Please open issues or pull requests.

## License

MIT License.
