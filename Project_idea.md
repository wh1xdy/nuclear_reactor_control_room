# Nuclear Reactor Control Room Simulation

## Overview

This repository contains a nuclear reactor control room simulator that provides simplified yet physics-informed models of three different types of reactors: a Pressurized Water Reactor (PWR), a Boiling Water Reactor (BWR), and an RBMK graphite moderated reactor. Each reactor is implemented as a standalone Python module that calculates reactor kinetics, thermal feedback, xenon poisoning and plant specific behaviour.

The simulator exposes a consistent API where each plant accepts control input (e.g. rod position, primary flow, recirc pump speed, turbine valve position, SCRAM) and returns a snapshot of the current state: power fraction, thermal/electric power, fuel and coolant temperatures, void fraction (for BWR and RBMK), reactivity, xenon inventory and other signals. This modular design allows different reactor types to be swapped without changing the UI.

A basic control room UI has been implemented using Pygame. It allows the user to select the reactor type (PWR/BWR/RBMK), adjust control rods, pumps and valves, and view live plant parameters. The UI is minimal but demonstrates how the physics engines can drive a real-time display.

## Goals

The goals of this project are:

* **Educational simulation** – provide an intuitive way to explore how different reactor types respond to control actions. The physics models aim to capture the qualitative behaviour of real plants while remaining computationally light enough to run in real time on a MacBook.
* **Multiple reactor support** – support at least three reactor designs: PWR, BWR and RBMK. Each should have realistic features such as delayed neutrons, Doppler feedback, void feedback and xenon poisoning.
* **Modular architecture** – decouple physics from UI. Physics engines are pure Python modules with no dependency on Pygame. This facilitates testing and future replacement of the UI.
* **Cross-platform packaging** – enable end users to run the simulator as a standalone desktop application. GitHub Actions are configured to build .app bundles and .dmg installers for macOS (Intel and ARM).
* **Extensibility** – provide a base for adding more detailed models, additional reactor types, improved control room interfaces (including web-based dashboards) and support for other operating systems.

## Structure

The repository is organised as follows:

* `physics_engine.py` – point kinetics and thermal-hydraulics model for a PWR. Includes xenon poisoning and a simplified steam generator.
* `physics_engine_bwr.py` – model of a BWR with boiling coolant, void feedback, and direct steam to turbine.
* `physics_engine_rbm.py` – model of an RBMK with graphite moderator and positive void coefficient.
* `run_control_room.py` – simple Pygame application that selects a reactor, accepts user input, and visualises the state.
* `pyproject.toml` – package metadata and dependencies.
* `.github/workflows/build-macos.yml` – GitHub Action for building universal macOS binaries (.app) and packaging them into .dmg installers.

## How to Run

1. Clone the repository.
2. Install dependencies: `pip install pygame numpy`.
3. (Optional) Install `iapws` for more accurate water/steam properties: `pip install iapws`.
4. Run the control room: `python run_control_room.py`.
5. Use keys to adjust controls (see on-screen help). Press `1` for PWR, `2` for BWR and `3` for RBMK.

To build the macOS application, tag a release (e.g., `git tag v0.1.0`) and push it. The GitHub Actions workflow will create .dmg packages for x86_64 and arm64 architectures.

## Future Work

* Enhance the UI with real-time charts, alarm handling, and panel layouts.
* Add more detailed thermal-hydraulic modeling and pressure dynamics.
* Provide additional reactor models (CANDU, AP1000, etc.).
* Package for Windows and Linux using similar workflows.

---

This file provides an overview and goals for the project separate from the existing README.md.
