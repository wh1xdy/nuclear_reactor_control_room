"""
physics_engine_rbm.py
---------------------

This module provides a lumped‑parameter physics engine for an
**RBMK** type nuclear power plant.  RBMK reactors are channel‑type
graphite moderated reactors cooled by light water.  The water
flows in pressure tubes through the graphite moderator and boils
partly within the core.  A notorious characteristic of RBMKs is
their **positive void coefficient**: when the coolant boils and
voids appear the moderation provided by water decreases and the
graphite alone over‑moderates the reactor, causing an increase
in reactivity.  This positive feedback contributed to the 1986
Chernobyl disaster.  Modern RBMK operation relies on careful
control rod design, negative temperature feedback in the fuel
and graphite, and strict procedural controls to mitigate this
instability.

The purpose of this implementation is educational.  It models
the essential dynamics of an RBMK while avoiding the complexity
of detailed thermal–hydraulic calculations.  The model is
composed of:

* A **six‑group point kinetics** model for the neutron
  population and delayed neutron precursors.  The same β_i and
  λ_i values as the PWR and BWR models are used【112925764221031†L1421-L1483】.
* **Three thermal nodes**: fuel, graphite moderator and
  coolant/steam mixture.  Heat flows from fuel to graphite,
  graphite to coolant and coolant to steam/condensor.
* A **void fraction estimator** based on the ratio of power to
  coolant flow.  Because the RBMK has a positive void
  coefficient, an increased void fraction raises reactivity.
* A **xenon/iodine poison model** identical to that used in
  the PWR and BWR implementations.

Each component exposes a ``step()`` method and the main plant
class orchestrates the integration.  Users control the reactor
with a control rod insertion fraction, a coolant pump speed and
a turbine valve opening.  The plant returns a snapshot of
observables after each integration call.

This module deliberately contains extensive docstrings and
comments explaining the rationale behind chosen values and the
qualitative behaviour.  It is longer than strictly necessary,
but the verbosity aids readability and understanding.  Feel free
to tune the parameter values in ``RBMKParams`` to explore
different dynamic responses.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

try:
    from iapws import IAPWS97
    HAVE_IAPWS = True
except ImportError:
    IAPWS97 = None  # type: ignore
    HAVE_IAPWS = False


@dataclass
class RBMKParams:
    """Physical constants and coefficients for the RBMK model.

    The values provided here are not tied to any specific real
    plant but instead represent reasonable magnitudes for a large
    graphite‑moderated channel reactor.  Most parameters mirror
    those found in ``BWRParams`` and ``PWRParams`` but with
    modifications to reflect RBMK characteristics.  See the
    accompanying documentation for further discussion.
    """

    beta: List[float] = field(default_factory=lambda: [
        2.15e-4,
        1.424e-3,
        1.274e-3,
        2.568e-3,
        7.48e-4,
        2.73e-4,
    ])
    lambda_d: List[float] = field(default_factory=lambda: [
        0.0124,
        0.0305,
        0.111,
        0.301,
        1.14,
        3.01,
    ])
    Lambda_prompt: float = 1.5e-4
    nominal_power: float = 3.2e9  # RBMKs are large (~1 GW electric)
    fuel_heat_capacity: float = 7.0e7
    graphite_heat_capacity: float = 1.5e8  # graphite moderator large mass
    coolant_heat_capacity: float = 5.0e7
    steam_heat_capacity: float = 1.2e8
    # Heat transfer coefficients [W/K]
    h_fuel_to_graphite: float = 6.0e7
    h_graphite_to_coolant: float = 4.5e7
    h_coolant_to_steam: float = 5.5e7
    h_steam_to_turbine: float = 4.0e7
    turbine_efficiency: float = 0.32
    rod_worth: float = -0.05  # control rods still negative
    fuel_temp_coeff: float = -1.0e-5  # negative doppler
    graphite_temp_coeff: float = -5.0e-6  # graphite also gives negative feedback
    void_coeff: float = +0.08  # positive void coefficient (≈ +8 pcm per % void)
    gamma_xe: float = 0.065
    lambda_I: float = 1.0 / (6.6 * 3600)
    lambda_Xe: float = 1.0 / (9.2 * 3600)
    xenon_burn_coeff: float = 0.08
    xenon_reactivity_coeff: float = 0.02
    nominal_fuel_temp: float = 600.0
    nominal_graphite_temp: float = 700.0
    nominal_void_fraction: float = 0.1
    nominal_steam_temp: float = 560.0
    internal_dt: float = 0.01


@dataclass
class RBMKControlInputs:
    """Control inputs for the RBMK reactor.

    Attributes
    ----------
    rod_position : float
        Fraction of control rod insertion [0–1].
    pump_speed : float
        Fraction of nominal coolant pump speed.  Higher pump speed
        reduces the void fraction by pushing more water through the
        channels.
    turbine_valve : float
        Opening of the turbine valve [0–1].
    scram : bool
        Emergency trip that forces rods fully inserted and
        latches them there.
    """
    rod_position: float = 0.0
    pump_speed: float = 1.0
    turbine_valve: float = 1.0
    scram: bool = False


@dataclass
class RBMKPlantSnapshot:
    """Snapshot of RBMK plant state returned to the UI."""
    time: float
    power_fraction: float
    power_watts: float
    fuel_temperature: float
    graphite_temperature: float
    coolant_temperature: float
    steam_temperature: float
    void_fraction: float
    reactivity: float
    xenon_inventory: float
    electric_power: float


class RBMKKinetics:
    """Point kinetics with six delayed neutron groups for the RBMK."""

    def __init__(self, params: RBMKParams):
        self.p = params
        self.n = 1.0
        self.c = [b / (self.p.Lambda_prompt * ld) for b, ld in zip(self.p.beta, self.p.lambda_d)]

    def reactivity_rhs(self, rho: float, n: float, c: List[float]) -> float:
        beta_sum = sum(self.p.beta)
        dn_dt = ((rho - beta_sum) / self.p.Lambda_prompt) * n
        dn_dt += sum(lambda_i * c_i for lambda_i, c_i in zip(self.p.lambda_d, c))
        return dn_dt

    def precursor_rhs(self, n: float, c: List[float]) -> List[float]:
        return [(beta_i / self.p.Lambda_prompt) * n - lambda_i * c_i
                for beta_i, lambda_i, c_i in zip(self.p.beta, self.p.lambda_d, c)]

    def step(self, dt: float, rho: float) -> None:
        # Second‑order Runge–Kutta
        n0 = self.n
        c0 = self.c[:]
        dn1 = self.reactivity_rhs(rho, n0, c0)
        dc1 = self.precursor_rhs(n0, c0)
        n_mid = n0 + 0.5 * dt * dn1
        c_mid = [ci + 0.5 * dt * dci for ci, dci in zip(c0, dc1)]
        dn2 = self.reactivity_rhs(rho, n_mid, c_mid)
        dc2 = self.precursor_rhs(n_mid, c_mid)
        self.n += dt * dn2
        for i in range(len(self.c)):
            self.c[i] += dt * dc2[i]


class XenonIodineRBMK:
    """Xenon/iodine dynamics for the RBMK (same as BWR/PWR)."""

    def __init__(self, params: RBMKParams):
        self.p = params
        self.I = 0.0
        self.Xe = 0.0

    def step(self, dt: float, power_fraction: float) -> float:
        dI_dt = self.p.gamma_xe * power_fraction - self.p.lambda_I * self.I
        dXe_dt = self.p.lambda_I * self.I - (self.p.lambda_Xe + self.p.xenon_burn_coeff * power_fraction) * self.Xe
        self.I += dI_dt * dt
        self.Xe += dXe_dt * dt
        return self.Xe


class RBMKThermal:
    """Thermal–hydraulic model for the RBMK reactor.

    The RBMK core is represented by four interacting nodes:
    fuel, graphite moderator, coolant mixture and steam.  Heat
    generated in the fuel flows to the graphite, from graphite to
    coolant and from coolant to steam.  A simple estimator
    computes the void fraction as a function of power and pump
    speed.  This class does not model individual channels but
    rather treats the core as a single averaged lump.
    """

    def __init__(self, params: RBMKParams):
        self.p = params
        self.T_fuel = params.nominal_fuel_temp
        self.T_graphite = params.nominal_graphite_temp
        self.T_coolant = params.nominal_steam_temp - 10.0
        self.T_steam = params.nominal_steam_temp
        self.alpha = params.nominal_void_fraction

    def saturation_temperature(self, pressure: float = 7.0e6) -> float:
        if HAVE_IAPWS:
            return IAPWS97(P=pressure / 1e6, x=0).T
        return 560.0 + 0.002 * (pressure - 7.0e6)

    def estimate_void_fraction(self, power_fraction: float, pump_speed: float) -> float:
        # Similar to BWR but with a lower nominal void and positive coefficient
        flow_ratio = max(pump_speed, 1e-3)
        alpha = self.p.nominal_void_fraction * (power_fraction / flow_ratio)
        return max(0.0, min(0.9, alpha))

    def step(self, dt: float, power_fraction: float, controls: RBMKControlInputs) -> None:
        # Heat generation
        Q_fuel = power_fraction * self.p.nominal_power
        # Heat flow from fuel to graphite
        dT_fg = self.T_fuel - self.T_graphite
        Q_fg = self.p.h_fuel_to_graphite * dT_fg
        # Heat flow from graphite to coolant
        dT_gw = self.T_graphite - self.T_coolant
        Q_gw = self.p.h_graphite_to_coolant * dT_gw
        # Heat flow from coolant to steam (boiling) depends on pump speed
        dT_ws = self.T_coolant - self.T_steam
        Q_ws = self.p.h_coolant_to_steam * controls.pump_speed * dT_ws
        # Heat removed by turbine
        dT_st = self.T_steam - (self.p.nominal_steam_temp - 50.0)
        Q_st = self.p.h_steam_to_turbine * controls.turbine_valve * dT_st
        # Update temperatures via energy balances
        self.T_fuel += dt * ((Q_fuel - Q_fg) / self.p.fuel_heat_capacity)
        self.T_graphite += dt * ((Q_fg - Q_gw) / self.p.graphite_heat_capacity)
        self.T_coolant += dt * ((Q_gw - Q_ws) / self.p.coolant_heat_capacity)
        self.T_steam += dt * ((Q_ws - Q_st) / self.p.steam_heat_capacity)
        # Update void fraction
        self.alpha = self.estimate_void_fraction(power_fraction, controls.pump_speed)


class RBMKPlant:
    """Top‑level class tying together the RBMK subsystems."""

    def __init__(self, params: Optional[RBMKParams] = None):
        self.p = params or RBMKParams()
        self.kinetics = RBMKKinetics(self.p)
        self.xenon = XenonIodineRBMK(self.p)
        self.thermal = RBMKThermal(self.p)
        self.time = 0.0
        self.scrammed = False

    def compute_reactivity(self, controls: RBMKControlInputs) -> float:
        rod_pos = 1.0 if self.scrammed or controls.scram else controls.rod_position
        if controls.scram:
            self.scrammed = True
        rho_rods = rod_pos * self.p.rod_worth
        rho_fuel = self.p.fuel_temp_coeff * (self.thermal.T_fuel - self.p.nominal_fuel_temp)
        rho_graphite = self.p.graphite_temp_coeff * (self.thermal.T_graphite - self.p.nominal_graphite_temp)
        rho_void = self.p.void_coeff * (self.thermal.alpha - self.p.nominal_void_fraction)
        rho_xe = -self.p.xenon_reactivity_coeff * self.xenon.Xe
        return rho_rods + rho_fuel + rho_graphite + rho_void + rho_xe

    def step(self, dt: float, controls: RBMKControlInputs) -> RBMKPlantSnapshot:
        n_steps = max(1, int(math.ceil(dt / self.p.internal_dt)))
        sub_dt = dt / n_steps
        for _ in range(n_steps):
            rho = self.compute_reactivity(controls)
            self.kinetics.step(sub_dt, rho)
            self.xenon.step(sub_dt, self.kinetics.n)
            self.thermal.step(sub_dt, self.kinetics.n, controls)
            self.time += sub_dt
        power_fraction = self.kinetics.n
        power_watts = power_fraction * self.p.nominal_power
        # Compute electric power from last Q_st
        dT_st = self.thermal.T_steam - (self.p.nominal_steam_temp - 50.0)
        Q_st = self.p.h_steam_to_turbine * controls.turbine_valve * dT_st
        electric_power = self.p.turbine_efficiency * Q_st
        snapshot = RBMKPlantSnapshot(
            time=self.time,
            power_fraction=power_fraction,
            power_watts=power_watts,
            fuel_temperature=self.thermal.T_fuel,
            graphite_temperature=self.thermal.T_graphite,
            coolant_temperature=self.thermal.T_coolant,
            steam_temperature=self.thermal.T_steam,
            void_fraction=self.thermal.alpha,
            reactivity=self.compute_reactivity(controls),
            xenon_inventory=self.xenon.Xe,
            electric_power=electric_power,
        )
        return snapshot

    def reset(self) -> None:
        self.__init__(self.p)


def main() -> None:
    """A simple test harness for the RBMK model.

    Running this module prints out the plant state for two
    successive changes in control rod position.  It demonstrates
    the positive void coefficient: withdrawing rods increases
    power which increases the void fraction, adding further
    reactivity and causing a faster rise than in the BWR.  The
    negative fuel and graphite temperature coefficients eventually
    stabilise the power.
    """
    plant = RBMKPlant()
    controls = RBMKControlInputs()
    # Start subcritical
    controls.rod_position = 0.9
    for i in range(300):
        snap = plant.step(0.1, controls)
        if i % 50 == 0:
            print(f"t={snap.time:.1f}s P={snap.power_fraction:.3f} alpha={snap.void_fraction:.2f}")
    # Reduce rod insertion to raise power
    controls.rod_position = 0.6
    for i in range(300):
        snap = plant.step(0.1, controls)
        if i % 50 == 0:
            print(f"t={snap.time:.1f}s P={snap.power_fraction:.3f} alpha={snap.void_fraction:.2f}")


if __name__ == "__main__":
    main()