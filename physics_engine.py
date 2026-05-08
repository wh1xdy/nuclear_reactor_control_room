"""
physics_engine.py
------------------

This module implements a simplified but reasonably realistic
physics engine for a pressurised water reactor (PWR) that can be
driven headless or integrated into a UI.  The goals of this
implementation are:

* Represent the dominant dynamic phenomena of a PWR without
  requiring a super‑computer.  A 2018 Intel MacBook should be
  able to run the simulation at real‑time or faster.
* Separate the core neutronics, thermal hydraulics and xenon
  poison dynamics into composable components.  Each component
  maintains its own internal state and exposes a ``step()``
  method that advances the state by a small time increment.
* Provide a clear API for control inputs and observable outputs so
  that a UI layer can send rod positions, pump speeds and
  turbine valve openings and receive updated snapshot data in
  return.

The model here is deliberately lumped–parameter.  Spatial
variations are ignored, and the reactor core is treated as a
single point.  Six delayed neutron groups capture the most
important neutronics dynamics.  Two thermal nodes (fuel and
coolant) capture the dominant heat transfer within the core.
A simple secondary side model approximates a steam generator and
turbine with a single energy inventory and a turbine valve.  A
xenon/iodine model gives the characteristic hours–long negative
reactivity swing after power changes.

Wherever possible the coefficients used are taken from openly
available literature and typical light water reactor values.  The
delayed neutron fractions and decay constants are those for
thermal fission of U‑235【112925764221031†L1421-L1483】.  The prompt neutron
lifetime is chosen in the middle of the typical range of
1e‑5–1e‑4 seconds【970620566235039†L37-L44】.  Temperature feedback
coefficients correspond to a few pcm per kelvin.  These
parameters can be adjusted in the ``PlantParams`` dataclass to
tune the behaviour.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

from decay_heat import DecayHeat

try:
    # If the iapws library is available, it can be used to compute
    # saturation temperatures and enthalpies for the steam generator.
    from iapws import IAPWS97
    HAVE_IAPWS = True
except ImportError:
    IAPWS97 = None  # type: ignore
    HAVE_IAPWS = False


@dataclass
class PlantParams:
    """Container for all tunable physical constants used by the plant.

    Most values come from typical PWR references and can be adjusted
    for a different reactor type.  The units are noted in the
    comments.  Changing these values will alter the time constants
    and static gains of the simulation.
    """

    # Six group delayed neutron fractions β_i for thermal fission of U‑235
    # Taken from M. Ragheb, *Point Reactor Kinetics*【112925764221031†L1421-L1483】.
    beta: List[float] = field(default_factory=lambda: [
        2.15e-4,
        1.424e-3,
        1.274e-3,
        2.568e-3,
        7.48e-4,
        2.73e-4,
    ])
    # Decay constants λ_i [1/s] for the six delayed neutron groups
    # Values correspond to the same reference【112925764221031†L1421-L1483】.
    lambda_d: List[float] = field(default_factory=lambda: [
        0.0124,
        0.0305,
        0.111,
        0.301,
        1.14,
        3.01,
    ])
    # Prompt neutron generation time Λ [s].  Typical values for light
    # water reactors are around 10‑4 s【970620566235039†L37-L44】.
    Lambda_prompt: float = 2.0e-5  # PWR thermal spectrum: ~20 µs
    # Nominal thermal power of the reactor [W].  A typical PWR is
    # around 3000 MW thermal.  This value sets the scale for all
    # energy flows.
    nominal_power: float = 3.0e9
    # Effective heat capacity of the fuel node [J/K].  This lumps
    # together the mass and specific heat of the fuel and cladding.  A
    # large value slows down fuel temperature changes.  The default
    # value gives a time constant on the order of tens of seconds.
    fuel_heat_capacity: float = 1.0e8
    # Effective heat capacity of the primary coolant node [J/K].
    coolant_heat_capacity: float = 4.0e7
    # Effective heat capacity of the secondary/steam generator node
    # [J/K].
    sg_heat_capacity: float = 1.0e8
    # Heat transfer coefficient between fuel and coolant [W/K].  This
    # constant determines how quickly heat flows from the fuel to the
    # coolant.  Chosen such that a fuel–coolant temperature difference
    # of ~30 K at 3 GW results in a sensible heat flux.
    h_fuel_to_coolant: float = 1.0e8
    # Heat transfer coefficient between coolant and steam generator
    # [W/K].  Like the above, this controls the core outlet
    # temperature and ultimately the secondary side temperature.
    h_coolant_to_sg: float = 8.0e7
    # Heat removal coefficient of the turbine/condensor system [W/K].
    # A higher value allows the secondary side to remove more heat
    # when the turbine valve is open.
    h_sg_to_turbine: float = 5.0e7
    # Thermal efficiency of the turbine/generator.  This converts
    # steam generator heat removal to electric power.  Typical
    # efficiency of a PWR is around 33 %.
    turbine_efficiency: float = 0.33
    # Reactivity worth of control rods when fully inserted
    # [Δk/k].  A negative value means inserting rods decreases
    # reactivity.  0.05 corresponds to 5000 pcm.  Adjust as needed.
    rod_worth: float = -0.05
    # Additional negative reactivity inserted on SCRAM (emergency boration
    # and fast rod insertion combined) [Δk/k].  Makes total SCRAM worth
    # consistent with BWR/RBMK models.
    scram_extra_worth: float = -0.12
    # Time for control rods to travel from their current position to
    # fully inserted on SCRAM [s].  PWR rods fall by gravity; typical
    # full-insertion time is about 2 s.
    rod_drop_tau: float = 2.0
    # Temperature coefficients of reactivity [Δk/k per K].
    # These values are expressed in units of Δk/k and can be
    # interpreted as pcm per kelvin × 1e-5.  Negative values
    # correspond to negative feedback: increasing temperature
    # decreases reactivity.  Typical PWR temperature coefficients are
    # a few pcm/K.  The defaults here correspond to −2 pcm/K for
    # fuel and −5 pcm/K for coolant.
    fuel_temp_coeff: float = -2e-5
    coolant_temp_coeff: float = -3e-4  # PWR moderator temp coeff: ~−20 to −50 pcm/K
    # Xenon dynamics parameters.  gamma_xe is the effective
    # production rate of xenon (including the branching ratio from
    # fission to I‑135, the yield of I‑135 and the conversion to Xe‑135).
    gamma_xe: float = 0.065  # fraction of fissions producing I‑135
    lambda_I: float = 0.6931 / (6.57 * 3600)   # I-135: λ = ln2/t½, t½ = 6.57 h
    lambda_Xe: float = 0.6931 / (9.17 * 3600)  # Xe-135: λ = ln2/t½, t½ = 9.17 h
    xenon_burn_coeff: float = 0.08  # effective burn rate per unit power
    xenon_reactivity_coeff: float = 0.02  # Δk/k per unit xenon inventory
    # Nominal temperatures [K].  These reference values define
    # zero reactivity feedback when the plant operates at design
    # conditions.
    nominal_fuel_temp: float = 600.0
    nominal_coolant_temp: float = 550.0
    nominal_sg_temp: float = 550.0
    # Internal integration time step [s].  Each call to
    # ``PWRPlant.step()`` can subdivide the supplied dt into
    # multiple substeps of this length to improve numerical
    # stability.  A smaller value gives better stability for
    # stiff systems like point kinetics at the cost of CPU
    # time.  The default of 10 ms works well for time steps
    # on the order of 50 ms–100 ms.
    internal_dt: float = 0.01   # thermal/xenon substep; kinetics uses its own adaptive splitting
    # External reactivity bias injected each step (used for boron, instructor overrides, etc.)
    external_reactivity: float = 0.0


@dataclass
class ControlInputs:
    """Control setpoints supplied to the plant at each time step.

    * ``rod_position`` should be between 0.0 and 1.0, where 0 means all
      rods withdrawn (maximum positive reactivity) and 1 means all
      rods inserted.  Reactivity worth is scaled linearly by
      ``PlantParams.rod_worth``.
    * ``primary_flow`` is a fraction of nominal flow rate in the
      primary loop.  It affects the heat transfer to the steam
      generator by scaling ``h_coolant_to_sg``.
    * ``secondary_flow`` is a fraction of nominal flow in the
      secondary side.  In this simplified model it is unused but
      reserved for future expansion.
    * ``turbine_valve`` controls the heat removed from the steam
      generator: 0 means the turbine is shut and no heat is removed;
      1 means maximum heat removal.
    * ``scram`` forces the control rods fully inserted in the next
      time step and latches them there.
    """

    rod_position: float = 0.0
    primary_flow: float = 1.0
    secondary_flow: float = 1.0
    turbine_valve: float = 1.0
    scram: bool = False


@dataclass
class PlantSnapshot:
    """Data returned to the UI after each plant step.

    All fields are in SI units unless otherwise noted.  Values are
    approximate but internally consistent.  Power values are
    presented both as a relative multiple of nominal power and in
    absolute watts.  Temperatures are in kelvin.  Reactivity is
    expressed in Δk/k units.
    """

    time: float
    power_fraction: float
    thermal_power: float
    electric_power: float
    fuel_temperature: float
    coolant_temperature: float
    sg_temperature: float
    reactivity: float
    xenon_inventory: float
    rod_position: float
    scrammed: bool


class PointKinetics:
    """Six group point kinetics model.

    Maintains the neutron population ``n`` (relative to nominal
    steady state) and the delayed neutron precursor concentrations
    ``C[i]`` for each group.  The ``step`` method updates these
    quantities given a reactivity input ρ.
    """

    def __init__(self, params: PlantParams) -> None:
        self.params = params
        self.n: float = 1.0
        self.C: List[float] = [
            params.beta[i] / (params.lambda_d[i] * params.Lambda_prompt)
            for i in range(len(params.beta))
        ]
        # Pre-allocated working arrays for RK2 — avoids per-substep list creation
        _n = len(params.beta)
        self._C0:         List[float] = [0.0] * _n
        self._Cmid:       List[float] = [0.0] * _n
        self._dC1:        List[float] = [0.0] * _n
        self._dC2:        List[float] = [0.0] * _n
        self._beta_total: float       = sum(params.beta)

    def _dn_dt(self, rho: float, n: float, C: List[float]) -> float:
        p = self.params
        beta_total = sum(p.beta)
        return ((rho - beta_total) / p.Lambda_prompt) * n + sum(l * c for l, c in zip(p.lambda_d, C))

    def _dC_dt(self, n: float, C: List[float]) -> List[float]:
        p = self.params
        return [(p.beta[i] / p.Lambda_prompt) * n - p.lambda_d[i] * C[i] for i in range(len(C))]

    def step(self, dt: float, rho: float) -> None:
        """Advance the kinetics state by dt seconds using RK2 (midpoint).

        Uses pre-allocated working arrays to avoid per-substep list creation.
        """
        if dt <= 0.0:
            return
        p = self.params
        lam = p.lambda_d
        bet = p.beta
        Lam = p.Lambda_prompt
        beta_total = self._beta_total
        C = self.C
        C0, Cmid, dC1, dC2 = self._C0, self._Cmid, self._dC1, self._dC2

        # Adaptive sub-step: dt_safe = 0.3 * n / |dn/dt|
        # Unrolled 6-group sums to avoid generator+sum overhead (hot path at 600×)
        _ps = lam[0]*C[0]+lam[1]*C[1]+lam[2]*C[2]+lam[3]*C[3]+lam[4]*C[4]+lam[5]*C[5]
        _pr = (rho - beta_total if rho > beta_total else beta_total - rho) / Lam
        _n0 = self.n if self.n > 1e-6 else 1e-6
        _dm = _pr * _n0 + (_ps if _ps > 0 else -_ps)
        # Skip adaptive splitting in deep shutdown: decay heat dominates, kinetics irrelevant
        if self.n < 1e-4:
            n_sub = 1
        else:
            dt_safe = 0.3 * _n0 / (_dm if _dm > 1e-10 else 1e-10)
            n_sub   = min(200, max(1, int(dt / (dt if dt < dt_safe else dt_safe))))
        sub_dt  = dt / n_sub

        # Hoist all loop-invariants out of the hot loop
        _rbl = (rho - beta_total) / Lam
        hdt  = 0.5 * sub_dt
        iLam = 1.0 / Lam
        # Precompute beta/Lambda for each group
        bL0=bet[0]*iLam; bL1=bet[1]*iLam; bL2=bet[2]*iLam
        bL3=bet[3]*iLam; bL4=bet[4]*iLam; bL5=bet[5]*iLam
        l0=lam[0]; l1=lam[1]; l2=lam[2]; l3=lam[3]; l4=lam[4]; l5=lam[5]

        for _ in range(n_sub):
            n0 = self.n
            C0[0]=C[0]; C0[1]=C[1]; C0[2]=C[2]; C0[3]=C[3]; C0[4]=C[4]; C0[5]=C[5]
            # Stage 1
            dn1 = _rbl*n0 + l0*C0[0]+l1*C0[1]+l2*C0[2]+l3*C0[3]+l4*C0[4]+l5*C0[5]
            dC1[0]=bL0*n0-l0*C0[0]; dC1[1]=bL1*n0-l1*C0[1]; dC1[2]=bL2*n0-l2*C0[2]
            dC1[3]=bL3*n0-l3*C0[3]; dC1[4]=bL4*n0-l4*C0[4]; dC1[5]=bL5*n0-l5*C0[5]
            n_mid = n0 + hdt * dn1
            Cmid[0]=C0[0]+hdt*dC1[0]; Cmid[1]=C0[1]+hdt*dC1[1]; Cmid[2]=C0[2]+hdt*dC1[2]
            Cmid[3]=C0[3]+hdt*dC1[3]; Cmid[4]=C0[4]+hdt*dC1[4]; Cmid[5]=C0[5]+hdt*dC1[5]
            # Stage 2
            dn2 = _rbl*n_mid + l0*Cmid[0]+l1*Cmid[1]+l2*Cmid[2]+l3*Cmid[3]+l4*Cmid[4]+l5*Cmid[5]
            dC2[0]=bL0*n_mid-l0*Cmid[0]; dC2[1]=bL1*n_mid-l1*Cmid[1]; dC2[2]=bL2*n_mid-l2*Cmid[2]
            dC2[3]=bL3*n_mid-l3*Cmid[3]; dC2[4]=bL4*n_mid-l4*Cmid[4]; dC2[5]=bL5*n_mid-l5*Cmid[5]
            self.n += sub_dt * dn2
            C[0]+=sub_dt*dC2[0]; C[1]+=sub_dt*dC2[1]; C[2]+=sub_dt*dC2[2]
            C[3]+=sub_dt*dC2[3]; C[4]+=sub_dt*dC2[4]; C[5]+=sub_dt*dC2[5]

        if self.n < 0.0:
            self.n = 0.0


class XenonIodine:
    """Simple xenon/iodine poison model.

    Two first‑order ODEs approximate the production and removal
    mechanisms for I‑135 and Xe‑135.  The fission yield term is
    scaled by the current neutron population n.  Xenon reactivity
    feedback is proportional to the xenon inventory.
    """

    def __init__(self, params: PlantParams) -> None:
        self.params = params
        # Start at equilibrium values for n=1.
        # Steady state xenon inventory X_ss = gamma / (λ_Xe + k_burn)
        self.X: float = params.gamma_xe / (params.lambda_Xe + params.xenon_burn_coeff)
        # Steady state iodine inventory I_ss = γ / λ_I
        self.I: float = params.gamma_xe / params.lambda_I

    def step(self, dt: float, n: float) -> None:
        p = self.params
        # Production of iodine from fission
        dI = p.gamma_xe * n - p.lambda_I * self.I
        # Production of xenon from decay of iodine minus its decay and
        # burn‑out under neutron flux
        dX = p.lambda_I * self.I - (p.lambda_Xe + p.xenon_burn_coeff * n) * self.X
        self.I += dI * dt
        self.X += dX * dt
        # Clamp to non‑negative values
        if self.I < 0.0:
            self.I = 0.0
        if self.X < 0.0:
            self.X = 0.0

    @property
    def reactivity(self) -> float:
        """Return the xenon reactivity in Δk/k units."""
        return -self.params.xenon_reactivity_coeff * self.X


class ThermalHydraulics:
    """Lumped thermal model for the primary and secondary systems.

    Maintains three temperatures: fuel, coolant and steam generator.
    Heat flows between nodes according to fixed coefficients.  The
    steam generator also converts heat removal to electric power via
    the turbine model.  The primary and secondary flow setpoints
    scale the corresponding heat transfer coefficients.
    """

    def __init__(self, params: PlantParams) -> None:
        self.params = params
        # Initialise temperatures to nominal operating conditions
        self.T_fuel: float = params.nominal_fuel_temp
        self.T_cool: float = params.nominal_coolant_temp
        self.T_sg: float = params.nominal_sg_temp
        # Track the most recent electric power output
        self.electric_power: float = params.nominal_power * params.turbine_efficiency

    def step(self, dt: float, P_th: float, ctrl: ControlInputs) -> float:
        """Advance the thermal state and return the electric power [W].

        :param dt: integration time step [s]
        :param P_th: heat generated in the core [W]
        :param ctrl: control inputs for flow and turbine valve
        :returns: electric power produced during this substep [W]
        """
        p = self.params
        # Scale heat transfer coefficients by flow fractions.  If flow
        # is zero, no heat is removed and temperatures diverge.
        h_fc = p.h_fuel_to_coolant * max(ctrl.primary_flow, 0.0)
        h_cs = p.h_coolant_to_sg * max(ctrl.primary_flow, 0.0)
        h_st = p.h_sg_to_turbine * max(ctrl.turbine_valve, 0.0)

        # Heat transfer from fuel to coolant
        Q_fc = h_fc * (self.T_fuel - self.T_cool)
        # Heat transfer from coolant to steam generator
        Q_cs = h_cs * (self.T_cool - self.T_sg)
        # Heat removed by turbine condenser from steam generator
        Q_st = h_st * (self.T_sg - p.nominal_coolant_temp)  # condensor sink at nominal coolant temp

        # Fuel temperature change
        dT_f = (P_th - Q_fc) / p.fuel_heat_capacity
        # Coolant temperature change
        dT_c = (Q_fc - Q_cs) / p.coolant_heat_capacity
        # Steam generator temperature change
        dT_sg = (Q_cs - Q_st) / p.sg_heat_capacity

        # Euler integration
        self.T_fuel += dT_f * dt
        self.T_cool += dT_c * dt
        self.T_sg += dT_sg * dt

        # Compute electric power.  We assume the turbine removes heat
        # Q_st from the secondary, of which a fraction becomes
        # electricity.  Negative Q_st would mean the secondary is
        # heating the reactor (unlikely); clamp to zero.
        Q_removed = max(Q_st, 0.0)
        elec = Q_removed * p.turbine_efficiency
        self.electric_power = elec
        return elec


class PWRPlant:
    """Top level class that composes kinetics, thermal and xenon models.

    The ``step`` method accepts a time step and control inputs and
    returns a ``PlantSnapshot``.  Internally, the method breaks
    the time step into smaller increments to improve numerical
    stability.  The ``time`` attribute accumulates the total
    simulation time.
    """

    def __init__(self, params: Optional[PlantParams] = None) -> None:
        self.params: PlantParams = params or PlantParams()
        self.kinetics = PointKinetics(self.params)
        self.xenon = XenonIodine(self.params)
        self.thermal = ThermalHydraulics(self.params)
        self.time: float = 0.0
        self.scrammed: bool = False
        # Effective rod position that moves gradually to 1.0 on SCRAM
        # instead of teleporting, giving realistic rod-drop dynamics.
        self.rod_pos_effective: float = 0.0
        # Four-group fission product decay heat model.
        self.decay_heat: DecayHeat = DecayHeat(initial_power=1.0)

    def _compute_reactivity(self, ctrl: ControlInputs) -> float:
        """Compute total reactivity (Δk/k) from control and feedback."""
        p = self.params
        # Rod contribution: fully withdrawn rods have 0 reactivity,
        # fully inserted rods have p.rod_worth reactivity.  If the
        # plant has been scrammed previously, force rod_position to 1.
        rod_pos = self.rod_pos_effective
        # S-curve (integrated cosine) rod worth: w(x) = 3x²-2x³
        w = 3.0 * rod_pos ** 2 - 2.0 * rod_pos ** 3
        rho_rods = w * p.rod_worth
        if self.scrammed:
            rho_rods += p.scram_extra_worth
        # Temperature feedback from fuel and coolant
        rho_temp = (
            p.fuel_temp_coeff * (self.thermal.T_fuel - p.nominal_fuel_temp)
            + p.coolant_temp_coeff * (self.thermal.T_cool - p.nominal_coolant_temp)
        )
        # Xenon poisoning
        rho_xe = self.xenon.reactivity
        # Sum contributions
        rho = rho_rods + rho_temp + rho_xe + p.external_reactivity
        # Clamp to prevent numerical blow‑up.  The maximum reactivity
        # insertion should be well below prompt critical (β).
        beta_total = sum(p.beta)
        rho = max(-0.9, min(beta_total * 1.5, rho))
        return rho

    def step(self, dt: float, ctrl: ControlInputs) -> PlantSnapshot:
        """Advance the plant model by dt seconds and return a snapshot."""
        p = self.params
        # Apply SCRAM logic: once scrammed, rods are latched in.
        if ctrl.scram:
            self.scrammed = True
        # Break the external step into smaller internal steps for stability
        n_steps = max(1, int(math.ceil(dt / p.internal_dt)))
        dt_sub = dt / n_steps
        for _ in range(n_steps):
            # Rod drop dynamics: on SCRAM move rod_pos_effective toward 1.0
            # linearly at rate 1/rod_drop_tau; otherwise track operator input.
            if self.scrammed:
                self.rod_pos_effective = min(
                    1.0, self.rod_pos_effective + dt_sub / p.rod_drop_tau
                )
            else:
                self.rod_pos_effective = max(0.0, min(1.0, ctrl.rod_position))
            # Compute reactivity for this substep
            rho = self._compute_reactivity(ctrl)
            # Step point kinetics
            self.kinetics.step(dt_sub, rho)
            # Step xenon model (uses current neutron population)
            self.xenon.step(dt_sub, self.kinetics.n)
            # Advance decay heat model; add to fission power for thermal load.
            decay_frac = self.decay_heat.step(dt_sub, self.kinetics.n)
            P_th = (self.kinetics.n + decay_frac) * p.nominal_power
            # Advance thermal hydraulics and get electric power
            self.thermal.step(dt_sub, P_th, ctrl)
            # Advance simulation time
            self.time += dt_sub
        # Compose snapshot
        snap = PlantSnapshot(
            time=self.time,
            power_fraction=self.kinetics.n,
            thermal_power=(self.kinetics.n + self.decay_heat.fraction) * p.nominal_power,
            electric_power=self.thermal.electric_power,
            fuel_temperature=self.thermal.T_fuel,
            coolant_temperature=self.thermal.T_cool,
            sg_temperature=self.thermal.T_sg,
            reactivity=self._compute_reactivity(ctrl),
            xenon_inventory=self.xenon.X,
            rod_position=self.rod_pos_effective,
            scrammed=self.scrammed,
        )
        return snap


def main() -> None:
    """Run a simple headless simulation as a smoke test.

    This function allows the module to be executed directly with
    ``python physics_engine.py``.  It simulates a reactor starting
    from steady state, performs a simple transient by inserting
    control rods, and prints a few key values.  This is useful for
    verifying that the model runs without requiring a UI.
    """
    params = PlantParams()
    plant = PWRPlant(params)
    ctrl = ControlInputs()
    sim_time = 0.0
    dt = 0.1  # 100 ms outer step
    # Let the plant stabilise for a second
    for _ in range(10):
        plant.step(dt, ctrl)
        sim_time += dt
    # Insert rods to 80 % inserted over 2 seconds
    for i in range(20):
        frac = i / 19
        ctrl.rod_position = 0.2 + 0.6 * frac
        snap = plant.step(dt, ctrl)
        sim_time += dt
        print(
            f"t={sim_time:5.1f}s, n={snap.power_fraction:6.3f}, "
            f"rho={snap.reactivity:+.4f}, Pth={snap.thermal_power/1e6:8.1f} MWt, "
            f"Tfuel={snap.fuel_temperature:6.1f} K, Xe={snap.xenon_inventory:.4f}"
        )
    # Hold rods
    for i in range(50):
        snap = plant.step(dt, ctrl)
        sim_time += dt
        if i % 10 == 0:
            print(
                f"t={sim_time:5.1f}s, n={snap.power_fraction:6.3f}, "
                f"rho={snap.reactivity:+.4f}, Pth={snap.thermal_power/1e6:8.1f} MWt, "
                f"Tfuel={snap.fuel_temperature:6.1f} K, Xe={snap.xenon_inventory:.4f}"
            )


if __name__ == "__main__":
    main()
