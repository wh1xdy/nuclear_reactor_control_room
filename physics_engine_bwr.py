"""
physics_engine_bwr.py
---------------------

This module implements a simplified but reasonably detailed
physics engine for a **boiling water reactor (BWR)**.  A BWR
differs from a pressurised water reactor in that the water
circulating through the core boils directly to produce steam
which then drives the turbine.  There is no separate steam
generator; instead the primary coolant water acts as both
moderator and working fluid.  In a BWR the **void fraction**
(the fraction of the core volume occupied by steam bubbles) has
a strong influence on reactivity.  When the power rises the
coolant boils more vigorously, increasing the void fraction and
decreasing the moderator density, which in turn introduces
negative reactivity and stabilises the system.  Conversely a
decrease in power reduces voids and adds positive reactivity.

The goals of this implementation mirror those of the PWR engine
provided elsewhere in this repository:

* Represent the dominant dynamic phenomena of a BWR without
  requiring a super‑computer.  A 2018 Intel MacBook should
  easily run the simulation at real time.
* Separate neutronics, thermal hydraulics and xenon poison
  dynamics into composable components.  Each component manages
  its own state and provides a ``step()`` method to advance
  by a small time increment.
* Provide a clear API for control inputs and observable
  outputs so that a UI layer can adjust control rods, recirculation
  pump speed and turbine valves and receive updated snapshot
  data.

The model is lumped‑parameter: spatial variations are ignored
and the reactor core is represented by a single point.  Six
groups of delayed neutron precursors capture the essential
neutronics dynamics using the same β and λ values as those for
thermal fission of U‑235【112925764221031†L1421-L1483】.  The prompt
neutron generation time is taken from typical light water
reactors【970620566235039†L37-L44】.  Two thermal nodes model the fuel
and the bulk coolant mixture, and a third node captures the
steam/water inventory leaving the core.  A simple turbine model
converts heat removal into electric power.  Xenon and iodine
dynamics are included using the same approach as in the PWR
implementation.

This module is intentionally verbose.  Extensive comments and
docstrings explain the purpose of each variable and the
approximate ranges of values.  If you are unfamiliar with BWR
operation you can use these comments as a learning resource.

"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

from decay_heat import DecayHeat

try:
    # We optionally support the IAPWS97 library for saturation
    # properties.  If it is available the steam node will use
    # accurate water/steam property tables.  Otherwise a simple
    # linear relationship between temperature and enthalpy is used.
    from iapws import IAPWS97
    HAVE_IAPWS = True
except ImportError:
    IAPWS97 = None  # type: ignore
    HAVE_IAPWS = False


@dataclass
class BWRParams:
    """Container for all tunable physical constants used in the BWR model.

    Many of these values are approximate and derived from typical
    light‑water BWR references.  Adjust them to tune the behaviour
    of the simulation.  All units are SI unless otherwise noted.
    Changing these values will alter the time constants and static
    gains of the system.  See the accompanying README for
    suggestions on how to choose reasonable numbers.

    Attributes
    ----------
    beta : List[float]
        Six group delayed neutron fractions β_i.  These are the
        same as those used in the PWR model【112925764221031†L1421-L1483】.
    lambda_d : List[float]
        Decay constants λ_i [1/s] for the delayed neutron groups.
    Lambda_prompt : float
        Prompt neutron generation time Λ [s].  This determines
        how quickly the reactor responds to changes in reactivity.
    nominal_power : float
        Nominal thermal power of the reactor [W].  A typical BWR
        operates at around 2.5–3.0 GW thermal.
    fuel_heat_capacity : float
        Effective heat capacity of the fuel node [J/K].
    coolant_heat_capacity : float
        Effective heat capacity of the coolant mixture [J/K].
    steam_heat_capacity : float
        Effective heat capacity of the steam node [J/K].  This
        represents the mass and specific heat of the mixture of
        steam and entrained water leaving the core.
    h_fuel_to_coolant : float
        Heat transfer coefficient between fuel and coolant [W/K].
    h_coolant_to_steam : float
        Heat transfer coefficient from coolant to steam [W/K].
        A larger value leads to a lower core outlet temperature.
    h_steam_to_turbine : float
        Heat removal coefficient of the turbine/condensor system [W/K].
    turbine_efficiency : float
        Efficiency of the turbine/generator (steam to electric).
    rod_worth : float
        Reactivity worth of control rods when fully inserted [Δk/k].
        Negative values mean rods remove reactivity.
    fuel_temp_coeff : float
        Temperature coefficient of reactivity for the fuel [Δk/k per K].
    void_coeff : float
        Void fraction coefficient of reactivity [Δk/k per fraction].
        This should be negative for a BWR as increased voids
        decrease moderation and thus reactivity.
    gamma_xe : float
        Effective production rate of xenon.
    lambda_I : float
        Decay constant of iodine‑135 [1/s].
    lambda_Xe : float
        Decay constant of xenon‑135 [1/s].
    xenon_burn_coeff : float
        Effective burn rate of xenon per unit power [1/s].
    xenon_reactivity_coeff : float
        Reactivity worth of xenon per unit inventory [Δk/k].
    nominal_fuel_temp : float
        Nominal fuel temperature [K] at which reactivity feedback is zero.
    nominal_void_fraction : float
        Nominal void fraction at which void reactivity feedback is zero.
    nominal_steam_temp : float
        Nominal steam temperature [K] for reference enthalpy.
    internal_dt : float
        Internal integration time step [s] used for sub‑stepping.
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
    Lambda_prompt: float = 4.0e-5  # BWR thermal spectrum: ~40 µs
    nominal_power: float = 2.7e9
    fuel_heat_capacity: float = 8.0e7
    coolant_heat_capacity: float = 5.0e7
    steam_heat_capacity: float = 1.0e8
    h_fuel_to_coolant: float = 8.0e7
    h_coolant_to_steam: float = 6.0e7
    h_steam_to_turbine: float = 4.0e7
    turbine_efficiency: float = 0.33
    rod_worth: float = -0.06  # rods have slightly more worth in BWRs
    scram_extra_worth: float = -0.12  # additional negative worth on full trip
    # Time for hydraulic control blades to fully insert on SCRAM [s].
    rod_drop_tau: float = 2.5
    fuel_temp_coeff: float = -1.5e-5
    void_coeff: float = -0.15  # BWR void reactivity: ~−100 to −200 pcm/% void = −0.10 to −0.20 per unit fraction
    gamma_xe: float = 0.065
    lambda_I: float = 0.6931 / (6.57 * 3600)   # I-135: λ = ln2/t½
    lambda_Xe: float = 0.6931 / (9.17 * 3600)  # Xe-135: λ = ln2/t½
    xenon_burn_coeff: float = 0.08
    xenon_reactivity_coeff: float = 0.02
    nominal_fuel_temp: float = 800.0  # BWR average fuel temp at full power ~800 K
    nominal_void_fraction: float = 0.4
    nominal_steam_temp: float = 540.0
    internal_dt: float = 0.01   # thermal/xenon substep; kinetics RK2 is sub-stepped adaptively


@dataclass
class BWRControlInputs:
    """Control setpoints supplied to the BWR at each time step.

    Parameters
    ----------
    rod_position : float
        Control rod insertion fraction [0–1].  0 means all rods
        withdrawn (maximum positive reactivity); 1 means fully
        inserted.  The linear scaling uses ``BWRParams.rod_worth``.
    recirc_pump : float
        Fraction of nominal recirculation pump speed.  Higher
        values decrease the void fraction by increasing core flow.
    turbine_valve : float
        Fractional opening of the turbine admission valve.  0
        means the turbine is isolated; 1 means fully open.
    scram : bool
        If True, rods are latched fully inserted on the next
        integration substep and remain there.
    """

    rod_position: float = 0.0
    recirc_pump: float = 1.0
    turbine_valve: float = 1.0
    scram: bool = False


@dataclass
class BWRPlantSnapshot:
    """Data returned to the UI after each integration step.

    All values are approximate but internally consistent.  Power
    values are presented both as a fraction of nominal and as
    absolute watts.  Temperatures are in kelvin.  Reactivity
    values are expressed in Δk/k units.  The void fraction is
    given as a number between 0 and 1.  Xenon inventory is
    included to give a sense of poison dynamics.
    """

    time: float
    power_fraction: float
    power_watts: float
    fuel_temperature: float
    coolant_temperature: float
    steam_temperature: float
    void_fraction: float
    reactivity: float
    xenon_inventory: float
    electric_power: float


class PointKinetics:
    """Six group point kinetics model.

    This class encapsulates the differential equations for a
    point‑kinetics reactor model with six delayed neutron groups.
    It maintains the neutron density (represented here as power
    fraction) and the precursor concentrations.  The ``step``
    method advances the state given a reactivity input and the
    current thermal power (used for xenon burnup).  Integration
    uses a simple explicit Euler or second‑order Runge–Kutta
    depending on the size of the timestep.
    """

    def __init__(self, params: BWRParams):
        self.p = params
        self.n = 1.0
        self.c = [b / (self.p.Lambda_prompt * ld) for b, ld in zip(self.p.beta, self.p.lambda_d)]
        # Pre-allocated working arrays for RK2
        self._c0:         List[float] = [0.0] * 6
        self._cmid:       List[float] = [0.0] * 6
        self._dc1:        List[float] = [0.0] * 6
        self._dc2:        List[float] = [0.0] * 6
        self._beta_total: float       = sum(self.p.beta)

    def reactivity_rhs(self, reactivity: float, n: float, c: List[float]) -> float:
        """Right‑hand side of the neutron balance equation.

        Parameters
        ----------
        reactivity : float
            Reactivity input Δk/k.
        n : float
            Current neutron population (power fraction).
        c : List[float]
            Delayed neutron precursor concentrations.

        Returns
        -------
        dn_dt : float
            Time derivative of neutron population.
        """
        beta_sum = sum(self.p.beta)
        # Point kinetics equation: dn/dt = ((ρ − β)/Λ) n + Σλ_i C_i
        dn_dt = ((reactivity - beta_sum) / self.p.Lambda_prompt) * n
        dn_dt += sum(lambda_i * c_i for lambda_i, c_i in zip(self.p.lambda_d, c))
        return dn_dt

    def precursor_rhs(self, n: float, c: List[float]) -> List[float]:
        """Right‑hand sides for the precursor concentrations.

        Parameters
        ----------
        n : float
            Current neutron population (power fraction).
        c : List[float]
            Current precursor concentrations.

        Returns
        -------
        dc_dt : List[float]
            Time derivatives of precursor concentrations.
        """
        dc_dt = []
        for beta_i, lambda_i, c_i in zip(self.p.beta, self.p.lambda_d, c):
            # dC_i/dt = β_i/Λ * n − λ_i C_i
            dc_dt.append((beta_i / self.p.Lambda_prompt) * n - lambda_i * c_i)
        return dc_dt

    def step(self, dt: float, reactivity: float) -> None:
        """Advance the point kinetics state by dt seconds using RK2 (midpoint).

        Adaptive sub-stepping (same as PWR) prevents instability at large dt.
        Uses pre-allocated working arrays to avoid per-substep list creation.
        """
        if dt <= 0.0:
            return
        p   = self.p
        lam = p.lambda_d
        bet = p.beta
        Lam = p.Lambda_prompt
        beta_total = self._beta_total
        c = self.c
        c0, cmid, dc1, dc2 = self._c0, self._cmid, self._dc1, self._dc2

        # Adaptive sub-step + unrolled 6-group sums (hot path at 600×)
        _ps = lam[0]*c[0]+lam[1]*c[1]+lam[2]*c[2]+lam[3]*c[3]+lam[4]*c[4]+lam[5]*c[5]
        _pr = (reactivity - beta_total if reactivity > beta_total else beta_total - reactivity) / Lam
        _n0 = self.n if self.n > 1e-6 else 1e-6
        _dm = _pr * _n0 + (_ps if _ps > 0 else -_ps)
        # Skip adaptive splitting in deep shutdown: decay heat dominates, kinetics irrelevant
        if self.n < 1e-4:
            n_sub = 1
        else:
            dt_safe = 0.3 * _n0 / (_dm if _dm > 1e-10 else 1e-10)
            n_sub   = min(200, max(1, int(dt / (dt if dt < dt_safe else dt_safe))))
        sub_dt  = dt / n_sub
        _rbl = (reactivity - beta_total) / Lam
        hdt  = 0.5 * sub_dt
        iLam = 1.0 / Lam
        bL0=bet[0]*iLam; bL1=bet[1]*iLam; bL2=bet[2]*iLam
        bL3=bet[3]*iLam; bL4=bet[4]*iLam; bL5=bet[5]*iLam
        l0=lam[0]; l1=lam[1]; l2=lam[2]; l3=lam[3]; l4=lam[4]; l5=lam[5]

        for _ in range(n_sub):
            n0 = self.n
            c0[0]=c[0]; c0[1]=c[1]; c0[2]=c[2]; c0[3]=c[3]; c0[4]=c[4]; c0[5]=c[5]
            dn1 = _rbl*n0 + l0*c0[0]+l1*c0[1]+l2*c0[2]+l3*c0[3]+l4*c0[4]+l5*c0[5]
            dc1[0]=bL0*n0-l0*c0[0]; dc1[1]=bL1*n0-l1*c0[1]; dc1[2]=bL2*n0-l2*c0[2]
            dc1[3]=bL3*n0-l3*c0[3]; dc1[4]=bL4*n0-l4*c0[4]; dc1[5]=bL5*n0-l5*c0[5]
            n_mid = n0 + hdt * dn1
            cmid[0]=c0[0]+hdt*dc1[0]; cmid[1]=c0[1]+hdt*dc1[1]; cmid[2]=c0[2]+hdt*dc1[2]
            cmid[3]=c0[3]+hdt*dc1[3]; cmid[4]=c0[4]+hdt*dc1[4]; cmid[5]=c0[5]+hdt*dc1[5]
            dn2 = _rbl*n_mid + l0*cmid[0]+l1*cmid[1]+l2*cmid[2]+l3*cmid[3]+l4*cmid[4]+l5*cmid[5]
            dc2[0]=bL0*n_mid-l0*cmid[0]; dc2[1]=bL1*n_mid-l1*cmid[1]; dc2[2]=bL2*n_mid-l2*cmid[2]
            dc2[3]=bL3*n_mid-l3*cmid[3]; dc2[4]=bL4*n_mid-l4*cmid[4]; dc2[5]=bL5*n_mid-l5*cmid[5]
            self.n += sub_dt * dn2
            c[0]+=sub_dt*dc2[0]; c[1]+=sub_dt*dc2[1]; c[2]+=sub_dt*dc2[2]
            c[3]+=sub_dt*dc2[3]; c[4]+=sub_dt*dc2[4]; c[5]+=sub_dt*dc2[5]

        if not math.isfinite(self.n) or self.n < 0.0:
            self.n = 0.0
        elif self.n > 5.0:
            self.n = 5.0


class XenonIodine:
    """Models the production and decay of iodine‑135 and xenon‑135.

    Xenon poisons the reactor because Xe‑135 has a very high
    neutron absorption cross‑section.  The production and decay
    rates are taken from typical nuclear data.  Xenon burns away
    when absorbing neutrons, so the burn coefficient multiplies
    the reactor power fraction.  This model returns the xenon
    inventory which is then converted into a reactivity term in
    the plant.
    """

    def __init__(self, params: BWRParams):
        self.p = params
        self.I = 0.0
        self.Xe = 0.0

    def step(self, dt: float, power_fraction: float) -> float:
        """Advance the iodine/xenon inventory by dt seconds.

        Parameters
        ----------
        dt : float
            Time increment [s].
        power_fraction : float
            Current reactor power as a fraction of nominal.

        Returns
        -------
        Xe : float
            Updated xenon inventory (dimensionless).
        """
        # Iodine production proportional to power
        dI_dt = self.p.gamma_xe * power_fraction - (self.p.lambda_I * self.I)
        dXe_dt = self.p.lambda_I * self.I - (self.p.lambda_Xe + self.p.xenon_burn_coeff * power_fraction) * self.Xe
        # Simple Euler integration; time constants are long (hours)
        self.I += dI_dt * dt
        self.Xe += dXe_dt * dt
        return self.Xe


class ThermalHydraulics:
    """Lumped thermal–hydraulic model for a BWR.

    This model keeps track of three variables: the fuel temperature,
    the bulk coolant mixture temperature and the steam temperature
    exiting the core.  It also estimates the void fraction as a
    function of coolant enthalpy and recirculation flow.  The
    ``step`` method updates these variables given a power input,
    control inputs and the current void fraction.  Although this
    representation is highly simplified compared with detailed
    two‑phase flow calculations, it captures the dominant
    qualitative behaviour: increasing power heats the coolant,
    raising the void fraction and thus decreasing reactivity.
    """

    def __init__(self, params: BWRParams):
        self.p = params
        # State variables
        self.T_fuel = params.nominal_fuel_temp
        self.T_coolant = params.nominal_steam_temp - 10.0  # slightly subcooled
        self.T_steam = params.nominal_steam_temp
        self.alpha = params.nominal_void_fraction

    def saturation_temperature(self, pressure: float = 7.0e6) -> float:
        """Return the saturation temperature for the given pressure [Pa].

        We assume a typical BWR operating pressure of around 7 MPa.
        If the ``iapws`` library is available this uses accurate
        steam tables; otherwise a linear approximation is used.
        """
        if HAVE_IAPWS:
            return IAPWS97(P=pressure / 1e6, x=0).T
        # Linear approximation: 7 MPa → 560 K; 6 MPa → 550 K
        return 540.0 + 0.002 * (pressure - 7.0e6)

    def mixture_density(self, T: float, alpha: float) -> float:
        """Approximate density [kg/m^3] of the two‑phase mixture.

        We linearly interpolate between saturated water density and
        saturated steam density at the given temperature.  At
        nominal conditions water at 560 K and 7 MPa has a density
        around 720 kg/m³ and saturated steam around 30 kg/m³.
        """
        rho_water = 720.0 - 0.3 * (T - 560.0)
        rho_steam = 30.0 + 0.05 * (T - 560.0)
        return (1 - alpha) * rho_water + alpha * rho_steam

    def estimate_void_fraction(self, power_fraction: float, recirc_pump: float) -> float:
        """Estimate the void fraction as a function of power and flow.

        A simple correlation is used: the void fraction scales with
        the ratio of power to recirculation flow.  At nominal power
        and flow the void fraction is ``p.nominal_void_fraction``.
        Doubling the power or halving the flow roughly doubles the
        excess void fraction.  The result is clamped between 0 and 1.
        """
        # Avoid division by zero; a tiny flow still produces large voids
        flow_ratio = max(recirc_pump, 1e-3)
        alpha = self.p.nominal_void_fraction * (power_fraction / flow_ratio)
        # Smooth clamp to [0, 0.95]
        return max(0.0, min(0.95, alpha))

    def step(self, dt: float, power_fraction: float, controls: BWRControlInputs) -> None:
        """Advance the thermal–hydraulic state by dt seconds.

        Parameters
        ----------
        dt : float
            Time step [s].
        power_fraction : float
            Current reactor power as a fraction of nominal.
        controls : BWRControlInputs
            Control inputs for this step.
        """
        # Compute heat generated in fuel [W]
        Q_fuel = power_fraction * self.p.nominal_power
        # Heat transferred from fuel to coolant [W]
        dT_fc = self.T_fuel - self.T_coolant
        Q_fc = self.p.h_fuel_to_coolant * dT_fc
        # Heat transferred from coolant to steam/water mixture [W]
        dT_cs = self.T_coolant - self.T_steam
        # Recirculation pump multiplies the transfer coefficient
        Q_cs = self.p.h_coolant_to_steam * controls.recirc_pump * dT_cs
        # Heat removed by turbine [W]
        # Opening the valve increases removal; a closed valve leads to heat
        # accumulation in the steam node causing temperature rise.
        dT_st = self.T_steam - (self.p.nominal_steam_temp - 40.0)  # condenser temp at 500 K
        Q_st = self.p.h_steam_to_turbine * controls.turbine_valve * dT_st
        # Update temperatures using energy balances
        self.T_fuel += dt * ((Q_fuel - Q_fc) / self.p.fuel_heat_capacity)
        self.T_coolant += dt * ((Q_fc - Q_cs) / self.p.coolant_heat_capacity)
        self.T_steam += dt * ((Q_cs - Q_st) / self.p.steam_heat_capacity)
        # Void fraction is managed externally by PlantSupervisor's mechanistic
        # void model.  When running standalone, alpha stays at its initial value
        # (nominal_void_fraction) or at whatever value was last set externally.
        # First-order nudge toward static estimate so standalone runs are stable.
        alpha_ss = self.estimate_void_fraction(power_fraction, controls.recirc_pump)
        self.alpha += dt * (alpha_ss - self.alpha) / 30.0
        self.alpha = max(0.0, min(0.95, self.alpha))


class BWRPlant:
    """Main class representing the boiling water reactor plant.

    This class composes the point kinetics, xenon/iodine and
    thermal–hydraulic models into a single plant.  The ``step``
    method performs sub‑stepped integration: the caller passes
    a larger dt (for example 0.1 s) and the plant subdivides it
    into multiples of ``p.internal_dt``.  At each substep the
    reactivity is calculated from control rod insertion, fuel
    temperature and void fraction, plus xenon poisoning.  The
    kinetic equations are advanced, xenon is updated and the
    thermal state is recalculated.  A ``BWRPlantSnapshot`` is
    returned to the caller containing all observable values.
    """

    def __init__(self, params: Optional[BWRParams] = None):
        self.p = params or BWRParams()
        self.kinetics = PointKinetics(self.p)
        self.xenon = XenonIodine(self.p)
        self.thermal = ThermalHydraulics(self.p)
        self.time = 0.0
        self.scrammed = False
        self.rod_pos_effective: float = 0.0
        self.decay_heat: DecayHeat = DecayHeat(initial_power=1.0)
        # Density wave: first-order lag (τ≈2s) on power→void coupling
        # Creates delayed void feedback → mild oscillatory tendency at low flow
        self._void_drive: float = self.p.nominal_void_fraction

    def compute_reactivity(self, controls: BWRControlInputs) -> float:
        """Compute the total reactivity [Δk/k] from controls and state.

        Rod insertion contributes a negative term scaled by
        ``rod_worth``.  Fuel temperature feedback and void
        feedback are relative to nominal values.  Xenon poison
        increases the negative reactivity proportional to xenon
        inventory.  The sum of these terms is returned.
        """
        rod_pos = self.rod_pos_effective
        # S-curve (integrated cosine) rod worth: w(x) = 3x²-2x³
        w = 3.0 * rod_pos ** 2 - 2.0 * rod_pos ** 3
        rho_rods = w * self.p.rod_worth
        if self.scrammed:
            rho_rods += self.p.scram_extra_worth
        rho_fuel = self.p.fuel_temp_coeff * (self.thermal.T_fuel - self.p.nominal_fuel_temp)
        rho_void = self.p.void_coeff * (self.thermal.alpha - self.p.nominal_void_fraction)
        rho_xe = -self.p.xenon_reactivity_coeff * self.xenon.Xe
        rho = rho_rods + rho_fuel + rho_void + rho_xe
        beta_total = sum(self.p.beta)
        rho = max(-0.9, min(beta_total * 1.5, rho))
        return rho

    def step(self, dt: float, controls: BWRControlInputs) -> BWRPlantSnapshot:
        """Advance the plant by dt seconds and return a snapshot.

        The integration subdivides dt into smaller increments to
        maintain stability of the stiff kinetics equations.  The
        number of internal steps is the ceiling of dt/internal_dt.
        """
        if controls.scram:
            self.scrammed = True
        n_steps = max(1, int(math.ceil(dt / self.p.internal_dt)))
        sub_dt = dt / n_steps
        for _ in range(n_steps):
            # Rod drop dynamics
            if self.scrammed:
                self.rod_pos_effective = min(
                    1.0, self.rod_pos_effective + sub_dt / self.p.rod_drop_tau
                )
            else:
                self.rod_pos_effective = max(0.0, min(1.0, controls.rod_position))
            # Density wave: update void_drive with 2s lag on kinetics.n→void
            tau_dw = 2.0
            void_target = self.p.nominal_void_fraction * self.kinetics.n / max(0.001, controls.recirc_pump)
            void_target = max(0.0, min(0.95, void_target))
            self._void_drive += sub_dt * (void_target - self._void_drive) / tau_dw
            # Nudge thermal alpha toward void_drive (slowly, so fast supervisor model can override)
            self.thermal.alpha += sub_dt * (self._void_drive - self.thermal.alpha) / 5.0
            self.thermal.alpha = max(0.0, min(0.95, self.thermal.alpha))
            # Compute reactivity from current state and control inputs
            rho = self.compute_reactivity(controls)
            # Advance kinetics
            self.kinetics.step(sub_dt, rho)
            # Update xenon inventory
            xe = self.xenon.step(sub_dt, self.kinetics.n)
            # Advance decay heat and compute total effective power
            decay_frac = self.decay_heat.step(sub_dt, self.kinetics.n)
            effective_power = self.kinetics.n + decay_frac
            # Advance thermal hydraulics
            self.thermal.step(sub_dt, effective_power, controls)
            # Advance time
            self.time += sub_dt
        # Compute power and electric output
        power_fraction = self.kinetics.n
        power_watts = (self.kinetics.n + self.decay_heat.fraction) * self.p.nominal_power
        # Electric power from turbine heat removal
        # Approximate heat removal as Q_st from last substep
        # dT_st was computed in thermal model using previous state; we repeat for reporting
        dT_st = self.thermal.T_steam - (self.p.nominal_steam_temp - 40.0)
        Q_st = self.p.h_steam_to_turbine * controls.turbine_valve * dT_st
        electric_power = self.p.turbine_efficiency * Q_st
        snapshot = BWRPlantSnapshot(
            time=self.time,
            power_fraction=power_fraction,
            power_watts=power_watts,
            fuel_temperature=self.thermal.T_fuel,
            coolant_temperature=self.thermal.T_coolant,
            steam_temperature=self.thermal.T_steam,
            void_fraction=self.thermal.alpha,
            reactivity=self.compute_reactivity(controls),
            xenon_inventory=self.xenon.Xe,
            electric_power=electric_power,
        )
        return snapshot

    def reset(self) -> None:
        """Reset the plant state to nominal conditions."""
        self.__init__(self.p)


def main() -> None:
    """Run a simple headless test of the BWR physics engine.

    This function demonstrates how to instantiate the plant and
    simulate a few tens of seconds without a UI.  It performs a
    step change in control rod position and prints the power
    response.  You can run this module as a script to verify
    correctness and stability.
    """
    plant = BWRPlant()
    controls = BWRControlInputs()
    # Start with rods at 0.8 insertion for subcriticality
    controls.rod_position = 0.8
    # Run for 30 seconds with this insertion
    for i in range(300):
        snapshot = plant.step(0.1, controls)
        if i % 50 == 0:
            print(f"t={snapshot.time:.1f}s P={snapshot.power_fraction:.3f} fuel={snapshot.fuel_temperature:.1f}K void={snapshot.void_fraction:.2f}")
    # Withdraw rods to 0.5 and see the power rise
    controls.rod_position = 0.5
    for i in range(300):
        snapshot = plant.step(0.1, controls)
        if i % 50 == 0:
            print(f"t={snapshot.time:.1f}s P={snapshot.power_fraction:.3f} fuel={snapshot.fuel_temperature:.1f}K void={snapshot.void_fraction:.2f}")


if __name__ == "__main__":
    main()