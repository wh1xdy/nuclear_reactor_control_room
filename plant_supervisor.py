from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Literal

from physics_engine import ControlInputs as PWRControlInputs, PWRPlant
from physics_engine_bwr import BWRControlInputs, BWRPlant
from physics_engine_rbm import RBMKControlInputs, RBMKPlant


ReactorType = Literal["PWR", "BWR", "RBMK"]


@dataclass
class SupervisorControls:
    rod_position: float = 0.0
    flow: float = 1.0
    turbine_valve: float = 1.0
    feedwater_valve: float = 0.7
    pressurizer_heater: float = 0.5
    startup_permit: bool = False
    turbine_trip: bool = False
    scram: bool = False
    fault_pump_degraded: bool = False
    fault_feedwater_loss: bool = False


@dataclass
class BalanceOfPlantState:
    pressure_mpa: float = 7.0
    steam_inventory: float = 1.0
    condenser_temp_k: float = 305.0
    feedwater_inventory: float = 1.0
    alarms: List[str] = field(default_factory=list)
    trips: List[str] = field(default_factory=list)
    unacked_alarms: List[str] = field(default_factory=list)
    acked_alarms: List[str] = field(default_factory=list)
    trip_latched: bool = False


@dataclass
class UnifiedSnapshot:
    reactor_type: ReactorType
    time: float
    power_fraction: float
    thermal_mw: float
    electric_mw: float
    reactivity: float
    fuel_temp_k: float
    coolant_temp_k: float
    steam_temp_k: float
    void_fraction: float
    pressure_mpa: float
    steam_inventory: float
    condenser_temp_k: float
    feedwater_inventory: float
    alarms: List[str]
    trips: List[str]
    unacked_alarms: List[str]
    trip_latched: bool


class PlantSupervisor:
    """Adds protection system, balance-of-plant and fault handling around core models."""

    def __init__(self, reactor_type: ReactorType = "PWR") -> None:
        self.reactor_type: ReactorType = reactor_type
        self.controls = SupervisorControls()
        self.bop = BalanceOfPlantState()
        self._init_plant()

    def _init_plant(self) -> None:
        if self.reactor_type == "PWR":
            self.plant = PWRPlant()
            self.ctrl = PWRControlInputs()
            self.bop.pressure_mpa = 15.5
        elif self.reactor_type == "BWR":
            self.plant = BWRPlant()
            self.ctrl = BWRControlInputs()
            self.bop.pressure_mpa = 7.0
        else:
            self.plant = RBMKPlant()
            self.ctrl = RBMKControlInputs()
            self.bop.pressure_mpa = 7.0

    def select(self, reactor_type: ReactorType) -> None:
        self.reactor_type = reactor_type
        self.controls = SupervisorControls()
        self.bop = BalanceOfPlantState()
        self._init_plant()

    def acknowledge_alarms(self) -> None:
        self.bop.acked_alarms = list(self.bop.alarms)
        self.bop.unacked_alarms = []

    def reset_trip_latch(self) -> bool:
        safe_pressure = self.bop.pressure_mpa < (15.9 if self.reactor_type == "PWR" else 7.5)
        safe_temp = True
        if hasattr(self.plant, "thermal") and hasattr(self.plant.thermal, "T_fuel"):
            safe_temp = self.plant.thermal.T_fuel < 900.0
        if safe_pressure and safe_temp and not self.controls.turbine_trip:
            self.bop.trip_latched = False
            self.controls.scram = False
            return True
        return False

    def _mechanistic_void_step(self, dt: float) -> None:
        if self.reactor_type not in ("BWR", "RBMK"):
            return
        thermal = self.plant.thermal
        flow = max(0.05, self.controls.flow)
        pressure_factor = max(0.4, min(1.4, 7.0 / max(5.0, self.bop.pressure_mpa)))
        coolant_superheat = max(0.0, thermal.T_coolant - (520.0 + 3.0 * (self.bop.pressure_mpa - 7.0)))
        boil_source = 0.015 * self.plant.kinetics.n * pressure_factor + 0.0008 * coolant_superheat
        collapse_sink = 0.04 * flow * (1.0 - 0.2 * pressure_factor)
        dalpha = boil_source * (1.0 - thermal.alpha) - collapse_sink * thermal.alpha
        thermal.alpha = max(0.0, min(0.95 if self.reactor_type == "BWR" else 0.9, thermal.alpha + dalpha * dt))

    def _update_bop(self, dt: float, thermal_power_w: float) -> None:
        c = self.controls
        if c.fault_feedwater_loss:
            c.feedwater_valve = min(c.feedwater_valve, 0.1)

        steam_gen = thermal_power_w / 3.0e9
        turbine_draw = c.turbine_valve * (0.0 if c.turbine_trip else 1.0)
        feed = c.feedwater_valve
        self.bop.steam_inventory += dt * (0.6 * steam_gen - 0.9 * turbine_draw)
        self.bop.steam_inventory += dt * 0.2 * feed
        self.bop.steam_inventory = max(0.0, min(2.0, self.bop.steam_inventory))

        pressure_nom = 15.5 if self.reactor_type == "PWR" else 7.0
        heater = c.pressurizer_heater if self.reactor_type == "PWR" else 0.2 * c.pressurizer_heater
        dpress = 0.12 * (self.bop.steam_inventory - 1.0) + 0.09 * heater - 0.08 * turbine_draw
        self.bop.pressure_mpa += dt * dpress
        self.bop.pressure_mpa += dt * 0.1 * (pressure_nom - self.bop.pressure_mpa)

        self.bop.condenser_temp_k += dt * (6.0 * turbine_draw - 0.06 * (self.bop.condenser_temp_k - 305.0))
        self.bop.feedwater_inventory += dt * (0.12 - 0.14 * feed)
        self.bop.feedwater_inventory = max(0.0, min(1.2, self.bop.feedwater_inventory))

    def _protection(self, snap: Dict[str, float]) -> None:
        alarms: List[str] = []
        trips: List[str] = []
        votes = 0

        if snap["fuel_temp_k"] > 1200:
            alarms.append("Fuel temperature high")
            votes += 1
        if self.bop.pressure_mpa > (16.2 if self.reactor_type == "PWR" else 7.8):
            alarms.append("Reactor pressure high")
            votes += 1
        if self.bop.feedwater_inventory < 0.2:
            alarms.append("Feedwater inventory low")
        if self.bop.condenser_temp_k > 330:
            alarms.append("Condenser temperature high")
        if self.controls.fault_pump_degraded:
            alarms.append("Primary pump degraded")

        if snap["fuel_temp_k"] > 1400:
            trips.append("AUTO SCRAM: high fuel temperature")
        if self.bop.pressure_mpa > (16.8 if self.reactor_type == "PWR" else 8.3):
            trips.append("AUTO SCRAM: high reactor pressure")
        if self.controls.turbine_trip:
            trips.append("Turbine trip active")

        if votes >= 2 and "AUTO SCRAM: 2oo2 process vote" not in trips:
            trips.append("AUTO SCRAM: 2oo2 process vote")

        if trips:
            self.controls.scram = True
            self.bop.trip_latched = True

        # if latch exists, keep SCRAM asserted until explicit reset
        if self.bop.trip_latched:
            self.controls.scram = True
            if not trips:
                trips.append("SCRAM latched")

        self.bop.alarms = alarms
        self.bop.trips = trips

        active_set = set(alarms)
        # Keep acknowledgements only for currently active alarms.
        self.bop.acked_alarms = [a for a in self.bop.acked_alarms if a in active_set]
        self.bop.unacked_alarms = [a for a in self.bop.unacked_alarms if a in active_set and a not in self.bop.acked_alarms]
        for a in alarms:
            if a not in self.bop.acked_alarms and a not in self.bop.unacked_alarms:
                self.bop.unacked_alarms.append(a)

    def step(self, dt: float) -> UnifiedSnapshot:
        c = self.controls
        flow = c.flow * (0.6 if c.fault_pump_degraded else 1.0)
        if not c.startup_permit:
            c.rod_position = max(c.rod_position, 0.85)

        if self.reactor_type == "PWR":
            self.ctrl.rod_position = c.rod_position
            self.ctrl.primary_flow = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            base = self.plant.step(dt, self.ctrl)
            thermal_power_w = base.thermal_power
            snap = {
                "time": base.time,
                "power_fraction": base.power_fraction,
                "thermal_mw": base.thermal_power / 1e6,
                "electric_mw": base.electric_power / 1e6,
                "reactivity": base.reactivity,
                "fuel_temp_k": base.fuel_temperature,
                "coolant_temp_k": base.coolant_temperature,
                "steam_temp_k": base.sg_temperature,
                "void_fraction": 0.0,
            }
        elif self.reactor_type == "BWR":
            self._mechanistic_void_step(dt)
            self.ctrl.rod_position = c.rod_position
            self.ctrl.recirc_pump = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            base = self.plant.step(dt, self.ctrl)
            thermal_power_w = base.power_watts
            snap = {
                "time": base.time,
                "power_fraction": base.power_fraction,
                "thermal_mw": base.power_watts / 1e6,
                "electric_mw": base.electric_power / 1e6,
                "reactivity": base.reactivity,
                "fuel_temp_k": base.fuel_temperature,
                "coolant_temp_k": base.coolant_temperature,
                "steam_temp_k": base.steam_temperature,
                "void_fraction": base.void_fraction,
            }
        else:
            self._mechanistic_void_step(dt)
            self.ctrl.rod_position = c.rod_position
            self.ctrl.pump_speed = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            base = self.plant.step(dt, self.ctrl)
            thermal_power_w = base.power_watts
            snap = {
                "time": base.time,
                "power_fraction": base.power_fraction,
                "thermal_mw": base.power_watts / 1e6,
                "electric_mw": base.electric_power / 1e6,
                "reactivity": base.reactivity,
                "fuel_temp_k": base.fuel_temperature,
                "coolant_temp_k": base.coolant_temperature,
                "steam_temp_k": base.steam_temperature,
                "void_fraction": base.void_fraction,
            }

        self._update_bop(dt, thermal_power_w)
        self._protection(snap)

        return UnifiedSnapshot(
            reactor_type=self.reactor_type,
            time=snap["time"],
            power_fraction=snap["power_fraction"],
            thermal_mw=snap["thermal_mw"],
            electric_mw=snap["electric_mw"],
            reactivity=snap["reactivity"],
            fuel_temp_k=snap["fuel_temp_k"],
            coolant_temp_k=snap["coolant_temp_k"],
            steam_temp_k=snap["steam_temp_k"],
            void_fraction=snap["void_fraction"],
            pressure_mpa=self.bop.pressure_mpa,
            steam_inventory=self.bop.steam_inventory,
            condenser_temp_k=self.bop.condenser_temp_k,
            feedwater_inventory=self.bop.feedwater_inventory,
            alarms=list(self.bop.alarms),
            trips=list(self.bop.trips),
            unacked_alarms=list(self.bop.unacked_alarms),
            trip_latched=self.bop.trip_latched,
        )
