from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Deque, Dict, List, Literal

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
    bypass_high_pressure_trip: bool = False
    bypass_high_fuel_trip: bool = False
    inhibit_auto_scram: bool = False


@dataclass
class SensorChannel:
    value: float
    bias: float = 0.0
    tau_s: float = 0.35
    failed_high: bool = False
    failed_low: bool = False
    frozen: bool = False

    def update(self, dt: float, true_value: float) -> float:
        if self.failed_high:
            self.value = true_value * 1.20 + 50.0
            return self.value
        if self.failed_low:
            self.value = true_value * 0.80 - 50.0
            return self.value
        if self.frozen:
            return self.value
        alpha = min(1.0, dt / max(self.tau_s, 1e-3))
        target = true_value + self.bias
        self.value += alpha * (target - self.value)
        return self.value


@dataclass
class Instrumentation:
    fuel_temp: List[SensorChannel] = field(default_factory=list)
    pressure: List[SensorChannel] = field(default_factory=list)


@dataclass
class AxialChannelState:
    quality: List[float]
    node_void: List[float]
    node_pressure_mpa: List[float]
    node_power: List[float]
    avg_void: float
    pressure_drop_mpa: float


@dataclass
class BalanceOfPlantState:
    pressure_mpa: float = 7.0
    steam_inventory: float = 1.0
    condenser_temp_k: float = 305.0
    feedwater_inventory: float = 1.0
    alarms: List[str] = field(default_factory=list)
    trips: List[str] = field(default_factory=list)
    diagnostics: List[str] = field(default_factory=list)
    unacked_alarms: List[str] = field(default_factory=list)
    acked_alarms: List[str] = field(default_factory=list)
    shelved_until: Dict[str, float] = field(default_factory=dict)
    trip_latched: bool = False
    trip_timers: Dict[str, float] = field(default_factory=lambda: {"fuel": 0.0, "pressure": 0.0})
    event_log: List[str] = field(default_factory=list)


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
    pressure_drop_mpa: float
    steam_inventory: float
    condenser_temp_k: float
    feedwater_inventory: float
    alarms: List[str]
    trips: List[str]
    diagnostics: List[str]
    unacked_alarms: List[str]
    trip_latched: bool
    checklist: Dict[str, bool]
    startup_checklist: Dict[str, bool]
    event_log: List[str]


class PlantSupervisor:
    def __init__(self, reactor_type: ReactorType = "PWR") -> None:
        self.reactor_type: ReactorType = reactor_type
        self.controls = SupervisorControls()
        self.bop = BalanceOfPlantState()
        self.channel_state = AxialChannelState([0.0] * 6, [0.0] * 6, [7.0] * 6, [1.0 / 6.0] * 6, 0.0, 0.0)
        self.instrumentation = Instrumentation()
        self.trend_power: Deque[float] = deque(maxlen=240)
        self.trend_pressure: Deque[float] = deque(maxlen=240)
        self.trend_fuel: Deque[float] = deque(maxlen=240)
        self._init_plant()

    def _init_plant(self) -> None:
        if self.reactor_type == "PWR":
            self.plant = PWRPlant()
            self.ctrl = PWRControlInputs()
            self.bop.pressure_mpa = 15.5
            self.channel_state = AxialChannelState([0.0] * 6, [0.0] * 6, [15.5] * 6, [1.0 / 6.0] * 6, 0.0, 0.0)
        elif self.reactor_type == "BWR":
            self.plant = BWRPlant()
            self.ctrl = BWRControlInputs()
            self.bop.pressure_mpa = 7.0
            self.channel_state = AxialChannelState([0.03] * 6, [0.03] * 6, [7.0] * 6, [1.0 / 6.0] * 6, 0.03, 0.0)
        else:
            self.plant = RBMKPlant()
            self.ctrl = RBMKControlInputs()
            self.bop.pressure_mpa = 7.0
            self.channel_state = AxialChannelState([0.02] * 6, [0.02] * 6, [7.0] * 6, [1.0 / 6.0] * 6, 0.02, 0.0)

        self.instrumentation = Instrumentation(
            fuel_temp=[
                SensorChannel(self.plant.thermal.T_fuel, tau_s=0.4),
                SensorChannel(self.plant.thermal.T_fuel, tau_s=0.5),
                SensorChannel(self.plant.thermal.T_fuel, tau_s=0.45),
            ],
            pressure=[
                SensorChannel(self.bop.pressure_mpa, tau_s=0.35),
                SensorChannel(self.bop.pressure_mpa, tau_s=0.4),
                SensorChannel(self.bop.pressure_mpa, tau_s=0.3),
            ],
        )

    def select(self, reactor_type: ReactorType) -> None:
        self.reactor_type = reactor_type
        self.controls = SupervisorControls()
        self.bop = BalanceOfPlantState()
        self.trend_power.clear()
        self.trend_pressure.clear()
        self.trend_fuel.clear()
        self._init_plant()

    def acknowledge_alarms(self) -> None:
        self.bop.acked_alarms = list(self.bop.alarms)
        self.bop.unacked_alarms = []

    def shelve_current_alarms(self, now: float, seconds: float = 20.0) -> None:
        for a in self.bop.alarms:
            self.bop.shelved_until[a] = now + seconds

    def reset_trip_latch(self) -> bool:
        checklist = self._build_reset_checklist()
        if all(checklist.values()):
            self.bop.trip_latched = False
            self.controls.scram = False
            return True
        return False

    def _build_reset_checklist(self) -> Dict[str, bool]:
        return {
            "pressure_ok": self.bop.pressure_mpa < (15.9 if self.reactor_type == "PWR" else 7.5),
            "fuel_temp_ok": self.plant.thermal.T_fuel < 900.0,
            "turbine_trip_clear": not self.controls.turbine_trip,
            "rods_inserted": self.controls.rod_position >= 0.95,
            "flow_available": self.controls.flow >= 0.8,
            "feedwater_ok": self.bop.feedwater_inventory >= 0.4,
            "no_inhibit": not self.controls.inhibit_auto_scram,
        }

    def _build_startup_checklist(self) -> Dict[str, bool]:
        return {
            "startup_permit": self.controls.startup_permit,
            "trip_not_latched": not self.bop.trip_latched,
            "feedwater_available": self.bop.feedwater_inventory > 0.4,
            "condenser_ready": self.bop.condenser_temp_k < 325.0,
            "pressure_within_band": abs(self.bop.pressure_mpa - (15.5 if self.reactor_type == "PWR" else 7.0)) < 1.5,
        }

    def _update_axial_channel(self, dt: float, power_fraction: float, flow: float) -> None:
        if self.reactor_type not in ("BWR", "RBMK"):
            self.channel_state.avg_void = 0.0
            self.channel_state.pressure_drop_mpa = 0.0
            return

        n = len(self.channel_state.quality)
        inlet_pressure = self.bop.pressure_mpa
        dp_total = 0.12 * power_fraction / max(flow, 0.1) + 0.015 * flow * flow
        self.channel_state.pressure_drop_mpa = max(0.0, dp_total)

        shape = []
        for i in range(n):
            z = (i + 0.5) / n
            shape.append(0.6 + 0.8 * (1.0 - (2.0 * z - 1.0) ** 2))
        ssum = sum(shape)
        self.channel_state.node_power = [power_fraction * s / ssum for s in shape]

        for i in range(n):
            z = (i + 0.5) / n
            local_p = inlet_pressure - z * self.channel_state.pressure_drop_mpa
            self.channel_state.node_pressure_mpa[i] = max(4.5, local_p)

            source = 0.16 * self.channel_state.node_power[i] / max(flow, 0.1)
            condense = 0.06 * self.channel_state.quality[i] * flow
            dquality = source - condense
            self.channel_state.quality[i] = max(0.0, min(1.0, self.channel_state.quality[i] + dt * dquality))

            slip_ratio = 1.0 + 2.0 * self.channel_state.quality[i]
            x = self.channel_state.quality[i]
            alpha = x / max(x + (1 - x) * slip_ratio * 0.05, 1e-4)
            self.channel_state.node_void[i] = max(0.0, min(0.97 if self.reactor_type == "BWR" else 0.92, alpha))

        self.channel_state.avg_void = sum(self.channel_state.node_void) / n
        self.plant.thermal.alpha = self.channel_state.avg_void

    def _update_bop(self, dt: float, thermal_power_w: float) -> None:
        c = self.controls
        if c.fault_feedwater_loss:
            c.feedwater_valve = min(c.feedwater_valve, 0.1)

        steam_gen = thermal_power_w / 3.0e9
        turbine_draw = c.turbine_valve * (0.0 if c.turbine_trip else 1.0)
        feed = c.feedwater_valve
        self.bop.steam_inventory += dt * (0.55 * steam_gen - 0.85 * turbine_draw)
        self.bop.steam_inventory += dt * 0.22 * feed
        self.bop.steam_inventory = max(0.0, min(2.0, self.bop.steam_inventory))

        pressure_nom = 15.5 if self.reactor_type == "PWR" else 7.0
        heater = c.pressurizer_heater if self.reactor_type == "PWR" else 0.2 * c.pressurizer_heater
        dpress = 0.15 * (self.bop.steam_inventory - 1.0) + 0.1 * heater - 0.08 * turbine_draw
        if self.reactor_type in ("BWR", "RBMK"):
            dpress -= 0.12 * self.channel_state.avg_void
        self.bop.pressure_mpa += dt * dpress
        self.bop.pressure_mpa += dt * 0.1 * (pressure_nom - self.bop.pressure_mpa)

        self.bop.condenser_temp_k += dt * (6.0 * turbine_draw - 0.06 * (self.bop.condenser_temp_k - 305.0))
        self.bop.feedwater_inventory += dt * (0.12 - 0.14 * feed)
        self.bop.feedwater_inventory = max(0.0, min(1.2, self.bop.feedwater_inventory))

    @staticmethod
    def _vote_2oo3(channels: List[float], threshold: float, high: bool = True) -> bool:
        votes = sum(1 for v in channels if (v > threshold if high else v < threshold))
        return votes >= 2

    def _protection(self, dt: float, true_fuel_temp: float, now: float) -> None:
        c = self.controls
        alarms: List[str] = []
        trips: List[str] = []
        diagnostics: List[str] = []

        meas_fuel = [ch.update(dt, true_fuel_temp) for ch in self.instrumentation.fuel_temp]
        meas_press = [ch.update(dt, self.bop.pressure_mpa) for ch in self.instrumentation.pressure]

        if max(meas_fuel) - min(meas_fuel) > 80.0:
            diagnostics.append("DIAG Fuel channel disagreement")
        if max(meas_press) - min(meas_press) > (0.8 if self.reactor_type == "PWR" else 0.35):
            diagnostics.append("DIAG Pressure channel disagreement")

        hi_fuel_alarm = self._vote_2oo3(meas_fuel, 1200.0, high=True)
        hi_press_alarm = self._vote_2oo3(meas_press, 16.2 if self.reactor_type == "PWR" else 7.8, high=True)

        if hi_fuel_alarm:
            alarms.append("HI Fuel temperature (2oo3)")
        if hi_press_alarm:
            alarms.append("HI Reactor pressure (2oo3)")
        if self.bop.feedwater_inventory < 0.2:
            alarms.append("MED Feedwater inventory low")
        if self.bop.condenser_temp_k > 330:
            alarms.append("MED Condenser temperature high")
        if c.fault_pump_degraded:
            alarms.append("MED Primary pump degraded")

        # apply alarm shelving
        alarms = [a for a in alarms if self.bop.shelved_until.get(a, 0.0) <= now]

        hi_hi_fuel_trip = self._vote_2oo3(meas_fuel, 1400.0, high=True)
        hi_hi_press_trip = self._vote_2oo3(meas_press, 16.8 if self.reactor_type == "PWR" else 8.3, high=True)

        self.bop.trip_timers["fuel"] = self.bop.trip_timers["fuel"] + dt if hi_hi_fuel_trip else max(0.0, self.bop.trip_timers["fuel"] - dt)
        self.bop.trip_timers["pressure"] = self.bop.trip_timers["pressure"] + dt if hi_hi_press_trip else max(0.0, self.bop.trip_timers["pressure"] - dt)

        if self.bop.trip_timers["fuel"] > 0.5 and not c.bypass_high_fuel_trip:
            trips.append("AUTO SCRAM: HH fuel temperature (2oo3 + delay)")
        if self.bop.trip_timers["pressure"] > 0.5 and not c.bypass_high_pressure_trip:
            trips.append("AUTO SCRAM: HH pressure (2oo3 + delay)")

        if hi_fuel_alarm and hi_press_alarm and not (c.bypass_high_fuel_trip or c.bypass_high_pressure_trip):
            trips.append("AUTO SCRAM: process 2oo2 vote")
        if c.turbine_trip:
            trips.append("Turbine trip active")

        if trips and not c.inhibit_auto_scram:
            c.scram = True
            self.bop.trip_latched = True

        if self.bop.trip_latched:
            c.scram = True
            if not trips:
                trips.append("SCRAM latched")

        self.bop.alarms = alarms
        self.bop.trips = trips
        self.bop.diagnostics = diagnostics

        active_set = set(alarms)
        self.bop.acked_alarms = [a for a in self.bop.acked_alarms if a in active_set]
        self.bop.unacked_alarms = [a for a in self.bop.unacked_alarms if a in active_set and a not in self.bop.acked_alarms]
        for a in alarms:
            if a not in self.bop.acked_alarms and a not in self.bop.unacked_alarms:
                self.bop.unacked_alarms.append(a)

        for a in alarms:
            msg = f"ALARM {a}"
            if not self.bop.event_log or self.bop.event_log[-1] != msg:
                self.bop.event_log.append(msg)
        for d in diagnostics:
            msg = f"DIAG {d}"
            if not self.bop.event_log or self.bop.event_log[-1] != msg:
                self.bop.event_log.append(msg)
        for tr in trips:
            msg = f"TRIP {tr}"
            if not self.bop.event_log or self.bop.event_log[-1] != msg:
                self.bop.event_log.append(msg)
        self.bop.event_log = self.bop.event_log[-120:]

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
            self._update_axial_channel(dt, base.power_fraction, flow)
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
            self.ctrl.rod_position = c.rod_position
            self.ctrl.recirc_pump = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            self._update_axial_channel(dt, self.plant.kinetics.n, flow)
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
                "void_fraction": self.channel_state.avg_void,
            }
        else:
            self.ctrl.rod_position = c.rod_position
            self.ctrl.pump_speed = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            self._update_axial_channel(dt, self.plant.kinetics.n, flow)
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
                "void_fraction": self.channel_state.avg_void,
            }

        self._update_bop(dt, thermal_power_w)
        self._protection(dt, snap["fuel_temp_k"], snap["time"])

        self.trend_power.append(snap["power_fraction"])
        self.trend_pressure.append(self.bop.pressure_mpa)
        self.trend_fuel.append(snap["fuel_temp_k"])

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
            pressure_drop_mpa=self.channel_state.pressure_drop_mpa,
            steam_inventory=self.bop.steam_inventory,
            condenser_temp_k=self.bop.condenser_temp_k,
            feedwater_inventory=self.bop.feedwater_inventory,
            alarms=list(self.bop.alarms),
            trips=list(self.bop.trips),
            diagnostics=list(self.bop.diagnostics),
            unacked_alarms=list(self.bop.unacked_alarms),
            trip_latched=self.bop.trip_latched,
            checklist=self._build_reset_checklist(),
            startup_checklist=self._build_startup_checklist(),
            event_log=list(self.bop.event_log),
        )
