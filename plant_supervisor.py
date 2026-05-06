from __future__ import annotations

import math
import pickle
import random
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Deque, Dict, List, Literal, Optional, Tuple

from physics_engine import ControlInputs as PWRControlInputs, PWRPlant
from physics_engine_bwr import BWRControlInputs, BWRPlant
from physics_engine_rbm import RBMKControlInputs, RBMKPlant


ReactorType = Literal["PWR", "BWR", "RBMK"]


@dataclass
class Alarm:
    id: str
    message: str
    priority: int   # 1=emergency/trip, 2=warning, 3=advisory
    time: float
    state: str = "unack"   # "unack" | "ack" | "clear"
    is_trip: bool = False
    first_out: bool = False


class AlarmManager:
    """Priority alarm system with acknowledgement, latching, and event log."""

    def __init__(self) -> None:
        self._active: Dict[str, Alarm] = {}
        self._first_out_id: Optional[str] = None
        self.log: Deque[Tuple[float, str, str]] = deque(maxlen=500)  # (time, id, event)

    def raise_alarm(self, id: str, priority: int, message: str, time: float, is_trip: bool = False) -> None:
        if id not in self._active:
            first_out = len(self._active) == 0
            alarm = Alarm(id=id, message=message, priority=priority,
                          time=time, state="unack", is_trip=is_trip, first_out=first_out)
            self._active[id] = alarm
            if first_out:
                self._first_out_id = id
            self.log.append((time, id, "raised"))

    def clear_alarm(self, id: str) -> None:
        if id in self._active:
            t = self._active[id].time
            del self._active[id]
            self.log.append((t, id, "cleared"))
            if self._first_out_id == id:
                self._first_out_id = None

    def acknowledge(self, id: str) -> None:
        if id in self._active and self._active[id].state == "unack":
            self._active[id].state = "ack"
            self.log.append((self._active[id].time, id, "acknowledged"))

    def acknowledge_all(self) -> None:
        for alarm in self._active.values():
            if alarm.state == "unack":
                alarm.state = "ack"

    def active_alarms(self) -> List[Alarm]:
        return sorted(self._active.values(), key=lambda a: a.priority)

    def clear_all(self) -> None:
        self._active.clear()
        self._first_out_id = None


@dataclass
class AutoController:
    """Simple proportional auto controller for one loop."""
    auto: bool = False
    setpoint: float = 0.0
    Kp: float = 1.0
    output_min: float = 0.0
    output_max: float = 1.0

    def compute(self, measurement: float) -> Optional[float]:
        if not self.auto:
            return None
        err = self.setpoint - measurement
        return min(self.output_max, max(self.output_min, 0.5 + self.Kp * err))


@dataclass
class SupervisorControls:
    rod_position: float = 0.0
    flow: float = 1.0
    turbine_valve: float = 1.0
    feedwater_valve: float = 0.7
    pressurizer_heater: float = 0.5
    boration_rate: float = 0.0
    dilution_rate: float = 0.0
    startup_permit: bool = False
    turbine_trip: bool = False
    scram: bool = False
    fault_pump_degraded: bool = False
    fault_feedwater_loss: bool = False
    fault_loca_break_area: float = 0.0   # m²; 0 = no LOCA


@dataclass
class BalanceOfPlantState:
    pressure_mpa: float = 7.0
    steam_inventory: float = 1.0
    condenser_temp_k: float = 305.0
    feedwater_inventory: float = 1.0
    omega_rcp: float = 1.0       # RCP fractional speed (coast-down)
    porv_open: bool = False       # pressurizer PORV status
    eccs_actuated: bool = False   # ECCS injection active
    alarms: List[str] = field(default_factory=list)
    trips: List[str] = field(default_factory=list)


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
    decay_heat_fraction: float = 0.0
    boron_ppm: float = 800.0
    alarm_objects: List[Alarm] = field(default_factory=list)
    scram_reset_message: str = ""
    porv_open: bool = False
    eccs_actuated: bool = False
    loca_area: float = 0.0
    channel_noise: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])


class PlantSupervisor:
    """Adds protection system, balance-of-plant and fault handling around core models."""

    def __init__(self, reactor_type: ReactorType = "PWR") -> None:
        self.reactor_type: ReactorType = reactor_type
        self.controls = SupervisorControls()
        self.bop = BalanceOfPlantState()
        self.alarm_mgr = AlarmManager()
        self._sim_time: float = 0.0
        self._scram_reset_message: str = ""
        self._boron_ppm: float = 800.0      # PWR boron concentration
        self.event_log: Deque[Tuple[float, str, str]] = deque(maxlen=500)
        # 2-of-3 protection channel noise (fractional)
        self._channel_noise: List[float] = [0.0, 0.0, 0.0]
        # Auto controllers
        self.auto_rod      = AutoController(setpoint=550.0, Kp=-0.006, output_min=0.0, output_max=1.0)
        self.auto_pressure = AutoController(setpoint=15.5,  Kp=0.25,  output_min=0.0, output_max=1.0)
        self._init_plant()

    def _init_plant(self) -> None:
        if self.reactor_type == "PWR":
            self.plant = PWRPlant()
            self.bop.pressure_mpa = 15.5
        elif self.reactor_type == "BWR":
            self.plant = BWRPlant()
            self.bop.pressure_mpa = 7.0
        else:
            self.plant = RBMKPlant()
            self.bop.pressure_mpa = 7.0
        self._init_ctrl()

    def _init_ctrl(self) -> None:
        if self.reactor_type == "PWR":
            self.ctrl = PWRControlInputs()
        elif self.reactor_type == "BWR":
            self.ctrl = BWRControlInputs()
        else:
            self.ctrl = RBMKControlInputs()

    def select(self, reactor_type: ReactorType) -> None:
        self.reactor_type = reactor_type
        self.controls = SupervisorControls()
        self.bop = BalanceOfPlantState()
        self.alarm_mgr.clear_all()
        self._sim_time = 0.0
        self._scram_reset_message = ""
        self._boron_ppm = 800.0
        self.auto_rod.auto = False
        self.auto_pressure.auto = False
        self._init_plant()

    def reset_scram(self) -> str:
        """Attempt to reset a latched SCRAM. Returns status message."""
        if not self.controls.scram:
            return "No SCRAM latched"
        active_trips = [a for a in self.alarm_mgr.active_alarms() if a.is_trip]
        if active_trips:
            return f"RESET BLOCKED: trip active — {active_trips[0].message}"
        if self.controls.rod_position < 0.95:
            return "RESET BLOCKED: rods must be >95% inserted"
        self.controls.scram = False
        self.plant.scrammed = False
        self.alarm_mgr.acknowledge_all()
        self.event_log.append((self._sim_time, "SCRAM_RESET", "SCRAM reset approved by operator"))
        return "SCRAM RESET APPROVED"

    def acknowledge_alarms(self) -> None:
        self.alarm_mgr.acknowledge_all()
        self.event_log.append((self._sim_time, "ACK_ALL", "Operator acknowledged all alarms"))

    def save_ic(self, path: str) -> None:
        state = {
            "reactor_type": self.reactor_type,
            "controls": self.controls,
            "bop": self.bop,
            "plant": self.plant,
            "sim_time": self._sim_time,
            "boron_ppm": self._boron_ppm,
        }
        with open(path, "wb") as f:
            pickle.dump(state, f)

    def load_ic(self, path: str) -> None:
        with open(path, "rb") as f:
            state = pickle.load(f)
        self.reactor_type = state["reactor_type"]
        self.controls = state["controls"]
        self.bop = state["bop"]
        self.plant = state["plant"]
        self._sim_time = state["sim_time"]
        self._boron_ppm = state.get("boron_ppm", 800.0)
        self.alarm_mgr.clear_all()
        self._init_ctrl()

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

        # RCP coast-down model
        if c.fault_pump_degraded:
            self.bop.omega_rcp += dt * (-self.bop.omega_rcp / 12.0)
        else:
            self.bop.omega_rcp += dt * (1.0 - self.bop.omega_rcp) / 30.0
        self.bop.omega_rcp = max(0.04, min(1.0, self.bop.omega_rcp))  # nat. circ floor 4%

        # PWR boron concentration (first-order mixing)
        if self.reactor_type == "PWR":
            V_primary = 300.0  # m³ (nominal PWR primary inventory)
            d_boron = (c.boration_rate * 2000.0 - c.dilution_rate * self._boron_ppm) / V_primary
            self._boron_ppm = max(0.0, min(4000.0, self._boron_ppm + d_boron * dt))

        steam_gen = thermal_power_w / 3.0e9
        turbine_draw = c.turbine_valve * (0.0 if c.turbine_trip else 1.0)
        feed = c.feedwater_valve
        self.bop.steam_inventory += dt * (0.6 * steam_gen - 0.9 * turbine_draw)
        self.bop.steam_inventory += dt * 0.2 * feed
        self.bop.steam_inventory = max(0.0, min(2.0, self.bop.steam_inventory))

        pressure_nom = 15.5 if self.reactor_type == "PWR" else 7.0

        if self.reactor_type == "PWR":
            # Pressurizer model: heater raises pressure, spray lowers it
            # PORV at 16.55 MPa, SRV at 17.2 MPa
            heater_dp = 0.12 * c.pressurizer_heater           # heater adds pressure
            spray_dp  = -0.10 * max(0.0, c.pressurizer_heater - 0.8) * 5.0  # spray ~above 80%
            steam_dp  = 0.08 * (self.bop.steam_inventory - 1.0)
            self.bop.porv_open = self.bop.pressure_mpa > 16.55
            porv_dp   = -0.40 if self.bop.porv_open else 0.0
            srv_dp    = -1.20 if self.bop.pressure_mpa > 17.2  else 0.0
            dpress = heater_dp + spray_dp + steam_dp + porv_dp + srv_dp
            self.bop.pressure_mpa += dt * dpress
            self.bop.pressure_mpa += dt * 0.06 * (pressure_nom - self.bop.pressure_mpa)

            # LOCA: primary break depressurizes and drains inventory
            A_break = c.fault_loca_break_area
            if A_break > 0 and self.bop.pressure_mpa > 0.1:
                dP_Pa = max(0.0, (self.bop.pressure_mpa - 0.1) * 1e6)
                dp_loca = 45.0 * A_break * math.sqrt(dP_Pa / 1e6)  # MPa/s scale
                self.bop.pressure_mpa -= dp_loca * dt
                self.bop.pressure_mpa = max(0.1, self.bop.pressure_mpa)

            # ECCS auto-actuates at 11.7 MPa
            if self.bop.pressure_mpa < 11.7 and A_break > 0:
                self.bop.eccs_actuated = True
            elif self.bop.pressure_mpa > 14.0:
                self.bop.eccs_actuated = False
        else:
            heater = 0.2 * c.pressurizer_heater
            dpress = 0.12 * (self.bop.steam_inventory - 1.0) + 0.09 * heater - 0.08 * turbine_draw
            self.bop.pressure_mpa += dt * dpress
            self.bop.pressure_mpa += dt * 0.1 * (pressure_nom - self.bop.pressure_mpa)

        self.bop.condenser_temp_k += dt * (6.0 * turbine_draw - 0.06 * (self.bop.condenser_temp_k - 305.0))
        self.bop.feedwater_inventory += dt * (0.12 - 0.14 * feed)
        self.bop.feedwater_inventory = max(0.0, min(1.2, self.bop.feedwater_inventory))

    def _protection(self, snap: Dict[str, float]) -> None:
        t = self._sim_time
        am = self.alarm_mgr
        noise = self._channel_noise  # per-channel measurement noise (fractional)

        def warn(id: str, msg: str, cond: bool) -> None:
            if cond:
                am.raise_alarm(id, 2, msg, t, is_trip=False)
            else:
                am.clear_alarm(id)

        def trip(id: str, msg: str, value: float, setpoint: float) -> None:
            """2-of-3 channel voting trip."""
            channels = [value * (1.0 + n) for n in noise]
            votes = sum(1 for ch in channels if ch > setpoint)
            cond = votes >= 2
            if cond:
                am.raise_alarm(id, 1, msg, t, is_trip=True)
                self.controls.scram = True
                self.event_log.append((t, "AUTO_SCRAM", msg))
                # Channel disagreement advisory
                if votes < 3:
                    am.raise_alarm(id + "_ch", 3, f"Protection channel split: {id}", t)
            else:
                am.clear_alarm(id)
                am.clear_alarm(id + "_ch")

        p_hi   = 16.2 if self.reactor_type == "PWR" else 7.8
        p_trip = 16.8 if self.reactor_type == "PWR" else 8.3

        warn("fuel_hi",    "Fuel temperature high",      snap["fuel_temp_k"] > 1200)
        warn("press_hi",   "Reactor pressure high",      self.bop.pressure_mpa > p_hi)
        warn("fw_lo",      "Feedwater inventory low",    self.bop.feedwater_inventory < 0.2)
        warn("cond_hi",    "Condenser temperature high", self.bop.condenser_temp_k > 330)
        warn("pump_deg",   "Primary pump degraded",      self.controls.fault_pump_degraded)
        warn("flux_hi",    "Neutron flux high",          snap["power_fraction"] > 1.1)
        warn("coolant_hi", "Coolant temperature high",   snap.get("coolant_temp_k", 0) > 610 and self.reactor_type == "PWR")
        warn("flow_lo",    "Primary flow low",           self.bop.omega_rcp < 0.87 and self.reactor_type == "PWR")
        warn("porv",       "PORV open — pressure relief",self.bop.porv_open)

        trip("fuel_trip",    "AUTO SCRAM: high fuel temperature",   snap["fuel_temp_k"],              1400.0)
        trip("press_trip",   "AUTO SCRAM: high reactor pressure",   self.bop.pressure_mpa,            p_trip)
        trip("flux_trip",    "AUTO SCRAM: high neutron flux",       snap["power_fraction"],            1.2)
        if self.reactor_type == "PWR":
            trip("coolant_trip", "AUTO SCRAM: high coolant temperature", snap.get("coolant_temp_k", 0), 616.0)
            if self.bop.omega_rcp < 0.87 and not self.controls.fault_pump_degraded and snap["power_fraction"] > 0.1:
                am.raise_alarm("flow_trip", 1, "AUTO SCRAM: low primary flow", t, is_trip=True)
                self.controls.scram = True
                self.event_log.append((t, "AUTO_SCRAM", "Low primary flow"))
            else:
                am.clear_alarm("flow_trip")

        # LOCA and ECCS alarms
        if self.controls.fault_loca_break_area > 0:
            am.raise_alarm("loca", 1, "LOCA: primary boundary breach", t, is_trip=True)
            self.controls.scram = True
        else:
            am.clear_alarm("loca")
        if self.bop.eccs_actuated:
            am.raise_alarm("eccs", 1, "ECCS INJECTION ACTIVE", t, is_trip=False)
        else:
            am.clear_alarm("eccs")

        if self.controls.turbine_trip:
            am.raise_alarm("tt", 2, "Turbine trip active", t, is_trip=False)
        else:
            am.clear_alarm("tt")

        active = self.alarm_mgr.active_alarms()
        self.bop.alarms = [a.message for a in active if not a.is_trip]
        self.bop.trips  = [a.message for a in active if a.is_trip]

    def step(self, dt: float) -> UnifiedSnapshot:
        self._sim_time += dt
        c = self.controls
        # Refresh 2-of-3 channel noise each step
        self._channel_noise = [random.gauss(0, 0.003) for _ in range(3)]
        # Apply auto controllers (override operator inputs for active loops)
        if self.reactor_type == "PWR":
            rod_out = self.auto_rod.compute(
                self.plant.thermal.T_cool if hasattr(self.plant, "thermal") else 550.0
            )
            if rod_out is not None:
                c.rod_position = rod_out
            press_out = self.auto_pressure.compute(self.bop.pressure_mpa)
            if press_out is not None:
                c.pressurizer_heater = press_out
        # Use RCP speed to scale flow (coast-down model); pump fault sets omega decay
        flow = c.flow * self.bop.omega_rcp
        if not c.startup_permit:
            c.rod_position = max(c.rod_position, 0.85)

        if self.reactor_type == "PWR":
            self.ctrl.rod_position = c.rod_position
            self.ctrl.primary_flow = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            # Inject boron reactivity: -8 pcm/ppm relative to 800 ppm nominal
            self.plant.params.external_reactivity = -8e-5 * (self._boron_ppm - 800.0)
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
                "decay_heat_fraction": self.plant.decay_heat.fraction,
            }
        elif self.reactor_type == "BWR":
            self.ctrl.rod_position = c.rod_position
            self.ctrl.recirc_pump = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            base = self.plant.step(dt, self.ctrl)
            self._mechanistic_void_step(dt)
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
                "void_fraction": self.plant.thermal.alpha,
                "decay_heat_fraction": self.plant.decay_heat.fraction,
            }
        else:
            self.ctrl.rod_position = c.rod_position
            self.ctrl.pump_speed = flow
            self.ctrl.turbine_valve = 0.0 if c.turbine_trip else c.turbine_valve
            self.ctrl.scram = c.scram
            base = self.plant.step(dt, self.ctrl)
            self._mechanistic_void_step(dt)
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
                "void_fraction": self.plant.thermal.alpha,
                "decay_heat_fraction": self.plant.decay_heat.fraction,
            }

        self._update_bop(dt, thermal_power_w)
        self._protection(snap)

        msg = self._scram_reset_message
        self._scram_reset_message = ""
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
            decay_heat_fraction=snap["decay_heat_fraction"],
            boron_ppm=self._boron_ppm,
            alarm_objects=self.alarm_mgr.active_alarms(),
            scram_reset_message=msg,
            porv_open=self.bop.porv_open,
            eccs_actuated=self.bop.eccs_actuated,
            loca_area=c.fault_loca_break_area,
            channel_noise=list(self._channel_noise),
        )
