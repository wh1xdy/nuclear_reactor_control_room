"""Validation harness with qualitative and benchmark-envelope checks."""

from __future__ import annotations

import json
from pathlib import Path

from plant_supervisor import PlantSupervisor

TIMES = [10, 20, 30, 40]
REPO_ROOT = Path(__file__).resolve().parent


def run_scram_transient(reactor: str) -> dict:
    sup = PlantSupervisor(reactor)  # type: ignore[arg-type]
    c = sup.controls
    c.startup_permit = True

    if reactor == "PWR":
        c.rod_position = 0.55
        c.flow = 1.0
        c.turbine_valve = 1.0
    elif reactor == "BWR":
        c.rod_position = 0.80
        c.flow = 1.0
        c.turbine_valve = 0.85
    else:
        c.rod_position = 0.70
        c.flow = 1.0
        c.turbine_valve = 0.80

    c.feedwater_valve = 0.7

    peak_power = 0.0
    peak_void = 0.0
    points: dict[int, dict[str, float]] = {}
    idx = 0

    for _ in range(500):
        s = sup.step(0.1)
        peak_power = max(peak_power, s.power_fraction)
        peak_void = max(peak_void, s.void_fraction)
        if idx < len(TIMES) and s.time >= TIMES[idx]:
            points[TIMES[idx]] = {
                "power": s.power_fraction,
                "pressure": s.pressure_mpa,
                "void": s.void_fraction,
            }
            idx += 1

    c.scram = True
    post_scram = []
    for _ in range(250):
        s = sup.step(0.1)
        post_scram.append(s.power_fraction)

    return {
        "peak": peak_power,
        "scram_tail": post_scram[-1],
        "scram_avg_last20": sum(post_scram[-20:]) / 20.0,
        "peak_void": peak_void,
        "points": points,
    }


def run_trip_transient() -> dict:
    sup = PlantSupervisor("PWR")
    c = sup.controls
    c.startup_permit = True
    c.rod_position = 0.4
    c.flow = 0.7
    c.turbine_valve = 1.0
    c.feedwater_valve = 0.2
    c.pressurizer_heater = 1.0
    c.turbine_trip = True

    observed_trip = False
    observed_scram = False
    observed_latch = False
    observed_unacked = 0
    observed_diag = False

    sup.instrumentation.pressure[2].failed_low = True

    for _ in range(800):
        s = sup.step(0.1)
        observed_trip = observed_trip or bool(s.trips)
        observed_scram = observed_scram or c.scram
        observed_latch = observed_latch or s.trip_latched
        observed_unacked = max(observed_unacked, len(s.unacked_alarms))
        observed_diag = observed_diag or bool(s.diagnostics)

    # Shelve current alarms, then re-step and ensure unacked can drop.
    sup.shelve_current_alarms(now=s.time, seconds=20.0)
    s_after_shelve = sup.step(0.1)
    shelve_reduced_unacked = len(s_after_shelve.unacked_alarms) <= observed_unacked

    sup.acknowledge_alarms()
    s_after_ack = sup.step(0.1)
    ack_cleared = len(s_after_ack.unacked_alarms) == 0

    c.turbine_trip = False
    c.pressurizer_heater = 0.0
    c.feedwater_valve = 0.9
    c.rod_position = 1.0
    c.flow = 1.0
    for _ in range(200):
        sup.step(0.1)

    latch_reset = sup.reset_trip_latch()
    c.inhibit_auto_scram = True
    s_after_reset = sup.step(0.1)

    return {
        "trip": observed_trip,
        "scram": observed_scram,
        "latch_seen": observed_latch,
        "peak_unacked": observed_unacked,
        "diag_seen": observed_diag,
        "shelve_reduced_unacked": shelve_reduced_unacked,
        "ack_cleared": ack_cleared,
        "reset_ok": latch_reset,
        "latched_after_reset": s_after_reset.trip_latched,
        "event_log_len": len(s_after_reset.event_log),
        "startup_checklist_present": bool(s_after_reset.startup_checklist),
    }


def assert_in_envelope(result: dict, env: dict, reactor: str) -> None:
    for i, t in enumerate(TIMES):
        p = result["points"][t]
        for key in ("power", "pressure", "void"):
            lo, hi = env[key][i]
            assert lo <= p[key] <= hi, f"{reactor} {key}@{t}s outside envelope: {p[key]} not in [{lo},{hi}]"


def assert_close_to_reference(result: dict, reference: dict, reactor: str) -> None:
    ref_by_time = {int(row["t"]): row for row in reference}
    for t in TIMES:
        p = result["points"][t]
        r = ref_by_time[t]

        # Keep tolerance loose enough for educational-model drift while still
        # catching larger regressions in trend shape.
        assert abs(p["power"] - r["power"]) <= 0.08, f"{reactor} power@{t}s diverges from reference: {p['power']} vs {r['power']}"
        assert abs(p["pressure"] - r["pressure"]) <= 0.8, f"{reactor} pressure@{t}s diverges from reference: {p['pressure']} vs {r['pressure']}"
        assert abs(p["void"] - r["void"]) <= 0.3, f"{reactor} void@{t}s diverges from reference: {p['void']} vs {r['void']}"


def main() -> None:
    envelopes = json.loads((REPO_ROOT / "data" / "benchmark_envelopes.json").read_text())
    references = json.loads((REPO_ROOT / "data" / "reference_transients.json").read_text())

    pwr = run_scram_transient("PWR")
    bwr = run_scram_transient("BWR")
    rbmk = run_scram_transient("RBMK")
    trips = run_trip_transient()

    assert pwr["peak"] < 10.0, f"PWR runaway detected: {pwr}"
    assert bwr["peak"] < 10.0, f"BWR runaway detected: {bwr}"
    assert rbmk["peak"] < 10.0, f"RBMK runaway detected: {rbmk}"

    assert pwr["scram_tail"] < 0.4 and pwr["scram_avg_last20"] < 0.45, f"PWR SCRAM ineffective: {pwr}"
    assert bwr["scram_tail"] < 0.5 and bwr["scram_avg_last20"] < 0.6, f"BWR SCRAM ineffective: {bwr}"
    assert rbmk["scram_tail"] < 0.7 and rbmk["scram_avg_last20"] < 0.8, f"RBMK SCRAM ineffective: {rbmk}"

    assert rbmk["peak_void"] >= 0.01, f"RBMK void unexpectedly low: {rbmk}"

    assert_in_envelope(pwr, envelopes["PWR"], "PWR")
    assert_in_envelope(bwr, envelopes["BWR"], "BWR")
    assert_in_envelope(rbmk, envelopes["RBMK"], "RBMK")
    assert_close_to_reference(pwr, references["PWR"], "PWR")
    assert_close_to_reference(bwr, references["BWR"], "BWR")
    assert_close_to_reference(rbmk, references["RBMK"], "RBMK")

    assert trips["trip"], f"Protection trips did not trigger: {trips}"
    assert trips["scram"], f"Trip did not latch SCRAM: {trips}"
    assert trips["latch_seen"], f"Trip latch was never observed: {trips}"
    assert trips["peak_unacked"] > 0, f"No unacked alarms observed: {trips}"
    assert trips["diag_seen"], f"No diagnostics observed under channel disagreement: {trips}"
    assert trips["shelve_reduced_unacked"], f"Alarm shelving path ineffective: {trips}"
    assert trips["ack_cleared"], f"Alarm acknowledge did not clear unacked set: {trips}"
    assert trips["reset_ok"], f"Latch reset permissive failed: {trips}"
    assert not trips["latched_after_reset"], f"Latch remained set after reset: {trips}"
    assert trips["event_log_len"] > 0, f"Event log did not record protection sequence: {trips}"
    assert trips["startup_checklist_present"], f"Startup checklist missing in snapshot: {trips}"

    print("Validation summary")
    print("PWR", pwr)
    print("BWR", bwr)
    print("RBMK", rbmk)
    print("TRIPS", trips)


if __name__ == "__main__":
    main()
