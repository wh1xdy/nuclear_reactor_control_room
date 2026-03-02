"""Validation harness with qualitative acceptance checks.

Step 2 focus: verify protection-system behavior with alarm acknowledgement
and SCRAM latch reset permissives in addition to transient stability.
"""

from plant_supervisor import PlantSupervisor


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
    else:  # RBMK
        c.rod_position = 0.70
        c.flow = 1.0
        c.turbine_valve = 0.80

    c.feedwater_valve = 0.7

    peak_power = 0.0
    peak_void = 0.0
    pre_steps = 200
    for _ in range(pre_steps):
        s = sup.step(0.1)
        peak_power = max(peak_power, s.power_fraction)
        peak_void = max(peak_void, s.void_fraction)

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
    for _ in range(800):
        s = sup.step(0.1)
        observed_trip = observed_trip or bool(s.trips)
        observed_scram = observed_scram or c.scram
        observed_latch = observed_latch or s.trip_latched
        observed_unacked = max(observed_unacked, len(s.unacked_alarms))

    # Acknowledge should clear unacked list while conditions may still alarm.
    sup.acknowledge_alarms()
    s_after_ack = sup.step(0.1)
    ack_cleared = len(s_after_ack.unacked_alarms) == 0

    # Remove turbine trip and reduce stress to satisfy reset permissives.
    c.turbine_trip = False
    c.pressurizer_heater = 0.0
    c.feedwater_valve = 0.9
    c.rod_position = 1.0
    c.flow = 1.0
    for _ in range(200):
        sup.step(0.1)

    latch_reset = sup.reset_trip_latch()
    s_after_reset = sup.step(0.1)

    return {
        "trip": observed_trip,
        "scram": observed_scram,
        "latch_seen": observed_latch,
        "peak_unacked": observed_unacked,
        "ack_cleared": ack_cleared,
        "reset_ok": latch_reset,
        "latched_after_reset": s_after_reset.trip_latched,
        "pressure_mpa": s_after_reset.pressure_mpa,
    }


def main() -> None:
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
    assert trips["trip"], f"Protection trips did not trigger: {trips}"
    assert trips["scram"], f"Trip did not latch SCRAM: {trips}"
    assert trips["latch_seen"], f"Trip latch was never observed: {trips}"
    assert trips["peak_unacked"] > 0, f"No unacked alarms observed: {trips}"
    assert trips["ack_cleared"], f"Alarm acknowledge did not clear unacked set: {trips}"
    assert trips["reset_ok"], f"Latch reset permissive failed: {trips}"
    assert not trips["latched_after_reset"], f"Latch remained set after reset: {trips}"

    print("Validation summary")
    print("PWR", pwr)
    print("BWR", bwr)
    print("RBMK", rbmk)
    print("TRIPS", trips)


if __name__ == "__main__":
    main()
