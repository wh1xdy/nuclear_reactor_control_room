"""Validation harness with qualitative acceptance checks.

Step 1 focus: ensure all reactor models remain numerically stable without
hard neutron-population clipping and that SCRAM behavior is effective for
representative operator conditions.
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
    for _ in range(800):
        s = sup.step(0.1)
        observed_trip = observed_trip or bool(s.trips)
        observed_scram = observed_scram or c.scram

    return {
        "trip": observed_trip,
        "scram": observed_scram,
        "pressure_mpa": s.pressure_mpa,
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

    print("Validation summary")
    print("PWR", pwr)
    print("BWR", bwr)
    print("RBMK", rbmk)
    print("TRIPS", trips)


if __name__ == "__main__":
    main()
