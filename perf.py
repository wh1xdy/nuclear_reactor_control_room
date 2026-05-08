"""
perf.py — physics-step benchmark against M5 Pro targets from CLAUDE.md

Targets:
  1× RT:   supervisor.step(dt) < 1 ms
  600× RT: ~600 substeps per frame must fit in < 6% of 16.6 ms → < 100 µs per step call
           (each step(9.6 s) fires internally ~9600 substeps; total wall budget ~96 ms/frame
            but we want the whole physics budget <6% of frame = <1 ms per step(dt=0.0167))

Usage:
  python perf.py
"""

import statistics
import time

from plant_supervisor import PlantSupervisor

REPS = 2000
WARMUP = 200


def bench(reactor_type: str, sim_dt: float) -> dict:
    sv = PlantSupervisor(reactor_type)
    sv.controls.startup_permit = True
    sv.controls.rod_position = 0.5
    # Warmup
    for _ in range(WARMUP):
        sv.step(sim_dt)
    # Timed
    times = []
    for _ in range(REPS):
        t0 = time.perf_counter()
        sv.step(sim_dt)
        times.append((time.perf_counter() - t0) * 1e6)  # µs
    return {
        "mean_us":   statistics.mean(times),
        "median_us": statistics.median(times),
        "p99_us":    sorted(times)[int(0.99 * REPS)],
        "max_us":    max(times),
    }


def main():
    print(f"{'Reactor':<8}  {'dt':>8}  {'mean':>8}  {'median':>8}  {'p99':>8}  {'max':>8}  target")
    print("-" * 72)

    for rt in ("PWR", "BWR", "RBMK"):
        # 1× real-time: dt = 1/60 s ≈ 16.7 ms
        r1 = bench(rt, 1 / 60)
        ok1 = "OK" if r1["p99_us"] < 1000 else "SLOW"
        print(f"{rt:<8}  {'1/60 s':>8}  {r1['mean_us']:>7.1f}µs  {r1['median_us']:>7.1f}µs  "
              f"{r1['p99_us']:>7.1f}µs  {r1['max_us']:>7.1f}µs  <1 ms {ok1}")

        # 600× real-time: dt = 600/60 = 10 s
        r6 = bench(rt, 600 / 60)
        ok6 = "OK" if r6["p99_us"] < 96_000 else "SLOW"
        print(f"{rt:<8}  {'600/60 s':>8}  {r6['mean_us']:>7.1f}µs  {r6['median_us']:>7.1f}µs  "
              f"{r6['p99_us']:>7.1f}µs  {r6['max_us']:>7.1f}µs  <96 ms {ok6}")
    print()
    print("Note: targets are for Apple M5 Pro (see CLAUDE.md). Times in µs.")


if __name__ == "__main__":
    main()
