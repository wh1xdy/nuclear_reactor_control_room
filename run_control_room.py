import time
import pygame

from plant_supervisor import PlantSupervisor


def clamp(x, lo, hi):
    return lo if x < lo else hi if x > hi else x


class Slider:
    def __init__(self, label, v=0.0, lo=0.0, hi=1.0):
        self.label = label
        self.v = float(v)
        self.lo = float(lo)
        self.hi = float(hi)

    def set(self, v):
        self.v = clamp(v, self.lo, self.hi)


def draw_trend(screen, x, y, w, h, values, color, y_min, y_max):
    pygame.draw.rect(screen, (25, 28, 35), (x, y, w, h), border_radius=6)
    pygame.draw.rect(screen, (60, 66, 80), (x, y, w, h), 1, border_radius=6)
    if len(values) < 2:
        return
    pts = []
    n = len(values)
    for i, v in enumerate(values):
        px = x + int(i * (w - 4) / max(1, n - 1)) + 2
        vv = clamp((v - y_min) / max(1e-6, y_max - y_min), 0.0, 1.0)
        py = y + h - 2 - int(vv * (h - 4))
        pts.append((px, py))
    pygame.draw.lines(screen, color, False, pts, 2)


def main():
    pygame.init()
    screen = pygame.display.set_mode((1450, 820))
    pygame.display.set_caption("Nuclear Reactor Control Room")
    font = pygame.font.SysFont("Menlo", 18)
    big = pygame.font.SysFont("Menlo", 24, bold=True)
    clock = pygame.time.Clock()

    supervisor = PlantSupervisor("PWR")
    controls = supervisor.controls

    rod = Slider("Rod position", 0.85)
    flow = Slider("Primary flow", 1.0)
    valve = Slider("Turbine valve", 1.0)
    feed = Slider("Feedwater valve", 0.7)
    press = Slider("Pressurizer heater", 0.5)

    last = time.time()

    def select(rt: str):
        supervisor.select(rt)
        nonlocal controls
        controls = supervisor.controls
        flow.label = "Primary flow" if rt == "PWR" else ("Recirc pump" if rt == "BWR" else "Pump speed")

    running = True
    while running:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                running = False
            elif ev.type == pygame.KEYDOWN:
                if ev.key == pygame.K_ESCAPE:
                    running = False
                elif ev.key == pygame.K_1:
                    select("PWR")
                elif ev.key == pygame.K_2:
                    select("BWR")
                elif ev.key == pygame.K_3:
                    select("RBMK")
                elif ev.key == pygame.K_SPACE:
                    controls.scram = True
                elif ev.key == pygame.K_r:
                    select(supervisor.reactor_type)
                elif ev.key == pygame.K_t:
                    controls.turbine_trip = not controls.turbine_trip
                elif ev.key == pygame.K_p:
                    controls.startup_permit = not controls.startup_permit
                elif ev.key == pygame.K_z:
                    controls.fault_pump_degraded = not controls.fault_pump_degraded
                elif ev.key == pygame.K_x:
                    controls.fault_feedwater_loss = not controls.fault_feedwater_loss
                elif ev.key == pygame.K_c:
                    supervisor.acknowledge_alarms()
                elif ev.key == pygame.K_l:
                    supervisor.reset_trip_latch()
                elif ev.key == pygame.K_b:
                    controls.bypass_high_pressure_trip = not controls.bypass_high_pressure_trip
                elif ev.key == pygame.K_g:
                    controls.bypass_high_fuel_trip = not controls.bypass_high_fuel_trip
                elif ev.key == pygame.K_i:
                    controls.inhibit_auto_scram = not controls.inhibit_auto_scram

        keys = pygame.key.get_pressed()
        step = 0.004
        if keys[pygame.K_LSHIFT] or keys[pygame.K_RSHIFT]:
            step = 0.02

        if keys[pygame.K_w]:
            rod.set(rod.v - step)
        if keys[pygame.K_s]:
            rod.set(rod.v + step)
        if keys[pygame.K_a]:
            flow.set(flow.v - step)
        if keys[pygame.K_d]:
            flow.set(flow.v + step)
        if keys[pygame.K_q]:
            valve.set(valve.v - step)
        if keys[pygame.K_e]:
            valve.set(valve.v + step)
        if keys[pygame.K_f]:
            feed.set(feed.v + step)
        if keys[pygame.K_v]:
            feed.set(feed.v - step)
        if keys[pygame.K_h]:
            press.set(press.v + step)
        if keys[pygame.K_n]:
            press.set(press.v - step)

        controls.rod_position = rod.v
        controls.flow = flow.v
        controls.turbine_valve = valve.v
        controls.feedwater_valve = feed.v
        controls.pressurizer_heater = press.v

        now = time.time()
        dt = clamp(now - last, 0.0, 0.2)
        last = now

        snapshot = supervisor.step(dt)

        screen.fill((10, 12, 16))
        header = big.render(
            f"Reactor: {snapshot.reactor_type}  1/2/3 switch  SPACE=SCRAM  T=turbine trip  P=startup permit",
            True,
            (220, 230, 240),
        )
        screen.blit(header, (20, 18))

        y = 70
        for sl in (rod, flow, valve, feed, press):
            txt = font.render(f"{sl.label}: {sl.v:.3f}", True, (200, 210, 220))
            screen.blit(txt, (20, y))
            bx, by, bw, bh = 20, y + 22, 360, 14
            pygame.draw.rect(screen, (35, 40, 52), (bx, by, bw, bh), border_radius=5)
            fill = int(bw * (sl.v - sl.lo) / (sl.hi - sl.lo))
            pygame.draw.rect(screen, (80, 160, 255), (bx, by, fill, bh), border_radius=5)
            y += 52

        status_lines = [
            f"Startup permit: {'ON' if controls.startup_permit else 'OFF'}",
            f"Turbine trip: {'ON' if controls.turbine_trip else 'OFF'}",
            f"Pump fault (Z): {'ON' if controls.fault_pump_degraded else 'OFF'}",
            f"Feedwater fault (X): {'ON' if controls.fault_feedwater_loss else 'OFF'}",
            f"Bypass pressure trip (B): {'ON' if controls.bypass_high_pressure_trip else 'OFF'}",
            f"Bypass fuel trip (G): {'ON' if controls.bypass_high_fuel_trip else 'OFF'}",
            f"Inhibit auto SCRAM (I): {'ON' if controls.inhibit_auto_scram else 'OFF'}",
            f"SCRAM latch: {'ON' if snapshot.trip_latched else 'OFF'}",
            f"Unacked alarms: {len(snapshot.unacked_alarms)} (C=ack)",
        ]
        for s in status_lines:
            screen.blit(font.render(s, True, (185, 200, 215)), (20, y))
            y += 22

        x0, y0 = 420, 70

        def line(label, value, unit=""):
            nonlocal y0
            t = font.render(f"{label:<24} {value} {unit}", True, (220, 230, 240))
            screen.blit(t, (x0, y0))
            y0 += 22

        line("Time", f"{snapshot.time:8.1f}", "s")
        line("Core power", f"{snapshot.power_fraction*100:8.2f}", "%")
        line("Thermal power", f"{snapshot.thermal_mw:8.1f}", "MWt")
        line("Electric power", f"{snapshot.electric_mw:8.1f}", "MWe")
        line("Reactivity", f"{snapshot.reactivity:+.5f}", "dk/k")
        line("Fuel temp", f"{snapshot.fuel_temp_k:8.1f}", "K")
        line("Coolant temp", f"{snapshot.coolant_temp_k:8.1f}", "K")
        line("Steam temp", f"{snapshot.steam_temp_k:8.1f}", "K")
        line("Void fraction", f"{snapshot.void_fraction:8.3f}")
        line("Pressure", f"{snapshot.pressure_mpa:8.3f}", "MPa")
        line("Channel dP", f"{snapshot.pressure_drop_mpa:8.3f}", "MPa")
        line("Steam inventory", f"{snapshot.steam_inventory:8.3f}")
        line("Feedwater inv.", f"{snapshot.feedwater_inventory:8.3f}")
        line("Condenser temp", f"{snapshot.condenser_temp_k:8.1f}", "K")

        # trends
        screen.blit(font.render("Trends (last 240 samples)", True, (210, 220, 230)), (790, 70))
        draw_trend(screen, 790, 95, 620, 110, list(supervisor.trend_power), (100, 200, 255), 0.0, 1.2)
        screen.blit(font.render("Power fraction", True, (170, 185, 205)), (795, 208))
        draw_trend(screen, 790, 230, 620, 110, list(supervisor.trend_pressure), (255, 180, 100), 6.0 if snapshot.reactor_type != "PWR" else 14.0, 17.5)
        screen.blit(font.render("Pressure MPa", True, (170, 185, 205)), (795, 343))
        draw_trend(screen, 790, 365, 620, 110, list(supervisor.trend_fuel), (255, 120, 120), 500.0, 1500.0)
        screen.blit(font.render("Fuel temperature K", True, (170, 185, 205)), (795, 478))

        # alarm priority / sequence
        y_alarm = 500
        screen.blit(font.render("Alarm sequence (priority)", True, (255, 200, 80)), (420, y_alarm))
        y_alarm += 26
        messages = [f"HI ! {a}" for a in snapshot.unacked_alarms] + snapshot.alarms + snapshot.trips
        if not messages:
            messages = ["None"]
        for m in messages[:12]:
            if m.startswith("HI") or "SCRAM" in m:
                col = (255, 90, 90)
            elif "MED" in m:
                col = (255, 210, 120)
            else:
                col = (205, 210, 220)
            screen.blit(font.render(f"- {m}", True, col), (420, y_alarm))
            y_alarm += 20

        # procedural checklist
        y_proc = 500
        screen.blit(font.render("Reset permissive checklist", True, (180, 220, 180)), (20, y_proc))
        y_proc += 24
        for k, ok in snapshot.checklist.items():
            mark = "[x]" if ok else "[ ]"
            col = (120, 220, 140) if ok else (220, 120, 120)
            screen.blit(font.render(f"{mark} {k}", True, col), (20, y_proc))
            y_proc += 20

        help1 = font.render(
            "W/S rods A/D flow Q/E turbine F/V feedwater H/N pressurizer C=ack L=reset B/G bypass I inhibit SHIFT=fast R=reset ESC=quit",
            True,
            (160, 170, 185),
        )
        screen.blit(help1, (20, 790))

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    main()
