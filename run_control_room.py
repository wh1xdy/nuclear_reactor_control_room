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


def main():
    pygame.init()
    screen = pygame.display.set_mode((1280, 760))
    pygame.display.set_caption("Nuclear Reactor Control Room")
    font = pygame.font.SysFont("Menlo", 18)
    big = pygame.font.SysFont("Menlo", 26, bold=True)
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

        y = 80
        for sl in (rod, flow, valve, feed, press):
            txt = font.render(f"{sl.label}: {sl.v:.3f}", True, (200, 210, 220))
            screen.blit(txt, (20, y))
            bx, by, bw, bh = 20, y + 24, 340, 16
            pygame.draw.rect(screen, (35, 40, 52), (bx, by, bw, bh), border_radius=6)
            fill = int(bw * (sl.v - sl.lo) / (sl.hi - sl.lo))
            pygame.draw.rect(screen, (80, 160, 255), (bx, by, fill, bh), border_radius=6)
            y += 62

        status_lines = [
            f"Startup permit: {'ON' if controls.startup_permit else 'OFF'}",
            f"Turbine trip: {'ON' if controls.turbine_trip else 'OFF'}",
            f"Pump fault (Z): {'ON' if controls.fault_pump_degraded else 'OFF'}",
            f"Feedwater fault (X): {'ON' if controls.fault_feedwater_loss else 'OFF'}",
        ]
        for s in status_lines:
            screen.blit(font.render(s, True, (185, 200, 215)), (20, y))
            y += 24

        x0, y0 = 420, 80

        def line(label, value, unit=""):
            nonlocal y0
            t = font.render(f"{label:<22} {value} {unit}", True, (220, 230, 240))
            screen.blit(t, (x0, y0))
            y0 += 24

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
        line("Steam inventory", f"{snapshot.steam_inventory:8.3f}")
        line("Feedwater inv.", f"{snapshot.feedwater_inventory:8.3f}")
        line("Condenser temp", f"{snapshot.condenser_temp_k:8.1f}", "K")

        y_alarm = 420
        alarm_title = font.render("Alarms / Trips", True, (255, 200, 80))
        screen.blit(alarm_title, (420, y_alarm))
        y_alarm += 28
        messages = snapshot.alarms + snapshot.trips
        if not messages:
            messages = ["None"]
        for m in messages[:10]:
            col = (255, 90, 90) if "SCRAM" in m or "trip" in m.lower() else (255, 220, 120)
            screen.blit(font.render(f"- {m}", True, col), (420, y_alarm))
            y_alarm += 22

        help1 = font.render(
            "W/S rods A/D flow Q/E turbine F/V feedwater H/N pressurizer  SHIFT=fast  R=reset ESC=quit",
            True,
            (160, 170, 185),
        )
        screen.blit(help1, (20, 724))

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    main()
