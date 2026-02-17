import time
import pygame

from physics_engine import PWRPlant, ControlInputs as PWRControls
from physics_engine_bwr import BWRPlant, BWRControlInputs
from physics_engine_rbm import RBMKPlant, RBMKControlInputs


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
    screen = pygame.display.set_mode((1100, 680))
    pygame.display.set_caption("Nuclear Reactor Control Room (MVP)")
    font = pygame.font.SysFont("Menlo", 18)
    big = pygame.font.SysFont("Menlo", 28, bold=True)
    clock = pygame.time.Clock()

    reactor_type = "PWR"
    plant = PWRPlant()
    controls = PWRControls()

    rod = Slider("Rod position", 0.0)
    flow = Slider("Primary flow", 1.0)
    valve = Slider("Turbine valve", 1.0)
    scrammed = False

    last = time.time()
    snapshot = None

    def select(rt: str):
        nonlocal reactor_type, plant, controls, scrammed
        reactor_type = rt
        scrammed = False
        if rt == "PWR":
            plant = PWRPlant()
            controls = PWRControls()
            flow.label = "Primary flow"
        elif rt == "BWR":
            plant = BWRPlant()
            controls = BWRControlInputs()
            flow.label = "Recirc pump"
        else:
            plant = RBMKPlant()
            controls = RBMKControlInputs()
            flow.label = "Pump speed"

    def apply_controls():
        nonlocal controls
        if reactor_type == "PWR":
            controls.rod_position = rod.v
            controls.primary_flow = flow.v
            controls.turbine_valve = valve.v
            controls.scram = scrammed
        elif reactor_type == "BWR":
            controls.rod_position = rod.v
            controls.recirc_pump = flow.v
            controls.turbine_valve = valve.v
            controls.scram = scrammed
        else:
            controls.rod_position = rod.v
            controls.pump_speed = flow.v
            controls.turbine_valve = valve.v
            controls.scram = scrammed

    select("PWR")

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
                    scrammed = True
                elif ev.key == pygame.K_r:
                    select(reactor_type)

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

        now = time.time()
        dt = clamp(now - last, 0.0, 0.2)
        last = now

        apply_controls()
        snapshot = plant.step(dt, controls)

        screen.fill((10, 12, 16))

        header = big.render(
            f"Reactor: {reactor_type}   (1=PWR, 2=BWR, 3=RBMK)   SPACE=SCRAM   R=Reset",
            True,
            (220, 230, 240),
        )
        screen.blit(header, (20, 18))

        y = 90
        for sl in (rod, flow, valve):
            txt = font.render(f"{sl.label}: {sl.v:.3f}", True, (200, 210, 220))
            screen.blit(txt, (20, y))
            bx, by, bw, bh = 20, y + 24, 320, 16
            pygame.draw.rect(screen, (35, 40, 52), (bx, by, bw, bh), border_radius=6)
            fill = int(bw * (sl.v - sl.lo) / (sl.hi - sl.lo))
            pygame.draw.rect(screen, (80, 160, 255), (bx, by, fill, bh), border_radius=6)
            y += 64

        scr = font.render(
            f"SCRAM: {'YES' if scrammed else 'no'}",
            True,
            (255, 80, 80) if scrammed else (200, 210, 220),
        )
        screen.blit(scr, (20, y))

        x0 = 380
        y0 = 90

        def line(label, value, unit=""):
            nonlocal y0
            t = font.render(f"{label:<18} {value} {unit}", True, (220, 230, 240))
            screen.blit(t, (x0, y0))
            y0 += 26

        if reactor_type == "PWR":
            line("Time", f"{snapshot.time:8.1f}", "s")
            line("Power", f"{snapshot.power_fraction*100:8.2f}", "%")
            line("Thermal", f"{snapshot.thermal_power/1e6:8.1f}", "MWt")
            line("Electric", f"{snapshot.electric_power/1e6:8.1f}", "MWe")
            line("Reactivity", f"{snapshot.reactivity:+.5f}", "dk/k")
            line("Xenon", f"{snapshot.xenon_inventory:8.4f}")
            line("Fuel T", f"{snapshot.fuel_temperature:8.1f}", "K")
            line("Coolant T", f"{snapshot.coolant_temperature:8.1f}", "K")
            line("SG T", f"{snapshot.sg_temperature:8.1f}", "K")
        elif reactor_type == "BWR":
            line("Time", f"{snapshot.time:8.1f}", "s")
            line("Power", f"{snapshot.power_fraction*100:8.2f}", "%")
            line("Thermal", f"{snapshot.power_watts/1e6:8.1f}", "MWt")
            line("Electric", f"{snapshot.electric_power/1e6:8.1f}", "MWe")
            line("Reactivity", f"{snapshot.reactivity:+.5f}", "dk/k")
            line("Xenon", f"{snapshot.xenon_inventory:8.4f}")
            line("Fuel T", f"{snapshot.fuel_temperature:8.1f}", "K")
            line("Coolant T", f"{snapshot.coolant_temperature:8.1f}", "K")
            line("Steam T", f"{snapshot.steam_temperature:8.1f}", "K")
            line("Void", f"{snapshot.void_fraction:8.3f}")
        else:
            line("Time", f"{snapshot.time:8.1f}", "s")
            line("Power", f"{snapshot.power_fraction*100:8.2f}", "%")
            line("Thermal", f"{snapshot.power_watts/1e6:8.1f}", "MWt")
            line("Electric", f"{snapshot.electric_power/1e6:8.1f}", "MWe")
            line("Reactivity", f"{snapshot.reactivity:+.5f}", "dk/k")
            line("Xenon", f"{snapshot.xenon_inventory:8.4f}")
            line("Fuel T", f"{snapshot.fuel_temperature:8.1f}", "K")
            line("Graphite", f"{snapshot.graphite_temperature:8.1f}", "K")
            line("Coolant T", f"{snapshot.coolant_temperature:8.1f}", "K")
            line("Steam T", f"{snapshot.steam_temperature:8.1f}", "K")
            line("Void", f"{snapshot.void_fraction:8.3f}")

        help1 = font.render(
            "Controls: W/S rods, A/D flow, Q/E valve, SHIFT=fast, SPACE=SCRAM, R=reset, ESC=quit",
            True,
            (160, 170, 185),
        )
        screen.blit(help1, (20, 640))

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    main()
