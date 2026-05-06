import csv
import os
import time
from collections import deque

import pygame

from plant_supervisor import PlantSupervisor

TREND_LEN = 600          # fast historian: ~10 s at 60 Hz
SLOW_HIST_LEN = 3600     # slow historian: 60 min at 1 Hz
TIME_WINDOWS = [("10 s", TREND_LEN), ("5 min", 300), ("30 min", 1800), ("60 min", SLOW_HIST_LEN)]
TIME_SPEEDS  = [1, 10, 60, 600]
W, H = 1400, 820

# ── Color palette ───────────────────────────────────────────────────────────────
C_BG       = (  8,  10,  14)
C_PANEL    = ( 12,  15,  22)
C_BORDER   = ( 45,  55,  75)
C_TEXT     = (210, 220, 235)
C_TEXT_DIM = (110, 125, 150)
C_TEXT_HDR = (160, 180, 210)
C_ACCENT   = ( 65, 130, 255)
C_GREEN    = ( 50, 200,  85)
C_YELLOW   = (255, 200,  55)
C_ORANGE   = (255, 140,  40)
C_RED      = (255,  55,  55)
C_SLIDER_BG = ( 28,  36,  52)
C_SLIDER_FG = ( 65, 130, 255)


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


def draw_panel(screen, x, y, w, h, title=None, font=None, border_color=None):
    bc = border_color or C_BORDER
    pygame.draw.rect(screen, C_PANEL, (x, y, w, h), border_radius=8)
    pygame.draw.rect(screen, bc,      (x, y, w, h), 1, border_radius=8)
    if title and font:
        lbl   = font.render(title, True, C_TEXT_HDR)
        sep_y = y + 7 + lbl.get_height() + 5
        screen.blit(lbl, (x + 14, y + 7))
        pygame.draw.line(screen, bc, (x + 8, sep_y), (x + w - 8, sep_y))
        return sep_y + 8
    return y + 10


def draw_led(screen, cx, cy, color, radius=7):
    pygame.draw.circle(screen, (18, 22, 32), (cx, cy), radius + 2)
    pygame.draw.circle(screen, color,        (cx, cy), radius)
    hi = tuple(min(255, c + 90) for c in color)
    pygame.draw.circle(screen, hi, (cx - radius // 3, cy - radius // 3), max(1, radius // 3))


def draw_trend(screen, x, y, w, h, values, color, y_min, y_max, label="", font=None):
    pygame.draw.rect(screen, (15, 18, 28), (x, y, w, h), border_radius=5)
    pygame.draw.rect(screen, C_BORDER,    (x, y, w, h), 1, border_radius=5)
    for i in range(1, 4):
        gy = y + int(h * i / 4)
        pygame.draw.line(screen, (28, 35, 52), (x + 2, gy), (x + w - 2, gy))
    if label and font:
        screen.blit(font.render(label, True, C_TEXT_DIM), (x + 5, y + 4))
        hi_s = f"{y_max:.3g}"
        lo_s = f"{y_min:.3g}"
        screen.blit(font.render(hi_s, True, (65, 80, 100)),
                    (x + w - font.size(hi_s)[0] - 4, y + 4))
        screen.blit(font.render(lo_s, True, (65, 80, 100)),
                    (x + w - font.size(lo_s)[0] - 4, y + h - font.get_height() - 2))
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
    if pts:
        pygame.draw.circle(screen, color, pts[-1], 3)


def main():
    pygame.init()
    screen = pygame.display.set_mode((W, H))
    pygame.display.set_caption("Nuclear Reactor Control Room")

    font_sm = pygame.font.SysFont("Menlo", 13)
    font    = pygame.font.SysFont("Menlo", 15)
    font_md = pygame.font.SysFont("Menlo", 17, bold=True)
    font_lg = pygame.font.SysFont("Menlo", 22, bold=True)
    font_xl = pygame.font.SysFont("Menlo", 32, bold=True)

    clock = pygame.time.Clock()

    supervisor = PlantSupervisor("PWR")
    controls   = supervisor.controls

    hist_power = deque([0.0]   * TREND_LEN, maxlen=TREND_LEN)
    hist_fuel  = deque([600.0] * TREND_LEN, maxlen=TREND_LEN)
    hist_react = deque([0.0]   * TREND_LEN, maxlen=TREND_LEN)
    hist_decay = deque([0.0]   * TREND_LEN, maxlen=TREND_LEN)

    slow_power = deque([0.0]   * SLOW_HIST_LEN, maxlen=SLOW_HIST_LEN)
    slow_fuel  = deque([600.0] * SLOW_HIST_LEN, maxlen=SLOW_HIST_LEN)
    slow_react = deque([0.0]   * SLOW_HIST_LEN, maxlen=SLOW_HIST_LEN)
    slow_decay = deque([0.0]   * SLOW_HIST_LEN, maxlen=SLOW_HIST_LEN)
    _slow_acc  = 0.0          # accumulate real time for 1-Hz slow sampling

    time_window_idx = 0       # index into TIME_WINDOWS
    time_speed_idx  = 0       # index into TIME_SPEEDS
    time_speed      = 1

    scram_msg       = ""      # feedback message for SCRAM reset
    scram_msg_timer = 0.0     # seconds to display it

    ic_path = os.path.join(os.path.dirname(__file__), "ic_quicksave.pkl")
    ic_msg       = ""
    ic_msg_timer = 0.0

    rod   = Slider("Rod position",       0.85)
    flow  = Slider("Primary flow",       1.0)
    valve = Slider("Turbine valve",      1.0)
    feed  = Slider("Feedwater valve",    0.7)
    press = Slider("Pressurizer heater", 0.5)
    bor   = Slider("Boration rate",      0.0)
    dil   = Slider("Dilution rate",      0.0)

    last = time.time()

    def select(rt: str):
        supervisor.select(rt)
        nonlocal controls
        controls   = supervisor.controls
        flow.label = ("Primary flow" if rt == "PWR"
                      else "Recirc pump" if rt == "BWR"
                      else "Pump speed")
        for h, v in [(hist_power, 0.0), (hist_fuel, 600.0), (hist_react, 0.0), (hist_decay, 0.0)]:
            h.clear(); h.extend([v] * TREND_LEN)
        for h, v in [(slow_power, 0.0), (slow_fuel, 600.0), (slow_react, 0.0), (slow_decay, 0.0)]:
            h.clear(); h.extend([v] * SLOW_HIST_LEN)
        rod.set(0.85); flow.set(1.0); valve.set(1.0); feed.set(0.7); press.set(0.5)
        bor.set(0.0); dil.set(0.0)

    running  = True
    snapshot = supervisor.step(0.0)

    while running:
        # ── Events ──────────────────────────────────────────────────────────
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                running = False
            elif ev.type == pygame.KEYDOWN:
                mods = pygame.key.get_mods()
                ctrl_held = bool(mods & (pygame.KMOD_LCTRL | pygame.KMOD_RCTRL))
                if ctrl_held and ev.key == pygame.K_s:
                    try:
                        supervisor.save_ic(ic_path)
                        ic_msg = "IC saved"; ic_msg_timer = 3.0
                    except Exception as e:
                        ic_msg = f"Save failed: {e}"; ic_msg_timer = 3.0
                elif ctrl_held and ev.key == pygame.K_l:
                    try:
                        supervisor.load_ic(ic_path)
                        controls = supervisor.controls
                        ic_msg = "IC loaded"; ic_msg_timer = 3.0
                    except Exception as e:
                        ic_msg = f"Load failed: {e}"; ic_msg_timer = 3.0
                elif ev.key == pygame.K_ESCAPE: running = False
                elif ev.key == pygame.K_1:      select("PWR")
                elif ev.key == pygame.K_2:      select("BWR")
                elif ev.key == pygame.K_3:      select("RBMK")
                elif ev.key == pygame.K_SPACE:  controls.scram = True
                elif ev.key == pygame.K_r:      select(supervisor.reactor_type)
                elif ev.key == pygame.K_t:      controls.turbine_trip          = not controls.turbine_trip
                elif ev.key == pygame.K_p:      controls.startup_permit        = not controls.startup_permit
                elif ev.key == pygame.K_z:      controls.fault_pump_degraded   = not controls.fault_pump_degraded
                elif ev.key == pygame.K_x:      controls.fault_feedwater_loss  = not controls.fault_feedwater_loss
                elif ev.key == pygame.K_l:
                    scram_msg = supervisor.reset_scram()
                    scram_msg_timer = 4.0
                elif ev.key == pygame.K_c:
                    supervisor.acknowledge_alarms()
                elif ev.key == pygame.K_LEFTBRACKET:
                    time_window_idx = (time_window_idx - 1) % len(TIME_WINDOWS)
                elif ev.key == pygame.K_RIGHTBRACKET:
                    time_window_idx = (time_window_idx + 1) % len(TIME_WINDOWS)
                elif ev.key == pygame.K_EQUALS or ev.key == pygame.K_PLUS:
                    time_speed_idx = min(time_speed_idx + 1, len(TIME_SPEEDS) - 1)
                    time_speed = TIME_SPEEDS[time_speed_idx]
                elif ev.key == pygame.K_MINUS:
                    time_speed_idx = max(time_speed_idx - 1, 0)
                    time_speed = TIME_SPEEDS[time_speed_idx]

        # ── Continuous keys ──────────────────────────────────────────────────
        keys = pygame.key.get_pressed()
        step = 0.02 if (keys[pygame.K_LSHIFT] or keys[pygame.K_RSHIFT]) else 0.004
        if keys[pygame.K_w]: rod.set(rod.v   - step)
        if keys[pygame.K_s]: rod.set(rod.v   + step)
        if keys[pygame.K_a]: flow.set(flow.v  - step)
        if keys[pygame.K_d]: flow.set(flow.v  + step)
        if keys[pygame.K_q]: valve.set(valve.v - step)
        if keys[pygame.K_e]: valve.set(valve.v + step)
        if keys[pygame.K_f]: feed.set(feed.v  + step)
        if keys[pygame.K_v]: feed.set(feed.v  - step)
        if keys[pygame.K_h]: press.set(press.v + step)
        if keys[pygame.K_n]: press.set(press.v - step)
        if keys[pygame.K_b]: bor.set(bor.v + step)
        if keys[pygame.K_g]: bor.set(bor.v - step)
        if keys[pygame.K_u]: dil.set(dil.v + step)
        if keys[pygame.K_i]: dil.set(dil.v - step)

        controls.rod_position       = rod.v
        controls.flow               = flow.v
        controls.turbine_valve      = valve.v
        controls.feedwater_valve    = feed.v
        controls.pressurizer_heater = press.v
        controls.boration_rate      = bor.v
        controls.dilution_rate      = dil.v

        # ── Physics step ─────────────────────────────────────────────────────
        now  = time.time()
        real_dt = clamp(now - last, 0.0, 0.2)
        last = now

        # Disable time acceleration when alarms are active
        effective_speed = time_speed if not snapshot.alarms and not snapshot.trips else 1
        dt = real_dt * effective_speed

        snapshot = supervisor.step(dt)
        hist_power.append(snapshot.power_fraction * 100.0)
        hist_fuel.append(snapshot.fuel_temp_k)
        hist_react.append(snapshot.reactivity)
        hist_decay.append(snapshot.decay_heat_fraction * 100.0)

        # Slow historian: sample at 1 Hz real time
        _slow_acc += real_dt
        if _slow_acc >= 1.0:
            _slow_acc -= 1.0
            slow_power.append(snapshot.power_fraction * 100.0)
            slow_fuel.append(snapshot.fuel_temp_k)
            slow_react.append(snapshot.reactivity)
            slow_decay.append(snapshot.decay_heat_fraction * 100.0)

        # Countdown timers for on-screen messages
        scram_msg_timer = max(0.0, scram_msg_timer - real_dt)
        ic_msg_timer    = max(0.0, ic_msg_timer    - real_dt)

        # ── Blink timing ─────────────────────────────────────────────────────
        t_ms       = pygame.time.get_ticks()
        blink_fast = (t_ms // 300) % 2 == 0
        blink_slow = (t_ms // 600) % 2 == 0
        has_alarms = bool(snapshot.alarms)
        has_trips  = bool(snapshot.trips)
        is_scram   = controls.scram

        # ── Background ───────────────────────────────────────────────────────
        screen.fill(C_BG)

        # ── Header ───────────────────────────────────────────────────────────
        HEADER_H = 52
        PAD      = 8

        if is_scram:
            scram_bg = (180, 15, 15) if blink_fast else (90, 8, 8)
            pygame.draw.rect(screen, scram_bg, (0, 0, W, HEADER_H))
            txt = font_xl.render("REACTOR SCRAM  --  REACTOR SCRAM  --  REACTOR SCRAM", True, (255, 255, 255))
            screen.blit(txt, (W // 2 - txt.get_width() // 2, HEADER_H // 2 - txt.get_height() // 2))
        else:
            pygame.draw.rect(screen, (12, 16, 26), (0, 0, W, HEADER_H))
            pygame.draw.line(screen, C_BORDER, (0, HEADER_H - 1), (W, HEADER_H - 1))

            # Reactor type tabs
            for i, (rt, base_col) in enumerate([
                ("PWR",  (50, 110, 210)),
                ("BWR",  (40, 160, 110)),
                ("RBMK", (200, 75,  50)),
            ]):
                active = snapshot.reactor_type == rt
                bx     = 16 + i * 115
                bg     = base_col if active else (25, 32, 46)
                pygame.draw.rect(screen, bg, (bx, 8, 105, 36), border_radius=6)
                if active:
                    bc2 = tuple(min(255, c + 70) for c in base_col)
                    pygame.draw.rect(screen, bc2, (bx, 8, 105, 36), 2, border_radius=6)
                tc = (255, 255, 255) if active else (80, 100, 130)
                t  = font_lg.render(rt, True, tc)
                screen.blit(t, (bx + 52 - t.get_width() // 2, 26 - t.get_height() // 2))

            # Time and speed (centered)
            speed_str = f"  {effective_speed}×" if effective_speed > 1 else ""
            t_surf = font_lg.render(f"T = {snapshot.time:8.1f} s{speed_str}", True, C_TEXT)
            screen.blit(t_surf, (W // 2 - t_surf.get_width() // 2, 14))

            # Alarm / trip badge (top-right)
            if has_trips:
                col = (220, 40, 40) if blink_fast else (100, 20, 20)
                pygame.draw.rect(screen, col, (W - 215, 8, 200, 36), border_radius=6)
                t = font_lg.render("TRIP ACTIVE", True, (255, 255, 255))
                screen.blit(t, (W - 215 + 100 - t.get_width() // 2, 26 - t.get_height() // 2))
            elif has_alarms:
                col = C_ORANGE if blink_slow else (100, 55, 0)
                pygame.draw.rect(screen, col, (W - 215, 8, 200, 36), border_radius=6)
                t = font_lg.render("ALARM", True, (255, 255, 255))
                screen.blit(t, (W - 215 + 100 - t.get_width() // 2, 26 - t.get_height() // 2))

        # ── Layout constants ─────────────────────────────────────────────────
        Y0     = HEADER_H + PAD
        BODY_H = H - HEADER_H - 28 - PAD    # 28 px reserved for help bar

        LP_X, LP_W = PAD,                    360
        CP_X, CP_W = LP_X + LP_W + PAD,     380
        RP_X, RP_W = CP_X + CP_W + PAD,     W - (CP_X + CP_W + PAD) - PAD

        # ── LEFT PANEL: Controls ─────────────────────────────────────────────
        cy = draw_panel(screen, LP_X, Y0, LP_W, BODY_H, "CONTROLS", font_sm)

        for sl in (rod, flow, valve, feed, press):
            lbl_s = font.render(sl.label, True, C_TEXT_DIM)
            val_s = font_md.render(f"{sl.v:.3f}", True, C_ACCENT)
            screen.blit(lbl_s, (LP_X + 14, cy))
            screen.blit(val_s, (LP_X + LP_W - val_s.get_width() - 14, cy))
            cy += max(lbl_s.get_height(), val_s.get_height()) + 4

            bx, bw, bh = LP_X + 14, LP_W - 28, 10
            pygame.draw.rect(screen, C_SLIDER_BG, (bx, cy, bw, bh), border_radius=4)
            fill = int(bw * (sl.v - sl.lo) / max(1e-6, sl.hi - sl.lo))
            if fill > 0:
                pygame.draw.rect(screen, C_SLIDER_FG, (bx, cy, fill, bh), border_radius=4)
            pygame.draw.rect(screen, (150, 190, 255), (bx + fill - 2, cy, 3, bh))
            cy += bh + 14

        pygame.draw.line(screen, C_BORDER, (LP_X + 10, cy), (LP_X + LP_W - 10, cy))
        cy += 12

        # Status LEDs
        screen.blit(font_sm.render("STATUS                    KEY", True, C_TEXT_DIM),
                    (LP_X + 14, cy))
        cy += 18

        for label, state, col_on, key in [
            ("Startup permit",  controls.startup_permit,       C_GREEN,  "P"),
            ("Turbine trip",    controls.turbine_trip,         C_ORANGE, "T"),
            ("Pump degraded",   controls.fault_pump_degraded,  C_RED,    "Z"),
            ("Feedwater fault", controls.fault_feedwater_loss, C_RED,    "X"),
        ]:
            led_col = col_on if state else (30, 38, 52)
            draw_led(screen, LP_X + 26, cy + 8, led_col)
            tc = col_on if state else C_TEXT_DIM
            screen.blit(font.render(label, True, tc), (LP_X + 44, cy))
            kt = font_sm.render(f"[{key}]", True, C_TEXT_DIM)
            screen.blit(kt, (LP_X + LP_W - kt.get_width() - 14, cy + 2))
            cy += 26

        pygame.draw.line(screen, C_BORDER, (LP_X + 10, cy + 4), (LP_X + LP_W - 10, cy + 4))
        cy += 16

        # SCRAM status
        if is_scram:
            scram_col = C_RED if blink_fast else C_ORANGE
            draw_led(screen, LP_X + 26, cy + 12, scram_col, radius=10)
            screen.blit(font_md.render("SCRAM ACTIVE",                True, scram_col),   (LP_X + 46, cy + 4))
            screen.blit(font_sm.render("L: reset  R: full restart",   True, C_TEXT_DIM), (LP_X + 46, cy + 26))
        else:
            draw_led(screen, LP_X + 26, cy + 12, C_GREEN, radius=10)
            screen.blit(font_md.render("REACTOR NORMAL",          True, C_GREEN),   (LP_X + 46, cy + 4))
            screen.blit(font_sm.render("[SPACE] Emergency SCRAM", True, C_TEXT_DIM), (LP_X + 46, cy + 26))
        cy += 48

        # SCRAM reset / IC feedback messages
        if scram_msg and scram_msg_timer > 0:
            ok = "APPROVED" in scram_msg
            mc = C_GREEN if ok else C_ORANGE
            screen.blit(font_sm.render(scram_msg, True, mc), (LP_X + 14, cy))
            cy += 18
        if ic_msg and ic_msg_timer > 0:
            screen.blit(font_sm.render(ic_msg, True, C_ACCENT), (LP_X + 14, cy))

        # ── CENTER PANEL: Plant readouts ─────────────────────────────────────
        ALARM_H   = 225
        READOUT_H = BODY_H - ALARM_H - PAD

        ry = draw_panel(screen, CP_X, Y0, CP_W, READOUT_H, "PLANT STATUS", font_sm)

        def rline(label, val_str, unit, col):
            nonlocal ry
            screen.blit(font.render(label, True, C_TEXT_DIM), (CP_X + 12, ry))
            vs  = font_md.render(val_str, True, col)
            uw  = font_sm.size(unit)[0]
            screen.blit(vs,  (CP_X + CP_W - uw - vs.get_width() - 18, ry - 1))
            screen.blit(font_sm.render(unit, True, C_TEXT_DIM), (CP_X + CP_W - uw - 10, ry + 2))
            ry += 24

        def rsep():
            nonlocal ry
            pygame.draw.line(screen, C_BORDER,
                             (CP_X + 10, ry + 2), (CP_X + CP_W - 10, ry + 2))
            ry += 12

        pf = snapshot.power_fraction * 100.0
        rline("Core power",    f"{pf:7.2f}",                     "%",
              C_RED if pf > 110 else C_YELLOW if pf > 90 else C_GREEN)
        rline("Thermal power", f"{snapshot.thermal_mw:7.1f}",    "MWt", C_TEXT)
        rline("Electric power",f"{snapshot.electric_mw:7.1f}",   "MWe", C_TEXT)
        dh = snapshot.decay_heat_fraction * 100.0
        rline("Decay heat",    f"{dh:7.2f}",                     "%",
              C_YELLOW if dh > 3.0 else C_TEXT)
        rsep()

        rk = snapshot.reactivity
        rline("Reactivity", f"{rk:+.5f}", "dk/k",
              C_RED if abs(rk) > 0.005 else C_YELLOW if abs(rk) > 0.001 else C_GREEN)
        rsep()

        ft = snapshot.fuel_temp_k
        rline("Fuel temp",    f"{ft:7.1f}", "K",
              C_RED if ft > 1400 else C_ORANGE if ft > 1200 else C_YELLOW if ft > 900 else C_GREEN)
        ct = snapshot.coolant_temp_k
        rline("Coolant temp", f"{ct:7.1f}", "K",
              C_RED if ct > 620 else C_YELLOW if ct > 580 else C_GREEN)
        rline("Steam temp",   f"{snapshot.steam_temp_k:7.1f}",    "K", C_TEXT)
        cond = snapshot.condenser_temp_k
        rline("Condenser",    f"{cond:7.1f}", "K",
              C_YELLOW if cond > 330 else C_TEXT)
        rsep()

        p_nom = 15.5 if snapshot.reactor_type == "PWR" else 7.0
        p     = snapshot.pressure_mpa
        rline("Pressure",      f"{p:7.3f}", "MPa",
              C_RED if p > p_nom * 1.08 else C_YELLOW if p > p_nom * 1.04 else C_GREEN)
        vf = snapshot.void_fraction
        rline("Void frac.",    f"{vf:7.3f}", "",
              C_YELLOW if vf > 0.3 else C_TEXT)
        si = snapshot.steam_inventory
        rline("Steam inv.",    f"{si:7.3f}", "",
              C_YELLOW if si > 1.5 or si < 0.3 else C_TEXT)
        fw = snapshot.feedwater_inventory
        rline("Feedwater inv.",f"{fw:7.3f}", "",
              C_RED if fw < 0.1 else C_YELLOW if fw < 0.25 else C_TEXT)
        if snapshot.reactor_type == "PWR":
            bp = snapshot.boron_ppm
            rline("Boron conc.",  f"{bp:7.1f}", "ppm",
                  C_YELLOW if bp < 100 else C_TEXT)
        rsep()

        rline("Sim time", f"{snapshot.time:7.1f}", "s", C_TEXT_DIM)

        # ── ALARM PANEL ──────────────────────────────────────────────────────
        alarm_y = Y0 + READOUT_H + PAD
        if has_trips:
            alarm_border = C_RED if blink_fast else (160, 30, 30)
        elif has_alarms:
            alarm_border = C_ORANGE
        else:
            alarm_border = C_BORDER

        ay = draw_panel(screen, CP_X, alarm_y, CP_W, ALARM_H,
                        "ALARMS & TRIPS", font_sm, border_color=alarm_border)

        messages = [(m, True) for m in snapshot.trips] + [(m, False) for m in snapshot.alarms]
        if not messages:
            ok = font_md.render("All systems nominal", True, C_GREEN)
            mid_x = CP_X + CP_W // 2
            mid_y = alarm_y + ALARM_H // 2
            draw_led(screen, mid_x - ok.get_width() // 2 - 14, mid_y, C_GREEN, radius=7)
            screen.blit(ok, (mid_x - ok.get_width() // 2, mid_y - ok.get_height() // 2))
        else:
            for m, is_trip in messages[:8]:
                if is_trip:
                    col    = C_RED if blink_fast else (180, 60, 60)
                    prefix = "TRIP"
                else:
                    col    = C_ORANGE
                    prefix = "WARN"
                screen.blit(font.render(f"[{prefix}]  {m}", True, col), (CP_X + 12, ay))
                ay += 22

        # ── RIGHT PANEL: Trends ──────────────────────────────────────────────
        tw_label, tw_len = TIME_WINDOWS[time_window_idx]
        trend_title = f"TRENDS  [{tw_label}]  [/]: window  [+/-]: speed"
        content_y = draw_panel(screen, RP_X, Y0, RP_W, BODY_H, trend_title, font_sm)

        GW  = RP_W - 20
        GX  = RP_X + 10
        GAP = 6
        GH  = (BODY_H - (content_y - Y0) - GAP * 3) // 4

        # Select data source based on window
        if tw_len <= TREND_LEN:
            p_data, f_data, r_data, d_data = hist_power, hist_fuel, hist_react, hist_decay
        else:
            n_slow = min(tw_len, SLOW_HIST_LEN)
            p_data = list(slow_power)[-n_slow:]
            f_data = list(slow_fuel)[-n_slow:]
            r_data = list(slow_react)[-n_slow:]
            d_data = list(slow_decay)[-n_slow:]

        for i, (hist, color, y_min, y_max, label) in enumerate([
            (p_data, C_GREEN,  0.0,   130.0, "Core power (%)"),
            (f_data, C_ORANGE, 500.0, 1600.0,"Fuel temp (K)"),
            (r_data, C_ACCENT, -0.01, 0.01,  "Reactivity (dk/k)"),
            (d_data, C_YELLOW, 0.0,   10.0,  "Decay heat (%)"),
        ]):
            gy = content_y + i * (GH + GAP)
            draw_trend(screen, GX, gy, GW, GH, list(hist),
                       color, y_min, y_max, label, font_sm)

        # ── Help bar ─────────────────────────────────────────────────────────
        help_y = H - 24
        pygame.draw.line(screen, C_BORDER, (0, help_y - 4), (W, help_y - 4))
        screen.blit(font_sm.render(
            "W/S:rods  A/D:flow  Q/E:turbine  F/V:feed  H/N:pzr  B/G:borate  SHIFT:fast  "
            "1/2/3:reactor  SPACE:SCRAM  L:reset  C:ack  [/]:window  +/-:speed  Ctrl+S:save  ESC:quit",
            True, C_TEXT_DIM,
        ), (10, help_y))

        pygame.display.flip()
        clock.tick(60)

    # Export event log on exit
    try:
        log_path = os.path.join(os.path.dirname(__file__), "event_log.csv")
        with open(log_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sim_time", "category", "message"])
            for row in supervisor.event_log:
                w.writerow(row)
    except Exception:
        pass
    pygame.quit()


if __name__ == "__main__":
    main()
