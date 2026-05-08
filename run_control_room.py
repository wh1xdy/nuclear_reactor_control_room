import csv
import math
import os
import time
from collections import deque

import pygame

from plant_supervisor import PlantSupervisor

TREND_LEN     = 600       # fast historian: ~10 s at 60 Hz
SLOW_HIST_LEN = 3600      # slow historian: 60 min at 1 Hz
TIME_WINDOWS  = [("10 s", TREND_LEN), ("5 min", 300), ("30 min", 1800), ("60 min", SLOW_HIST_LEN)]
TIME_SPEEDS   = [1, 10, 60, 600]
W, H          = 1920, 1080

# ── Color palette ────────────────────────────────────────────────────────────
C_BG        = (  8,  10,  14)
C_PANEL     = ( 12,  15,  22)
C_BORDER    = ( 45,  55,  75)
C_TEXT      = (210, 220, 235)
C_TEXT_DIM  = (110, 125, 150)
C_TEXT_HDR  = (160, 180, 210)
C_ACCENT    = ( 65, 130, 255)
C_GREEN     = ( 50, 200,  85)
C_YELLOW    = (255, 200,  55)
C_ORANGE    = (255, 140,  40)
C_RED       = (255,  55,  55)
C_SLIDER_BG = ( 28,  36,  52)
C_SLIDER_FG = ( 65, 130, 255)
C_BLUE_PIPE = ( 40, 100, 200)
C_CYAN_PIPE = ( 50, 180, 210)
C_GRAY_PIPE = (140, 150, 165)


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


# ── Generic drawing helpers ──────────────────────────────────────────────────

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
        for val_s, vy in [(f"{y_max:.3g}", y + 4), (f"{y_min:.3g}", y + h - font.get_height() - 2)]:
            screen.blit(font.render(val_s, True, (65, 80, 100)),
                        (x + w - font.size(val_s)[0] - 4, vy))
    if len(values) < 2:
        return
    n = len(values)
    pts = []
    for i, v in enumerate(values):
        px = x + int(i * (w - 4) / max(1, n - 1)) + 2
        vv = clamp((v - y_min) / max(1e-6, y_max - y_min), 0.0, 1.0)
        pts.append((px, y + h - 2 - int(vv * (h - 4))))
    pygame.draw.lines(screen, color, False, pts, 2)
    if pts:
        pygame.draw.circle(screen, color, pts[-1], 3)


def draw_bar_meter(screen, x, y, w, h, value, lo, hi, font,
                   label="", unit="", bands=None, trip_lo=None, trip_hi=None):
    """Vertical bar meter with optional color bands and trip lines."""
    pygame.draw.rect(screen, (10, 13, 20), (x, y, w, h), border_radius=4)
    pygame.draw.rect(screen, C_BORDER,    (x, y, w, h), 1, border_radius=4)
    frac = clamp((value - lo) / max(1e-9, hi - lo), 0.0, 1.0)
    bar_h = int(frac * (h - 4))
    bar_y = y + h - 2 - bar_h
    if bands:
        col = C_GREEN
        for band_lo, band_hi, band_col in bands:
            if value >= band_lo:
                col = band_col
    else:
        col = C_GREEN
    if bar_h > 0:
        pygame.draw.rect(screen, col, (x + 2, bar_y, w - 4, bar_h), border_radius=2)
    if trip_hi is not None:
        ty = y + h - 2 - int(clamp((trip_hi - lo) / max(1e-9, hi - lo), 0, 1) * (h - 4))
        pygame.draw.line(screen, C_RED, (x, ty), (x + w, ty), 2)
    if trip_lo is not None:
        ty = y + h - 2 - int(clamp((trip_lo - lo) / max(1e-9, hi - lo), 0, 1) * (h - 4))
        pygame.draw.line(screen, C_ORANGE, (x, ty), (x + w, ty), 2)
    if label and font:
        lbl = font.render(label, True, C_TEXT_DIM)
        screen.blit(lbl, (x + w // 2 - lbl.get_width() // 2, y - lbl.get_height() - 2))
    if font:
        val_s = f"{value:.1f}" if abs(value) < 1000 else f"{value:.0f}"
        vs = font.render(val_s, True, C_TEXT)
        screen.blit(vs, (x + w // 2 - vs.get_width() // 2, y + h + 2))
        if unit:
            us = font.render(unit, True, C_TEXT_DIM)
            screen.blit(us, (x + w // 2 - us.get_width() // 2, y + h + 2 + vs.get_height()))


def draw_dcs_box(screen, x, y, w, h, label, value_str, col, font_lbl, font_val):
    """DCS-style digital readout box."""
    pygame.draw.rect(screen, (6, 9, 14), (x, y, w, h), border_radius=4)
    pygame.draw.rect(screen, C_BORDER,   (x, y, w, h), 1, border_radius=4)
    lbl = font_lbl.render(label, True, C_TEXT_DIM)
    val = font_val.render(value_str, True, col)
    screen.blit(lbl, (x + 6, y + 4))
    screen.blit(val, (x + w - val.get_width() - 6, y + h - val.get_height() - 4))


# ── P&ID schematics ──────────────────────────────────────────────────────────

def _pipe(screen, pts, color, w=3):
    if len(pts) >= 2:
        pygame.draw.lines(screen, color, False, pts, w)


def _box(screen, x, y, w, h, label, font, col=C_BORDER, bg=(18, 24, 36), alarm=False, blink=False):
    bg_use = (60, 10, 10) if alarm and blink else (30, 5, 5) if alarm else bg
    pygame.draw.rect(screen, bg_use, (x, y, w, h), border_radius=5)
    pygame.draw.rect(screen, col if not alarm else C_RED, (x, y, w, h), 2, border_radius=5)
    if label and font:
        lbl = font.render(label, True, C_TEXT if not alarm else C_RED)
        screen.blit(lbl, (x + w // 2 - lbl.get_width() // 2, y + h // 2 - lbl.get_height() // 2))


def draw_pid_pwr(screen, x, y, w, h, snap, font, blink):
    """Simplified PWR P&ID schematic."""
    mx, my = x + w // 2, y + h // 2
    fuel_alarm = snap.fuel_temp_k > 1200
    press_alarm = snap.pressure_mpa > 16.5
    frac = clamp(snap.power_fraction, 0, 1)
    flow_col = C_BLUE_PIPE if snap.coolant_temp_k < 580 else C_CYAN_PIPE

    # Reactor vessel
    rv_x, rv_y, rv_w, rv_h = x + 80, y + h // 2 - 60, 90, 120
    _box(screen, rv_x, rv_y, rv_w, rv_h, "REACTOR", font, alarm=fuel_alarm, blink=blink)
    # Glow inside vessel proportional to power
    glow_h = int(rv_h * frac * 0.8)
    glow_col = (int(50 + 180 * frac), int(80 * (1 - frac)), 0)
    if glow_h > 0:
        pygame.draw.rect(screen, glow_col, (rv_x + 6, rv_y + rv_h - glow_h - 6, rv_w - 12, glow_h), border_radius=3)

    # Steam generator
    sg_x, sg_y, sg_w, sg_h = x + 300, y + 60, 70, 130
    _box(screen, sg_x, sg_y, sg_w, sg_h, "SG", font, col=(60, 100, 180))

    # Pressurizer
    pz_x, pz_y, pz_w, pz_h = x + 80, y + 30, 70, 60
    _box(screen, pz_x, pz_y, pz_w, pz_h, "PZR", font,
         col=C_RED if press_alarm else (80, 60, 160), alarm=press_alarm, blink=blink)

    # Turbine / condenser / feedwater pump
    tb_x, tb_y, tb_w, tb_h = x + 460, y + 60, 70, 60
    _box(screen, tb_x, tb_y, tb_w, tb_h, "TURB", font, col=(60, 140, 80))

    cd_x, cd_y, cd_w, cd_h = x + 460, y + 180, 70, 55
    _box(screen, cd_x, cd_y, cd_w, cd_h, "COND", font, col=(50, 90, 120))

    fw_x, fw_y, fw_w, fw_h = x + 300, y + 240, 70, 50
    _box(screen, fw_x, fw_y, fw_w, fw_h, "FW PMP", font, col=(60, 100, 160))

    # RCP
    rcp_x, rcp_y, rcp_w, rcp_h = x + 80, y + h - 80, 70, 50
    _box(screen, rcp_x, rcp_y, rcp_w, rcp_h, "RCP", font, col=(70, 110, 180))

    # Primary loop pipes (hot leg: vessel top-right → SG left)
    hot_col  = (200, 80, 40) if snap.coolant_temp_k > 560 else flow_col
    cold_col = C_BLUE_PIPE
    _pipe(screen, [(rv_x + rv_w, rv_y + 30), (sg_x, sg_y + 35)], hot_col, 4)
    # Cold leg: SG bottom → RCP → vessel bottom
    _pipe(screen, [(sg_x, sg_y + sg_h - 30), (sg_x - 20, sg_y + sg_h + 20),
                   (rcp_x + rcp_w, rcp_y + 25), (rv_x + rv_w, rv_y + rv_h - 30)], cold_col, 4)
    # Pressurizer surge line
    _pipe(screen, [(pz_x + pz_w // 2, pz_y + pz_h), (rv_x + rv_w // 2, rv_y)], (130, 80, 200), 2)

    # Secondary loop (steam side)
    steam_col = C_GRAY_PIPE if snap.turbine_valve < 0.1 else C_CYAN_PIPE
    _pipe(screen, [(sg_x + sg_w, sg_y + 20), (tb_x, tb_y + 30)], steam_col, 3)
    # Turbine → condenser
    _pipe(screen, [(tb_x + tb_w // 2, tb_y + tb_h), (tb_x + tb_w // 2, cd_y)], C_GRAY_PIPE, 3)
    # Condenser → FW pump → SG
    _pipe(screen, [(cd_x, cd_y + cd_h // 2), (fw_x + fw_w, fw_y + fw_h // 2),
                   (fw_x, fw_y + fw_h // 2), (sg_x + sg_w // 2, sg_y + sg_h)], C_BLUE_PIPE, 3)

    # Live value overlays
    def small_val(sx, sy, txt, col=C_TEXT):
        lbl = font.render(txt, True, col)
        screen.blit(lbl, (sx, sy))

    small_val(rv_x + rv_w + 4,  rv_y + 10,   f"{snap.power_fraction*100:.0f}%",
              C_RED if frac > 1.1 else C_GREEN)
    small_val(rv_x + rv_w + 4,  rv_y + 28,   f"{snap.fuel_temp_k:.0f}K",
              C_RED if snap.fuel_temp_k > 1200 else C_TEXT)
    small_val(pz_x + pz_w + 4,  pz_y + 10,   f"{snap.pressure_mpa:.2f}MPa",
              C_RED if press_alarm else C_TEXT)
    small_val(sg_x + sg_w + 4,  sg_y + 10,   f"{snap.steam_temp_k:.0f}K", C_TEXT)
    small_val(cd_x + cd_w + 4,  cd_y + 10,   f"{snap.condenser_temp_k:.0f}K", C_TEXT)


def draw_pid_bwr(screen, x, y, w, h, snap, font, blink):
    """Simplified BWR P&ID schematic."""
    frac = clamp(snap.power_fraction, 0, 1)
    fuel_alarm = snap.fuel_temp_k > 1200

    # Reactor pressure vessel (tall, contains steam space)
    rv_x, rv_y, rv_w, rv_h = x + 70, y + 40, 100, 200
    _box(screen, rv_x, rv_y, rv_w, rv_h, "", font, alarm=fuel_alarm, blink=blink)
    font.render("REACTOR", True, C_TEXT)
    lbl = font.render("REACTOR", True, C_TEXT)
    screen.blit(lbl, (rv_x + rv_w // 2 - lbl.get_width() // 2, rv_y + 6))
    # Two-phase region
    two_phase_h = int(rv_h * 0.4 * frac)
    pygame.draw.rect(screen, C_CYAN_PIPE, (rv_x + 4, rv_y + rv_h - two_phase_h - 4, rv_w - 8, two_phase_h), border_radius=3)

    # Steam line to turbine
    tb_x, tb_y, tb_w, tb_h = x + 290, y + 40, 70, 60
    _box(screen, tb_x, tb_y, tb_w, tb_h, "TURB", font, col=(60, 140, 80))
    steam_col = C_GRAY_PIPE if snap.turbine_valve < 0.1 else C_CYAN_PIPE
    _pipe(screen, [(rv_x + rv_w, rv_y + 20), (tb_x, tb_y + 30)], steam_col, 4)

    # Condenser
    cd_x, cd_y, cd_w, cd_h = x + 290, y + 160, 70, 55
    _box(screen, cd_x, cd_y, cd_w, cd_h, "COND", font, col=(50, 90, 120))
    _pipe(screen, [(tb_x + tb_w // 2, tb_y + tb_h), (tb_x + tb_w // 2, cd_y)], C_GRAY_PIPE, 3)

    # Feedwater
    _pipe(screen, [(cd_x, cd_y + cd_h // 2), (rv_x + rv_w, rv_y + rv_h - 40)], C_BLUE_PIPE, 3)

    # Recirculation loop
    rp_x, rp_y, rp_w, rp_h = x + 70, y + h - 60, 70, 44
    _box(screen, rp_x, rp_y, rp_w, rp_h, "RECIRC", font, col=(70, 110, 180))
    _pipe(screen, [(rv_x + rv_w // 2, rv_y + rv_h), (rv_x + rv_w // 2, rp_y),
                   (rp_x + rp_w, rp_y + rp_h // 2), (rv_x, rv_y + rv_h - 60)], C_BLUE_PIPE, 3)

    # Void fraction bar
    vf = snap.void_fraction
    vbar_x = rv_x - 28
    vbar_h = int(rv_h * vf)
    pygame.draw.rect(screen, C_BORDER, (vbar_x, rv_y, 18, rv_h), border_radius=3)
    if vbar_h > 0:
        pygame.draw.rect(screen, C_CYAN_PIPE, (vbar_x, rv_y + rv_h - vbar_h, 18, vbar_h), border_radius=2)
    vf_lbl = font.render(f"α={vf:.2f}", True, C_CYAN_PIPE)
    screen.blit(vf_lbl, (vbar_x - 4, rv_y + rv_h + 4))

    def small_val(sx, sy, txt, col=C_TEXT):
        screen.blit(font.render(txt, True, col), (sx, sy))

    small_val(rv_x + rv_w + 4, rv_y + 10, f"{snap.power_fraction*100:.0f}%",
              C_RED if frac > 1.1 else C_GREEN)
    small_val(rv_x + rv_w + 4, rv_y + 28, f"{snap.pressure_mpa:.2f}MPa", C_TEXT)
    small_val(cd_x + cd_w + 4, cd_y + 10, f"{snap.condenser_temp_k:.0f}K", C_TEXT)


def draw_pid_rbmk(screen, x, y, w, h, snap, font, blink):
    """Simplified RBMK P&ID schematic (channel-tube style)."""
    frac = clamp(snap.power_fraction, 0, 1)
    fuel_alarm = snap.fuel_temp_k > 1200

    # Graphite block representation
    gx, gy, gw, gh = x + 50, y + 30, 160, 200
    pygame.draw.rect(screen, (30, 30, 20), (gx, gy, gw, gh), border_radius=6)
    pygame.draw.rect(screen, (80, 80, 50), (gx, gy, gw, gh), 2, border_radius=6)
    lbl = font.render("GRAPHITE", True, (120, 120, 80))
    screen.blit(lbl, (gx + gw // 2 - lbl.get_width() // 2, gy + 6))

    # Fuel channels (3 representative)
    for i, ch_x in enumerate([gx + 30, gx + 70, gx + 110]):
        ch_col = (min(255, 50 + int(200 * frac)), 30, 10) if frac > 0.1 else (50, 60, 80)
        pygame.draw.rect(screen, ch_col, (ch_x, gy + 20, 22, gh - 40), border_radius=3)
        void_h = int((gh - 40) * snap.void_fraction)
        if void_h > 0:
            pygame.draw.rect(screen, C_CYAN_PIPE, (ch_x + 2, gy + 20, 18, void_h), border_radius=2)

    # Header drum
    hd_x, hd_y, hd_w, hd_h = x + 270, y + 30, 80, 50
    _box(screen, hd_x, hd_y, hd_w, hd_h, "DRUM", font, col=(70, 100, 160))

    # Steam pipes to turbine
    tb_x, tb_y, tb_w, tb_h = x + 420, y + 30, 70, 60
    _box(screen, tb_x, tb_y, tb_w, tb_h, "TURB", font, col=(60, 140, 80))
    steam_col = C_GRAY_PIPE if snap.turbine_valve < 0.1 else C_CYAN_PIPE
    _pipe(screen, [(gx + gw, gy + 40), (hd_x, hd_y + 25)], steam_col, 3)
    _pipe(screen, [(hd_x + hd_w, hd_y + 25), (tb_x, tb_y + 30)], steam_col, 3)

    # Condenser / feedwater
    cd_x, cd_y, cd_w, cd_h = x + 420, y + 160, 70, 55
    _box(screen, cd_x, cd_y, cd_w, cd_h, "COND", font, col=(50, 90, 120))
    _pipe(screen, [(tb_x + tb_w // 2, tb_y + tb_h), (tb_x + tb_w // 2, cd_y)], C_GRAY_PIPE, 3)
    _pipe(screen, [(cd_x, cd_y + cd_h // 2), (gx, gy + gh - 40)], C_BLUE_PIPE, 3)

    # Pump
    pp_x, pp_y, pp_w, pp_h = x + 50, y + h - 70, 70, 44
    _box(screen, pp_x, pp_y, pp_w, pp_h, "PUMP", font, col=(70, 110, 180))
    _pipe(screen, [(gx + gw // 2, gy + gh), (gx + gw // 2, pp_y), (pp_x + pp_w, pp_y + pp_h // 2),
                   (gx, gy + gh - 60)], C_BLUE_PIPE, 3)

    def small_val(sx, sy, txt, col=C_TEXT):
        screen.blit(font.render(txt, True, col), (sx, sy))

    small_val(gx + gw + 4, gy + 10, f"{snap.power_fraction*100:.0f}%",
              C_RED if frac > 1.1 else C_GREEN)
    small_val(gx + gw + 4, gy + 28, f"void={snap.void_fraction:.2f}", C_CYAN_PIPE)
    small_val(gx + gw + 4, gy + 46, f"{snap.fuel_temp_k:.0f}K",
              C_RED if fuel_alarm else C_TEXT)


# ── Tab screens ──────────────────────────────────────────────────────────────

TABS = [
    ("F1", "OVERVIEW"),
    ("F2", "PRIMARY"),
    ("F3", "SECONDARY"),
    ("F4", "ALARMS"),
    ("F5", "I&C"),
    ("F6", "REACTIVITY"),
]


def draw_tab_bar(screen, active_tab, font, y, w):
    tab_w = w // len(TABS)
    for i, (fkey, name) in enumerate(TABS):
        bx = i * tab_w
        active = i == active_tab
        bg = (35, 55, 90) if active else (16, 20, 30)
        bc = C_ACCENT if active else C_BORDER
        pygame.draw.rect(screen, bg, (bx, y, tab_w, 34))
        pygame.draw.rect(screen, bc, (bx, y, tab_w, 34), 1 if not active else 2)
        lbl = font.render(f"{fkey}  {name}", True, C_TEXT if active else C_TEXT_DIM)
        screen.blit(lbl, (bx + tab_w // 2 - lbl.get_width() // 2, y + 17 - lbl.get_height() // 2))


def draw_screen_overview(screen, snap, fonts, layout, blink_fast, blink_slow, controls,
                         scram_msg, scram_msg_timer, ic_msg, ic_msg_timer, tw_label,
                         p_data, f_data, r_data, d_data):
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy_body, cw, rp_x, rp_w, body_h = layout

    # Center: P&ID schematic (top) + readouts (bottom)
    pid_h = 310
    readout_h = body_h - pid_h - 8

    pid_y = cy_body
    ro_y  = pid_y + pid_h + 8

    draw_panel(screen, cx, pid_y, cw, pid_h, "PLANT OVERVIEW", font_sm)
    rt = snap.reactor_type
    pid_draw = {"PWR": draw_pid_pwr, "BWR": draw_pid_bwr, "RBMK": draw_pid_rbmk}[rt]
    pid_draw(screen, cx + 8, pid_y + 28, cw - 16, pid_h - 36, snap, font_sm, blink_fast)

    # Right: bar meters + trends
    bm_h = 200
    trend_y = cy_body + bm_h + 8
    trend_h = body_h - bm_h - 8

    draw_panel(screen, rp_x, cy_body, rp_w, bm_h, "INSTRUMENTS", font_sm)

    bm_area_x = rp_x + 12
    bm_area_w = rp_w - 24
    n_bars     = 4
    bar_w      = (bm_area_w - (n_bars - 1) * 8) // n_bars
    bar_h_px   = bm_h - 56

    pf = snap.power_fraction * 100
    pnom = 15.5 if snap.reactor_type == "PWR" else 7.0

    bar_defs = [
        ("Power", "%",   pf,               0, 130,
         [(90, 110, C_YELLOW), (110, 130, C_RED)], None, 120.0),
        ("Pressure", "MPa", snap.pressure_mpa, 0, pnom * 1.15,
         [(pnom * 1.04, pnom * 1.15, C_RED)], None, pnom * 1.08),
        ("Fuel T", "K",   snap.fuel_temp_k, 400, 1800,
         [(900, 1200, C_YELLOW), (1200, 1800, C_RED)], None, 1400.0),
        ("Void", "",      snap.void_fraction, 0, 1.0,
         [(0.4, 0.7, C_YELLOW), (0.7, 1.0, C_RED)], None, None),
    ]
    for i, (lbl, unit, val, lo, hi, bands, tlo, thi) in enumerate(bar_defs):
        bx = bm_area_x + i * (bar_w + 8)
        draw_bar_meter(screen, bx, cy_body + 32, bar_w, bar_h_px, val, lo, hi,
                       font_sm, label=lbl, unit=unit, bands=bands, trip_lo=tlo, trip_hi=thi)

    # Readout text panel (below P&ID)
    ry = draw_panel(screen, cx, ro_y, cw, readout_h, "PLANT STATUS", font_sm)

    def rline(label, val_str, unit, col):
        nonlocal ry
        screen.blit(font.render(label, True, C_TEXT_DIM), (cx + 12, ry))
        vs  = font_md.render(val_str, True, col)
        uw  = font_sm.size(unit)[0]
        screen.blit(vs,  (cx + cw - uw - vs.get_width() - 18, ry - 1))
        screen.blit(font_sm.render(unit, True, C_TEXT_DIM), (cx + cw - uw - 10, ry + 2))
        ry += 22

    def rsep():
        nonlocal ry
        pygame.draw.line(screen, C_BORDER, (cx + 10, ry + 2), (cx + cw - 10, ry + 2))
        ry += 10

    rline("Core power",    f"{pf:7.2f}",                     "%",
          C_RED if pf > 110 else C_YELLOW if pf > 90 else C_GREEN)
    rline("Thermal power", f"{snap.thermal_mw:7.1f}",        "MWt", C_TEXT)
    rline("Electric power",f"{snap.electric_mw:7.1f}",       "MWe", C_TEXT)
    dh = snap.decay_heat_fraction * 100.0
    rline("Decay heat",    f"{dh:7.2f}",                     "%",
          C_YELLOW if dh > 3.0 else C_TEXT)
    rsep()
    rk = snap.reactivity
    rline("Reactivity",    f"{rk:+.5f}",                     "dk/k",
          C_RED if abs(rk) > 0.005 else C_YELLOW if abs(rk) > 0.001 else C_GREEN)
    rsep()
    ft = snap.fuel_temp_k
    rline("Fuel temp",     f"{ft:7.1f}",                     "K",
          C_RED if ft > 1400 else C_ORANGE if ft > 1200 else C_YELLOW if ft > 900 else C_GREEN)
    ct = snap.coolant_temp_k
    rline("Coolant temp",  f"{ct:7.1f}",                     "K",
          C_RED if ct > 620 else C_YELLOW if ct > 580 else C_GREEN)
    rline("Steam temp",    f"{snap.steam_temp_k:7.1f}",      "K", C_TEXT)
    rline("Condenser",     f"{snap.condenser_temp_k:7.1f}",  "K",
          C_YELLOW if snap.condenser_temp_k > 330 else C_TEXT)
    rsep()
    p = snap.pressure_mpa
    rline("Pressure",      f"{p:7.3f}",                      "MPa",
          C_RED if p > pnom * 1.08 else C_YELLOW if p > pnom * 1.04 else C_GREEN)
    rline("Void frac.",    f"{snap.void_fraction:7.3f}",     "",
          C_YELLOW if snap.void_fraction > 0.3 else C_TEXT)
    rline("Steam inv.",    f"{snap.steam_inventory:7.3f}",   "",
          C_YELLOW if snap.steam_inventory > 1.5 or snap.steam_inventory < 0.3 else C_TEXT)
    fw = snap.feedwater_inventory
    rline("Feedwater",     f"{fw:7.3f}",                     "",
          C_RED if fw < 0.1 else C_YELLOW if fw < 0.25 else C_TEXT)
    if snap.reactor_type == "PWR":
        bp = snap.boron_ppm
        rline("Boron",     f"{bp:7.1f}",                     "ppm",
              C_YELLOW if bp < 100 else C_TEXT)
    rsep()
    rline("Sim time",      f"{snap.time:7.1f}",              "s", C_TEXT_DIM)

    # Trends (right panel)
    tw_label_str = f"TRENDS  [{tw_label}]"
    ty = draw_panel(screen, rp_x, trend_y, rp_w, trend_h, tw_label_str, font_sm)
    GW = rp_w - 20
    GX = rp_x + 10
    GAP = 5
    GH  = (trend_h - (ty - trend_y) - GAP * 3) // 4
    for i, (hist, color, ylo, yhi, lbl2) in enumerate([
        (p_data, C_GREEN,  0.0,   130.0, "Core power (%)"),
        (f_data, C_ORANGE, 400.0, 1800.0,"Fuel temp (K)"),
        (r_data, C_ACCENT, -0.01, 0.01,  "Reactivity (dk/k)"),
        (d_data, C_YELLOW, 0.0,   10.0,  "Decay heat (%)"),
    ]):
        gy = ty + i * (GH + GAP)
        draw_trend(screen, GX, gy, GW, GH, list(hist), color, ylo, yhi, lbl2, font_sm)


def draw_screen_primary(screen, snap, supervisor, fonts, layout, blink_fast):
    """F2: Primary system — pressurizer, RCP, pressure detail."""
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy, cw, rp_x, rp_w, body_h = layout
    pnom = 15.5 if snap.reactor_type == "PWR" else 7.0

    cy2 = draw_panel(screen, cx, cy, cw + rp_w + 8, body_h, "PRIMARY SYSTEM", font_sm)

    # Large bar meters
    bm_y  = cy2 + 20
    bm_h2 = body_h - (cy2 - cy) - 140
    bar_w2 = 60
    bar_gap = 30

    defs = [
        ("Pressure", "MPa", snap.pressure_mpa, 0, pnom * 1.2,
         [(pnom * 1.04, pnom * 1.1, C_YELLOW), (pnom * 1.1, pnom * 1.2, C_RED)], None, pnom * 1.1),
        ("Fuel T", "K", snap.fuel_temp_k, 300, 1800,
         [(900, 1200, C_YELLOW), (1200, 1800, C_RED)], None, 1400.0),
        ("Coolant T", "K", snap.coolant_temp_k, 300, 700,
         [(580, 620, C_YELLOW), (620, 700, C_RED)], None, 620.0),
        ("Flow", "%", supervisor.bop.omega_rcp * 100, 0, 110,
         [(0, 87, C_ORANGE)], 87.0, None),
        ("Steam inv.", "", snap.steam_inventory, 0, 2.0,
         [(0, 0.3, C_RED), (1.5, 2.0, C_YELLOW)], 0.3, None),
    ]
    bx_start = cx + 40
    for i, (lbl, unit, val, lo, hi, bands, tlo, thi) in enumerate(defs):
        bx = bx_start + i * (bar_w2 + bar_gap)
        draw_bar_meter(screen, bx, bm_y, bar_w2, bm_h2, val, lo, hi,
                       font_sm, label=lbl, unit=unit, bands=bands, trip_lo=tlo, trip_hi=thi)

    # DCS readout grid
    box_y = cy2 + 20 + bm_h2 + 40
    porv_str = "OPEN" if snap.porv_open else "SHUT"
    eccs_str = "ACTIVE" if snap.eccs_actuated else "STANDBY"
    box_defs = [
        ("Core power", f"{snap.power_fraction*100:.2f} %",  C_GREEN if snap.power_fraction < 1.1 else C_RED),
        ("Pressure",   f"{snap.pressure_mpa:.3f} MPa",      C_GREEN if snap.pressure_mpa < pnom*1.04 else C_YELLOW),
        ("Fuel temp",  f"{snap.fuel_temp_k:.1f} K",         C_GREEN if snap.fuel_temp_k < 900 else C_ORANGE),
        ("Coolant T",  f"{snap.coolant_temp_k:.1f} K",      C_GREEN if snap.coolant_temp_k < 580 else C_YELLOW),
        ("RCP speed",  f"{supervisor.bop.omega_rcp*100:.1f} %",
         C_GREEN if supervisor.bop.omega_rcp > 0.90 else C_ORANGE),
        ("PORV",       porv_str,  C_ORANGE if snap.porv_open else C_GREEN),
        ("ECCS",       eccs_str,  C_RED if snap.eccs_actuated else C_TEXT_DIM),
        ("Reactivity", f"{snap.reactivity:+.5f}",           C_GREEN if abs(snap.reactivity) < 0.001 else C_YELLOW),
    ]
    bw3, bh3 = 182, 52  # 8 boxes × (182+10) = 1536, fits in cx+cw+rp_w=1546 with 10px margin
    for i, (lbl3, val3, col3) in enumerate(box_defs):
        bx3 = cx + 20 + i * (bw3 + 10)
        draw_dcs_box(screen, bx3, box_y, bw3, bh3, lbl3, val3, col3, font_sm, font_md)


def draw_screen_secondary(screen, snap, fonts, layout, blink_fast):
    """F3: Secondary system — SG, feedwater, turbine."""
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy, cw, rp_x, rp_w, body_h = layout

    cy2 = draw_panel(screen, cx, cy, cw + rp_w + 8, body_h, "SECONDARY SYSTEM", font_sm)

    bar_defs = [
        ("Steam temp", "K",   snap.steam_temp_k,       400, 700,
         [(600, 650, C_YELLOW), (650, 700, C_RED)], None, None),
        ("Condenser",  "K",   snap.condenser_temp_k,   280, 380,
         [(330, 360, C_YELLOW), (360, 380, C_RED)], None, 340.0),
        ("Steam inv.", "",    snap.steam_inventory,    0, 2.0,
         [(0, 0.3, C_RED), (1.5, 2.0, C_YELLOW)], 0.3, None),
        ("FW inv.",    "",    snap.feedwater_inventory, 0, 1.2,
         [(0, 0.1, C_RED), (0.1, 0.25, C_YELLOW)], 0.1, None),
        ("Turb valve", "%",   snap.electric_mw / max(1, snap.thermal_mw) * 100, 0, 105,
         [], None, None),
        ("Elec power", "MWe", snap.electric_mw,        0, 1500,
         [], None, None),
    ]
    bm_h2 = body_h - (cy2 - cy) - 140
    bar_w2 = 65
    bar_gap = 24
    bm_y = cy2 + 20
    for i, (lbl, unit, val, lo, hi, bands, tlo, thi) in enumerate(bar_defs):
        bx = cx + 40 + i * (bar_w2 + bar_gap)
        draw_bar_meter(screen, bx, bm_y, bar_w2, bm_h2, val, lo, hi,
                       font_sm, label=lbl, unit=unit, bands=bands, trip_lo=tlo, trip_hi=thi)

    box_y = cy2 + 20 + bm_h2 + 40
    box_defs = [
        ("Thermal power", f"{snap.thermal_mw:.1f} MWt",       C_TEXT),
        ("Electric power",f"{snap.electric_mw:.1f} MWe",      C_TEXT),
        ("Steam temp",    f"{snap.steam_temp_k:.1f} K",        C_TEXT),
        ("Condenser",     f"{snap.condenser_temp_k:.1f} K",
         C_YELLOW if snap.condenser_temp_k > 330 else C_TEXT),
        ("Feedwater",     f"{snap.feedwater_inventory:.3f}",
         C_RED if snap.feedwater_inventory < 0.1 else C_TEXT),
        ("Decay heat",    f"{snap.decay_heat_fraction*100:.2f} %",
         C_YELLOW if snap.decay_heat_fraction > 0.03 else C_TEXT),
    ]
    bw3, bh3 = 200, 52
    for i, (lbl3, val3, col3) in enumerate(box_defs):
        bx3 = cx + 20 + i * (bw3 + 12)
        if bx3 + bw3 < cx + cw + rp_w:
            draw_dcs_box(screen, bx3, box_y, bw3, bh3, lbl3, val3, col3, font_sm, font_md)


def draw_screen_alarms(screen, snap, supervisor, fonts, layout, blink_fast, blink_slow):
    """F4: Alarm annunciator grid + event sequence recorder."""
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy, cw, rp_x, rp_w, body_h = layout
    full_w = cw + rp_w + 8

    alarm_panel_h = body_h // 2 - 4
    log_panel_h   = body_h - alarm_panel_h - 8

    ay = draw_panel(screen, cx, cy, full_w, alarm_panel_h, "ALARM ANNUNCIATOR", font_sm)

    all_alarms = snap.alarm_objects
    TILE_W, TILE_H, TILE_GAP = 220, 44, 6
    cols = (full_w - 20) // (TILE_W + TILE_GAP)
    for idx, alarm in enumerate(all_alarms):
        col_i = idx % cols
        row_i = idx // cols
        tx = cx + 10 + col_i * (TILE_W + TILE_GAP)
        ty = ay + row_i * (TILE_H + TILE_GAP)
        if ty + TILE_H > cy + alarm_panel_h:
            break
        if alarm.is_trip:
            bg = (100, 10, 10) if not blink_fast else (180, 20, 20)
            bc = C_RED
        elif alarm.state == "unack":
            bg = (80, 50, 0) if not blink_slow else (140, 90, 0)
            bc = C_ORANGE
        else:
            bg = (40, 50, 30)
            bc = C_GREEN
        pygame.draw.rect(screen, bg, (tx, ty, TILE_W, TILE_H), border_radius=4)
        pygame.draw.rect(screen, bc, (tx, ty, TILE_W, TILE_H), 1, border_radius=4)
        tag = "TRIP" if alarm.is_trip else "WARN" if alarm.state == "unack" else "ACK"
        tag_col = C_RED if alarm.is_trip else C_ORANGE if alarm.state == "unack" else C_GREEN
        screen.blit(font_sm.render(tag, True, tag_col), (tx + 5, ty + 4))
        if alarm.first_out:
            fo = font_sm.render("FO", True, (255, 220, 0))
            screen.blit(fo, (tx + TILE_W - fo.get_width() - 5, ty + 4))
        msg_lbl = font_sm.render(alarm.message[:28], True, C_TEXT)
        screen.blit(msg_lbl, (tx + 5, ty + TILE_H - msg_lbl.get_height() - 4))

    if not all_alarms:
        ok_lbl = font_md.render("All systems nominal", True, C_GREEN)
        mid_x = cx + full_w // 2
        mid_y = cy + alarm_panel_h // 2
        draw_led(screen, mid_x - ok_lbl.get_width() // 2 - 16, mid_y, C_GREEN, 7)
        screen.blit(ok_lbl, (mid_x - ok_lbl.get_width() // 2, mid_y - ok_lbl.get_height() // 2))

    # Event sequence recorder
    log_y = cy + alarm_panel_h + 8
    ly = draw_panel(screen, cx, log_y, full_w, log_panel_h, "EVENT SEQUENCE RECORDER", font_sm)
    log_entries = list(supervisor.event_log)[-20:]
    for entry in reversed(log_entries):
        sim_t, cat, msg = entry
        line = font_sm.render(f"T={sim_t:8.1f}s  [{cat:12s}]  {msg}", True, C_TEXT_DIM)
        screen.blit(line, (cx + 10, ly))
        ly += font_sm.get_height() + 2
        if ly > log_y + log_panel_h - 10:
            break


def draw_screen_ic(screen, snap, fonts, layout, blink_fast):
    """F5: I&C / protection channels — real 2-of-3 voting display."""
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy, cw, rp_x, rp_w, body_h = layout
    full_w = cw + rp_w + 8
    cy2 = draw_panel(screen, cx, cy, full_w, body_h, "I&C / PROTECTION CHANNELS  (2-of-3 voting)", font_sm)

    pnom = 15.5 if snap.reactor_type == "PWR" else 7.0
    params = [
        ("Neutron flux",      snap.power_fraction * 100, "%",   120.0,  130.0),
        ("Fuel temperature",  snap.fuel_temp_k,          "K",   1400.0, None),
        ("Reactor pressure",  snap.pressure_mpa,         "MPa", pnom * 1.08, pnom * 1.1),
        ("Coolant temp",      snap.coolant_temp_k,       "K",   610.0,  616.0),
    ]
    noise = snap.channel_noise  # live noise from supervisor
    row_h = 65
    for ri, (param_name, value, unit, warn_sp, trip_sp) in enumerate(params):
        py = cy2 + ri * (row_h + 8)
        screen.blit(font_md.render(param_name, True, C_TEXT_HDR), (cx + 10, py + 4))
        val_str = f"{value:.3f} {unit}"
        warn_ok = value < warn_sp if warn_sp else True
        col = C_GREEN if warn_ok else (C_YELLOW if not trip_sp or value < trip_sp else C_RED)
        screen.blit(font_md.render(val_str, True, col), (cx + 280, py + 4))

        # Real 2-of-3 channel display using live noise
        ch_vals = [value * (1.0 + n) for n in noise]
        votes = sum(1 for v in ch_vals if trip_sp and v > trip_sp)
        vote_col = C_RED if votes >= 2 else C_YELLOW if votes == 1 else C_TEXT_DIM
        screen.blit(font_sm.render(f"Votes: {votes}/3", True, vote_col), (cx + 480, py + 4))
        for ch_i, ch_val in enumerate(ch_vals):
            ch_x = cx + 560 + ch_i * 160
            ch_str = f"Ch {chr(65+ch_i)}: {ch_val:.3f}"
            ch_col = C_GREEN
            if warn_sp and ch_val > warn_sp:
                ch_col = C_YELLOW
            if trip_sp and ch_val > trip_sp:
                ch_col = C_RED
            screen.blit(font.render(ch_str, True, ch_col), (ch_x, py + 4))

        setpt_str = (f"Warn: {warn_sp:.2f} {unit}  Trip: {trip_sp:.2f} {unit}"
                     if trip_sp else f"Warn: {warn_sp:.2f} {unit}")
        screen.blit(font_sm.render(setpt_str, True, C_TEXT_DIM), (cx + 10, py + 30))
        pygame.draw.line(screen, C_BORDER, (cx + 10, py + row_h + 2), (cx + full_w - 10, py + row_h + 2))


def draw_screen_reactivity(screen, snap, fonts, layout, p_data, f_data):
    """F6: Reactivity budget breakdown + xenon/iodine trends."""
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy, cw, rp_x, rp_w, body_h = layout
    full_w = cw + rp_w + 8

    cy2 = draw_panel(screen, cx, cy, full_w, body_h, "REACTIVITY BUDGET & POISON TRENDS", font_sm)

    # Stacked horizontal bar chart (reactivity components, qualitative)
    rho_total = snap.reactivity
    pf = snap.power_fraction
    pnom = 15.5 if snap.reactor_type == "PWR" else 7.0

    # Approximate components from visible state
    rod_pos = 0.5  # unknown here, so approximate
    rho_decay = snap.decay_heat_fraction * 0.001

    components = [
        ("Total reactivity",  rho_total,  C_ACCENT),
        ("Xenon (approx)",    -0.02 * (1 - pf),   C_YELLOW),
        ("Decay heat effect", rho_decay,  C_ORANGE),
    ]

    bar_y = cy2 + 20
    bar_max_w = full_w - 80
    zero_x = cx + full_w // 2

    screen.blit(font_md.render("Reactivity Components (qualitative)", True, C_TEXT_HDR), (cx + 10, bar_y))
    bar_y += 28
    pygame.draw.line(screen, C_TEXT_DIM, (zero_x, bar_y), (zero_x, bar_y + len(components) * 40 + 10))

    scale = bar_max_w / 2 / max(0.05, max(abs(v) for _, v, _ in components))
    for name, val, col in components:
        blen = int(abs(val) * scale)
        bx = zero_x if val >= 0 else zero_x - blen
        pygame.draw.rect(screen, col, (bx, bar_y + 5, blen, 24), border_radius=3)
        screen.blit(font_sm.render(f"{name}: {val:+.5f}", True, C_TEXT), (cx + 10, bar_y + 8))
        bar_y += 40

    # Large trends: power + xenon proxy
    trend_y = bar_y + 20
    trend_h  = body_h - (trend_y - cy)
    TW = (full_w - 24) // 2
    draw_trend(screen, cx + 8, trend_y, TW, trend_h, list(p_data),
               C_GREEN, 0.0, 130.0, "Core power (%)", font_sm)
    draw_trend(screen, cx + TW + 16, trend_y, TW, trend_h, list(f_data),
               C_ORANGE, 400.0, 1800.0, "Fuel temperature (K)", font_sm)


# ── Main loop ────────────────────────────────────────────────────────────────

def main():
    pygame.init()
    screen = pygame.display.set_mode((W, H))
    pygame.display.set_caption("Nuclear Reactor Control Room")

    font_sm = pygame.font.SysFont("Menlo", 13)
    font    = pygame.font.SysFont("Menlo", 15)
    font_md = pygame.font.SysFont("Menlo", 17, bold=True)
    font_lg = pygame.font.SysFont("Menlo", 22, bold=True)
    font_xl = pygame.font.SysFont("Menlo", 32, bold=True)
    fonts   = (font_sm, font, font_md, font_lg, font_xl)

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
    _slow_acc  = 0.0

    time_window_idx = 0
    time_speed_idx  = 0
    time_speed      = 1
    active_tab      = 0   # 0-5 for F1-F6

    scram_msg       = ""
    scram_msg_timer = 0.0
    ic_path         = os.path.join(os.path.dirname(__file__), "ic_quicksave.pkl")
    ic_msg          = ""
    ic_msg_timer    = 0.0

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
    effective_speed = 1

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
                elif ev.key == pygame.K_ESCAPE:    running = False
                elif ev.key == pygame.K_1:          select("PWR")
                elif ev.key == pygame.K_2:          select("BWR")
                elif ev.key == pygame.K_3:          select("RBMK")
                elif ev.key == pygame.K_SPACE:      controls.scram = True
                elif ev.key == pygame.K_r:          select(supervisor.reactor_type)
                elif ev.key == pygame.K_t:          controls.turbine_trip         = not controls.turbine_trip
                elif ev.key == pygame.K_p:          controls.startup_permit       = not controls.startup_permit
                elif ev.key == pygame.K_z:          controls.fault_pump_degraded  = not controls.fault_pump_degraded
                elif ev.key == pygame.K_x:          controls.fault_feedwater_loss = not controls.fault_feedwater_loss
                elif ev.key == pygame.K_i:
                    # Cycle LOCA break area: 0 → 0.001 → 0.01 → 0
                    LOCA_AREAS = [0.0, 0.001, 0.01]
                    cur = controls.fault_loca_break_area
                    idx_next = (LOCA_AREAS.index(min(LOCA_AREAS, key=lambda a: abs(a - cur))) + 1) % len(LOCA_AREAS)
                    controls.fault_loca_break_area = LOCA_AREAS[idx_next]
                    supervisor.event_log.append((snapshot.time, "INSTRUCTOR", f"LOCA break area set to {LOCA_AREAS[idx_next]} m²"))
                elif ev.key == pygame.K_o:
                    if supervisor.reactor_type == "PWR":
                        supervisor.auto_rod.auto = not supervisor.auto_rod.auto
                        supervisor.event_log.append((snapshot.time, "AUTO_CTRL", f"Auto rod: {'ON' if supervisor.auto_rod.auto else 'OFF'}"))
                elif ev.key == pygame.K_m:
                    if supervisor.reactor_type == "PWR":
                        supervisor.auto_pressure.auto = not supervisor.auto_pressure.auto
                        supervisor.event_log.append((snapshot.time, "AUTO_CTRL", f"Auto pressure: {'ON' if supervisor.auto_pressure.auto else 'OFF'}"))
                elif ev.key == pygame.K_l:
                    scram_msg = supervisor.reset_scram()
                    scram_msg_timer = 4.0
                elif ev.key == pygame.K_c:          supervisor.acknowledge_alarms()
                elif ev.key == pygame.K_LEFTBRACKET:
                    time_window_idx = (time_window_idx - 1) % len(TIME_WINDOWS)
                elif ev.key == pygame.K_RIGHTBRACKET:
                    time_window_idx = (time_window_idx + 1) % len(TIME_WINDOWS)
                elif ev.key in (pygame.K_EQUALS, pygame.K_PLUS):
                    time_speed_idx = min(time_speed_idx + 1, len(TIME_SPEEDS) - 1)
                    time_speed = TIME_SPEEDS[time_speed_idx]
                elif ev.key == pygame.K_MINUS:
                    time_speed_idx = max(time_speed_idx - 1, 0)
                    time_speed = TIME_SPEEDS[time_speed_idx]
                elif ev.key == pygame.K_F1:  active_tab = 0
                elif ev.key == pygame.K_F2:  active_tab = 1
                elif ev.key == pygame.K_F3:  active_tab = 2
                elif ev.key == pygame.K_F4:  active_tab = 3
                elif ev.key == pygame.K_F5:  active_tab = 4
                elif ev.key == pygame.K_F6:  active_tab = 5

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
        if keys[pygame.K_j]: dil.set(dil.v - step)

        controls.rod_position       = rod.v
        controls.flow               = flow.v
        controls.turbine_valve      = valve.v
        controls.feedwater_valve    = feed.v
        controls.pressurizer_heater = press.v
        controls.boration_rate      = bor.v
        controls.dilution_rate      = dil.v

        # ── Physics step ─────────────────────────────────────────────────────
        now     = time.time()
        real_dt = clamp(now - last, 0.0, 0.2)
        last    = now

        effective_speed = time_speed if not snapshot.alarms and not snapshot.trips else 1
        dt = real_dt * effective_speed

        snapshot = supervisor.step(dt)

        # Sync sliders to auto controller outputs so display tracks correctly
        if supervisor.auto_rod.auto and supervisor.reactor_type == "PWR":
            rod.set(controls.rod_position)
        if supervisor.auto_pressure.auto and supervisor.reactor_type == "PWR":
            press.set(controls.pressurizer_heater)
        hist_power.append(snapshot.power_fraction * 100.0)
        hist_fuel.append(snapshot.fuel_temp_k)
        hist_react.append(snapshot.reactivity)
        hist_decay.append(snapshot.decay_heat_fraction * 100.0)

        _slow_acc += real_dt
        if _slow_acc >= 1.0:
            _slow_acc -= 1.0
            slow_power.append(snapshot.power_fraction * 100.0)
            slow_fuel.append(snapshot.fuel_temp_k)
            slow_react.append(snapshot.reactivity)
            slow_decay.append(snapshot.decay_heat_fraction * 100.0)

        scram_msg_timer = max(0.0, scram_msg_timer - real_dt)
        ic_msg_timer    = max(0.0, ic_msg_timer    - real_dt)

        # ── Render ───────────────────────────────────────────────────────────
        t_ms       = pygame.time.get_ticks()
        blink_fast = (t_ms // 300) % 2 == 0
        blink_slow = (t_ms // 600) % 2 == 0
        has_alarms = bool(snapshot.alarms)
        has_trips  = bool(snapshot.trips)
        is_scram   = controls.scram

        screen.fill(C_BG)

        # ── Header ───────────────────────────────────────────────────────────
        HEADER_H = 52
        TAB_H    = 36
        PAD      = 8
        HELP_H   = 24

        if is_scram:
            scram_bg = (180, 15, 15) if blink_fast else (90, 8, 8)
            pygame.draw.rect(screen, scram_bg, (0, 0, W, HEADER_H))
            txt = font_xl.render("⚠  REACTOR SCRAM  ──  REACTOR SCRAM  ──  REACTOR SCRAM  ⚠", True, (255, 255, 255))
            screen.blit(txt, (W // 2 - txt.get_width() // 2, HEADER_H // 2 - txt.get_height() // 2))
        else:
            pygame.draw.rect(screen, (12, 16, 26), (0, 0, W, HEADER_H))
            pygame.draw.line(screen, C_BORDER, (0, HEADER_H - 1), (W, HEADER_H - 1))
            for i, (rt, base_col) in enumerate([("PWR", (50, 110, 210)), ("BWR", (40, 160, 110)), ("RBMK", (200, 75, 50))]):
                active = snapshot.reactor_type == rt
                bx = 16 + i * 120
                bg = base_col if active else (25, 32, 46)
                pygame.draw.rect(screen, bg, (bx, 8, 110, 36), border_radius=6)
                if active:
                    bc2 = tuple(min(255, c + 70) for c in base_col)
                    pygame.draw.rect(screen, bc2, (bx, 8, 110, 36), 2, border_radius=6)
                t = font_lg.render(rt, True, (255, 255, 255) if active else (80, 100, 130))
                screen.blit(t, (bx + 55 - t.get_width() // 2, 26 - t.get_height() // 2))

            speed_str = f"  {effective_speed}×" if effective_speed > 1 else ""
            t_surf = font_lg.render(f"T = {snapshot.time:8.1f} s{speed_str}", True, C_TEXT)
            screen.blit(t_surf, (W // 2 - t_surf.get_width() // 2, 14))

            if has_trips:
                col = (220, 40, 40) if blink_fast else (100, 20, 20)
                pygame.draw.rect(screen, col, (W - 220, 8, 210, 36), border_radius=6)
                t = font_lg.render("TRIP ACTIVE", True, (255, 255, 255))
                screen.blit(t, (W - 220 + 105 - t.get_width() // 2, 26 - t.get_height() // 2))
            elif has_alarms:
                col = C_ORANGE if blink_slow else (100, 55, 0)
                pygame.draw.rect(screen, col, (W - 220, 8, 210, 36), border_radius=6)
                t = font_lg.render("ALARM", True, (255, 255, 255))
                screen.blit(t, (W - 220 + 105 - t.get_width() // 2, 26 - t.get_height() // 2))

        # ── Tab bar ──────────────────────────────────────────────────────────
        draw_tab_bar(screen, active_tab, font, HEADER_H, W)

        # ── Body layout ──────────────────────────────────────────────────────
        Y0     = HEADER_H + TAB_H + PAD
        BODY_H = H - HEADER_H - TAB_H - HELP_H - PAD * 2

        LP_X, LP_W = PAD, 340
        CP_X = LP_X + LP_W + PAD
        CP_W = 560
        RP_X = CP_X + CP_W + PAD
        RP_W = W - RP_X - PAD

        layout = (CP_X, Y0, CP_W, RP_X, RP_W, BODY_H)

        # ── Left panel (always visible): Controls ────────────────────────────
        cy = draw_panel(screen, LP_X, Y0, LP_W, BODY_H, "CONTROLS", font_sm)

        sliders_for_rt = [rod, flow, valve, feed, press]
        if snapshot.reactor_type == "PWR":
            sliders_for_rt += [bor]

        for sl in sliders_for_rt:
            lbl_s = font.render(sl.label, True, C_TEXT_DIM)
            val_s = font_md.render(f"{sl.v:.3f}", True, C_ACCENT)
            screen.blit(lbl_s, (LP_X + 14, cy))
            screen.blit(val_s, (LP_X + LP_W - val_s.get_width() - 14, cy))
            cy += max(lbl_s.get_height(), val_s.get_height()) + 3
            bx, bw, bh = LP_X + 14, LP_W - 28, 9
            pygame.draw.rect(screen, C_SLIDER_BG, (bx, cy, bw, bh), border_radius=4)
            fill = int(bw * (sl.v - sl.lo) / max(1e-6, sl.hi - sl.lo))
            if fill > 0:
                pygame.draw.rect(screen, C_SLIDER_FG, (bx, cy, fill, bh), border_radius=4)
            pygame.draw.rect(screen, (150, 190, 255), (bx + fill - 2, cy, 3, bh))
            cy += bh + 10

        pygame.draw.line(screen, C_BORDER, (LP_X + 10, cy), (LP_X + LP_W - 10, cy))
        cy += 10

        for label, state, col_on, key in [
            ("Startup permit",  controls.startup_permit,       C_GREEN,  "P"),
            ("Turbine trip",    controls.turbine_trip,         C_ORANGE, "T"),
            ("Pump degraded",   controls.fault_pump_degraded,  C_RED,    "Z"),
            ("Feedwater fault", controls.fault_feedwater_loss, C_RED,    "X"),
            ("LOCA fault",      controls.fault_loca_break_area > 0, C_RED, "I"),
        ]:
            led_col = col_on if state else (30, 38, 52)
            draw_led(screen, LP_X + 22, cy + 8, led_col)
            tc = col_on if state else C_TEXT_DIM
            screen.blit(font.render(label, True, tc), (LP_X + 40, cy))
            kt = font_sm.render(f"[{key}]", True, C_TEXT_DIM)
            screen.blit(kt, (LP_X + LP_W - kt.get_width() - 12, cy + 2))
            cy += 24

        pygame.draw.line(screen, C_BORDER, (LP_X + 10, cy + 2), (LP_X + LP_W - 10, cy + 2))
        cy += 12

        if is_scram:
            scram_col = C_RED if blink_fast else C_ORANGE
            draw_led(screen, LP_X + 22, cy + 10, scram_col, radius=10)
            screen.blit(font_md.render("SCRAM ACTIVE",              True, scram_col),   (LP_X + 44, cy + 2))
            screen.blit(font_sm.render("L: reset  R: full restart", True, C_TEXT_DIM),  (LP_X + 44, cy + 22))
        else:
            draw_led(screen, LP_X + 22, cy + 10, C_GREEN, radius=10)
            screen.blit(font_md.render("REACTOR NORMAL",          True, C_GREEN),   (LP_X + 44, cy + 2))
            screen.blit(font_sm.render("[SPACE] Emergency SCRAM", True, C_TEXT_DIM), (LP_X + 44, cy + 22))
        cy += 44

        if scram_msg and scram_msg_timer > 0:
            mc = C_GREEN if "APPROVED" in scram_msg else C_ORANGE
            screen.blit(font_sm.render(scram_msg[:38], True, mc), (LP_X + 12, cy))
            cy += 18
        if ic_msg and ic_msg_timer > 0:
            screen.blit(font_sm.render(ic_msg, True, C_ACCENT), (LP_X + 12, cy))
            cy += 18

        # PORV / ECCS status (PWR only)
        if snapshot.reactor_type == "PWR":
            if snapshot.porv_open:
                draw_led(screen, LP_X + 22, cy + 8, C_ORANGE if blink_fast else (100, 55, 0))
                screen.blit(font_sm.render("PORV OPEN", True, C_ORANGE), (LP_X + 40, cy + 2))
                cy += 20
            if snapshot.eccs_actuated:
                draw_led(screen, LP_X + 22, cy + 8, C_RED if blink_fast else (80, 0, 0))
                screen.blit(font_sm.render("ECCS ACTUATED", True, C_RED), (LP_X + 40, cy + 2))
                cy += 20
            if snapshot.loca_area > 0:
                screen.blit(font_sm.render(f"LOCA {snapshot.loca_area*10000:.0f}cm²",
                                           True, C_RED), (LP_X + 12, cy))
                cy += 18

        # Auto controller indicators (PWR only)
        if snapshot.reactor_type == "PWR":
            auto_rod_col  = C_ACCENT if supervisor.auto_rod.auto      else C_TEXT_DIM
            auto_pres_col = C_ACCENT if supervisor.auto_pressure.auto  else C_TEXT_DIM
            screen.blit(font_sm.render(f"[O] Auto rod:  {'AUTO' if supervisor.auto_rod.auto else 'MAN '}",
                                       True, auto_rod_col), (LP_X + 12, cy))
            cy += 16
            screen.blit(font_sm.render(f"[M] Auto pzr:  {'AUTO' if supervisor.auto_pressure.auto else 'MAN '}",
                                       True, auto_pres_col), (LP_X + 12, cy))

        # ── Trend data selection ──────────────────────────────────────────────
        tw_label, tw_len = TIME_WINDOWS[time_window_idx]
        if tw_len <= TREND_LEN:
            p_data, f_data, r_data, d_data = hist_power, hist_fuel, hist_react, hist_decay
        else:
            n_slow = min(tw_len, SLOW_HIST_LEN)
            p_data = list(slow_power)[-n_slow:]
            f_data = list(slow_fuel)[-n_slow:]
            r_data = list(slow_react)[-n_slow:]
            d_data = list(slow_decay)[-n_slow:]

        # ── Tab content ───────────────────────────────────────────────────────
        if active_tab == 0:
            draw_screen_overview(screen, snapshot, fonts, layout, blink_fast, blink_slow,
                                 controls, scram_msg, scram_msg_timer, ic_msg, ic_msg_timer,
                                 tw_label, p_data, f_data, r_data, d_data)
        elif active_tab == 1:
            draw_screen_primary(screen, snapshot, supervisor, fonts, layout, blink_fast)
        elif active_tab == 2:
            draw_screen_secondary(screen, snapshot, fonts, layout, blink_fast)
        elif active_tab == 3:
            draw_screen_alarms(screen, snapshot, supervisor, fonts, layout, blink_fast, blink_slow)
        elif active_tab == 4:
            draw_screen_ic(screen, snapshot, fonts, layout, blink_fast)
        elif active_tab == 5:
            draw_screen_reactivity(screen, snapshot, fonts, layout, p_data, f_data)

        # ── Help bar ─────────────────────────────────────────────────────────
        help_y = H - HELP_H
        pygame.draw.line(screen, C_BORDER, (0, help_y - 3), (W, help_y - 3))
        screen.blit(font_sm.render(
            "W/S:rods  A/D:flow  Q/E:turbine  F/V:feed  H/N:pzr  B/G:borate  SHIFT:fast  "
            "O:auto-rod  M:auto-pzr  I:LOCA  1/2/3:reactor  SPACE:SCRAM  L:reset  C:ack  "
            "[/]:window  +/-:speed  Ctrl+S:save  F1-F6:tabs  ESC:quit",
            True, C_TEXT_DIM,
        ), (10, help_y))

        pygame.display.flip()
        clock.tick(60)

    # Export event log
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
