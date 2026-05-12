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
# Logical design resolution — all coordinates in the code use these units.
# SCALE (set in main() after detecting the display) maps them to native pixels.
W, H  = 1920, 1080
SCALE = 1.0


def _s(n: int) -> int:
    """Scale a design-space pixel value to native display pixels."""
    return max(1, int(n * SCALE))

# ── Color palette (ISA DCS standard) ─────────────────────────────────────────
C_BG        = ( 10,  12,  14)   # very dark background
C_PANEL     = ( 16,  19,  22)   # panel fill
C_PANEL_HDR = ( 20,  24,  28)   # title strip
C_BORDER    = ( 44,  50,  56)   # borders
C_SEP       = ( 30,  34,  38)   # separator lines

C_TEXT      = (215, 218, 215)   # normal values — neutral white
C_TEXT_DIM  = ( 85,  94,  88)   # labels / units — muted
C_TEXT_HDR  = (135, 145, 140)   # panel header text

# ISA status colors
C_NORMAL    = ( 55, 185,  75)   # green — within limits
C_CAUTION   = (205, 160,  18)   # amber — approaching limit
C_WARNING   = (195,  95,  18)   # orange — exceeded
C_ALARM     = (205,  38,  38)   # red — trip / emergency
C_ACCENT    = ( 48, 125, 195)   # blue — active / selected element

# Backward-compat aliases used in draw functions
C_GREEN     = C_NORMAL
C_YELLOW    = C_CAUTION
C_ORANGE    = C_WARNING
C_RED       = C_ALARM

# P&ID fluid colors
C_BLUE_PIPE = ( 38,  95, 185)   # subcooled liquid
C_CYAN_PIPE = ( 28, 145, 185)   # two-phase / boiling
C_GRAY_PIPE = (105, 118, 128)   # steam / vapor
C_HOT_LEG   = (185,  65,  28)   # hot primary coolant

# Controls
C_SLIDER_BG = ( 28,  32,  36)
C_SLIDER_FG = ( 48, 125, 195)   # blue fill (matches C_ACCENT)
C_BEZEL     = ( 22,  24,  26)


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

_arc_bg_cache: dict = {}   # keyed by radius → pre-rendered Surface

def draw_arc_gauge(screen, cx, cy, radius, value, lo, hi, font,
                   label="", unit="", bands=None, trip_hi=None):
    """Circular arc gauge — 270° sweep. Background ring is cached per radius."""
    START_DEG = 225
    SWEEP_DEG = 270
    frac = clamp((value - lo) / max(1e-9, hi - lo), 0.0, 1.0)

    r_outer = radius
    r_inner = int(radius * 0.68)
    r_tick  = int(radius * 0.82)
    r_tip   = int(radius * 0.62)
    n_segs  = min(90, max(40, radius))   # proportional but capped

    col = C_GREEN
    if bands:
        for blo, bhi, bcol in bands:
            if value >= blo:
                col = bcol

    # ── Background ring (cached) ──
    if radius not in _arc_bg_cache:
        sz = radius * 2 + 4
        bg_surf = pygame.Surface((sz, sz), pygame.SRCALPHA)
        rcx = rcy = radius + 2
        for i in range(n_segs):
            a0 = math.radians(START_DEG - i * SWEEP_DEG / n_segs)
            a1 = math.radians(START_DEG - (i + 1) * SWEEP_DEG / n_segs)
            sf = i / n_segs
            sc = (20, 50, 22) if sf < 0.67 else (50, 48, 10) if sf < 0.85 else (55, 16, 10)
            pts = [
                (rcx + int(r_outer * math.cos(a0)), rcy - int(r_outer * math.sin(a0))),
                (rcx + int(r_outer * math.cos(a1)), rcy - int(r_outer * math.sin(a1))),
                (rcx + int(r_inner * math.cos(a1)), rcy - int(r_inner * math.sin(a1))),
                (rcx + int(r_inner * math.cos(a0)), rcy - int(r_inner * math.sin(a0))),
            ]
            pygame.draw.polygon(bg_surf, sc, pts)
        _arc_bg_cache[radius] = bg_surf
    screen.blit(_arc_bg_cache[radius], (cx - radius - 2, cy - radius - 2))

    # ── Filled arc (live value) ──
    n_fill = max(1, int(frac * n_segs))
    for i in range(n_fill):
        a0 = math.radians(START_DEG - i * SWEEP_DEG / n_segs)
        a1 = math.radians(START_DEG - (i + 1) * SWEEP_DEG / n_segs)
        sf = i / n_segs
        fc = C_GREEN if sf < 0.67 else C_YELLOW if sf < 0.85 else C_RED
        pts = [
            (cx + int(r_outer * math.cos(a0)), cy - int(r_outer * math.sin(a0))),
            (cx + int(r_outer * math.cos(a1)), cy - int(r_outer * math.sin(a1))),
            (cx + int(r_inner * math.cos(a1)), cy - int(r_inner * math.sin(a1))),
            (cx + int(r_inner * math.cos(a0)), cy - int(r_inner * math.sin(a0))),
        ]
        pygame.draw.polygon(screen, fc, pts)

    # ── Tick marks ──
    for ti in range(6):
        ta = math.radians(START_DEG - ti * SWEEP_DEG / 5)
        x0 = cx + int(r_tick * math.cos(ta));  y0 = cy - int(r_tick * math.sin(ta))
        x1 = cx + int(r_outer * math.cos(ta)); y1 = cy - int(r_outer * math.sin(ta))
        pygame.draw.line(screen, (90, 112, 145), (x0, y0), (x1, y1), _s(2))

    # ── Trip line ──
    if trip_hi is not None:
        tf = clamp((trip_hi - lo) / max(1e-9, hi - lo), 0, 1)
        ta = math.radians(START_DEG - tf * SWEEP_DEG)
        x0 = cx + int(r_inner * math.cos(ta)); y0 = cy - int(r_inner * math.sin(ta))
        x1 = cx + int(r_outer * math.cos(ta)); y1 = cy - int(r_outer * math.sin(ta))
        pygame.draw.line(screen, C_RED, (x0, y0), (x1, y1), _s(3))

    # ── Needle ──
    needle_a = math.radians(START_DEG - frac * SWEEP_DEG)
    nx = cx + int(r_tip * math.cos(needle_a))
    ny = cy - int(r_tip * math.sin(needle_a))
    pygame.draw.line(screen, (240, 242, 255), (cx, cy), (nx, ny), _s(2))
    pygame.draw.circle(screen, (40, 50, 72), (cx, cy), _s(7))
    pygame.draw.circle(screen, (170, 185, 215), (cx, cy), _s(5))

    # ── Digital readout ──
    if font:
        val_s = f"{value:.1f}" if abs(value) < 100 else f"{value:.0f}"
        vs = font.render(val_s, True, col)
        screen.blit(vs, (cx - vs.get_width() // 2, cy + int(radius * 0.22)))
        if unit:
            us = font.render(unit, True, C_TEXT_DIM)
            screen.blit(us, (cx - us.get_width() // 2, cy + int(radius * 0.22) + vs.get_height()))
        if label:
            ls = font.render(label, True, C_TEXT_DIM)
            screen.blit(ls, (cx - ls.get_width() // 2, cy - radius - ls.get_height() - _s(2)))


def draw_valve(screen, cx, cy, radius, open_frac, color=C_BORDER):
    """Butterfly valve symbol — two triangles forming a bowtie, gap proportional to opening."""
    r = max(2, radius)
    gap = int(r * 0.3 * open_frac)
    pygame.draw.circle(screen, (14, 18, 28), (cx, cy), r + _s(2))
    pygame.draw.circle(screen, color, (cx, cy), r + _s(2), 1)
    # Left plate
    pygame.draw.polygon(screen, color, [(cx - r, cy), (cx - gap, cy - r//2), (cx - gap, cy + r//2)])
    # Right plate
    pygame.draw.polygon(screen, color, [(cx + r, cy), (cx + gap, cy - r//2), (cx + gap, cy + r//2)])


def draw_pump_symbol(screen, cx, cy, radius, running, color=C_BORDER):
    """IEC pump symbol — circle with a chord indicating impeller."""
    r = max(3, radius)
    bg = (18, 26, 42) if running else (12, 14, 18)
    pygame.draw.circle(screen, bg, (cx, cy), r)
    pygame.draw.circle(screen, color, (cx, cy), r, _s(2))
    # Impeller line (chord at 45°)
    dx = int(r * 0.7)
    pygame.draw.line(screen, color, (cx, cy), (cx + dx, cy - dx), _s(2))
    if running:
        # Small rotation dots
        for a in range(0, 360, 90):
            ax = cx + int((r - _s(3)) * math.cos(math.radians(a)))
            ay = cy - int((r - _s(3)) * math.sin(math.radians(a)))
            pygame.draw.circle(screen, color, (ax, ay), _s(2))


_scanline_surf: pygame.Surface | None = None

def draw_scanlines(screen, w, h, alpha=18):
    global _scanline_surf
    if _scanline_surf is None or _scanline_surf.get_size() != (w, h):
        _scanline_surf = pygame.Surface((w, h), pygame.SRCALPHA)
        _scanline_surf.fill((0, 0, 0, 0))
        for row in range(0, h, 2):
            pygame.draw.line(_scanline_surf, (0, 0, 0, alpha), (0, row), (w, row))
    screen.blit(_scanline_surf, (0, 0))


def draw_grid_bg(screen, x, y, w, h, spacing=40):
    """Subtle dot-grid background for the body area."""
    sp = _s(spacing)
    col = (22, 26, 30)
    for gx in range(x, x + w, sp):
        pygame.draw.line(screen, col, (gx, y), (gx, y + h))
    for gy in range(y, y + h, sp):
        pygame.draw.line(screen, col, (x, gy), (x + w, gy))


def draw_corner_brackets(screen, rect, color, arm=10, thickness=1):
    x, y, w, h = rect.x, rect.y, rect.width, rect.height
    a = _s(arm)
    for bx, by, dx, dy in [(x, y, 1, 1), (x+w, y, -1, 1), (x, y+h, 1, -1), (x+w, y+h, -1, -1)]:
        pygame.draw.line(screen, color, (bx, by), (bx + dx*a, by), thickness)
        pygame.draw.line(screen, color, (bx, by), (bx, by + dy*a), thickness)


def draw_panel(screen, x, y, w, h, title=None, font=None, border_color=None):
    bc = border_color or C_BORDER
    r  = _s(3)
    pygame.draw.rect(screen, C_PANEL, (x, y, w, h), border_radius=r)
    pygame.draw.rect(screen, bc, (x, y, w, h), 1, border_radius=r)
    # Corner brackets in accent color
    draw_corner_brackets(screen, pygame.Rect(x, y, w, h), C_ACCENT, arm=8, thickness=1)
    if title and font:
        lbl     = font.render(title, True, C_TEXT_HDR)
        title_h = lbl.get_height()
        pad     = _s(5)
        strip_h = title_h + _s(10)
        pygame.draw.rect(screen, C_PANEL_HDR, (x + 1, y + 1, w - 2, strip_h), border_radius=r)
        sep_y = y + pad + title_h + _s(5)
        screen.blit(lbl, (x + _s(10), y + pad))
        # Accent underline on header
        pygame.draw.line(screen, C_ACCENT, (x + 1, sep_y), (x + w - 1, sep_y))
        return sep_y + _s(6)
    return y + _s(8)


def draw_led(screen, cx, cy, color, radius=7):
    r = _s(radius)
    # Glow halo
    glow = pygame.Surface((r*6, r*6), pygame.SRCALPHA)
    gr, gg, gb = color
    pygame.draw.circle(glow, (gr, gg, gb, 40), (r*3, r*3), r*3)
    pygame.draw.circle(glow, (gr, gg, gb, 90), (r*3, r*3), r*2)
    screen.blit(glow, (cx - r*3, cy - r*3))
    pygame.draw.circle(screen, (18, 22, 32), (cx, cy), r + _s(1))
    pygame.draw.circle(screen, color, (cx, cy), r)
    hi = tuple(min(255, c + 100) for c in color)
    pygame.draw.circle(screen, hi, (cx - r // 3, cy - r // 3), max(1, r // 3))


def draw_trend(screen, x, y, w, h, values, color, y_min, y_max, label="", font=None):
    r = _s(3)
    pygame.draw.rect(screen, (8, 10, 12), (x, y, w, h), border_radius=r)
    pygame.draw.rect(screen, C_BORDER,   (x, y, w, h), 1, border_radius=r)
    # Grid lines (4 horizontal, 3 vertical)
    for i in range(1, 4):
        gy = y + int(h * i / 4)
        pygame.draw.line(screen, (22, 26, 30), (x + _s(2), gy), (x + w - _s(2), gy))
    for i in range(1, 4):
        gx = x + int(w * i / 4)
        pygame.draw.line(screen, (18, 22, 26), (gx, y + _s(2)), (gx, y + h - _s(2)))
    if label and font:
        screen.blit(font.render(label, True, C_TEXT_DIM), (x + _s(5), y + _s(4)))
        for val_s, vy in [(f"{y_max:.3g}", y + _s(4)), (f"{y_min:.3g}", y + h - font.get_height() - _s(2))]:
            screen.blit(font.render(val_s, True, (60, 75, 100)),
                        (x + w - font.size(val_s)[0] - _s(4), vy))
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
    """Vertical bar meter — industrial style with tick marks and color bands."""
    r4 = _s(4)
    # Outer bezel
    pygame.draw.rect(screen, (16, 20, 30), (x - _s(2), y - _s(2), w + _s(4), h + _s(4)), border_radius=r4)
    # Inner face (darker)
    pygame.draw.rect(screen, (8, 11, 18), (x, y, w, h), border_radius=r4)
    pygame.draw.rect(screen, C_BORDER,   (x, y, w, h), 1, border_radius=r4)

    frac  = clamp((value - lo) / max(1e-9, hi - lo), 0.0, 1.0)
    bar_h = int(frac * (h - _s(4)))
    bar_y = y + h - _s(2) - bar_h

    if bands:
        col = C_GREEN
        for band_lo, band_hi, band_col in bands:
            if value >= band_lo:
                col = band_col
    else:
        col = C_GREEN

    if bar_h > 0:
        pygame.draw.rect(screen, col, (x + _s(2), bar_y, w - _s(4), bar_h), border_radius=_s(2))
        # Bright top edge on the fill bar
        pygame.draw.line(screen, tuple(min(255, c + 60) for c in col),
                         (x + _s(2), bar_y), (x + w - _s(2), bar_y))

    # Tick marks (5 major, left side)
    for ti in range(6):
        ty2 = y + int((h - _s(4)) * (1.0 - ti / 5.0)) + _s(2)
        tick_w = _s(5) if ti % 5 == 0 else _s(3)
        pygame.draw.line(screen, (70, 85, 110), (x, ty2), (x + tick_w, ty2))

    if trip_hi is not None:
        ty2 = y + h - _s(2) - int(clamp((trip_hi - lo) / max(1e-9, hi - lo), 0, 1) * (h - _s(4)))
        pygame.draw.line(screen, C_RED,    (x - _s(2), ty2), (x + w + _s(2), ty2), _s(2))
    if trip_lo is not None:
        ty2 = y + h - _s(2) - int(clamp((trip_lo - lo) / max(1e-9, hi - lo), 0, 1) * (h - _s(4)))
        pygame.draw.line(screen, C_ORANGE, (x - _s(2), ty2), (x + w + _s(2), ty2), _s(2))

    if label and font:
        lbl = font.render(label, True, C_TEXT_DIM)
        screen.blit(lbl, (x + w // 2 - lbl.get_width() // 2, y - lbl.get_height() - _s(2)))
    if font:
        val_s = f"{value:.1f}" if abs(value) < 1000 else f"{value:.0f}"
        vs = font.render(val_s, True, C_TEXT)
        screen.blit(vs, (x + w // 2 - vs.get_width() // 2, y + h + _s(2)))
        if unit:
            us = font.render(unit, True, C_TEXT_DIM)
            screen.blit(us, (x + w // 2 - us.get_width() // 2, y + h + _s(2) + vs.get_height()))


def draw_dcs_box(screen, x, y, w, h, label, value_str, col, font_lbl, font_val):
    """DCS-style digital readout box — industrial Foxboro/Yokogawa look."""
    r = _s(2)
    # Outer bezel
    pygame.draw.rect(screen, (20, 22, 26), (x, y, w, h), border_radius=r)
    # Inner recessed face
    inset = _s(2)
    pygame.draw.rect(screen, (8, 10, 12), (x + inset, y + inset, w - inset*2, h - inset*2))
    # Border
    pygame.draw.rect(screen, C_BORDER, (x, y, w, h), 1, border_radius=r)
    # Tag label strip
    strip_h = font_lbl.get_height() + _s(4)
    pygame.draw.rect(screen, (10, 12, 15), (x + inset, y + inset, w - inset*2, strip_h))
    lbl = font_lbl.render(label, True, C_TEXT_DIM)
    screen.blit(lbl, (x + _s(6), y + _s(3)))
    # Value area
    val_area_y = y + inset + strip_h
    val_area_h = h - inset - strip_h - inset
    pygame.draw.rect(screen, (6, 8, 10), (x + inset, val_area_y, w - inset*2, val_area_h))
    val = font_val.render(value_str, True, col)
    screen.blit(val, (x + w - val.get_width() - _s(6), y + h - val.get_height() - _s(4)))


# ── P&ID schematics ──────────────────────────────────────────────────────────

# Animated flow dots: list of [segment_id, position 0-1]
_flow_dots: list = []
_flow_initialized = False

def _init_flow_dots(n_segs: int):
    global _flow_dots, _flow_initialized
    _flow_dots = [[i % n_segs, i / n_segs] for i in range(n_segs * 3)]
    _flow_initialized = True

def _update_flow_dots(n_segs: int, flow_rate: float, dt: float):
    speed = flow_rate * 0.4 * dt
    for dot in _flow_dots:
        dot[1] += speed
        if dot[1] >= 1.0:
            dot[1] -= 1.0

def _draw_dots_on_segment(screen, p0, p1, seg_id: int, color, n_segs: int):
    """Draw animated flow dots along the line from p0 to p1."""
    dx = p1[0] - p0[0]; dy = p1[1] - p0[1]
    length = math.hypot(dx, dy)
    if length < 1:
        return
    for dot in _flow_dots:
        if dot[0] % n_segs != seg_id % n_segs:
            continue
        t = dot[1]
        px = int(p0[0] + dx * t)
        py = int(p0[1] + dy * t)
        r = max(1, _s(3))
        bright = tuple(min(255, c + 80) for c in color)
        pygame.draw.circle(screen, bright, (px, py), r)


def _pipe(screen, pts, color, w=3):
    if len(pts) >= 2:
        pygame.draw.lines(screen, color, False, pts, max(1, int(w * SCALE)))


def _box(screen, x, y, w, h, label, font, col=C_BORDER, bg=(18, 24, 36), alarm=False, blink=False):
    bg_use = (60, 10, 10) if alarm and blink else (30, 5, 5) if alarm else bg
    r = _s(5)
    pygame.draw.rect(screen, bg_use, (x, y, w, h), border_radius=r)
    # Beveled inner highlight for depth
    pygame.draw.line(screen, tuple(min(255, c + 30) for c in bg_use),
                     (x + r, y + 1), (x + w - r, y + 1))
    pygame.draw.rect(screen, col if not alarm else C_RED, (x, y, w, h), _s(2), border_radius=r)
    if label and font:
        lbl = font.render(label, True, C_TEXT if not alarm else C_RED)
        screen.blit(lbl, (x + w // 2 - lbl.get_width() // 2, y + h // 2 - lbl.get_height() // 2))


def draw_pid_pwr(screen, x, y, w, h, snap, font, blink, flow_rate=1.0):
    """Simplified PWR P&ID schematic — left: primary loop, right: secondary loop."""
    fuel_alarm  = snap.fuel_temp_k > 1200
    press_alarm = snap.pressure_mpa > 16.5
    frac     = clamp(snap.power_fraction, 0, 1)

    # Primary loop components — left column
    rv_x  = x + _s(24);  rv_y = y + _s(50);  rv_w = _s(100); rv_h = _s(160)
    pz_x  = x + _s(30);  pz_y = y + _s(8);   pz_w = _s(70);  pz_h = _s(40)
    rcp_x = x + _s(210); rcp_y = y + h - _s(55); rcp_w = _s(70); rcp_h = _s(44)
    sg_x  = x + _s(210); sg_y = y + _s(30);   sg_w = _s(70);  sg_h = _s(160)

    # Secondary loop components — right column
    tb_x  = x + _s(360); tb_y = y + _s(30);  tb_w = _s(65);  tb_h = _s(55)
    cd_x  = x + _s(360); cd_y = y + h - _s(80); cd_w = _s(65); cd_h = _s(50)

    # ── Primary loop boxes ──
    _box(screen, rv_x, rv_y, rv_w, rv_h, "", font, alarm=fuel_alarm, blink=blink)
    # Fuel glow (bottom portion)
    glow_h = int((rv_h - _s(30)) * frac)
    if glow_h > 0:
        glow_col = (min(255, 60 + int(180*frac)), int(80*(1-frac)), 0)
        pygame.draw.rect(screen, glow_col,
                         (rv_x + _s(5), rv_y + rv_h - glow_h - _s(5), rv_w - _s(10), glow_h),
                         border_radius=_s(2))
    # Fuel rods (leave top 20px clear for label)
    rod_top = rv_y + _s(22)
    for ci in range(3):
        cx2 = rv_x + _s(12 + ci * 26)
        ch_col = (min(255, 60 + int(180*frac)), 30, 10) if frac > 0.05 else (40, 50, 70)
        pygame.draw.rect(screen, ch_col, (cx2, rod_top, _s(12), rv_h - _s(30)), border_radius=_s(2))
    # Label drawn last (on top of everything inside box)
    lbl_rv = font.render("REACTOR", True, C_RED if fuel_alarm else C_TEXT_DIM)
    screen.blit(lbl_rv, (rv_x + rv_w // 2 - lbl_rv.get_width() // 2, rv_y + _s(4)))

    _box(screen, pz_x, pz_y, pz_w, pz_h, "PZR", font,
         col=C_RED if press_alarm else (80, 60, 160), alarm=press_alarm, blink=blink)
    _box(screen, rcp_x, rcp_y, rcp_w, rcp_h, "RCP", font, col=(70, 110, 180))
    _box(screen, sg_x, sg_y, sg_w, sg_h, "S/G", font, col=(60, 100, 180))

    # ── Secondary loop boxes ──
    tb_x  = x + _s(370); tb_y = y + _s(30);  tb_w = _s(65);  tb_h = _s(55)
    cd_x  = x + _s(370); cd_y = y + h - _s(80); cd_w = _s(65); cd_h = _s(50)
    _box(screen, tb_x, tb_y, tb_w, tb_h, "TURB", font, col=(60, 140, 80))
    _box(screen, cd_x, cd_y, cd_w, cd_h, "COND", font, col=(50, 90, 120))

    # ── Primary loop pipes ──
    hot_col  = (200, 80, 40) if snap.coolant_temp_k > 560 else C_CYAN_PIPE
    cold_col = C_BLUE_PIPE
    mid_x = (rv_x + rv_w + sg_x) // 2   # horizontal run between reactor and SG
    # Hot leg: reactor right → across → down into SG top
    _pipe(screen, [
        (rv_x + rv_w, rv_y + _s(30)),
        (mid_x, rv_y + _s(30)),
        (mid_x, sg_y + _s(15)),
        (sg_x, sg_y + _s(15)),
    ], hot_col, 4)
    # Cold leg: SG bottom → down → across → up into reactor bottom
    cold_y = rcp_y + rcp_h // 2
    _pipe(screen, [
        (sg_x, sg_y + sg_h - _s(15)),
        (mid_x, sg_y + sg_h - _s(15)),
        (mid_x, cold_y),
        (rcp_x + rcp_w, cold_y),
    ], cold_col, 4)
    _pipe(screen, [
        (rcp_x, cold_y),
        (rv_x + rv_w // 2, cold_y),
        (rv_x + rv_w // 2, rv_y + rv_h),
    ], cold_col, 4)
    # Pressurizer surge line
    _pipe(screen, [(pz_x + pz_w // 2, pz_y + pz_h), (rv_x + rv_w // 2, rv_y + _s(15))], (130, 80, 200), 2)

    # ── Secondary loop pipes ──
    steam_col = C_GRAY_PIPE if snap.turbine_valve < 0.1 else C_CYAN_PIPE
    # Main steam: SG right → turbine left
    _pipe(screen, [
        (sg_x + sg_w, sg_y + _s(20)),
        (tb_x, tb_y + tb_h // 2),
    ], steam_col, 3)
    draw_valve(screen, sg_x + sg_w + _s(20), sg_y + _s(20), _s(7), snap.turbine_valve, steam_col)
    # Exhaust: turbine bottom → condenser top
    _pipe(screen, [
        (tb_x + tb_w // 2, tb_y + tb_h),
        (tb_x + tb_w // 2, cd_y),
    ], C_GRAY_PIPE, 3)
    # Feedwater: condenser left → across bottom → up → SG bottom
    fw_y_run = cd_y + cd_h // 2
    _pipe(screen, [
        (cd_x, fw_y_run),
        (sg_x + sg_w, fw_y_run),
        (sg_x + sg_w, sg_y + sg_h - _s(15)),
    ], C_BLUE_PIPE, 3)
    # Feedwater pump symbol on feedwater line
    fw_mid_x = (cd_x + sg_x + sg_w) // 2
    draw_pump_symbol(screen, fw_mid_x, fw_y_run, _s(10), snap.feedwater_inventory > 0.15, C_BLUE_PIPE)
    # RCP pump symbol
    draw_pump_symbol(screen, rcp_x + rcp_w // 2, rcp_y + rcp_h // 2, _s(11),
                     snap.flow > 0.1, (70, 110, 200))

    # ── Animated flow dots on primary loop ──
    mid_x2 = (rv_x + rv_w + sg_x) // 2
    cold_y2 = rcp_y + rcp_h // 2
    N = 8
    _draw_dots_on_segment(screen, (rv_x + rv_w, rv_y + _s(30)),      (mid_x2, rv_y + _s(30)),        0, hot_col,  N)
    _draw_dots_on_segment(screen, (mid_x2, rv_y + _s(30)),           (mid_x2, sg_y + _s(15)),         1, hot_col,  N)
    _draw_dots_on_segment(screen, (mid_x2, sg_y + _s(15)),           (sg_x, sg_y + _s(15)),           2, hot_col,  N)
    _draw_dots_on_segment(screen, (sg_x, sg_y + sg_h - _s(15)),      (mid_x2, sg_y + sg_h - _s(15)), 3, cold_col, N)
    _draw_dots_on_segment(screen, (mid_x2, sg_y + sg_h - _s(15)),    (mid_x2, cold_y2),               4, cold_col, N)
    _draw_dots_on_segment(screen, (mid_x2, cold_y2),                  (rcp_x + rcp_w, cold_y2),       5, cold_col, N)
    _draw_dots_on_segment(screen, (rcp_x, cold_y2),                   (rv_x + rv_w // 2, cold_y2),    6, cold_col, N)
    _draw_dots_on_segment(screen, (rv_x + rv_w // 2, cold_y2),        (rv_x + rv_w // 2, rv_y + rv_h),7, cold_col, N)

    # ── Annotations ──
    def small_val(sx, sy, txt, col=C_TEXT):
        screen.blit(font.render(txt, True, col), (sx, sy))
    fh = font.get_height()
    small_val(rv_x - _s(2), rv_y + _s(4),
              f"{snap.power_fraction*100:.0f}%", C_RED if frac > 1.1 else C_GREEN)
    small_val(rv_x - _s(2), rv_y + _s(4) + fh,
              f"{snap.fuel_temp_k:.0f}K", C_RED if snap.fuel_temp_k > 1200 else C_TEXT_DIM)
    small_val(pz_x + pz_w + _s(3), pz_y + _s(4),
              f"{snap.pressure_mpa:.2f}MPa", C_RED if press_alarm else C_TEXT_DIM)
    small_val(tb_x + tb_w + _s(3), tb_y + _s(4), f"{snap.steam_temp_k:.0f}K", C_TEXT_DIM)
    small_val(cd_x + cd_w + _s(3), cd_y + _s(4), f"{snap.condenser_temp_k:.0f}K", C_TEXT_DIM)


def draw_pid_bwr(screen, x, y, w, h, snap, font, blink):
    """Simplified BWR P&ID schematic."""
    frac       = clamp(snap.power_fraction, 0, 1)
    fuel_alarm = snap.fuel_temp_k > 1200

    rv_x = x + _s(70);  rv_y = y + _s(40);  rv_w = _s(100); rv_h = _s(200)
    tb_x = x + _s(290); tb_y = y + _s(40);  tb_w = _s(70);  tb_h = _s(60)
    cd_x = x + _s(290); cd_y = y + _s(160); cd_w = _s(70);  cd_h = _s(55)
    rp_x = x + _s(70);  rp_y = y + h - _s(60); rp_w = _s(70); rp_h = _s(44)

    _box(screen, rv_x, rv_y, rv_w, rv_h, "", font, alarm=fuel_alarm, blink=blink)
    lbl = font.render("REACTOR", True, C_TEXT)
    screen.blit(lbl, (rv_x + rv_w // 2 - lbl.get_width() // 2, rv_y + _s(6)))
    # Two-phase region + boiling glow
    two_phase_h = int(rv_h * 0.4 * frac)
    if two_phase_h > 0:
        pygame.draw.rect(screen, C_CYAN_PIPE,
                         (rv_x + _s(4), rv_y + rv_h - two_phase_h - _s(4), rv_w - _s(8), two_phase_h),
                         border_radius=_s(3))
    # Fuel channels
    for ci in range(4):
        cx2 = rv_x + _s(10 + ci * 20)
        ch_col = (min(255, 50 + int(180*frac)), 30, 10) if frac > 0.05 else (40, 50, 70)
        pygame.draw.rect(screen, ch_col, (cx2, rv_y + two_phase_h + _s(20), _s(10), rv_h - two_phase_h - _s(40)),
                         border_radius=_s(2))

    _box(screen, tb_x, tb_y, tb_w, tb_h, "TURB", font, col=(60, 140, 80))
    steam_col = C_GRAY_PIPE if snap.turbine_valve < 0.1 else C_CYAN_PIPE
    _pipe(screen, [(rv_x + rv_w, rv_y + _s(20)), (tb_x, tb_y + _s(30))], steam_col, 4)

    _box(screen, cd_x, cd_y, cd_w, cd_h, "COND", font, col=(50, 90, 120))
    _pipe(screen, [(tb_x + tb_w // 2, tb_y + tb_h), (tb_x + tb_w // 2, cd_y)], C_GRAY_PIPE, 3)
    # Feedwater: condenser left → down → across bottom → up → reactor right
    fw_bot = cd_y + cd_h + _s(14)
    _pipe(screen, [
        (cd_x, cd_y + cd_h // 2),
        (rv_x + rv_w + _s(10), cd_y + cd_h // 2),
        (rv_x + rv_w + _s(10), rv_y + rv_h - _s(40)),
        (rv_x + rv_w, rv_y + rv_h - _s(40)),
    ], C_BLUE_PIPE, 3)
    draw_pump_symbol(screen, (cd_x + rv_x + rv_w + _s(10)) // 2, cd_y + cd_h // 2,
                     _s(9), snap.feedwater_inventory > 0.15, C_BLUE_PIPE)

    _box(screen, rp_x, rp_y, rp_w, rp_h, "RECIRC", font, col=(70, 110, 180))
    _pipe(screen, [(rv_x + rv_w // 2, rv_y + rv_h), (rv_x + rv_w // 2, rp_y),
                   (rp_x + rp_w, rp_y + rp_h // 2), (rv_x, rv_y + rv_h - _s(60))], C_BLUE_PIPE, 3)
    # Recirc pump symbol
    draw_pump_symbol(screen, rp_x + rp_w // 2, rp_y + rp_h // 2, _s(10), snap.flow > 0.1, (70, 110, 200))

    vf    = snap.void_fraction
    vbw   = _s(18)
    vbar_x = rv_x - _s(28)
    vbar_h = int(rv_h * vf)
    pygame.draw.rect(screen, C_BORDER,   (vbar_x, rv_y, vbw, rv_h), border_radius=_s(3))
    if vbar_h > 0:
        pygame.draw.rect(screen, C_CYAN_PIPE, (vbar_x, rv_y + rv_h - vbar_h, vbw, vbar_h), border_radius=_s(2))
    screen.blit(font.render(f"α={vf:.2f}", True, C_CYAN_PIPE), (vbar_x - _s(4), rv_y + rv_h + _s(4)))

    def small_val(sx, sy, txt, col=C_TEXT):
        screen.blit(font.render(txt, True, col), (sx, sy))
    fh = font.get_height()
    small_val(rv_x + rv_w + _s(4), rv_y + _s(10),
              f"{snap.power_fraction*100:.0f}%", C_RED if frac > 1.1 else C_GREEN)
    small_val(rv_x + rv_w + _s(4), rv_y + _s(10) + fh,
              f"{snap.pressure_mpa:.2f}MPa", C_TEXT)
    small_val(cd_x + cd_w + _s(4), cd_y + _s(10), f"{snap.condenser_temp_k:.0f}K", C_TEXT)


def draw_pid_rbmk(screen, x, y, w, h, snap, font, blink):
    """Simplified RBMK P&ID schematic (channel-tube style)."""
    frac       = clamp(snap.power_fraction, 0, 1)
    fuel_alarm = snap.fuel_temp_k > 1200

    gx = x + _s(50);  gy = y + _s(30);  gw = _s(160); gh = _s(200)
    hd_x = x + _s(270); hd_y = y + _s(30); hd_w = _s(80); hd_h = _s(50)
    tb_x = x + _s(420); tb_y = y + _s(30); tb_w = _s(70); tb_h = _s(60)
    cd_x = x + _s(420); cd_y = y + _s(160);cd_w = _s(70); cd_h = _s(55)
    pp_x = x + _s(50);  pp_y = y + h - _s(70); pp_w = _s(70); pp_h = _s(44)

    # Graphite moderator block
    pygame.draw.rect(screen, (30, 30, 20), (gx, gy, gw, gh), border_radius=_s(6))
    pygame.draw.rect(screen, (80, 80, 50), (gx, gy, gw, gh), _s(2), border_radius=_s(6))
    glbl = font.render("GRAPHITE", True, (120, 120, 80))
    screen.blit(glbl, (gx + gw // 2 - glbl.get_width() // 2, gy + _s(6)))

    # Fuel channels (3 representative)
    for i, ch_rel in enumerate([_s(30), _s(70), _s(110)]):
        ch_x  = gx + ch_rel
        ch_col = (min(255, 50 + int(200 * frac)), 30, 10) if frac > 0.1 else (50, 60, 80)
        pygame.draw.rect(screen, ch_col, (ch_x, gy + _s(20), _s(22), gh - _s(40)), border_radius=_s(3))
        void_h = int((gh - _s(40)) * snap.void_fraction)
        if void_h > 0:
            pygame.draw.rect(screen, C_CYAN_PIPE, (ch_x + _s(2), gy + _s(20), _s(18), void_h), border_radius=_s(2))

    _box(screen, hd_x, hd_y, hd_w, hd_h, "DRUM", font, col=(70, 100, 160))
    _box(screen, tb_x, tb_y, tb_w, tb_h, "TURB", font, col=(60, 140, 80))
    steam_col = C_GRAY_PIPE if snap.turbine_valve < 0.1 else C_CYAN_PIPE
    _pipe(screen, [(gx + gw, gy + _s(40)), (hd_x, hd_y + _s(25))], steam_col, 3)
    _pipe(screen, [(hd_x + hd_w, hd_y + _s(25)), (tb_x, tb_y + _s(30))], steam_col, 3)

    _box(screen, cd_x, cd_y, cd_w, cd_h, "COND", font, col=(50, 90, 120))
    _pipe(screen, [(tb_x + tb_w // 2, tb_y + tb_h), (tb_x + tb_w // 2, cd_y)], C_GRAY_PIPE, 3)
    _pipe(screen, [(cd_x, cd_y + cd_h // 2), (gx, gy + gh - _s(40))], C_BLUE_PIPE, 3)

    _box(screen, pp_x, pp_y, pp_w, pp_h, "PUMP", font, col=(70, 110, 180))
    _pipe(screen, [(gx + gw // 2, gy + gh), (gx + gw // 2, pp_y),
                   (pp_x + pp_w, pp_y + pp_h // 2), (gx, gy + gh - _s(60))], C_BLUE_PIPE, 3)
    draw_pump_symbol(screen, pp_x + pp_w // 2, pp_y + pp_h // 2, _s(11), snap.flow > 0.1, (70, 110, 200))

    def small_val(sx, sy, txt, col=C_TEXT):
        screen.blit(font.render(txt, True, col), (sx, sy))
    fh = font.get_height()
    small_val(gx + gw + _s(4), gy + _s(10),
              f"{snap.power_fraction*100:.0f}%", C_RED if frac > 1.1 else C_GREEN)
    small_val(gx + gw + _s(4), gy + _s(10) + fh,
              f"void={snap.void_fraction:.2f}", C_CYAN_PIPE)
    small_val(gx + gw + _s(4), gy + _s(10) + fh * 2,
              f"{snap.fuel_temp_k:.0f}K", C_RED if fuel_alarm else C_TEXT)


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
    tab_h = _s(34)
    tab_w = w // len(TABS)
    rects = []
    # Full background strip
    pygame.draw.rect(screen, (12, 14, 17), (0, y, w, tab_h))
    pygame.draw.line(screen, C_BORDER, (0, y + tab_h - 1), (w, y + tab_h - 1))
    for i, (fkey, name) in enumerate(TABS):
        bx = i * tab_w
        active = i == active_tab
        rect = pygame.Rect(bx, y, tab_w, tab_h)
        rects.append(rect)
        if active:
            pygame.draw.rect(screen, (22, 28, 36), rect)
            # Left/right borders only for active tab
            pygame.draw.line(screen, C_BORDER, (bx, y), (bx, y + tab_h))
            pygame.draw.line(screen, C_BORDER, (bx + tab_w - 1, y), (bx + tab_w - 1, y + tab_h))
            # Bold accent underline
            pygame.draw.rect(screen, C_ACCENT, (bx + 1, y + tab_h - _s(3), tab_w - 2, _s(3)))
        else:
            # Subtle top line on hover hint
            pygame.draw.line(screen, C_BORDER, (bx, y), (bx + tab_w, y))
        fkey_s = font.render(fkey, True, C_ACCENT if active else C_TEXT_DIM)
        name_s = font.render(f"  {name}", True, C_TEXT if active else C_TEXT_DIM)
        total_w = fkey_s.get_width() + name_s.get_width()
        tx = bx + tab_w // 2 - total_w // 2
        ty = y + tab_h // 2 - fkey_s.get_height() // 2
        screen.blit(fkey_s, (tx, ty))
        screen.blit(name_s, (tx + fkey_s.get_width(), ty))
    return rects


def draw_screen_overview(screen, snap, fonts, layout, blink_fast, blink_slow, controls,
                         scram_msg, scram_msg_timer, ic_msg, ic_msg_timer, tw_label,
                         p_data, f_data, r_data, d_data):
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy_body, cw, rp_x, rp_w, body_h = layout

    # Center: P&ID schematic (top) + readouts (bottom)
    pid_h     = _s(310)
    readout_h = body_h - pid_h - _s(8)

    pid_y = cy_body
    ro_y  = pid_y + pid_h + _s(8)

    draw_panel(screen, cx, pid_y, cw, pid_h, "PLANT OVERVIEW", font_sm)
    rt = snap.reactor_type
    if rt == "PWR":
        draw_pid_pwr(screen, cx + _s(8), pid_y + _s(28), cw - _s(16), pid_h - _s(36), snap, font_sm, blink_fast, snap.flow)
    elif rt == "BWR":
        draw_pid_bwr(screen, cx + _s(8), pid_y + _s(28), cw - _s(16), pid_h - _s(36), snap, font_sm, blink_fast)
    else:
        draw_pid_rbmk(screen, cx + _s(8), pid_y + _s(28), cw - _s(16), pid_h - _s(36), snap, font_sm, blink_fast)

    # Right: bar meters + trends
    bm_h = _s(200)
    trend_y = cy_body + bm_h + _s(8)
    trend_h = body_h - bm_h - _s(8)

    draw_panel(screen, rp_x, cy_body, rp_w, bm_h, "INSTRUMENTS", font_sm)

    bm_area_x = rp_x + _s(12)
    bm_area_w = rp_w - _s(24)
    n_bars     = 4
    bar_w      = (bm_area_w - (n_bars - 1) * _s(8)) // n_bars
    bar_h_px   = bm_h - _s(56)

    pf = snap.power_fraction * 100
    pnom = 15.5 if snap.reactor_type == "PWR" else 7.0

    # Power: large circular arc gauge (left slot)
    gauge_r = min(bar_h_px // 2, bar_w) - _s(4)
    gauge_cx = bm_area_x + gauge_r + _s(4)
    gauge_cy = cy_body + _s(32) + bar_h_px // 2
    draw_arc_gauge(screen, gauge_cx, gauge_cy, gauge_r, pf, 0, 130, font_sm,
                   label="Power", unit="%",
                   bands=[(90, 110, C_YELLOW), (110, 130, C_RED)],
                   trip_hi=120.0)

    # Remaining bar meters (pressure, fuel temp, void) shifted right
    bar_defs = [
        ("Pressure", "MPa", snap.pressure_mpa, 0, pnom * 1.15,
         [(pnom * 1.04, pnom * 1.15, C_RED)], None, pnom * 1.08),
        ("Fuel T", "K",   snap.fuel_temp_k, 400, 1800,
         [(900, 1200, C_YELLOW), (1200, 1800, C_RED)], None, 1400.0),
        ("Void", "",      snap.void_fraction, 0, 1.0,
         [(0.4, 0.7, C_YELLOW), (0.7, 1.0, C_RED)], None, None),
    ]
    gauge_slot_w = gauge_r * 2 + _s(16)
    bar_start_x  = bm_area_x + gauge_slot_w
    for i, (lbl, unit, val, lo, hi, bands, tlo, thi) in enumerate(bar_defs):
        bx = bar_start_x + i * (bar_w + _s(8))
        draw_bar_meter(screen, bx, cy_body + _s(32), bar_w, bar_h_px, val, lo, hi,
                       font_sm, label=lbl, unit=unit, bands=bands, trip_lo=tlo, trip_hi=thi)

    # Readout text panel (below P&ID)
    ry = draw_panel(screen, cx, ro_y, cw, readout_h, "PLANT STATUS", font_sm)

    def rline(label, val_str, unit, col):
        nonlocal ry
        screen.blit(font.render(label, True, C_TEXT_DIM), (cx + _s(12), ry))
        vs  = font_md.render(val_str, True, col)
        uw  = font_sm.size(unit)[0]
        screen.blit(vs,  (cx + cw - uw - vs.get_width() - _s(18), ry - 1))
        screen.blit(font_sm.render(unit, True, C_TEXT_DIM), (cx + cw - uw - _s(10), ry + _s(2)))
        ry += _s(26)

    def rsep():
        nonlocal ry
        pygame.draw.line(screen, C_BORDER, (cx + _s(10), ry + _s(2)), (cx + cw - _s(10), ry + _s(2)))
        ry += _s(10)

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
    GW = rp_w - _s(20)
    GX = rp_x + _s(10)
    GAP = _s(5)
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
    bm_y   = cy2 + _s(20)
    bm_h2  = body_h - (cy2 - cy) - _s(140)
    bar_w2 = _s(60)
    bar_gap = _s(30)

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
    bx_start = cx + _s(40)
    for i, (lbl, unit, val, lo, hi, bands, tlo, thi) in enumerate(defs):
        bx = bx_start + i * (bar_w2 + bar_gap)
        draw_bar_meter(screen, bx, bm_y, bar_w2, bm_h2, val, lo, hi,
                       font_sm, label=lbl, unit=unit, bands=bands, trip_lo=tlo, trip_hi=thi)

    # DCS readout grid — width computed dynamically so boxes always fit
    box_y   = cy2 + _s(20) + bm_h2 + _s(40)
    n_boxes = 8
    box_gap = _s(10)
    bw3 = (cw + rp_w + _s(8) - _s(20) - box_gap * (n_boxes - 1)) // n_boxes
    bh3 = _s(52)
    porv_str = "OPEN" if snap.porv_open else "SHUT"
    eccs_str = "ACTIVE" if snap.eccs_actuated else "STANDBY"
    box_defs = [
        ("Core power", f"{snap.power_fraction*100:.2f} %",  C_GREEN if snap.power_fraction < 1.1 else C_RED),
        ("Pressure",   f"{snap.pressure_mpa:.3f} MPa",      C_GREEN if snap.pressure_mpa < pnom*1.04 else C_YELLOW),
        ("Fuel temp",  f"{snap.fuel_temp_k:.1f} K",         C_GREEN if snap.fuel_temp_k < 1100 else C_ORANGE),
        ("Coolant T",  f"{snap.coolant_temp_k:.1f} K",      C_GREEN if snap.coolant_temp_k < 580 else C_YELLOW),
        ("RCP speed",  f"{supervisor.bop.omega_rcp*100:.1f} %",
         C_GREEN if supervisor.bop.omega_rcp > 0.90 else C_ORANGE),
        ("PORV",       porv_str,  C_ORANGE if snap.porv_open else C_GREEN),
        ("ECCS",       eccs_str,  C_RED if snap.eccs_actuated else C_TEXT_DIM),
        ("Reactivity", f"{snap.reactivity:+.5f}",           C_GREEN if abs(snap.reactivity) < 0.001 else C_YELLOW),
    ]
    for i, (lbl3, val3, col3) in enumerate(box_defs):
        bx3 = cx + _s(10) + i * (bw3 + box_gap)
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
    bm_h2  = body_h - (cy2 - cy) - _s(140)
    bar_w2 = _s(65)
    bar_gap = _s(24)
    bm_y   = cy2 + _s(20)
    for i, (lbl, unit, val, lo, hi, bands, tlo, thi) in enumerate(bar_defs):
        bx = cx + _s(40) + i * (bar_w2 + bar_gap)
        draw_bar_meter(screen, bx, bm_y, bar_w2, bm_h2, val, lo, hi,
                       font_sm, label=lbl, unit=unit, bands=bands, trip_lo=tlo, trip_hi=thi)

    box_y   = cy2 + _s(20) + bm_h2 + _s(40)
    n_boxes2 = 6
    box_gap2 = _s(12)
    bw3 = (cw + rp_w + _s(8) - _s(20) - box_gap2 * (n_boxes2 - 1)) // n_boxes2
    bh3 = _s(52)
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
    for i, (lbl3, val3, col3) in enumerate(box_defs):
        bx3 = cx + _s(10) + i * (bw3 + box_gap2)
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
    TILE_W = _s(220); TILE_H = _s(44); TILE_GAP = _s(6)
    cols = max(1, (full_w - _s(20)) // (TILE_W + TILE_GAP))
    for idx, alarm in enumerate(all_alarms):
        col_i = idx % cols
        row_i = idx // cols
        tx = cx + _s(10) + col_i * (TILE_W + TILE_GAP)
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
        r4 = _s(4)
        pygame.draw.rect(screen, bg, (tx, ty, TILE_W, TILE_H), border_radius=r4)
        pygame.draw.rect(screen, bc, (tx, ty, TILE_W, TILE_H), 1, border_radius=r4)
        tag = "TRIP" if alarm.is_trip else "WARN" if alarm.state == "unack" else "ACK"
        tag_col = C_RED if alarm.is_trip else C_ORANGE if alarm.state == "unack" else C_GREEN
        screen.blit(font_sm.render(tag, True, tag_col), (tx + _s(5), ty + _s(4)))
        if alarm.first_out:
            fo = font_sm.render("FO", True, (255, 220, 0))
            screen.blit(fo, (tx + TILE_W - fo.get_width() - _s(5), ty + _s(4)))
        msg_lbl = font_sm.render(alarm.message[:28], True, C_TEXT)
        screen.blit(msg_lbl, (tx + _s(5), ty + TILE_H - msg_lbl.get_height() - _s(4)))

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
    noise  = snap.channel_noise
    row_h  = _s(65)
    col_name_w  = _s(280)
    col_val_w   = _s(200)
    col_vote_w  = _s(100)
    col_ch_w    = _s(160)
    for ri, (param_name, value, unit, warn_sp, trip_sp) in enumerate(params):
        py = cy2 + ri * (row_h + _s(8))
        screen.blit(font_md.render(param_name, True, C_TEXT_HDR), (cx + _s(10), py + _s(4)))
        val_str = f"{value:.3f} {unit}"
        warn_ok = value < warn_sp if warn_sp else True
        col = C_GREEN if warn_ok else (C_YELLOW if not trip_sp or value < trip_sp else C_RED)
        screen.blit(font_md.render(val_str, True, col), (cx + col_name_w, py + _s(4)))

        ch_vals = [value * (1.0 + n) for n in noise]
        votes = sum(1 for v in ch_vals if trip_sp and v > trip_sp)
        vote_col = C_RED if votes >= 2 else C_YELLOW if votes == 1 else C_TEXT_DIM
        screen.blit(font_sm.render(f"Votes: {votes}/3", True, vote_col),
                    (cx + col_name_w + col_val_w, py + _s(4)))
        for ch_i, ch_val in enumerate(ch_vals):
            ch_x = cx + col_name_w + col_val_w + col_vote_w + ch_i * col_ch_w
            ch_str = f"Ch {chr(65+ch_i)}: {ch_val:.3f}"
            ch_col = C_GREEN
            if warn_sp and ch_val > warn_sp:
                ch_col = C_YELLOW
            if trip_sp and ch_val > trip_sp:
                ch_col = C_RED
            screen.blit(font.render(ch_str, True, ch_col), (ch_x, py + _s(4)))

        setpt_str = (f"Warn: {warn_sp:.2f} {unit}  Trip: {trip_sp:.2f} {unit}"
                     if trip_sp else f"Warn: {warn_sp:.2f} {unit}")
        screen.blit(font_sm.render(setpt_str, True, C_TEXT_DIM), (cx + _s(10), py + _s(30)))
        pygame.draw.line(screen, C_BORDER,
                         (cx + _s(10), py + row_h + _s(2)),
                         (cx + full_w - _s(10), py + row_h + _s(2)))


def draw_screen_reactivity(screen, snap, fonts, layout, p_data, f_data):
    """F6: Reactivity budget breakdown + xenon/iodine trends."""
    font_sm, font, font_md, font_lg, font_xl = fonts
    cx, cy, cw, rp_x, rp_w, body_h = layout
    full_w = cw + rp_w + 8

    cy2 = draw_panel(screen, cx, cy, full_w, body_h, "REACTIVITY BUDGET & POISON TRENDS", font_sm)

    # ── RBMK ORM / SKALA panel ───────────────────────────────────────────────
    if snap.reactor_type == "RBMK":
        orm_live  = snap.orm
        orm_skala = snap.skala_orm
        age_min   = snap.skala_age_s / 60.0

        if orm_live < 15:   orm_col = C_RED
        elif orm_live < 26: orm_col = C_YELLOW
        else:               orm_col = C_GREEN

        orm_panel_y = cy2
        orm_panel_h = _s(90)
        pw = full_w - _s(16)

        pygame.draw.rect(screen, (10, 14, 22), (cx + _s(8), orm_panel_y, pw, orm_panel_h), border_radius=_s(2))
        pygame.draw.rect(screen, (40, 60, 100),(cx + _s(8), orm_panel_y, pw, _s(18)), border_radius=_s(2))
        pygame.draw.rect(screen, (40, 60, 100),(cx + _s(8), orm_panel_y, pw, orm_panel_h), 1, border_radius=_s(2))
        screen.blit(font_sm.render("СКАЛА — OPERATIONAL REACTIVITY MARGIN  (ORM / ОРМ)", True, (180, 200, 240)),
                    (cx + _s(16), orm_panel_y + _s(3)))

        # Live ORM gauge box
        gx2 = cx + _s(16);  gy2 = orm_panel_y + _s(22);  gw2 = _s(200);  gh2 = _s(62)
        pygame.draw.rect(screen, (6, 8, 14), (gx2, gy2, gw2, gh2), border_radius=_s(2))
        pygame.draw.rect(screen, orm_col,    (gx2, gy2, gw2, gh2), _s(2), border_radius=_s(2))
        screen.blit(font_sm.render("LIVE  ORM", True, C_TEXT_DIM), (gx2 + _s(6), gy2 + _s(4)))
        val_surf = font_xl.render(f"{orm_live:5.1f}", True, orm_col)
        screen.blit(val_surf, (gx2 + gw2 - val_surf.get_width() - _s(8), gy2 + _s(22)))
        screen.blit(font_sm.render("rod equiv.", True, C_TEXT_DIM), (gx2 + _s(6), gy2 + gh2 - _s(16)))

        # SKALA printout box
        sx2 = cx + _s(226);  sy2 = orm_panel_y + _s(22);  sw2 = _s(210);  sH2 = _s(62)
        age_border = C_RED if age_min > 30 else C_YELLOW if age_min > 10 else C_BORDER
        pygame.draw.rect(screen, (6, 8, 14), (sx2, sy2, sw2, sH2), border_radius=_s(2))
        pygame.draw.rect(screen, age_border, (sx2, sy2, sw2, sH2), _s(2), border_radius=_s(2))
        screen.blit(font_sm.render("SKALA PRINTOUT", True, C_TEXT_DIM), (sx2 + _s(6), sy2 + _s(4)))
        sk_surf = font_lg.render(f"{orm_skala:5.1f}", True, C_TEXT)
        screen.blit(sk_surf, (sx2 + sw2 - sk_surf.get_width() - _s(8), sy2 + _s(22)))
        age_col = C_RED if age_min > 30 else C_YELLOW if age_min > 10 else C_TEXT_DIM
        screen.blit(font_sm.render(f"age: {age_min:.0f} min", True, age_col), (sx2 + _s(6), sy2 + sH2 - _s(16)))

        # ORM scale bar
        bx0 = cx + _s(446);  by0 = orm_panel_y + _s(22);  bw0 = _s(460);  bh0 = _s(62)
        pygame.draw.rect(screen, (6, 8, 14), (bx0, by0, bw0, bh0), border_radius=_s(2))
        pygame.draw.rect(screen, C_BORDER,   (bx0, by0, bw0, bh0), 1, border_radius=_s(2))
        screen.blit(font_sm.render("ORM SCALE  (0 – 211 rods)", True, C_TEXT_DIM), (bx0 + _s(6), by0 + _s(4)))
        bar_inner_x = bx0 + _s(8);  bar_inner_y = by0 + _s(22)
        bar_inner_w = bw0 - _s(16); bar_inner_h = _s(18)
        pygame.draw.rect(screen, (20, 25, 35), (bar_inner_x, bar_inner_y, bar_inner_w, bar_inner_h))
        zone_scale = bar_inner_w / 211.0
        z15 = int(15 * zone_scale);  z26 = int(26 * zone_scale)
        pygame.draw.rect(screen, (80, 0, 0),  (bar_inner_x,       bar_inner_y, z15,             bar_inner_h))
        pygame.draw.rect(screen, (80, 60, 0), (bar_inner_x + z15, bar_inner_y, z26 - z15,       bar_inner_h))
        pygame.draw.rect(screen, (0, 60, 20), (bar_inner_x + z26, bar_inner_y, bar_inner_w-z26, bar_inner_h))
        needle_x = bar_inner_x + int(min(orm_live, 211) * zone_scale)
        pygame.draw.rect(screen, orm_col, (needle_x - _s(2), bar_inner_y - _s(2), _s(4), bar_inner_h + _s(4)))
        sk_needle_x = bar_inner_x + int(min(orm_skala, 211) * zone_scale)
        for dy in range(0, bar_inner_h, _s(4)):
            pygame.draw.line(screen, C_TEXT_DIM, (sk_needle_x, bar_inner_y + dy),
                             (sk_needle_x, min(bar_inner_y + dy + _s(2), bar_inner_y + bar_inner_h)))
        for lim, lbl_s, lim_col in [(15, "15", C_RED), (26, "26", C_YELLOW)]:
            lx = bar_inner_x + int(lim * zone_scale)
            pygame.draw.line(screen, lim_col, (lx, bar_inner_y - _s(3)), (lx, bar_inner_y + bar_inner_h + _s(3)), 1)
            lt = font_sm.render(lbl_s, True, lim_col)
            screen.blit(lt, (lx - lt.get_width() // 2, bar_inner_y + bar_inner_h + _s(5)))

        ref_x = cx + _s(916)
        screen.blit(font_sm.render("LIMITS:", True, C_TEXT_DIM), (ref_x, orm_panel_y + _s(22)))
        screen.blit(font_sm.render("< 26: min normal ops (authorisation req.)", True, C_YELLOW),  (ref_x, orm_panel_y + _s(38)))
        screen.blit(font_sm.render("< 15: ABSOLUTE MINIMUM — SHUTDOWN",         True, C_RED),     (ref_x, orm_panel_y + _s(54)))
        screen.blit(font_sm.render("  ~8: Chernobyl Unit 4, 01:22 26 Apr 1986", True, C_TEXT_DIM),(ref_x, orm_panel_y + _s(70)))

        cy2 += orm_panel_h + _s(8)

    # Stacked horizontal bar chart (reactivity components, qualitative)
    rho_total = snap.reactivity
    pf = snap.power_fraction
    pnom = 15.5 if snap.reactor_type == "PWR" else 7.0

    rho_decay = snap.decay_heat_fraction * 0.001
    # Approximate rod reactivity from rod_position via S-curve
    rp = snap.rod_position
    rho_rods_approx = -(3*rp**2 - 2*rp**3) * 0.05  # −rod_worth × w(rod_pos)

    components = [
        ("Total reactivity",  rho_total,         C_ACCENT),
        ("Rods (approx)",     rho_rods_approx,   C_CYAN_PIPE),
        ("Xenon (approx)",    -0.02 * (1 - pf),  C_YELLOW),
        ("Decay heat effect", rho_decay,          C_ORANGE),
    ]

    bar_y = cy2 + _s(20)
    bar_max_w = full_w - _s(80)
    zero_x = cx + full_w // 2

    screen.blit(font_md.render("Reactivity Components (qualitative)", True, C_TEXT_HDR), (cx + _s(10), bar_y))
    bar_y += _s(28)
    pygame.draw.line(screen, C_TEXT_DIM, (zero_x, bar_y), (zero_x, bar_y + len(components) * _s(38) + _s(10)))

    scale = bar_max_w / 2 / max(0.05, max(abs(v) for _, v, _ in components))
    for name, val, col in components:
        screen.blit(font_sm.render(f"{name}:  {val:+.6f}", True, C_TEXT), (cx + _s(10), bar_y))
        bar_y += _s(16)
        blen = int(abs(val) * scale)
        bx = zero_x if val >= 0 else zero_x - blen
        pygame.draw.rect(screen, col, (bx, bar_y, max(_s(2), blen), _s(14)), border_radius=_s(2))
        bar_y += _s(22)

    trend_y = bar_y + _s(20)
    trend_h  = body_h - (trend_y - cy)
    TW = (full_w - _s(24)) // 2
    draw_trend(screen, cx + _s(8), trend_y, TW, trend_h, list(p_data),
               C_GREEN, 0.0, 130.0, "Core power (%)", font_sm)
    draw_trend(screen, cx + TW + _s(16), trend_y, TW, trend_h, list(f_data),
               C_ORANGE, 400.0, 1800.0, "Fuel temperature (K)", font_sm)


# ── Main loop ────────────────────────────────────────────────────────────────

def main():
    global SCALE, W, H
    pygame.init()

    # On macOS Retina, request native pixel resolution via ALLOW_HIGHDPI.
    # screen.get_size() then returns physical pixels (~2× logical), so
    # SCALE = actual_W / 1920 ≈ 1.78 on a 27" 5K Retina — all _s() calls
    # produce sharp native-pixel coordinates automatically.
    info  = pygame.display.Info()
    log_w = max(1280, info.current_w - 20)
    log_h = max(720,  info.current_h - 60)
    _HIGHDPI = getattr(pygame, "ALLOW_HIGHDPI", 0x00002000)
    screen = pygame.display.set_mode((log_w, log_h), pygame.RESIZABLE | _HIGHDPI)
    pygame.display.set_caption("Nuclear Reactor Control Room")

    W, H  = screen.get_size()
    SCALE = W / 1920.0   # on Retina ≈ 1.78; on standard display ≈ 0.89

    def _make_fonts():
        return (
            pygame.font.SysFont("Menlo", _s(13)),
            pygame.font.SysFont("Menlo", _s(15)),
            pygame.font.SysFont("Menlo", _s(17), bold=True),
            pygame.font.SysFont("Menlo", _s(22), bold=True),
            pygame.font.SysFont("Menlo", _s(32), bold=True),
        )

    font_sm, font, font_md, font_lg, font_xl = _make_fonts()
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

    rod   = Slider("Rod position",       0.0)
    flow  = Slider("Primary flow",       1.0)
    valve = Slider("Turbine valve",      1.0)
    feed  = Slider("Feedwater valve",    0.7)
    press = Slider("Pressurizer heater", 0.5)
    bor   = Slider("Boration rate",      0.0)
    dil   = Slider("Dilution rate",      0.0)

    last = time.time()

    # ── Mouse interaction state ──────────────────────────────────────────────
    slider_dragging  = None   # Slider object currently being dragged
    slider_rects     = []     # [(pygame.Rect, Slider)]  — updated each frame
    tab_rects_prev   = []     # [pygame.Rect]            — from previous frame
    btn_rects_prev   = []     # [(pygame.Rect, str)]     — reactor selector buttons
    toggle_rects_prev= []     # [(pygame.Rect, callable)]

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
        rod.set(0.0); flow.set(1.0); valve.set(1.0); feed.set(0.7); press.set(0.5)
        bor.set(0.0); dil.set(0.0)

    running  = True
    snapshot = supervisor.step(0.0)
    effective_speed = 1

    while running:
        # ── Events ──────────────────────────────────────────────────────────
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                running = False
            elif ev.type == pygame.VIDEORESIZE:
                new_w = max(960, ev.w)
                new_h = max(540, ev.h)
                screen = pygame.display.set_mode((new_w, new_h), pygame.RESIZABLE)
                W, H = screen.get_size()
                SCALE = W / 1920.0
                _arc_bg_cache.clear()          # invalidate gauge ring cache
                if hasattr(draw_scanlines, '_surf'):
                    del draw_scanlines._surf   # invalidate scanline cache
                font_sm, font, font_md, font_lg, font_xl = _make_fonts()
                fonts = (font_sm, font, font_md, font_lg, font_xl)
            elif ev.type == pygame.MOUSEBUTTONDOWN and ev.button == 1:
                mx, my = ev.pos
                # tabs
                for i, rect in enumerate(tab_rects_prev):
                    if rect.collidepoint(mx, my):
                        active_tab = i
                # reactor selector buttons
                for rect, rt in btn_rects_prev:
                    if rect.collidepoint(mx, my):
                        select(rt)
                # sliders
                for rect, sl in slider_rects:
                    if rect.collidepoint(mx, my):
                        slider_dragging = sl
                        frac = clamp((mx - rect.x) / max(1, rect.w), 0.0, 1.0)
                        sl.set(sl.lo + frac * (sl.hi - sl.lo))
                # toggles
                for rect, action in toggle_rects_prev:
                    if rect.collidepoint(mx, my):
                        action()
            elif ev.type == pygame.MOUSEMOTION:
                if slider_dragging is not None and ev.buttons[0]:
                    mx = ev.pos[0]
                    for rect, sl in slider_rects:
                        if sl is slider_dragging:
                            frac = clamp((mx - rect.x) / max(1, rect.w), 0.0, 1.0)
                            sl.set(sl.lo + frac * (sl.hi - sl.lo))
            elif ev.type == pygame.MOUSEBUTTONUP and ev.button == 1:
                slider_dragging = None
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

        # ── Scaled layout constants (recomputed each frame; SCALE fixed after init) ──
        HEADER_H = _s(56)
        TAB_H    = _s(34)
        PAD      = _s(8)
        HELP_H   = _s(22)
        BTN_H    = _s(42);  BTN_Y = _s(7);  BTN_W = _s(162);  BTN_GAP = _s(170)

        if is_scram:
            scram_bg = (180, 15, 15) if blink_fast else (90, 8, 8)
            pygame.draw.rect(screen, scram_bg, (0, 0, W, HEADER_H))
            txt = font_xl.render("⚠  REACTOR SCRAM  ──  REACTOR SCRAM  ──  REACTOR SCRAM  ⚠", True, (255, 255, 255))
            screen.blit(txt, (W // 2 - txt.get_width() // 2, HEADER_H // 2 - txt.get_height() // 2))
        else:
            # DCS header bar
            pygame.draw.rect(screen, (8, 12, 20), (0, 0, W, HEADER_H))
            pygame.draw.rect(screen, (30, 50, 90), (0, HEADER_H - _s(2), W, _s(2)))

            # Reactor selector buttons (left)
            rt_info = [("PWR",  "PRESSURISED WATER",  (40, 100, 200)),
                       ("BWR",  "BOILING WATER",       (30, 150, 100)),
                       ("RBMK", "RBMK-1000",           (190, 65, 40))]
            btn_rects_prev = []
            for i, (rt, rt_long, base_col) in enumerate(rt_info):
                active = snapshot.reactor_type == rt
                bx = _s(12) + i * BTN_GAP
                btn_rects_prev.append((pygame.Rect(bx, BTN_Y, BTN_W, BTN_H), rt))
                bg = tuple(int(c * 0.6) for c in base_col) if active else (18, 24, 36)
                pygame.draw.rect(screen, bg,       (bx, BTN_Y, BTN_W, BTN_H), border_radius=_s(3))
                pygame.draw.rect(screen, base_col if active else (40, 52, 72),
                                 (bx, BTN_Y, BTN_W, BTN_H), 1, border_radius=_s(3))
                label_t = font_sm.render(rt_long, True, (200, 210, 225) if active else (55, 70, 90))
                num_t   = font_lg.render(rt,       True, (255, 255, 255) if active else (60, 80, 110))
                screen.blit(num_t,   (bx + _s(6), BTN_Y + _s(4)))
                screen.blit(label_t, (bx + _s(6) + num_t.get_width() + _s(6), BTN_Y + _s(10)))

            # Facility ID (centre)
            unit_names = {"PWR": "UNIT 1", "BWR": "UNIT 2", "RBMK": "UNIT 4"}
            fac_names  = {"PWR": "PRESSURISED WATER REACTOR PLANT", "BWR": "BOILING WATER REACTOR PLANT",
                          "RBMK": "V.I. LENIN NUCLEAR POWER STATION"}
            unit_t = font_lg.render(unit_names[snapshot.reactor_type], True, (180, 200, 230))
            fac_t  = font_sm.render(fac_names[snapshot.reactor_type],  True, (90, 110, 140))
            cx_hdr = W // 2
            screen.blit(unit_t, (cx_hdr - unit_t.get_width() // 2, BTN_Y))
            screen.blit(fac_t,  (cx_hdr - fac_t.get_width()  // 2, BTN_Y + unit_t.get_height() + _s(2)))

            # Sim time clock (HH:MM:SS) + speed multiplier
            st = int(snapshot.time)
            hh, rem = divmod(st, 3600)
            mm, ss  = divmod(rem, 60)
            speed_str = f"  ×{effective_speed}" if effective_speed > 1 else ""
            clock_str = f"{hh:02d}:{mm:02d}:{ss:02d}{speed_str}"
            t_surf  = font_lg.render(clock_str, True, (140, 165, 210))
            t_label = font_sm.render("SIM TIME", True, C_TEXT_DIM)
            t_x = cx_hdr + _s(300)
            screen.blit(t_label, (t_x, BTN_Y + _s(2)))
            screen.blit(t_surf,  (t_x, BTN_Y + t_label.get_height() + _s(2)))

            # Alarm/trip indicator (right)
            ind_w = _s(220);  ind_x = W - ind_w - _s(8);  ind_cx = ind_x + ind_w // 2
            if has_trips:
                col = (220, 40, 40) if blink_fast else (100, 20, 20)
                pygame.draw.rect(screen, col, (ind_x, BTN_Y, ind_w, BTN_H), border_radius=_s(3))
                t = font_lg.render("▶  TRIP ACTIVE", True, (255, 255, 255))
                screen.blit(t, (ind_cx - t.get_width() // 2, BTN_Y + (BTN_H - t.get_height()) // 2))
            elif has_alarms:
                col = C_ORANGE if blink_slow else (100, 55, 0)
                pygame.draw.rect(screen, col, (ind_x, BTN_Y, ind_w, BTN_H), border_radius=_s(3))
                t = font_lg.render("▶  ALARM", True, (255, 255, 255))
                screen.blit(t, (ind_cx - t.get_width() // 2, BTN_Y + (BTN_H - t.get_height()) // 2))
            else:
                pygame.draw.rect(screen, (14, 36, 14), (ind_x, BTN_Y, ind_w, BTN_H), border_radius=_s(3))
                pygame.draw.rect(screen, (30, 80, 30), (ind_x, BTN_Y, ind_w, BTN_H), 1, border_radius=_s(3))
                t = font_md.render("ALL SYSTEMS NORMAL", True, (50, 200, 80))
                screen.blit(t, (ind_cx - t.get_width() // 2, BTN_Y + (BTN_H - t.get_height()) // 2))

        # ── Reactor mode badge (between alarm indicator and right edge) ──────
        if not is_scram:
            pf_now = snapshot.power_fraction * 100
            if controls.scram or pf_now < 1.0:
                mode_str, mode_col = "SHUTDOWN", C_TEXT_DIM
            elif pf_now < 20:
                mode_str, mode_col = "STARTUP", C_YELLOW
            elif pf_now > 95:
                mode_str, mode_col = "FULL POWER", C_GREEN
            else:
                mode_str, mode_col = "POWER OPS", C_ACCENT
            mode_surf = font_sm.render(mode_str, True, mode_col)
            mode_x = W - _s(8) - mode_surf.get_width()
            mode_y = BTN_Y + BTN_H - mode_surf.get_height() - _s(3)
            screen.blit(mode_surf, (mode_x, mode_y))

        # ── Tab bar ──────────────────────────────────────────────────────────
        tab_rects_prev = draw_tab_bar(screen, active_tab, font, HEADER_H, W)

        # ── Body layout ──────────────────────────────────────────────────────
        Y0     = HEADER_H + TAB_H + PAD
        BODY_H = H - HEADER_H - TAB_H - HELP_H - PAD * 2

        LP_X, LP_W = PAD, _s(340)
        CP_X = LP_X + LP_W + PAD
        CP_W = _s(560)
        RP_X = CP_X + CP_W + PAD
        RP_W = W - RP_X - PAD

        layout = (CP_X, Y0, CP_W, RP_X, RP_W, BODY_H)

        # ── Subtle grid background on body area ──────────────────────────────
        draw_grid_bg(screen, 0, Y0, W, BODY_H, spacing=48)

        # ── Animated flow dots update ─────────────────────────────────────────
        if not _flow_initialized:
            _init_flow_dots(8)
        _update_flow_dots(8, snapshot.flow, dt * time_speed)

        # ── Left panel (always visible): Controls ────────────────────────────
        cy = draw_panel(screen, LP_X, Y0, LP_W, BODY_H, "CONTROLS", font_sm)

        sliders_for_rt = [rod, flow, valve, feed, press]
        if snapshot.reactor_type == "PWR":
            sliders_for_rt += [bor]

        SL_PAD = _s(10); SL_H = _s(16); SL_HIT_H = _s(32); LED_X = _s(22); LED_TXT = _s(40)
        new_slider_rects = []
        # Effective rod position may differ from commanded (startup permit clamp, SCRAM)
        effective_rod = snapshot.rod_position
        rod_clamped = abs(effective_rod - rod.v) > 0.01 and not controls.scram

        def _sl_unit(sl, v):
            if sl is rod:   return f"{v*100:.0f} %"
            if sl is flow:  return f"{v*100:.0f} %"
            if sl is valve: return f"{v*100:.0f} %"
            if sl is feed:  return f"{v*100:.0f} %"
            if sl is press: return f"{v*100:.0f} %"
            if sl is bor:   return f"{v*100:.0f} %"
            if sl is dil:   return f"{v*100:.0f} %"
            return f"{v*100:.0f} %"

        for sl in sliders_for_rt:
            is_dragging = sl is slider_dragging
            clamped = rod_clamped and sl is rod
            lbl_s = font_sm.render(sl.label, True, C_TEXT_HDR if is_dragging else
                                   C_ORANGE if clamped else C_TEXT_DIM)
            display_v = effective_rod if sl is rod else sl.v
            val_col = C_ORANGE if clamped else (C_TEXT if is_dragging else C_ACCENT)
            val_s = font_sm.render(_sl_unit(sl, display_v), True, val_col)
            screen.blit(lbl_s, (LP_X + SL_PAD, cy))
            screen.blit(val_s, (LP_X + LP_W - val_s.get_width() - SL_PAD, cy))
            cy += max(lbl_s.get_height(), val_s.get_height()) + _s(3)
            bx, bw = LP_X + SL_PAD, LP_W - SL_PAD * 2
            # Larger invisible hit rect centered on track (easier to click)
            hit_rect = pygame.Rect(bx, cy - (SL_HIT_H - SL_H)//2, bw, SL_HIT_H)
            new_slider_rects.append((hit_rect, sl))
            sl_rect = pygame.Rect(bx, cy, bw, SL_H)
            track_col = (35, 42, 48) if is_dragging else C_SLIDER_BG
            pygame.draw.rect(screen, track_col, sl_rect, border_radius=_s(4))
            fill_v = display_v if sl is rod else sl.v
            fill = int(bw * (fill_v - sl.lo) / max(1e-6, sl.hi - sl.lo))
            if fill > 0:
                fill_col = C_ORANGE if clamped else (C_TEXT_HDR if is_dragging else C_SLIDER_FG)
                pygame.draw.rect(screen, fill_col, (bx, cy, fill, SL_H), border_radius=_s(4))
            # thumb
            thumb_x = bx + fill - _s(4)
            thumb_col = (200, 140, 40) if clamped else ((235, 240, 245) if is_dragging else (130, 150, 170))
            pygame.draw.rect(screen, thumb_col,
                             (thumb_x, cy - _s(2), _s(8), SL_H + _s(4)), border_radius=_s(3))
            if clamped:
                clamp_lbl = font_sm.render("CLAMPED-PERMIT", True, C_ORANGE)
                screen.blit(clamp_lbl, (bx, cy + SL_H + _s(2)))
            cy += SL_H + (_s(22) if clamped else _s(14))
        slider_rects = new_slider_rects

        pygame.draw.line(screen, C_BORDER, (LP_X + _s(10), cy), (LP_X + LP_W - _s(10), cy))
        cy += _s(10)

        new_toggle_rects = []
        def _toggle_startup():   controls.startup_permit       = not controls.startup_permit
        def _toggle_turbtrip():  controls.turbine_trip         = not controls.turbine_trip
        def _toggle_pumpdeg():   controls.fault_pump_degraded  = not controls.fault_pump_degraded
        def _toggle_fwloss():    controls.fault_feedwater_loss = not controls.fault_feedwater_loss
        def _toggle_loca():
            LOCA_AREAS = [0.0, 0.001, 0.01]
            cur = controls.fault_loca_break_area
            idx_next = (LOCA_AREAS.index(min(LOCA_AREAS, key=lambda a: abs(a - cur))) + 1) % len(LOCA_AREAS)
            controls.fault_loca_break_area = LOCA_AREAS[idx_next]
        TGL_H = _s(28)
        for label, state, col_on, key, action in [
            ("Startup permit",  controls.startup_permit,            C_GREEN,  "P", _toggle_startup),
            ("Turbine trip",    controls.turbine_trip,              C_ORANGE, "T", _toggle_turbtrip),
            ("Pump degraded",   controls.fault_pump_degraded,       C_RED,    "Z", _toggle_pumpdeg),
            ("Feedwater fault", controls.fault_feedwater_loss,      C_RED,    "X", _toggle_fwloss),
            ("LOCA fault",      controls.fault_loca_break_area > 0, C_RED,    "I", _toggle_loca),
        ]:
            btn_rect = pygame.Rect(LP_X + _s(4), cy, LP_W - _s(8), TGL_H)
            new_toggle_rects.append((btn_rect, action))
            if state:
                r, g, b = col_on
                bg_col  = (int(r*0.15), int(g*0.15), int(b*0.15))
                brd_col = col_on
            else:
                bg_col  = C_PANEL
                brd_col = C_BORDER
            pygame.draw.rect(screen, bg_col,  btn_rect, border_radius=_s(3))
            pygame.draw.rect(screen, brd_col, btn_rect, 1, border_radius=_s(3))
            led_col = col_on if state else (45, 55, 65)
            draw_led(screen, LP_X + _s(18), cy + TGL_H//2, led_col)
            tc = col_on if state else C_TEXT_DIM
            screen.blit(font_sm.render(label, True, tc), (LP_X + _s(34), cy + (TGL_H - font_sm.get_height())//2))
            kt = font_sm.render(f"[{key}]", True, C_TEXT_DIM)
            screen.blit(kt, (LP_X + LP_W - kt.get_width() - _s(16), cy + (TGL_H - kt.get_height())//2))
            cy += TGL_H + _s(3)

        pygame.draw.line(screen, C_BORDER, (LP_X + _s(10), cy + _s(2)), (LP_X + LP_W - _s(10), cy + _s(2)))
        cy += _s(12)

        scram_led_x = LP_X + LED_X;  scram_txt_x = LP_X + _s(44)
        scram_row_h = _s(42)
        scram_btn_rect = pygame.Rect(LP_X + _s(4), cy, LP_W - _s(8), scram_row_h)
        if is_scram:
            scram_col = C_RED if blink_fast else C_ORANGE
            pygame.draw.rect(screen, (40, 8, 8), scram_btn_rect, border_radius=_s(3))
            pygame.draw.rect(screen, scram_col, scram_btn_rect, 1, border_radius=_s(3))
            draw_led(screen, scram_led_x, cy + _s(12), scram_col, radius=8)
            screen.blit(font_md.render("SCRAM ACTIVE",   True, scram_col),   (scram_txt_x, cy + _s(2)))
            screen.blit(font_sm.render("click: reset  R: restart", True, C_TEXT_DIM), (scram_txt_x, cy + _s(24)))
            def _do_scram_reset():
                nonlocal scram_msg, scram_msg_timer
                scram_msg = supervisor.reset_scram()
                scram_msg_timer = 4.0
            new_toggle_rects.append((scram_btn_rect, _do_scram_reset))
        else:
            scram_fire_rect = pygame.Rect(LP_X + _s(4), cy, LP_W - _s(8), scram_row_h)
            new_toggle_rects.append((scram_fire_rect, lambda: controls.__setattr__('scram', True)))
            pygame.draw.rect(screen, (10, 28, 10), scram_btn_rect, border_radius=_s(3))
            pygame.draw.rect(screen, (30, 80, 30), scram_btn_rect, 1, border_radius=_s(3))
            draw_led(screen, scram_led_x, cy + _s(12), C_GREEN, radius=8)
            screen.blit(font_md.render("REACTOR NORMAL",  True, C_GREEN),    (scram_txt_x, cy + _s(2)))
            screen.blit(font_sm.render("[SPACE] Emergency SCRAM", True, C_TEXT_DIM), (scram_txt_x, cy + _s(24)))
        toggle_rects_prev = new_toggle_rects
        cy += scram_row_h + _s(4)

        if scram_msg and scram_msg_timer > 0:
            mc = C_GREEN if "APPROVED" in scram_msg else C_ORANGE
            screen.blit(font_sm.render(scram_msg[:38], True, mc), (LP_X + _s(12), cy))
            cy += _s(18)
        if ic_msg and ic_msg_timer > 0:
            screen.blit(font_sm.render(ic_msg, True, C_ACCENT), (LP_X + _s(12), cy))
            cy += _s(18)

        # PORV / ECCS status (PWR only)
        if snapshot.reactor_type == "PWR":
            if snapshot.porv_open:
                draw_led(screen, LP_X + LED_X, cy + _s(8), C_ORANGE if blink_fast else (100, 55, 0))
                screen.blit(font_sm.render("PORV OPEN", True, C_ORANGE), (LP_X + LED_TXT, cy + _s(2)))
                cy += _s(20)
            if snapshot.eccs_actuated:
                draw_led(screen, LP_X + LED_X, cy + _s(8), C_RED if blink_fast else (80, 0, 0))
                screen.blit(font_sm.render("ECCS ACTUATED", True, C_RED), (LP_X + LED_TXT, cy + _s(2)))
                cy += _s(20)
            if snapshot.loca_area > 0:
                screen.blit(font_sm.render(f"LOCA {snapshot.loca_area*10000:.0f}cm²",
                                           True, C_RED), (LP_X + _s(12), cy))
                cy += _s(18)

        # Auto controller indicators (PWR only)
        if snapshot.reactor_type == "PWR":
            auto_rod_col  = C_ACCENT if supervisor.auto_rod.auto      else C_TEXT_DIM
            auto_pres_col = C_ACCENT if supervisor.auto_pressure.auto  else C_TEXT_DIM
            screen.blit(font_sm.render(f"[O] Auto rod:  {'AUTO' if supervisor.auto_rod.auto else 'MAN '}",
                                       True, auto_rod_col), (LP_X + _s(12), cy))
            cy += _s(16)
            screen.blit(font_sm.render(f"[M] Auto pzr:  {'AUTO' if supervisor.auto_pressure.auto else 'MAN '}",
                                       True, auto_pres_col), (LP_X + _s(12), cy))

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

        # ── CRT scanline overlay (applied last, over everything) ──────────────
        draw_scanlines(screen, W, H, alpha=12)

        # ── Mouse cursor feedback ─────────────────────────────────────────────
        mx, my = pygame.mouse.get_pos()
        _hover = (
            any(r.collidepoint(mx, my) for r in tab_rects_prev) or
            any(r.collidepoint(mx, my) for r, _ in btn_rects_prev) or
            any(r.collidepoint(mx, my) for r, _ in slider_rects) or
            any(r.collidepoint(mx, my) for r, _ in toggle_rects_prev) or
            scram_btn_rect.collidepoint(mx, my)
        )
        pygame.mouse.set_cursor(
            pygame.SYSTEM_CURSOR_HAND if _hover else pygame.SYSTEM_CURSOR_ARROW
        )

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
