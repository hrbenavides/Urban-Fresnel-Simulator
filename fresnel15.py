#!/usr/bin/env python3
"""
Simulador de elipsoide de Fresnel para enlaces de radio urbano.
Universidad Mayor de San Andrés — Facultad de Ingeniería
Univ. Henry Rafael Benavides Gutierrez

Versión con cálculo de potencia (FSPL + difracción + ganancias de antena,
sensibilidad del receptor y margen de seguridad por entorno),
optimizada, limitada a 10 km, con selección de unidades de ganancia: dB, dBi, dBd.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
import math

C_LIGHT = 3e8  # m/s
MAX_DIST_KM = 10.0   # límite para no considerar curvatura terrestre

# ══════════════════════════════════════════════════════════════════════════════
#  PALETA  (Blender Dark)
# ══════════════════════════════════════════════════════════════════════════════
_BG  = '#0d0d0d'
_PNL = '#161616'
_BDR = '#333340'
_GR1 = '#2a2a2a'
_GR2 = '#1e1e1e'
_TXT = '#e8e8e8'
_MUT = '#7a8a9a'
_ACC = '#5ba3ff'
_LOS = '#ffaa00'
_FK  = '#00ddff'
_FB  = '#ff4444'
_ANT = '#ffe055'

_PA = ['#bbd4ff', '#7aa4e8', '#5888cc', '#1a3a8a']
_PB = ['#bbffcc', '#66dd99', '#44bb77', '#166644']

_POBS = [
    ['#ffdc99', '#e8a84a', '#cc8830', '#995520'],
    ['#ffbbee', '#ee77cc', '#cc55aa', '#881144'],
    ['#dbbeff', '#aa66ee', '#8844cc', '#550088'],
    ['#99eeff', '#44ccee', '#22aacc', '#085868'],
    ['#ffaaaa', '#ee5555', '#cc3333', '#881111'],
    ['#aaffaa', '#55ee55', '#33cc33', '#117711'],
]

_OBS_HCOL = ['#ff9955', '#ff66cc', '#cc66ff', '#33eeff', '#ff5555', '#66ff66']
_OBS_NAMES = ['C', 'D', 'E', 'F', 'G', 'H']
_OBS_DEF = [
    dict(h=70, d=30,  off=10,   bw=15, rot=0),
    dict(h=55, d=50,  off=-8,  bw=12, rot=0),
    dict(h=65, d=65,  off=10,  bw=10, rot=0),
    dict(h=50, d=75,  off=-5,  bw=15, rot=0),
    dict(h=45, d=82,  off=8,   bw=12, rot=0),
    dict(h=52, d=88,  off=-3,  bw=10, rot=0),
]

MAX_OBS = 6

# ══════════════════════════════════════════════════════════════════════════════
#  GEOMETRÍA 2-D / 3-D (optimizada)
# ══════════════════════════════════════════════════════════════════════════════

def _verts(cx, cy, w, d, h, rot_deg=0.0):
    hw, hd = w / 2.0, d / 2.0
    offsets = [(-hw, -hd), (hw, -hd), (hw, hd), (-hw, hd)]
    th = math.radians(float(rot_deg))
    cth, sth = math.cos(th), math.sin(th)
    v_bot, v_top = [], []
    for dx, dy in offsets:
        rx = cth * dx - sth * dy
        ry = sth * dx + cth * dy
        v_bot.append([cx + rx, cy + ry, 0.0])
        v_top.append([cx + rx, cy + ry, float(h)])
    v = v_bot + v_top
    return {
        'top': [v[4], v[5], v[6], v[7]],
        'frt': [v[0], v[1], v[5], v[4]],
        'bk':  [v[2], v[3], v[7], v[6]],
        'lft': [v[0], v[3], v[7], v[4]],
        'rgt': [v[1], v[2], v[6], v[5]],
    }

def _point_in_poly(x, y, poly):
    inside = False
    n = len(poly)
    for i in range(n):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % n]
        if ((y0 > y) != (y1 > y)) and (
                x < (x1 - x0) * (y - y0) / (y1 - y0 + 1e-12) + x0):
            inside = not inside
    return inside

def _dist_point_to_segment_2d(px, py, x0, y0, x1, y1):
    dx, dy = x1 - x0, y1 - y0
    if dx == 0 and dy == 0:
        return math.hypot(px - x0, py - y0)
    t = ((px - x0) * dx + (py - y0) * dy) / (dx * dx + dy * dy)
    t = max(0.0, min(1.0, t))
    return math.hypot(px - (x0 + t * dx), py - (y0 + t * dy))

def _dist_point_to_polygon_2d(px, py, poly):
    if _point_in_poly(px, py, poly):
        return 0.0
    return min(
        _dist_point_to_segment_2d(px, py,
                                   poly[i][0], poly[i][1],
                                   poly[(i + 1) % len(poly)][0],
                                   poly[(i + 1) % len(poly)][1])
        for i in range(len(poly))
    )

def _closest_point_on_segment_2d(px, py, x0, y0, x1, y1):
    dx, dy = x1 - x0, y1 - y0
    if dx == 0 and dy == 0:
        return x0, y0
    t = max(0.0, min(1.0, ((px - x0) * dx + (py - y0) * dy) / (dx * dx + dy * dy)))
    return x0 + t * dx, y0 + t * dy

def _closest_point_on_polygon_2d(px, py, poly):
    best, best_d = None, float('inf')
    n = len(poly)
    for i in range(n):
        cx2, cy2 = _closest_point_on_segment_2d(
            px, py, poly[i][0], poly[i][1],
            poly[(i + 1) % n][0], poly[(i + 1) % n][1])
        d = math.hypot(px - cx2, py - cy2)
        if d < best_d:
            best_d = d
            best = (cx2, cy2)
    return best

def closest_point_to_prism_3d(px, py, pz, poly, hC):
    inside_xy = _point_in_poly(px, py, poly)
    cpx, cpy = (px, py) if inside_xy else _closest_point_on_polygon_2d(px, py, poly)
    cpz = min(max(pz, 0.0), float(hC))
    return np.array([cpx, cpy, cpz], float)

def dist_point_to_prism_3d(px, py, pz, poly, hC):
    inside_xy = _point_in_poly(px, py, poly)
    h2d = 0.0 if inside_xy else _dist_point_to_polygon_2d(px, py, poly)
    dz = max(0.0, pz - hC) if pz > hC else (max(0.0, -pz) if pz < 0.0 else 0.0)
    return dz if inside_xy else math.hypot(h2d, dz)

def find_closest_los_to_prism(A, B, poly, hC, n_coarse=200):
    d_3D = float(np.linalg.norm(B - A))
    if d_3D < 1e-9:
        return 0.5, 0.0, 0.0, 0.0, 0.0, (A + B) / 2.0

    ts = np.linspace(0.0, 1.0, n_coarse)
    min_d, t_star = float('inf'), 0.5
    for t in ts:
        P = A + t * (B - A)
        d = dist_point_to_prism_3d(P[0], P[1], P[2], poly, hC)
        if d < min_d:
            min_d, t_star = d, t

    step = 1.0 / n_coarse
    lo, hi = max(0.0, t_star - 2 * step), min(1.0, t_star + 2 * step)
    phi = (math.sqrt(5) - 1) / 2
    for _ in range(40):
        if hi - lo < 1e-9:
            break
        t1, t2 = hi - phi * (hi - lo), lo + phi * (hi - lo)
        P1, P2 = A + t1 * (B - A), A + t2 * (B - A)
        d1 = dist_point_to_prism_3d(P1[0], P1[1], P1[2], poly, hC)
        d2 = dist_point_to_prism_3d(P2[0], P2[1], P2[2], poly, hC)
        if d1 < d2:
            hi = t2
        else:
            lo = t1
    t_star = (lo + hi) / 2.0

    los_pt   = A + t_star * (B - A)
    min_dist = dist_point_to_prism_3d(los_pt[0], los_pt[1], los_pt[2], poly, hC)
    return t_star, d_3D, t_star * d_3D, (1.0 - t_star) * d_3D, min_dist, los_pt

def roof_antenna_position(cx, cy, w, rot_deg, choice):
    hw = w / 2.0
    offsets = [(-hw, -hw), (hw, -hw), (hw, hw), (-hw, hw)]
    ch = choice.lower()
    if ch.startswith('cent') or ch == 'center':
        return float(cx), float(cy)
    try:
        idx = int(ch[3:]) - 1 if ch.startswith('esq') else int(ch) - 1
    except Exception:
        idx = 0
    dx, dy = offsets[int(idx) % 4]
    th = math.radians(float(rot_deg))
    return float(cx + math.cos(th) * dx - math.sin(th) * dy), \
           float(cy + math.sin(th) * dx + math.cos(th) * dy)

# ══════════════════════════════════════════════════════════════════════════════
#  RENDERIZADO (optimizado)
# ══════════════════════════════════════════════════════════════════════════════

def draw_building(ax, cx, cy, w, d, h, pal, alpha=0.92, label=None, show_labels=True):
    rot_deg = 0.0
    if isinstance(d, dict):
        rot_deg = d.get('rot', 0.0)
        d = d.get('d', w)
    f = _verts(cx, cy, w, d, h, rot_deg=rot_deg)
    for name, fc in [('top', pal[0]), ('frt', pal[1]), ('bk', pal[2]),
                     ('lft', pal[2]), ('rgt', pal[1])]:
        ax.add_collection3d(Poly3DCollection(
            [f[name]], alpha=alpha, facecolor=fc,
            edgecolor=pal[3], linewidth=0.45, zorder=3))
    if label and show_labels:
        pad = max(3, h * 0.05)
        ax.text(cx, cy, h + pad, label,
                color='white', fontsize=12, fontweight='bold',
                ha='center', va='bottom', zorder=10)

def draw_antenna(ax, x, y, h, extra=8):
    ax.plot([x, x], [y, y], [h, h + extra], color=_ANT, lw=2.5, zorder=8)
    ax.scatter([x], [y], [h + extra], color=_ANT, s=60, depthshade=False, zorder=9)

def draw_cota_3d(ax, p1, p2, text, offset, color, fs=8, show_labels=True):
    p1, p2, off = np.array(p1, float), np.array(p2, float), np.array(offset, float)
    v1, v2 = p1 + off, p2 + off
    ax.plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]], color=color, lw=1.2, zorder=12)
    ax.plot([p1[0], v1[0]], [p1[1], v1[1]], [p1[2], v1[2]], color=color, lw=0.7, ls=':', alpha=0.6, zorder=12)
    ax.plot([p2[0], v2[0]], [p2[1], v2[1]], [p2[2], v2[2]], color=color, lw=0.7, ls=':', alpha=0.6, zorder=12)
    ax.scatter([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]], color=color, s=10, marker='o', depthshade=False, zorder=13)
    if show_labels:
        mid = (v1 + v2) / 2.0
        ax.text(mid[0], mid[1], mid[2], text,
                color='white', fontsize=fs, fontweight='bold',
                ha='center', va='center', zorder=15,
                bbox=dict(facecolor=_BG, edgecolor='none', alpha=0.80, pad=2))

def draw_grid(ax, cx, cy, d_ab):
    ext, spc, spc2 = max(d_ab * 2.0, 500), 50, 25
    x0 = np.floor((cx - ext / 2) / spc) * spc
    x1 = np.ceil( (cx + ext / 2) / spc) * spc
    y0 = np.floor((cy - ext / 2) / spc) * spc
    y1 = np.ceil( (cy + ext / 2) / spc) * spc
    for x in np.arange(x0, x1 + spc, spc):
        a = max(0.02, 0.30 * (1.0 - (abs(x - cx) / (ext / 2 + 1.0)) ** 0.55))
        ax.plot([x, x], [y0, y1], [0, 0], color=_GR1, alpha=a, lw=0.65, zorder=0)
    for y in np.arange(y0, y1 + spc, spc):
        a = max(0.02, 0.30 * (1.0 - (abs(y - cy) / (ext / 2 + 1.0)) ** 0.55))
        ax.plot([x0, x1], [y, y], [0, 0], color=_GR1, alpha=a, lw=0.65, zorder=0)
    half = ext * 0.38
    sx0, sx1, sy0, sy1 = cx - half, cx + half, cy - half, cy + half
    for x in np.arange(np.floor(sx0 / spc2) * spc2, sx1 + spc2, spc2):
        if abs(round(x) % spc) > 1:
            a = max(0.01, 0.11 * (1.0 - abs(x - cx) / (half + 1.0)))
            ax.plot([x, x], [sy0, sy1], [0, 0], color=_GR2, alpha=a, lw=0.3, zorder=0)
    for y in np.arange(np.floor(sy0 / spc2) * spc2, sy1 + spc2, spc2):
        if abs(round(y) % spc) > 1:
            a = max(0.01, 0.11 * (1.0 - abs(y - cy) / (half + 1.0)))
            ax.plot([sx0, sx1], [y, y], [0, 0], color=_GR2, alpha=a, lw=0.3, zorder=0)
    ax.plot([x0, x1], [0, 0], [0, 0], color='#cc4444', alpha=0.5, lw=1.8, zorder=1)
    ax.plot([0,  0],  [y0, y1], [0, 0], color='#4444cc', alpha=0.5, lw=1.8, zorder=1)

def draw_fresnel_wire(ax, A, B, lam, zf, color, nu=30, nv=14):
    A, B = np.array(A, float), np.array(B, float)
    d = float(np.linalg.norm(B - A))
    if d < 1e-9:
        return
    u = (B - A) / d
    ref = np.array([1, 0, 0], float) if abs(u[0]) < 0.85 else np.array([0, 1, 0], float)
    v_ = np.cross(ref, u); v_ /= np.linalg.norm(v_)
    w_ = np.cross(u, v_)
    t  = np.linspace(0.0, d, nu)
    th = np.linspace(0.0, 2 * np.pi, nv, endpoint=False)
    r  = np.where((t <= 0) | (t >= d), 0.0,
                  zf * np.sqrt(np.maximum(lam * t * (d - t) / d, 0.0)))
    Xw = np.zeros((nu, nv)); Yw = np.zeros((nu, nv)); Zw = np.zeros((nu, nv))
    for j in range(nu):
        c, rj = A + t[j] * u, r[j]
        for i in range(nv):
            pt = c + rj * (np.cos(th[i]) * v_ + np.sin(th[i]) * w_)
            Xw[j, i], Yw[j, i], Zw[j, i] = pt
    n_m = min(nv, 10); step = max(1, nv // n_m)
    for i in range(0, nv, step):
        ax.plot(Xw[:, i], Yw[:, i], Zw[:, i], color=color, alpha=0.42, lw=0.95, zorder=5)
    n_r = 16; step = max(1, nu // n_r)
    for j in range(0, nu, step):
        xi = np.append(Xw[j, :], Xw[j, 0])
        yi = np.append(Yw[j, :], Yw[j, 0])
        zi = np.append(Zw[j, :], Zw[j, 0])
        ring_a = 0.18 + 0.62 * np.sin(np.pi * (t[j] / d if d > 0 else 0.0))
        ax.plot(xi, yi, zi, color=color, alpha=ring_a, lw=0.72, zorder=5)

# ══════════════════════════════════════════════════════════════════════════════
#  FUNCIONES DE ATENUACIÓN POR DIFRACCIÓN (UIT-R P.526)
# ══════════════════════════════════════════════════════════════════════════════

def diffraction_loss(nu):
    if nu <= -0.7:
        return 0.0
    term = math.sqrt((nu - 0.1)**2 + 1) + nu - 0.1
    if term <= 0:
        return 0.0
    return 6.9 + 20.0 * math.log10(term)

def knife_edge_attenuation(clearance, r_fresnel, d1, d2, wavelength):
    if clearance >= r_fresnel:
        return 0.0
    H = r_fresnel - clearance
    if H <= 0:
        return 0.0
    if d1 <= 0 or d2 <= 0:
        return 0.0
    nu = H * math.sqrt(2.0 * (d1 + d2) / (wavelength * d1 * d2))
    return diffraction_loss(nu)

# ══════════════════════════════════════════════════════════════════════════════
#  APLICACIÓN PRINCIPAL
# ══════════════════════════════════════════════════════════════════════════════

class FresnelApp:
    def __init__(self, root):
        self.root = root
        root.title("Simulador Elipsoide de Fresnel — Radio Enlace Urbano v16")
        root.geometry("1440x880")
        root.configure(bg=_BG)

        self._vars = {}
        self._show_ell = tk.BooleanVar(value=True)
        self._xray_mode = tk.BooleanVar(value=True)
        self._show_labels = tk.BooleanVar(value=True)
        self._vars['uniform_scale'] = tk.BooleanVar(value=False)

        self._n_obs_var = tk.IntVar(value=1)
        self._obs_vars = []
        self._obs_container = None
        self._right_canvas = None
        self._right_inner_id = None

        self._setup_style()
        self._build_layout()
        self.root.bind('<Return>', self.simular)
        self.root.bind('<KP_Enter>', self.simular)
        self._drag_info = {}
        self._update_cam()
        self.simular()

    # ── Estilo ttk ─────────────────────────────────────────────────────────────

    def _setup_style(self):
        s = ttk.Style()
        s.theme_use('clam')
        for w in ('TFrame', 'TLabel', 'TCheckbutton'):
            try:
                s.configure(w, background=_PNL, foreground=_TXT, font=('Consolas', 9))
            except Exception:
                pass
        s.configure('TEntry', fieldbackground='#0a0a0a', foreground=_TXT,
                    insertcolor=_TXT, bordercolor=_BDR)
        s.configure('TScale', background=_PNL, troughcolor='#2a2a2a')
        s.configure('TButton', background='#262626', foreground=_TXT,
                    relief='flat', padding=(6, 4))
        s.map('TButton', background=[('active', '#353535')])

    # ── Layout principal ───────────────────────────────────────────────────────

    def _build_layout(self):
        self.root.columnconfigure(0, weight=0, minsize=350)
        self.root.columnconfigure(1, weight=1)
        self.root.columnconfigure(2, weight=0, minsize=350)
        self.root.rowconfigure(0, weight=1)

        left   = tk.Frame(self.root, bg=_PNL, width=350)
        center = tk.Frame(self.root, bg=_BG)
        right  = tk.Frame(self.root, bg=_PNL, width=350)

        left  .grid(row=0, column=0, sticky='nswe', padx=(8, 4), pady=8)
        center.grid(row=0, column=1, sticky='nswe', padx=4,      pady=8)
        right .grid(row=0, column=2, sticky='nswe', padx=(4, 8), pady=8)

        left .grid_propagate(False)
        right.grid_propagate(False)

        self._build_panel(left)
        self._build_canvas(center)
        self._build_side_panel(right)

    # ── Helpers de formulario ─────────────────────────────────────────────────

    def _section(self, parent, title, color):
        tk.Label(parent, text=f"▸ {title}", bg=_PNL, fg=color,
                 font=('Consolas', 9, 'bold')).pack(anchor='w', padx=10, pady=(10, 0))
        tk.Frame(parent, bg=color, height=1).pack(fill=tk.X, padx=10, pady=(1, 4))
        body = tk.Frame(parent, bg=_PNL)
        body.pack(fill=tk.X, padx=16, pady=2)
        body.columnconfigure(0, weight=1)
        return body

    def _field(self, body, label, key, default, row):
        tk.Label(body, text=label, bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=row, column=0, sticky='w', pady=2)
        var = tk.StringVar(value=str(default))
        self._vars[key] = var
        tk.Entry(body, textvariable=var, bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
                 font=('Consolas', 10), width=7, relief='flat',
                 highlightbackground=_BDR, highlightcolor=_ACC, highlightthickness=1
                 ).grid(row=row, column=1, sticky='e', pady=2, padx=(4, 0))

    def _dropdown(self, body, label, key, options, default, row):
        tk.Label(body, text=label, bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=row, column=0, sticky='w', pady=2)
        var = tk.StringVar(value=str(default))
        self._vars[key] = var
        om = ttk.OptionMenu(body, var, default, *options,
                            command=lambda _=None: self.simular())
        om.config(width=8)
        om.grid(row=row, column=1, sticky='e', pady=2, padx=(4, 0))

    def _obs_field(self, body, label, var, row):
        tk.Label(body, text=label, bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=row, column=0, sticky='w', pady=2)
        tk.Entry(body, textvariable=var, bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
                 font=('Consolas', 10), width=7, relief='flat',
                 highlightbackground=_BDR, highlightcolor=_ACC, highlightthickness=1
                 ).grid(row=row, column=1, sticky='e', pady=2, padx=(4, 0))

    # ── Panel izquierdo (scrollable) ──────────────────────────────────────────

    def _build_panel(self, parent):
        vscroll = tk.Scrollbar(parent, orient='vertical', bg='#1e1e1e',
                               troughcolor=_BDR, activebackground='#444', width=10)
        vscroll.pack(side=tk.RIGHT, fill=tk.Y)

        self._left_canvas = tk.Canvas(parent, bg=_PNL,
                                       yscrollcommand=vscroll.set,
                                       highlightthickness=0)
        self._left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vscroll.config(command=self._left_canvas.yview)

        inner = tk.Frame(self._left_canvas, bg=_PNL)
        self._left_inner_id = self._left_canvas.create_window(
            (0, 0), window=inner, anchor='nw')

        inner.bind('<Configure>',
                   lambda e: self._left_canvas.configure(
                       scrollregion=self._left_canvas.bbox('all')))
        self._left_canvas.bind('<Configure>',
                               lambda e: self._left_canvas.itemconfig(
                                   self._left_inner_id, width=e.width))

        def _pw(ev):
            if   getattr(ev, 'num', None) == 4: self._left_canvas.yview_scroll(-1, 'units')
            elif getattr(ev, 'num', None) == 5: self._left_canvas.yview_scroll( 1, 'units')
            elif getattr(ev, 'delta', 0):
                self._left_canvas.yview_scroll(int(-ev.delta / 120), 'units')

        def _bind_pw(e):
            self._left_canvas.bind_all('<MouseWheel>', _pw)
            self._left_canvas.bind_all('<Button-4>',   _pw)
            self._left_canvas.bind_all('<Button-5>',   _pw)

        def _unbind_pw(e):
            self._left_canvas.unbind_all('<MouseWheel>')
            self._left_canvas.unbind_all('<Button-4>')
            self._left_canvas.unbind_all('<Button-5>')

        inner.bind('<Enter>', _bind_pw)
        inner.bind('<Leave>', _unbind_pw)

        self._build_panel_content(inner)

    def _build_panel_content(self, parent):
        tk.Label(parent, text="📡  FRESNEL", bg=_PNL, fg=_ACC,
                 font=('Consolas', 14, 'bold')).pack(anchor='w', padx=10, pady=(16, 2))
        tk.Label(parent, text="Radio Enlace Urbano · 3D · Factibilidad · hasta 10 km",
                 bg=_PNL, fg=_MUT, font=('Consolas', 8)
                 ).pack(anchor='w', padx=10, pady=(0, 6))
        tk.Frame(parent, bg=_BDR, height=1).pack(fill=tk.X, padx=8)

        # ── 1. ENLACE ─────────────────────────────────────────────────────────
        b = self._section(parent, "ENLACE", '#ffcc44')
        self._field(b, "Distancia A ↔ B (m)", 'd_AB', 100, 0)
        self._field(b, "Frecuencia",           'freq', 2.4,  1)
        self._dropdown(b, "Unidad", 'freq_unit', ['Hz', 'kHz', 'MHz', 'GHz'], 'GHz', 2)

        # ── 2. ZONA DE FRESNEL ────────────────────────────────────────────────
        tk.Label(parent, text="▸ ZONA DE FRESNEL", bg=_PNL, fg='#cc88ff',
                 font=('Consolas', 9, 'bold')
                 ).pack(anchor='w', padx=10, pady=(10, 0))
        tk.Frame(parent, bg='#cc88ff', height=1).pack(fill=tk.X, padx=10, pady=(1, 4))

        self._vars['zone_60']            = tk.BooleanVar(value=True)
        self._vars['zone_80']            = tk.BooleanVar(value=False)
        self._vars['zone_100']           = tk.BooleanVar(value=False)
        self._vars['zone_other_enabled'] = tk.BooleanVar(value=False)
        self._vars['zone_other']         = tk.StringVar(value="0.90")

        zf = tk.Frame(parent, bg=_PNL)
        zf.pack(fill=tk.X, padx=16)
        tk.Label(zf, text="Elipsoides de Fresnel:", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9)).pack(anchor='w')
        row = tk.Frame(zf, bg=_PNL); row.pack(fill=tk.X, pady=(4, 2))
        for txt, key, col in [("60%", 'zone_60', '#22aaee'),
                               ("80%", 'zone_80', '#9955ee'),
                               ("100%",'zone_100','#eeaa22')]:
            tk.Checkbutton(row, text=txt, variable=self._vars[key],
                           command=self.simular, bg=_PNL, fg=_TXT,
                           selectcolor=col, activebackground=_PNL,
                           activeforeground=_ACC,
                           font=('Consolas', 9, 'bold')
                           ).pack(side=tk.LEFT, padx=(0, 6))
        other_row = tk.Frame(zf, bg=_PNL); other_row.pack(fill=tk.X, pady=(2, 0))
        tk.Checkbutton(other_row, text="Otro",
                       variable=self._vars['zone_other_enabled'],
                       command=self._toggle_zone_other, bg=_PNL, fg=_TXT,
                       selectcolor='#cccccc', activebackground=_PNL,
                       activeforeground=_ACC,
                       font=('Consolas', 9, 'bold')
                       ).pack(side=tk.LEFT)
        self._zone_other_entry = tk.Entry(
            other_row, textvariable=self._vars['zone_other'],
            bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
            font=('Consolas', 9), width=6, relief='flat',
            highlightbackground=_BDR, highlightcolor=_ACC,
            highlightthickness=1, state='disabled')
        self._zone_other_entry.pack(side=tk.LEFT, padx=(6, 0))
        self._zone_other_entry.bind('<Return>', self.simular)
        tk.Label(other_row, text="valor 0.01-1.00", bg=_PNL, fg=_MUT,
                 font=('Consolas', 7)).pack(side=tk.LEFT, padx=(8, 0))

        # ── 3. POTENCIA, GANANCIAS Y MARGEN DE SEGURIDAD ─────────────────────
        b = self._section(parent, "POTENCIA Y GANANCIAS", '#88bbff')

        # Fila 0: Potencia Tx
        tk.Label(b, text="Potencia Tx", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=0, column=0, sticky='w', pady=2)
        self._vars['tx_power_val'] = tk.StringVar(value="20")
        tk.Entry(b, textvariable=self._vars['tx_power_val'],
                 bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
                 font=('Consolas', 10), width=7, relief='flat',
                 highlightbackground=_BDR, highlightcolor=_ACC, highlightthickness=1
                 ).grid(row=0, column=1, sticky='e', pady=2, padx=(4, 0))
        self._vars['tx_power_unit'] = tk.StringVar(value="dBm")
        tpu = ttk.OptionMenu(b, self._vars['tx_power_unit'], "dBm", "dBm", "W", "kW",
                             command=lambda _: self.simular())
        tpu.config(width=7)
        tpu.grid(row=0, column=2, sticky='e', pady=2, padx=(4, 0))

        # Fila 1: Ganancia Tx
        tk.Label(b, text="Ganancia Tx", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=1, column=0, sticky='w', pady=2)
        self._vars['gtx'] = tk.StringVar(value="0")
        tk.Entry(b, textvariable=self._vars['gtx'],
                 bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
                 font=('Consolas', 10), width=7, relief='flat',
                 highlightbackground=_BDR, highlightcolor=_ACC, highlightthickness=1
                 ).grid(row=1, column=1, sticky='e', pady=2, padx=(4, 0))
        self._vars['gtx_unit'] = tk.StringVar(value="dBi")
        gtu = ttk.OptionMenu(b, self._vars['gtx_unit'], "dBi", "dB", "dBi", "dBd",
                             command=lambda _: self.simular())
        gtu.config(width=5)
        gtu.grid(row=1, column=2, sticky='e', pady=2, padx=(4, 0))

        # Fila 2: Ganancia Rx
        tk.Label(b, text="Ganancia Rx", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=2, column=0, sticky='w', pady=2)
        self._vars['grx'] = tk.StringVar(value="0")
        tk.Entry(b, textvariable=self._vars['grx'],
                 bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
                 font=('Consolas', 10), width=7, relief='flat',
                 highlightbackground=_BDR, highlightcolor=_ACC, highlightthickness=1
                 ).grid(row=2, column=1, sticky='e', pady=2, padx=(4, 0))
        self._vars['grx_unit'] = tk.StringVar(value="dBi")
        gru = ttk.OptionMenu(b, self._vars['grx_unit'], "dBi", "dB", "dBi", "dBd",
                             command=lambda _: self.simular())
        gru.config(width=5)
        gru.grid(row=2, column=2, sticky='e', pady=2, padx=(4, 0))

        # Fila 3: Margen de seguridad (pérdidas por entorno)
        tk.Label(b, text="Margen de seguridad (dB)", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=3, column=0, sticky='w', pady=2)
        self._vars['margin'] = tk.StringVar(value="0")
        tk.Entry(b, textvariable=self._vars['margin'],
                 bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
                 font=('Consolas', 10), width=7, relief='flat',
                 highlightbackground=_BDR, highlightcolor=_ACC, highlightthickness=1
                 ).grid(row=3, column=1, sticky='e', pady=2, padx=(4, 0))
        tk.Label(b, text="(0 = sin pérdidas adicionales)", bg=_PNL, fg=_MUT,
                 font=('Consolas', 7), anchor='w'
                 ).grid(row=3, column=2, sticky='w', pady=2, padx=(4, 0))

        # ── 4. SENSIBILIDAD DEL RECEPTOR (Pmin) ──────────────────────────────
        s = self._section(parent, "SENSIBILIDAD DEL RECEPTOR", '#ff9977')
        tk.Label(s, text="P mínima (dBm/W/µW/nW)", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), anchor='w'
                 ).grid(row=0, column=0, sticky='w', pady=2)
        self._vars['rx_sens_val'] = tk.StringVar(value="-90")
        tk.Entry(s, textvariable=self._vars['rx_sens_val'],
                 bg='#0a0a0a', fg=_TXT, insertbackground=_TXT,
                 font=('Consolas', 10), width=7, relief='flat',
                 highlightbackground=_BDR, highlightcolor=_ACC, highlightthickness=1
                 ).grid(row=0, column=1, sticky='e', pady=2, padx=(4, 0))
        self._vars['rx_sens_unit'] = tk.StringVar(value="dBm")
        rsu = ttk.OptionMenu(s, self._vars['rx_sens_unit'], "dBm",
                             "dBm", "W", "kW", "µW", "nW",
                             command=lambda _: self.simular())
        rsu.config(width=7)
        rsu.grid(row=0, column=2, sticky='e', pady=2, padx=(4, 0))

        # ── 5. EDIFICIO A ─────────────────────────────────────────────────────
        b = self._section(parent, "PUNTO DE TRANSMISIÓN A  —  Transmisor", '#6699ff')
        self._field(b, "Altura edificio A (m)",  'hA',      60, 0)
        self._field(b, "Altura mástil Tx (m)",    'hA_ant',  10, 1)
        self._field(b, "Rotación edificio A (°)", 'rotA',     0, 2)
        self._field(b, "Ancho edificio A (m)",    'bwA',     15, 3)
        self._dropdown(b, "Posición antena A",
                       'antA_pos', ['Centro', 'Esq1', 'Esq2', 'Esq3', 'Esq4'], 'Centro', 4)

        # ── 6. EDIFICIO B ─────────────────────────────────────────────────────
        b = self._section(parent, "PUNTO DE TRANSMISIÓN B  —  Receptor", '#44cc77')
        self._field(b, "Altura edificio B (m)",  'hB',       50, 0)
        self._field(b, "Altura mástil Rx (m)",    'hB_ant',   5, 1)
        self._field(b, "Rotación edificio B (°)", 'rotB',     0, 2)
        self._field(b, "Ancho edificio B (m)",    'bwB',     15, 3)
        self._dropdown(b, "Posición antena B",
                       'antB_pos', ['Centro', 'Esq1', 'Esq2', 'Esq3', 'Esq4'], 'Centro', 4)

        # ── 7. OBSTÁCULOS múltiples ───────────────────────────────────────────
        tk.Label(parent, text="▸ OBSTÁCULOS DE OBSTRUCCIÓN", bg=_PNL,
                 fg='#ffaa55', font=('Consolas', 9, 'bold')
                 ).pack(anchor='w', padx=10, pady=(12, 0))
        tk.Frame(parent, bg='#ffaa55', height=1).pack(fill=tk.X, padx=10, pady=(1, 6))

        nrow = tk.Frame(parent, bg=_PNL)
        nrow.pack(fill=tk.X, padx=16, pady=(0, 6))
        tk.Label(nrow, text="Cantidad de obstáculos:", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9)).pack(side=tk.LEFT)
        btn_minus = tk.Button(nrow, text="−", bg='#1e1e1e', fg=_TXT,
                              font=('Consolas', 10, 'bold'), relief='flat',
                              width=2, command=self._dec_obs)
        btn_minus.pack(side=tk.LEFT, padx=(8, 2))
        self._n_obs_label = tk.Label(nrow, textvariable=self._n_obs_var,
                                      bg='#0a0a0a', fg='#ffe055',
                                      font=('Consolas', 11, 'bold'), width=2,
                                      relief='flat')
        self._n_obs_label.pack(side=tk.LEFT, padx=2)
        btn_plus = tk.Button(nrow, text="+", bg='#1e1e1e', fg=_TXT,
                             font=('Consolas', 10, 'bold'), relief='flat',
                             width=2, command=self._inc_obs)
        btn_plus.pack(side=tk.LEFT, padx=(2, 0))
        tk.Label(nrow, text=f"(máx {MAX_OBS})", bg=_PNL, fg=_MUT,
                 font=('Consolas', 7)).pack(side=tk.LEFT, padx=(8, 0))

        self._obs_container = tk.Frame(parent, bg=_PNL)
        self._obs_container.pack(fill=tk.X)

        self._rebuild_obs_panels()

        # ── Botón principal ───────────────────────────────────────────────────
        tk.Button(parent, text="▶   SIMULAR",
                  bg=_ACC, fg='white', font=('Consolas', 12, 'bold'),
                  relief='flat', pady=10, cursor='hand2',
                  command=self.simular
                  ).pack(fill=tk.X, padx=10, pady=(10, 4))

    # ── Gestión de obstáculos dinámicos ───────────────────────────────────────

    def _inc_obs(self):
        n = self._n_obs_var.get()
        if n < MAX_OBS:
            self._n_obs_var.set(n + 1)
            self._rebuild_obs_panels()
            self.simular()

    def _dec_obs(self):
        n = self._n_obs_var.get()
        if n > 1:
            self._n_obs_var.set(n - 1)
            self._rebuild_obs_panels()
            self.simular()

    def _rebuild_obs_panels(self):
        n = self._n_obs_var.get()
        for child in self._obs_container.winfo_children():
            child.destroy()

        while len(self._obs_vars) < n:
            i = len(self._obs_vars)
            d = _OBS_DEF[min(i, MAX_OBS - 1)]
            self._obs_vars.append({
                'h':   tk.StringVar(value=str(d['h'])),
                'd':   tk.StringVar(value=str(d['d'])),
                'off': tk.StringVar(value=str(d['off'])),
                'bw':  tk.StringVar(value=str(d['bw'])),
                'rot': tk.StringVar(value=str(d['rot'])),
            })

        for i in range(n):
            name = _OBS_NAMES[min(i, MAX_OBS - 1)]
            hcol = _OBS_HCOL[min(i, MAX_OBS - 1)]
            v    = self._obs_vars[i]

            tk.Label(self._obs_container,
                     text=f"▸ OBSTÁCULO {name}", bg=_PNL, fg=hcol,
                     font=('Consolas', 9, 'bold')
                     ).pack(anchor='w', padx=10, pady=(10, 0))
            tk.Frame(self._obs_container, bg=hcol, height=1
                     ).pack(fill=tk.X, padx=10, pady=(1, 4))

            body = tk.Frame(self._obs_container, bg=_PNL)
            body.pack(fill=tk.X, padx=16, pady=2)
            body.columnconfigure(0, weight=1)

            for row_i, (lbl, var_key) in enumerate([
                (f"Altura edif. {name} (m)", 'h'),
                (f"Dist. A → {name} (m)",    'd'),
                ("Offset lateral (m)",        'off'),
                (f"Ancho edif. {name} (m)",   'bw'),
                ("Rotación (°)",              'rot'),
            ]):
                self._obs_field(body, lbl, v[var_key], row_i)

    # ── Panel derecho ─────────────────────────────────────────────────────────

    def _card(self, parent, border_color, pady=(4, 4)):
        outer = tk.Frame(parent, bg=border_color)
        outer.pack(fill=tk.X, padx=8, pady=pady)
        inner = tk.Frame(outer, bg='#0e0e12', padx=6, pady=4)
        inner.pack(fill=tk.X, padx=2)
        return inner

    def _info_row(self, parent, label, key, fg_val=None):
        row = tk.Frame(parent, bg='#0e0e12')
        row.pack(fill=tk.X, pady=1)
        tk.Label(row, text=label, bg='#0e0e12', fg=_MUT,
                 font=('Consolas', 8), anchor='w').pack(side=tk.LEFT)
        val = tk.Label(row, text="—", bg='#0e0e12',
                       fg=fg_val if fg_val else _TXT,
                       font=('Consolas', 9, 'bold'))
        val.pack(side=tk.RIGHT)
        self._info[key] = val

    def _build_side_panel(self, parent):
        vscroll = tk.Scrollbar(parent, orient='vertical', bg='#1e1e1e',
                               troughcolor=_BDR, activebackground='#444', width=10)
        vscroll.pack(side=tk.RIGHT, fill=tk.Y)

        self._right_canvas = tk.Canvas(parent, bg=_PNL,
                                       yscrollcommand=vscroll.set,
                                       highlightthickness=0)
        self._right_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vscroll.config(command=self._right_canvas.yview)

        inner = tk.Frame(self._right_canvas, bg=_PNL)
        self._right_inner_id = self._right_canvas.create_window(
            (0, 0), window=inner, anchor='nw')

        inner.bind('<Configure>',
                   lambda e: self._right_canvas.configure(
                       scrollregion=self._right_canvas.bbox('all')))
        self._right_canvas.bind('<Configure>',
                                lambda e: self._right_canvas.itemconfig(
                                    self._right_inner_id, width=e.width))

        def _pw_side(ev):
            if getattr(ev, 'num', None) == 4:
                self._right_canvas.yview_scroll(-1, 'units')
            elif getattr(ev, 'num', None) == 5:
                self._right_canvas.yview_scroll(1, 'units')
            elif getattr(ev, 'delta', 0):
                self._right_canvas.yview_scroll(int(-ev.delta / 120), 'units')

        def _bind_side(e):
            self._right_canvas.bind_all('<MouseWheel>', _pw_side)
            self._right_canvas.bind_all('<Button-4>', _pw_side)
            self._right_canvas.bind_all('<Button-5>', _pw_side)

        def _unbind_side(e):
            self._right_canvas.unbind_all('<MouseWheel>')
            self._right_canvas.unbind_all('<Button-4>')
            self._right_canvas.unbind_all('<Button-5>')

        inner.bind('<Enter>', _bind_side)
        inner.bind('<Leave>', _unbind_side)

        # ----------------------------------------------
        # Contenido del panel derecho dentro de `inner`
        # ----------------------------------------------
        tk.Label(inner, text="📊  ANÁLISIS DEL ENLACE", bg=_PNL, fg=_ACC,
                 font=('Consolas', 13, 'bold')).pack(anchor='w', padx=10, pady=(14, 1))
        tk.Label(inner, text="Geometría · Potencia · Margen",
                 bg=_PNL, fg=_MUT, font=('Consolas', 8)
                 ).pack(anchor='w', padx=10, pady=(0, 5))
        tk.Frame(inner, bg=_BDR, height=1).pack(fill=tk.X, padx=8)

        self._stlbl = tk.Label(inner, text="", bg='#0e0e12', fg=_TXT,
                               font=('Consolas', 13, 'bold'),
                               wraplength=285, justify='center',
                               relief='flat', bd=0, padx=10, pady=12)
        self._stlbl.pack(fill=tk.X, padx=8, pady=(8, 4))
        self._stlbl_border = tk.Frame(inner, bg=_BDR, height=2)
        self._stlbl_border.pack(fill=tk.X, padx=8, pady=(0, 4))

        self._info = {}

        # Card 1: ENLACE RF
        tk.Label(inner, text="▸ ENLACE RF", bg=_PNL, fg='#88aaff',
                 font=('Consolas', 9, 'bold')).pack(anchor='w', padx=10, pady=(6, 0))
        c1 = self._card(inner, '#2a3a60')
        self._info_row(c1, 'λ longitud de onda', 'lam')
        self._info_row(c1, 'Distancia LOS', 'd3d')
        self._info_row(c1, 'r₁ máx. Fresnel', 'r1')

        # Card 2: OBSTÁCULO CRÍTICO
        tk.Label(inner, text="▸ OBSTÁCULO CRÍTICO", bg=_PNL, fg='#ffaa55',
                 font=('Consolas', 9, 'bold')).pack(anchor='w', padx=10, pady=(6, 0))
        c2 = self._card(inner, '#604020')
        self._info_row(c2, 'd₁* desde Tx', 'd1c')
        self._info_row(c2, 'd₂* hacia Rx', 'd2c')
        self._info_row(c2, 'Distancia mínima', 'min')
        self._info_row(c2, 'Radios de zona', 'rz')
        row_clr = tk.Frame(c2, bg='#0e0e12'); row_clr.pack(fill=tk.X, pady=2)
        tk.Label(row_clr, text='Despeje al 60%', bg='#0e0e12', fg=_MUT,
                 font=('Consolas', 8), anchor='w').pack(side=tk.LEFT)
        val_clr = tk.Label(row_clr, text="—", bg='#0e0e12', fg=_TXT,
                           font=('Consolas', 10, 'bold'))
        val_clr.pack(side=tk.RIGHT)
        self._info['clr'] = val_clr

        # Card 3: POTENCIA Y MARGEN
        tk.Label(inner, text="▸ POTENCIA Y MARGEN", bg=_PNL, fg='#88ddff',
                 font=('Consolas', 9, 'bold')).pack(anchor='w', padx=10, pady=(6, 0))
        c3 = self._card(inner, '#204858', pady=(4, 6))
        self._info_row(c3, 'FSPL (espacio libre)', 'fspl')
        self._info_row(c3, 'Atenuación difracción', 'adif', fg_val='#ffaa55')
        self._info_row(c3, 'Ganancia total (Tx+Rx)', 'gtotal')
        tk.Frame(c3, bg='#2a3a48', height=1).pack(fill=tk.X, pady=(4, 3))

        # Prx en dBm
        prx_row = tk.Frame(c3, bg='#0e0e12'); prx_row.pack(fill=tk.X, pady=2)
        tk.Label(prx_row, text='Prx (dBm)', bg='#0e0e12', fg='#aaddff',
                 font=('Consolas', 9, 'bold'), anchor='w').pack(side=tk.LEFT)
        val_prx_dbm = tk.Label(prx_row, text="—", bg='#0e0e12', fg=_TXT,
                               font=('Consolas', 12, 'bold'))
        val_prx_dbm.pack(side=tk.RIGHT)
        self._info['prx_dbm'] = val_prx_dbm

        # Prx en vatios
        prx_w_row = tk.Frame(c3, bg='#0e0e12'); prx_w_row.pack(fill=tk.X, pady=2)
        tk.Label(prx_w_row, text='Prx (W)', bg='#0e0e12', fg='#aaddff',
                 font=('Consolas', 9, 'bold'), anchor='w').pack(side=tk.LEFT)
        val_prx_w = tk.Label(prx_w_row, text="—", bg='#0e0e12', fg=_TXT,
                             font=('Consolas', 10))
        val_prx_w.pack(side=tk.RIGHT)
        self._info['prx_w'] = val_prx_w

        # Margen de enlace (calculado)
        margin_row = tk.Frame(c3, bg='#0e0e12'); margin_row.pack(fill=tk.X, pady=2)
        tk.Label(margin_row, text='Margen de enlace', bg='#0e0e12', fg='#ffdd99',
                 font=('Consolas', 9, 'bold'), anchor='w').pack(side=tk.LEFT)
        val_margin = tk.Label(margin_row, text="—", bg='#0e0e12', fg=_TXT,
                              font=('Consolas', 12, 'bold'))
        val_margin.pack(side=tk.RIGHT)
        self._info['margin'] = val_margin

        # Sensibilidad configurada
        sens_row = tk.Frame(c3, bg='#0e0e12'); sens_row.pack(fill=tk.X, pady=2)
        tk.Label(sens_row, text='Sensibilidad (Pmin)', bg='#0e0e12', fg=_MUT,
                 font=('Consolas', 8), anchor='w').pack(side=tk.LEFT)
        val_sens = tk.Label(sens_row, text="—", bg='#0e0e12', fg=_TXT,
                            font=('Consolas', 9))
        val_sens.pack(side=tk.RIGHT)
        self._info['sens'] = val_sens

        # Margen de seguridad configurado
        margin_set_row = tk.Frame(c3, bg='#0e0e12'); margin_set_row.pack(fill=tk.X, pady=2)
        tk.Label(margin_set_row, text='Margen de seguridad configurado', bg='#0e0e12', fg=_MUT,
                 font=('Consolas', 8), anchor='w').pack(side=tk.LEFT)
        val_margin_set = tk.Label(margin_set_row, text="—", bg='#0e0e12', fg=_TXT,
                                  font=('Consolas', 9))
        val_margin_set.pack(side=tk.RIGHT)
        self._info['margin_set'] = val_margin_set

        # Tabla de obstáculos
        tk.Label(inner, text="▸ TABLA DE OBSTÁCULOS", bg=_PNL, fg='#bbbbcc',
                 font=('Consolas', 9, 'bold')).pack(anchor='w', padx=10, pady=(8, 0))
        self._obs_table_frame = tk.Frame(inner, bg='#0e0e12')
        self._obs_table_frame.pack(fill=tk.X, padx=8, pady=(2, 4))

        # Cámara
        self._vars['zoom'] = tk.DoubleVar(value=1.5)
        tk.Frame(inner, bg=_BDR, height=1).pack(fill=tk.X, padx=8, pady=(6, 4))
        tk.Label(inner, text="▸ CÁMARA  &  VISTA", bg=_PNL, fg=_MUT,
                 font=('Consolas', 9, 'bold')).pack(anchor='w', padx=10)
        cf = tk.Frame(inner, bg=_PNL); cf.pack(fill=tk.X, padx=10, pady=2)

        for lbl, key, dflt, lo, hi in [
            ("Azimut :", 'azim', 45, 0, 360),
            ("Elevac. :", 'elev', 30, 0, 90),
        ]:
            rw = tk.Frame(cf, bg=_PNL); rw.pack(fill=tk.X, pady=2)
            tk.Label(rw, text=lbl, bg=_PNL, fg=_TXT,
                     font=('Consolas', 9), width=9).pack(side=tk.LEFT)
            var = tk.IntVar(value=dflt)
            self._vars[key] = var
            tk.Scale(rw, from_=lo, to=hi, orient=tk.HORIZONTAL,
                     variable=var, command=self._update_cam,
                     bg=_PNL, fg=_TXT, troughcolor='#2a2a3a',
                     highlightthickness=0, activebackground='#555',
                     length=120, font=('Consolas', 7)
                     ).pack(side=tk.LEFT, expand=True)

        zr = tk.Frame(cf, bg=_PNL); zr.pack(fill=tk.X, pady=2)
        tk.Label(zr, text="Zoom:", bg=_PNL, fg=_TXT,
                 font=('Consolas', 9), width=9).pack(side=tk.LEFT)
        tk.Button(zr, text="–", command=lambda: self._adjust_zoom(0.85),
                  bg='#222230', fg=_TXT, font=('Consolas', 9, 'bold'),
                  relief='flat', width=2, padx=4).pack(side=tk.LEFT, padx=(0, 4))
        tk.Button(zr, text="+", command=lambda: self._adjust_zoom(1.15),
                  bg='#222230', fg=_TXT, font=('Consolas', 9, 'bold'),
                  relief='flat', width=2, padx=4).pack(side=tk.LEFT)
        self._zoom_label = tk.Label(zr, text="150%", bg=_PNL, fg=_ACC,
                                    font=('Consolas', 9, 'bold'))
        self._zoom_label.pack(side=tk.LEFT, padx=(10, 0))

        # Checkboxes
        tk.Frame(inner, bg=_BDR, height=1).pack(fill=tk.X, padx=8, pady=(6, 2))
        for txt, var, sel in [
            ("Mostrar Elipsoide de Fresnel", self._show_ell, '#003366'),
            ("Modo Rayos X  (translúcido)", self._xray_mode, '#221133'),
            ("Escala 3D proporcional", self._vars['uniform_scale'], '#222222'),
            ("Mostrar etiquetas de texto", self._show_labels, '#222222'),
        ]:
            tk.Checkbutton(inner, text=txt, variable=var,
                           command=self.simular, bg=_PNL, fg=_TXT,
                           selectcolor=sel, activebackground=_PNL,
                           activeforeground=_ACC, font=('Consolas', 9)
                           ).pack(anchor='w', pady=2, padx=12)

        # Botones
        btn_row = tk.Frame(inner, bg=_PNL)
        btn_row.pack(fill=tk.X, padx=10, pady=(8, 6))
        tk.Button(btn_row, text="↺ Reset cámara",
                  bg='#222230', fg=_MUT, font=('Consolas', 9),
                  relief='flat', command=self._reset_cam, padx=6, pady=3
                  ).pack(side=tk.LEFT, padx=(0, 6))
        tk.Button(btn_row, text="⊕ Centrar antenas",
                  bg='#222230', fg=_MUT, font=('Consolas', 9),
                  relief='flat', command=self._center_antennas, padx=6, pady=3
                  ).pack(side=tk.LEFT)

    def _update_obs_table(self, results):
        for child in self._obs_table_frame.winfo_children():
            child.destroy()

        BG2 = '#0e0e12'
        colors = ['#111118', '#15151c']

        cols = [
            ('Obs',      50, 'w'),
            ('d_mín',    70, 'e'),
            ('r₆₀',      70, 'e'),
            ('Despeje',  80, 'e'),
            ('Aten.dB',  70, 'e'),
            ('',         30, 'e'),
        ]

        hdr = tk.Frame(self._obs_table_frame, bg=BG2)
        hdr.pack(fill=tk.X, padx=4, pady=(3, 0))
        for i, (txt, w, anchor) in enumerate(cols):
            lbl = tk.Label(hdr, text=txt, bg=BG2, fg='#aabbcc',
                           font=('Consolas', 8, 'bold'))
            lbl.grid(row=0, column=i, sticky=anchor, padx=2, pady=2)
            hdr.columnconfigure(i, minsize=w, weight=0)

        sep = tk.Frame(self._obs_table_frame, bg='#334455', height=1)
        sep.pack(fill=tk.X, padx=4, pady=2)

        for idx, r in enumerate(results):
            ok = r['clr'] > 0.0
            cc = '#44ee77' if ok else '#ff5555'
            bg_r = colors[idx % 2] if ok else '#180808'
            row = tk.Frame(self._obs_table_frame, bg=bg_r)
            row.pack(fill=tk.X, padx=4, pady=1)
            for i, (_, w, anchor) in enumerate(cols):
                row.columnconfigure(i, minsize=w, weight=0)

            tk.Label(row, text=r['name'], bg=bg_r, fg=r['hcol'],
                     font=('Consolas', 9, 'bold'), anchor='w'
                     ).grid(row=0, column=0, sticky='w', padx=2, pady=2)
            tk.Label(row, text=f"{r['dist_min']:.2f}m", bg=bg_r, fg=_TXT,
                     font=('Consolas', 8), anchor='e'
                     ).grid(row=0, column=1, sticky='e', padx=2, pady=2)
            tk.Label(row, text=f"{r['r60']:.2f}m", bg=bg_r, fg=_TXT,
                     font=('Consolas', 8), anchor='e'
                     ).grid(row=0, column=2, sticky='e', padx=2, pady=2)
            tk.Label(row, text=f"{r['clr']:+.2f}m", bg=bg_r, fg=cc,
                     font=('Consolas', 9, 'bold'), anchor='e'
                     ).grid(row=0, column=3, sticky='e', padx=2, pady=2)
            att_fg = '#ffbb55' if r['adif'] > 0.5 else _MUT
            tk.Label(row, text=f"{r['adif']:.1f}", bg=bg_r, fg=att_fg,
                     font=('Consolas', 9, 'bold'), anchor='e'
                     ).grid(row=0, column=4, sticky='e', padx=2, pady=2)
            tk.Label(row, text='✅' if ok else '❌', bg=bg_r, fg=cc,
                     font=('Consolas', 9), anchor='e'
                     ).grid(row=0, column=5, sticky='e', padx=2, pady=2)

    # ── Canvas matplotlib ─────────────────────────────────────────────────────

    def _build_canvas(self, parent):
        self.fig = Figure(figsize=(12.0, 9.0), dpi=96)
        self.fig.patch.set_facecolor(_BG)
        self.ax  = self.fig.add_subplot(111, projection='3d')
        self.ax.set_facecolor(_BG)
        self.fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
        self.cvs = FigureCanvasTkAgg(self.fig, master=parent)
        w = self.cvs.get_tk_widget()
        w.pack(fill=tk.BOTH, expand=True)
        w.bind('<MouseWheel>',     self._on_mousewheel)
        w.bind('<Button-4>',       self._on_mousewheel)
        w.bind('<Button-5>',       self._on_mousewheel)
        w.bind('<ButtonPress-1>',  self._on_mouse_press)
        w.bind('<B1-Motion>',      self._on_mouse_drag)
        w.bind('<ButtonRelease-1>',self._on_mouse_release)
        w.bind('<ButtonPress-2>',  lambda e: 'break')
        w.bind('<B2-Motion>',      lambda e: 'break')
        w.bind('<ButtonPress-3>',  lambda e: 'break')
        w.bind('<B3-Motion>',      lambda e: 'break')

    # ── Cámara e interacción ──────────────────────────────────────────────────

    def _update_cam(self, _=None):
        if hasattr(self, 'ax'):
            self.ax.view_init(elev=self._vars['elev'].get(),
                              azim=self._vars['azim'].get())
            self.cvs.draw_idle()

    def _reset_cam(self):
        self._vars['azim'].set(45)
        self._vars['elev'].set(30)
        self._update_cam()

    def _center_antennas(self):
        self._vars['antA_pos'].set('Centro')
        self._vars['antB_pos'].set('Centro')
        self.simular()

    def _toggle_zone_other(self):
        active = self._vars['zone_other_enabled'].get()
        self._zone_other_entry.config(state='normal' if active else 'disabled')
        if active:
            self._zone_other_entry.focus_set()
        self.simular()

    def _adjust_zoom(self, factor):
        z = float(np.clip(self._vars['zoom'].get() * factor, 0.4, 2.5))
        self._vars['zoom'].set(z)
        if hasattr(self, '_zoom_label'):
            self._zoom_label.config(text=f"{int(z * 100)}%")
        self.simular()

    def _on_mousewheel(self, event):
        if getattr(event, 'delta', 0):
            factor = 1.15 if event.delta > 0 else 0.85
        elif getattr(event, 'num', None) in (4, 5):
            factor = 1.15 if event.num == 4 else 0.85
        else:
            return
        self._adjust_zoom(factor)
        return "break"

    def _on_mouse_press(self, event):
        self._drag_info = {
            'x': event.x, 'y': event.y,
            'azim': self._vars['azim'].get(),
            'elev': self._vars['elev'].get(),
        }
        return 'break'

    def _on_mouse_release(self, event):
        self._drag_info.clear()
        return 'break'

    def _on_mouse_drag(self, event):
        if not self._drag_info:
            return 'break'
        dx = event.x - self._drag_info['x']
        dy = event.y - self._drag_info['y']
        self._vars['azim'].set(int(round((self._drag_info['azim'] - dx * 0.35) % 360)))
        self._vars['elev'].set(int(round(float(np.clip(
            self._drag_info['elev'] + dy * 0.22, 0, 90)))))
        self._update_cam()
        return 'break'

    # ── Lectura y validación de parámetros ───────────────────────────────────

    def _params(self):
        try:
            p = {}
            for k, v in self._vars.items():
                if k in ('azim', 'elev'):
                    p[k] = int(v.get())
                elif k in ('zone_60', 'zone_80', 'zone_100', 'zone_other_enabled',
                           'uniform_scale'):
                    p[k] = bool(v.get())
                elif k == 'zone_other':
                    p[k] = float(v.get())
                elif k in ('antA_pos', 'antB_pos', 'freq_unit',
                           'tx_power_unit', 'rx_sens_unit',
                           'gtx_unit', 'grx_unit'):
                    p[k] = str(v.get())
                elif k == 'zoom':
                    p[k] = float(v.get())
                else:
                    p[k] = float(v.get())

            # Verificar que las claves esenciales existen
            required = ('freq', 'freq_unit', 'd_AB', 'hA', 'hA_ant', 'hB', 'hB_ant',
                        'bwA', 'bwB', 'rotA', 'rotB',
                        'tx_power_val', 'tx_power_unit',
                        'gtx', 'grx', 'gtx_unit', 'grx_unit',
                        'rx_sens_val', 'rx_sens_unit', 'margin')
            for r in required:
                if r not in p:
                    raise ValueError(f"Falta el parámetro {r} en la configuración.")

            factor_freq = {'Hz': 1.0, 'kHz': 1e3, 'MHz': 1e6, 'GHz': 1e9}.get(p['freq_unit'], 1e9)
            p['freq_hz'] = p['freq'] * factor_freq

            # Convertir potencia Tx a dBm
            tx_val = p['tx_power_val']
            tx_unit = p['tx_power_unit']
            if tx_unit == 'dBm':
                p['tx_power_dBm'] = tx_val
            elif tx_unit == 'W':
                p['tx_power_dBm'] = 10.0 * math.log10(tx_val * 1000.0)
            elif tx_unit == 'kW':
                p['tx_power_dBm'] = 10.0 * math.log10(tx_val * 1e6)
            else:
                p['tx_power_dBm'] = tx_val

            # Ganancias de antena: dB y dBi se tratan igual, dBd se convierte sumando 2.15
            gtx_val = p['gtx']
            gtx_unit = p['gtx_unit']
            if gtx_unit == 'dBd':
                gtx_dBi = gtx_val + 2.15
            else:
                gtx_dBi = gtx_val
            p['gtx'] = gtx_dBi

            grx_val = p['grx']
            grx_unit = p['grx_unit']
            if grx_unit == 'dBd':
                grx_dBi = grx_val + 2.15
            else:
                grx_dBi = grx_val
            p['grx'] = grx_dBi

            # Convertir sensibilidad a dBm
            sens_val = p['rx_sens_val']
            sens_unit = p['rx_sens_unit']
            if sens_unit == 'dBm':
                p['rx_sens_dBm'] = sens_val
            elif sens_unit == 'W':
                p['rx_sens_dBm'] = 10.0 * math.log10(sens_val * 1000.0)
            elif sens_unit == 'kW':
                p['rx_sens_dBm'] = 10.0 * math.log10(sens_val * 1e6)
            elif sens_unit == 'µW':
                p['rx_sens_dBm'] = 10.0 * math.log10(sens_val * 1e-3)
            elif sens_unit == 'nW':
                p['rx_sens_dBm'] = 10.0 * math.log10(sens_val * 1e-6)
            else:
                p['rx_sens_dBm'] = sens_val

            # Margen de seguridad
            p['margin'] = p['margin']
            if p['margin'] < 0:
                raise ValueError("El margen de seguridad debe ser >= 0 dB.")

            # Validar límite de distancia (10 km)
            if p['d_AB'] > MAX_DIST_KM * 1000:
                raise ValueError(f"La distancia máxima permitida es {MAX_DIST_KM} km (no se considera curvatura terrestre).")

            if p['hA'] < 0 or p['hB'] < 0 or p['hA_ant'] < 0 or p['hB_ant'] < 0:
                raise ValueError('Las alturas deben ser positivas.')
            if p['zone_other_enabled'] and not (0.01 <= p['zone_other'] <= 1.0):
                raise ValueError('Zona "Otro" debe estar entre 0.01 y 1.00.')
            if not any((p['zone_60'], p['zone_80'], p['zone_100'],
                        p['zone_other_enabled'])):
                raise ValueError('Seleccione al menos una zona de Fresnel.')
            if p['bwA'] <= 0 or p['bwB'] <= 0:
                raise ValueError("Ancho de edificio debe ser positivo para A y B.")
            return p
        except (ValueError, AssertionError) as e:
            messagebox.showerror("Parámetro inválido", str(e))
            return None

    # ══════════════════════════════════════════════════════════════════════════
    #  NÚCLEO DE SIMULACIÓN
    # ══════════════════════════════════════════════════════════════════════════

    def simular(self, event=None):
        p = self._params()
        if not p:
            return "break" if event is not None else None

        lam  = C_LIGHT / p['freq_hz']
        d_AB = p['d_AB']
        bwA  = p['bwA']
        bwB  = p['bwB']

        # ── Antenas ──────────────────────────────────────────────────────────
        ax_bx, ax_by = roof_antenna_position(
            0.0,  0.0,  bwA, p.get('rotA', 0.0), p.get('antA_pos', 'Centro'))
        bx_bx, bx_by = roof_antenna_position(
            d_AB, 0.0,  bwB, p.get('rotB', 0.0), p.get('antB_pos', 'Centro'))

        A = np.array([ax_bx, ax_by, p['hA'] + p['hA_ant']], float)
        B = np.array([bx_bx, bx_by, p['hB'] + p['hB_ant']], float)
        d_3D = float(np.linalg.norm(B - A))

        # FSPL (Free Space Path Loss) en dB
        fspl = 20.0 * math.log10(d_3D) + 20.0 * math.log10(p['freq_hz']) - 147.55
        if fspl < 0:
            fspl = 0.0

        r1_max = math.sqrt(lam * (d_3D / 2.0) ** 2 / (d_3D + 1e-12))

        # ── Zonas seleccionadas ───────────────────────────────────────────────
        zones = []
        if p['zone_60']:           zones.append(0.60)
        if p['zone_80']:           zones.append(0.80)
        if p['zone_100']:          zones.append(1.00)
        if p['zone_other_enabled']: zones.append(p['zone_other'])
        zones = sorted(set(zones))

        # ── Análisis por cada obstáculo ───────────────────────────────────────
        n_obs    = self._n_obs_var.get()
        obs_data = []

        for i in range(n_obs):
            v = self._obs_vars[i]
            try:
                h_i   = float(v['h'].get())
                d_i   = float(v['d'].get())
                off_i = float(v['off'].get())
                bw_i  = float(v['bw'].get())
                rot_i = float(v['rot'].get())
            except ValueError:
                messagebox.showerror("Parámetro inválido",
                                     f"Obstáculo {_OBS_NAMES[i]}: valor no numérico.")
                return

            if bw_i <= 0:
                messagebox.showerror("Parámetro inválido",
                                     f"Obstáculo {_OBS_NAMES[i]}: ancho debe ser > 0.")
                return
            if not (0 < d_i < d_AB):
                messagebox.showerror("Parámetro inválido",
                                     f"Obstáculo {_OBS_NAMES[i]}: "
                                     f"distancia debe estar entre 0 y d_AB={d_AB}.")
                return

            f_i    = _verts(d_i, off_i, bw_i, bw_i, h_i, rot_deg=rot_i)
            poly_i = [(vv[0], vv[1]) for vv in f_i['top']]

            t_i, _, d1c_i, d2c_i, dm_i, lp_i = find_closest_los_to_prism(
                A, B, poly_i, h_i)

            r1_i  = (math.sqrt(lam * d1c_i * d2c_i / (d_3D + 1e-12))
                     if d1c_i > 0 and d2c_i > 0 else 0.0)
            r_ref = 0.60 * r1_i
            clearance = dm_i - r_ref
            adif = 0.0
            if clearance < 0:
                H = r_ref - dm_i
                if d1c_i > 0 and d2c_i > 0:
                    nu = H * math.sqrt(2.0 * (d1c_i + d2c_i) / (lam * d1c_i * d2c_i))
                    adif = diffraction_loss(nu)
                else:
                    adif = 0.0

            obs_data.append({
                'name': _OBS_NAMES[min(i, MAX_OBS - 1)],
                'pal':  _POBS[min(i, MAX_OBS - 1)],
                'hcol': _OBS_HCOL[min(i, MAX_OBS - 1)],
                'cx': d_i, 'cy': off_i, 'h': h_i, 'bw': bw_i, 'rot': rot_i,
                'poly': poly_i,
                't_star': t_i, 'd1c': d1c_i, 'd2c': d2c_i,
                'dist_min': dm_i, 'los_pt': lp_i,
                'r1': r1_i, 'r60': 0.60 * r1_i, 'clr': clearance,
                'adif': adif,
                'penetration': dm_i < 1e-6,
                'r_visual': [zf * r1_i for zf in zones],
            })

        # ── Obstáculo más crítico ─────────────────────────────────────────────
        worst = max(obs_data, key=lambda x: (x['adif'], -x['clr']))
        clr_60      = worst['clr']
        dist_min    = worst['dist_min']
        ok_geo      = clr_60 > 0.0
        penetration = worst['penetration']
        d1c         = worst['d1c']
        d2c         = worst['d2c']
        r_visual    = worst['r_visual']
        total_adif  = worst['adif']

        # ── Cálculo de potencia con ganancias y margen de seguridad ───────────
        tx_power_dBm = p['tx_power_dBm']
        gtx = p['gtx']
        grx = p['grx']
        gtotal = gtx + grx
        margin_setting = p['margin']
        prx_dBm = tx_power_dBm + gtotal - fspl - total_adif - margin_setting

        # Sensibilidad
        rx_sens_dBm = p['rx_sens_dBm']
        margin = prx_dBm - rx_sens_dBm
        ok_pwr = margin > 0.0

        # Veredicto combinado (geometría y potencia)
        ok = ok_geo and ok_pwr

        # ── Actualizar panel derecho ──────────────────────────────────────────
        cc_geo = '#44cc66' if ok_geo else _FB
        cc_pwr = '#44cc66' if ok_pwr else _FB
        cc_global = '#44cc66' if ok else _FB

        self._info['lam'].config(text=f"{lam * 100:.2f} cm")
        self._info['r1' ].config(text=f"{r1_max:.3f} m")
        self._info['d3d'].config(text=f"{d_3D:.3f} m")
        self._info['d1c'].config(text=f"{d1c:.3f} m")
        self._info['d2c'].config(text=f"{d2c:.3f} m")
        self._info['min'].config(text=f"{dist_min:.3f} m")
        self._info['clr'].config(text=f"{clr_60:+.3f} m", fg=cc_geo)

        if r_visual:
            self._info['rz'].config(
                text='\n'.join(f"r{i+1}: {rv:.3f} m" for i, rv in enumerate(r_visual)))
        else:
            self._info['rz'].config(text="—")

        self._info['fspl'].config(text=f"{fspl:.2f} dB")
        if total_adif > 0.1:
            self._info['adif'].config(text=f"{total_adif:.2f} dB", fg='#ffbb55')
        else:
            self._info['adif'].config(text=f"0.00 dB", fg=_MUT)
        self._info['gtotal'].config(text=f"{gtotal:+.1f} dBi")

        # Potencia recibida
        if prx_dBm > rx_sens_dBm:
            prx_color = '#44ff66'
        elif prx_dBm > rx_sens_dBm - 10:
            prx_color = '#ffcc44'
        else:
            prx_color = '#ff6666'
        self._info['prx_dbm'].config(text=f"{prx_dBm:.1f} dBm", fg=prx_color)

        # Prx en vatios (formato legible)
        prx_watts = 10.0 ** ((prx_dBm - 30) / 10.0)
        if prx_watts >= 1.0:
            prx_w_str = f"{prx_watts:.3f} W"
        elif prx_watts >= 1e-3:
            prx_w_str = f"{prx_watts*1000:.2f} mW"
        elif prx_watts >= 1e-6:
            prx_w_str = f"{prx_watts*1e6:.2f} µW"
        else:
            prx_w_str = f"{prx_watts*1e9:.2f} nW"
        self._info['prx_w'].config(text=prx_w_str, fg=prx_color)

        # Margen de enlace
        margin_color = '#44ff66' if margin > 0 else '#ff6666'
        self._info['margin'].config(text=f"{margin:+.1f} dB", fg=margin_color)

        # Sensibilidad mostrada
        sens_display = f"{rx_sens_dBm:.1f} dBm"
        self._info['sens'].config(text=sens_display)

        # Margen de seguridad configurado
        self._info['margin_set'].config(text=f"{margin_setting:.1f} dB")

        # ── Veredicto combinado ───────────────────────────────────────────────
        if ok:
            status_txt = (f"✅ ENLACE FACTIBLE\n"
                          f"Geometría: despeje {clr_60:+.2f} m\n"
                          f"Potencia: margen {margin:+.1f} dB\n"
                          f"Prx: {prx_dBm:.1f} dBm")
        else:
            if not ok_geo and not ok_pwr:
                status_txt = (f"❌ NO FACTIBLE\n"
                              f"Falla geometría y potencia.\n"
                              f"Despeje: {clr_60:+.2f} m\n"
                              f"Margen: {margin:+.1f} dB")
            elif not ok_geo:
                status_txt = (f"❌ NO FACTIBLE (geometría)\n"
                              f"Obstáculo {worst['name']} invade zona 60%\n"
                              f"Faltan {abs(clr_60):.2f} m")
            else:  # no ok_pwr
                status_txt = (f"❌ NO FACTIBLE (potencia)\n"
                              f"Prx ({prx_dBm:.1f} dBm) por debajo de sensibilidad\n"
                              f"Faltan {abs(margin):.1f} dB")

        self._stlbl.config(text=status_txt, fg=cc_global)
        self._stlbl_border.config(bg=cc_global)

        # Tabla por obstáculo
        self._update_obs_table([
            {'name': o['name'], 'hcol': o['hcol'],
             'dist_min': o['dist_min'], 'r60': o['r60'], 'clr': o['clr'],
             'adif': o['adif']}
            for o in obs_data
        ])

        # ── Empaquetar datos para el renderizado ──────────────────────────────
        p['ax_base'] = (ax_bx, ax_by)
        p['bx_base'] = (bx_bx, bx_by)
        p['obs_data'] = obs_data
        p['worst_idx'] = obs_data.index(worst)

        fcol = _FK if ok_geo else _FB
        self._render(p, A, B, lam, fcol, r_visual, zones,
                     ok_geo, self._show_labels.get())

    # ══════════════════════════════════════════════════════════════════════════
    #  RENDERIZADO 3-D (optimizado)
    # ══════════════════════════════════════════════════════════════════════════

    def _render(self, p, A, B, lam, fcol, r_visual, zones,
                ok, show_labels=True):
        ax = self.ax
        ax.cla()
        ax.set_facecolor(_BG)
        self.fig.patch.set_facecolor(_BG)

        d_AB      = p['d_AB']
        cx        = d_AB / 2.0
        bwA       = p['bwA']
        bwB       = p['bwB']
        obs_data  = p['obs_data']
        worst_idx = p['worst_idx']
        is_xray   = self._xray_mode.get()
        alpha_ed  = 0.20 if is_xray else 0.85
        max_bw    = max(bwA, bwB, *(o['bw'] for o in obs_data))

        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pane.fill = False; pane.set_edgecolor('none')
        ax.grid(False)
        ax.xaxis.line.set_linewidth(0)
        ax.yaxis.line.set_linewidth(0)
        ax.zaxis.line.set_linewidth(0)

        draw_grid(ax, cx, 0, d_AB)

        draw_building(ax, 0.0,  0.0, bwA, {'d': bwA, 'rot': p.get('rotA', 0.0)},
                      p['hA'], _PA, alpha=alpha_ed, label='A', show_labels=show_labels)
        draw_building(ax, d_AB, 0.0, bwB, {'d': bwB, 'rot': p.get('rotB', 0.0)},
                      p['hB'], _PB, alpha=alpha_ed, label='B', show_labels=show_labels)

        for i, o in enumerate(obs_data):
            is_worst = (i == worst_idx)
            a = (alpha_ed * 1.15 if is_worst else alpha_ed)
            a = min(a, 0.95)
            draw_building(ax, o['cx'], o['cy'], o['bw'],
                          {'d': o['bw'], 'rot': o['rot']},
                          o['h'], o['pal'], alpha=a,
                          label=o['name'], show_labels=show_labels)

        ax_base = p.get('ax_base', (0.0, 0.0))
        bx_base = p.get('bx_base', (d_AB, 0.0))
        draw_antenna(ax, ax_base[0], ax_base[1], p['hA'], extra=p['hA_ant'])
        draw_antenna(ax, bx_base[0], bx_base[1], p['hB'], extra=p['hB_ant'])

        draw_cota_3d(ax, [0, bwA*0.7, 0], [0, bwA*0.7, p['hA']],
                     f"Edif. A: {p['hA']:.0f}m",
                     [0, bwA*0.6, 0], '#5588ff', show_labels=show_labels)
        draw_cota_3d(ax, [0, bwA*0.7, p['hA']],
                     [0, bwA*0.7, p['hA']+p['hA_ant']],
                     f"Mástil Tx: +{p['hA_ant']:.0f}m",
                     [0, bwA*0.6, 0], _ANT, show_labels=show_labels)
        draw_cota_3d(ax, [d_AB, bwB*0.7, 0], [d_AB, bwB*0.7, p['hB']],
                     f"Edif. B: {p['hB']:.0f}m",
                     [0, bwB*0.6, 0], '#44bb66', show_labels=show_labels)
        draw_cota_3d(ax, [d_AB, bwB*0.7, p['hB']],
                     [d_AB, bwB*0.7, p['hB']+p['hB_ant']],
                     f"Mástil Rx: +{p['hB_ant']:.0f}m",
                     [0, bwB*0.6, 0], _ANT, show_labels=show_labels)

        for o in obs_data:
            bw = o['bw']
            draw_cota_3d(ax,
                         [o['cx'], o['cy'] - bw*0.7, 0],
                         [o['cx'], o['cy'] - bw*0.7, o['h']],
                         f"Obs. {o['name']}: {o['h']:.0f}m",
                         [0, -bw*0.6, 0], o['hcol'], show_labels=show_labels)

        ax.plot([A[0], B[0]], [A[1], B[1]], [A[2], B[2]],
                '--', color=_LOS, lw=2.1, alpha=0.88, zorder=6)

        d_3D = float(np.linalg.norm(B - A))
        if d_3D > 1e-9:
            u_los = (B - A) / d_3D
            ref   = (np.array([1, 0, 0], float) if abs(u_los[0]) < 0.85
                     else np.array([0, 1, 0], float))
            v_    = np.cross(ref, u_los); v_ /= np.linalg.norm(v_)
            w_    = np.cross(u_los, v_)
        else:
            u_los = v_ = w_ = np.array([0, 0, 1], float)

        for i, o in enumerate(obs_data):
            poly_i = o['poly']
            los_pt = o['los_pt']
            is_worst = (i == worst_idx)

            if is_worst and o['dist_min'] > 1e-6:
                closest_prism = closest_point_to_prism_3d(
                    los_pt[0], los_pt[1], los_pt[2], poly_i, o['h'])
                ax.plot([los_pt[0], closest_prism[0]],
                        [los_pt[1], closest_prism[1]],
                        [los_pt[2], closest_prism[2]],
                        ':', color='#ffcc00', lw=1.8, alpha=0.9, zorder=20)
                ax.scatter([closest_prism[0]], [closest_prism[1]], [closest_prism[2]],
                           color='#ffcc00', s=28, zorder=21)
                if show_labels:
                    mid = (los_pt + closest_prism) / 2.0
                    ax.text(mid[0]+1, mid[1]+1, mid[2]+5,
                            f"δ={o['dist_min']:.2f} m",
                            color='#ffcc88', fontsize=8, ha='center', va='bottom',
                            zorder=22,
                            bbox=dict(facecolor=_BG, edgecolor='none', alpha=0.80, pad=1))

            scan_t   = np.linspace(0.0, 1.0, 300)
            pts_in   = []
            for t in scan_t:
                P = A + t * (B - A)
                if _point_in_poly(P[0], P[1], poly_i) and 0.0 <= P[2] <= o['h']:
                    pts_in.append(P)
            if pts_in:
                pts_in = np.array(pts_in)
                ax.plot(pts_in[:, 0], pts_in[:, 1], pts_in[:, 2],
                        color='#ff2222', lw=5.0, zorder=25)

            r_max_vis = max(o['r_visual']) if o['r_visual'] else 0.0
            if o['dist_min'] < r_max_vis and r_max_vis > 0:
                angles = np.linspace(0, 2 * np.pi, 120, endpoint=False)
                hits   = []
                for ang in angles:
                    pt = los_pt + r_max_vis * (math.cos(ang) * v_ +
                                               math.sin(ang) * w_)
                    if _point_in_poly(pt[0], pt[1], poly_i) and 0.0 <= pt[2] <= o['h']:
                        hits.append(pt)
                if hits:
                    hx = [h[0] for h in hits]
                    hy = [h[1] for h in hits]
                    hz = [h[2] for h in hits]
                    ax.plot(hx, hy, hz, color='#ff3333', marker='.', ms=5,
                            ls='None', zorder=30)
                    ax.plot(hx, hy, [o['h']] * len(hx),
                            color='#ff1111', lw=3.5, zorder=31)

        if self._show_ell.get():
            zone_colors = {0.60: fcol, 0.80: '#a96cff', 1.00: '#ffcc66'}
            for z in zones:
                draw_fresnel_wire(ax, A, B, lam, z,
                                  zone_colors.get(z, '#ffffff'), nu=24, nv=12)

        off_abs    = max(abs(o['cy']) for o in obs_data)
        dim_y      = -(max(16, max_bw * 1.7))
        grid_color = '#bfbfbf'

        ax.plot([0.0, 0.0],   [dim_y, dim_y+8], [0, 0], color=grid_color, lw=1.2, alpha=0.95, zorder=4)
        ax.plot([d_AB, d_AB], [dim_y, dim_y+8], [0, 0], color=grid_color, lw=1.2, alpha=0.95, zorder=4)
        ax.plot([0.0, d_AB],  [dim_y, dim_y],   [0, 0], color=grid_color, lw=1.6, alpha=0.95, zorder=4)
        ax.scatter([0.0, d_AB], [dim_y, dim_y], [0, 0], color=grid_color, s=18, zorder=5)
        if show_labels:
            ax.text(d_AB/2, dim_y-6, 0, f"d_AB = {d_AB:.0f} m",
                    color='#ffffff', fontsize=8, ha='center', va='top', zorder=5,
                    bbox=dict(facecolor=_BG, edgecolor=grid_color,
                              boxstyle='round,pad=0.22', alpha=0.88))

        for i, o in enumerate(obs_data):
            dy_cota = dim_y + 14 + i * 14
            hcol    = o['hcol']
            ax.plot([0, o['cx']], [dy_cota, dy_cota], [0, 0],
                    color=hcol, lw=1.4, alpha=0.85, zorder=4)
            ax.scatter([0, o['cx']], [dy_cota, dy_cota], [0, 0],
                       color=hcol, s=14, zorder=5)
            if show_labels:
                ax.text(o['cx']/2, dy_cota+6, 0,
                        f"d{o['name']}={o['cx']:.0f}m",
                        color='#ffffff', fontsize=7, ha='center', va='bottom', zorder=5,
                        bbox=dict(facecolor=_BG, edgecolor=hcol,
                                  boxstyle='round,pad=0.18', alpha=0.85))

        ax.set_xlabel(''); ax.set_ylabel(''); ax.set_zlabel(''); ax.set_title('')
        ax.tick_params(axis='both', colors='#3e3e3e', labelsize=6)
        ax.tick_params(axis='z',    colors='#3e3e3e', labelsize=6)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])

        zmax = max(A[2], B[2], *(o['h'] for o in obs_data)) + 28
        ymg  = max(45, off_abs + max_bw * 0.6 + 18)
        zoom = self._vars['zoom'].get()

        if self._vars['uniform_scale'].get():
            axis_range = max(d_AB + 60.0, ymg * 2.0, zmax)
            half = axis_range / 2.0 / zoom
            ax.set_xlim(cx - half, cx + half)
            ax.set_ylim(-half, half)
            ax.set_zlim(0, axis_range / zoom)
        else:
            half_x = (d_AB / 2.0 + 30) / zoom
            ax.set_xlim(cx - half_x, cx + half_x)
            ax.set_ylim(-ymg / zoom, ymg / zoom)
            ax.set_zlim(0, zmax / zoom)

        handles = [
            mpatches.Patch(facecolor=_PA[0], edgecolor=_PA[3], label='Edificio A (Tx)'),
            mpatches.Patch(facecolor=_PB[0], edgecolor=_PB[3], label='Edificio B (Rx)'),
        ]
        for o in obs_data:
            handles.append(mpatches.Patch(
                facecolor=o['pal'][0], edgecolor=o['pal'][3],
                label=f"Obstáculo {o['name']} "
                      f"({'✅' if o['clr']>0 else '❌'} {o['clr']:+.1f}m)"))
        handles.append(Line2D([0], [0], color=_LOS, ls='--', lw=1.8, label='LOS'))

        if self._show_ell.get():
            zone_colors = {0.60: fcol, 0.80: '#a96cff', 1.00: '#ffcc66'}
            zone_labels = {0.60: 'Fresnel 60%', 0.80: 'Fresnel 80%', 1.00: 'Fresnel 100%'}
            for z in zones:
                if z in zone_labels:
                    handles.append(Line2D([0], [0], color=zone_colors.get(z, '#fff'),
                                          lw=2.0, label=zone_labels[z]))

        ax.legend(handles=handles, loc='upper left', fontsize=7,
                  facecolor=_PNL, edgecolor=_BDR, labelcolor=_TXT, framealpha=0.92)

        self.cvs.draw()

# ══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    root = tk.Tk()
    app  = FresnelApp(root)
    root.mainloop()