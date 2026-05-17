#!/usr/bin/env python3
"""
Generate Paper Figure 1: Observational constraints from Matthew 2:9
  Panel (a): All-sky polar chart (altazimuth) showing the guidance-phase
             azimuth corridor [190°, 210°] and the stopping-point star symbol.
  Panel (b): Azimuth corridor diagram from Bethlehem — road bearing 203°,
             corridor sweep 190°–210°, with Jerusalem marked to the NNE.

Output: figure_paper1_sky_geometry.png  (600 DPI)

Aaron Adair / 2026
"""
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

OUTDIR = os.path.dirname(os.path.abspath(__file__)) + os.sep

# ── figure layout ────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 5.8))
ax1 = fig.add_subplot(121, projection='polar')
ax2 = fig.add_subplot(122)

# ── colour constants ─────────────────────────────────────────────────────────
C_CORRIDOR = '#2ca02c'   # green
C_STAR     = '#1f77b4'   # blue  (guidance track)
C_STOP     = '#d62728'   # red   (stopping symbol)
C_ROAD     = '#8c564b'   # brown

# ════════════════════════════════════════════════════════════════════════════
# Panel (a): polar sky chart
#   theta = azimuth (N=top, clockwise)
#   r     = (90 - altitude) / 90   →  zenith=0, horizon=1
# ════════════════════════════════════════════════════════════════════════════
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1)           # clockwise azimuth
ax1.set_rlim(0, 1)

# altitude circles at 30° and 60°
for alt_deg in [30, 60]:
    r_circ = 1.0 - alt_deg / 90.0
    th = np.linspace(0, 2 * np.pi, 360)
    ax1.plot(th, np.full_like(th, r_circ), '-', color='grey',
             lw=0.6, alpha=0.4, zorder=1)

# horizon
th = np.linspace(0, 2 * np.pi, 360)
ax1.plot(th, np.ones_like(th), 'k-', lw=1.2, zorder=2)

# shaded corridor [190°, 210°]
az_lo, az_hi = np.radians(190), np.radians(210)
az_fill = np.linspace(az_lo, az_hi, 120)
r_fill  = np.ones(120)
th_patch = np.concatenate([[0.0], az_fill, az_fill[::-1], [0.0]])
r_patch  = np.concatenate([[0.0], r_fill,  np.zeros(120), [0.0]])
ax1.fill(th_patch, r_patch, color=C_CORRIDOR, alpha=0.22, zorder=3)

# corridor boundary lines
for az_b in [np.radians(190), np.radians(210)]:
    ax1.plot([az_b, az_b], [0, 1], '--', color=C_CORRIDOR,
             lw=1.3, alpha=0.7, zorder=4)

# central bearing 203°
ax1.plot([np.radians(203), np.radians(203)], [0, 1],
         '-', color=C_CORRIDOR, lw=1.8, alpha=0.9, zorder=4)

# guidance track: arrow from alt≈52° down to alt≈22° along az=201°
az_tr = np.radians(201)
r_start = 1.0 - 52.0 / 90.0    # ≈ 0.422
r_end   = 1.0 - 22.0 / 90.0    # ≈ 0.756
ax1.annotate('',
             xy=(az_tr, r_end), xytext=(az_tr, r_start),
             arrowprops=dict(
                 arrowstyle='->', color=C_STAR,
                 lw=2.2, mutation_scale=18),
             zorder=6)
ax1.plot(az_tr, r_start, 'o', color=C_STAR, ms=7, zorder=7)

# label for guidance track
ax1.text(np.radians(168), 0.62,
         'Guidance phase\n(προῆγεν)', ha='center', va='center',
         fontsize=8.0, color=C_STAR,
         bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                   edgecolor=C_STAR, alpha=0.9), zorder=8)

# stopping point: star symbol near horizon along 203°
r_stop = 1.0 - 8.0 / 90.0     # ≈ 0.911
ax1.plot(np.radians(203), r_stop, '*', color=C_STOP,
         ms=18, mec='darkred', mew=0.8, zorder=9)
ax1.text(np.radians(228), r_stop + 0.06,
         'Stopping\n(ἐστάθη)', ha='center', va='center',
         fontsize=8.0, color=C_STOP,
         bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                   edgecolor=C_STOP, alpha=0.9), zorder=8)

# altitude labels on radial axis
ax1.set_yticks([0.0, 1.0 - 60/90, 1.0 - 30/90, 1.0])
ax1.set_yticklabels(['90°', '60°', '30°', '0°\n(horiz.)'],
                    fontsize=8, color='grey')

# cardinal direction labels
for az_deg, label in [(0,'N'),(90,'E'),(180,'S'),(270,'W')]:
    ax1.text(np.radians(az_deg), 1.17, label,
             ha='center', va='center', fontsize=10, fontweight='bold')

# corridor annotation on plot
ax1.text(np.radians(203), 0.30,
         '[190°–210°]\ncorridor', ha='center', va='center',
         fontsize=7.5, color=C_CORRIDOR,
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8),
         zorder=8)

ax1.set_title('(a)', fontsize=12, fontweight='bold', pad=16, loc='left')

# ════════════════════════════════════════════════════════════════════════════
# Panel (b): azimuth corridor diagram from Bethlehem (cartesian, top-down)
#   North = up, East = right
#   x = sin(az), y = cos(az)   (standard bearing convention)
# ════════════════════════════════════════════════════════════════════════════
ax2.set_aspect('equal')
ax2.set_xlim(-1.5, 1.5)
ax2.set_ylim(-1.5, 1.5)
ax2.axis('off')

# horizon circle
th_c = np.linspace(0, 2 * np.pi, 360)
ax2.plot(np.cos(th_c), np.sin(th_c), 'k-', lw=1.5)

# shaded corridor wedge  [190°, 210°]
az_wedge = np.linspace(np.radians(190), np.radians(210), 120)
xw = np.concatenate([[0], np.sin(az_wedge), [0]])
yw = np.concatenate([[0], np.cos(az_wedge), [0]])
ax2.fill(xw, yw, color=C_CORRIDOR, alpha=0.25)

# corridor boundary lines
for az_b in [190, 210]:
    xb, yb = np.sin(np.radians(az_b)), np.cos(np.radians(az_b))
    ax2.plot([0, xb], [0, yb], '--', color=C_CORRIDOR, lw=1.5, alpha=0.75)
    # degree labels at the horizon
    ax2.text(xb * 1.15, yb * 1.15, f'{az_b}°',
             ha='center', va='center', fontsize=8.5, color=C_CORRIDOR)

# central bearing (straight-line road direction, 203°)
xc, yc = np.sin(np.radians(203)), np.cos(np.radians(203))
ax2.annotate('', xy=(xc, yc), xytext=(0, 0),
             arrowprops=dict(arrowstyle='->', color=C_ROAD,
                             lw=2.0, mutation_scale=14))
ax2.text(xc * 0.58, yc * 0.58 - 0.13,
         '203°\n(straight-\nline bearing)',
         ha='center', va='top', fontsize=8.0, color=C_ROAD,
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85))

# Bethlehem observer at centre
ax2.plot(0, 0, 'ko', ms=9, zorder=5)
ax2.text(0.10, -0.13, 'Bethlehem\n(observer)',
         ha='left', va='top', fontsize=8.5, fontweight='bold')

# Jerusalem: NNE of Bethlehem, bearing ~23°, distance ≈ 0.58 (normalised)
j_az = 23.0
jx, jy = np.sin(np.radians(j_az)) * 0.58, np.cos(np.radians(j_az)) * 0.58
ax2.plot(jx, jy, 's', color='steelblue', ms=9, zorder=5)
ax2.text(jx + 0.09, jy + 0.04, 'Jerusalem',
         ha='left', va='bottom', fontsize=8.5)

# road path: schematic curve from Bethlehem toward Jerusalem,
# deviating around a ridge to sweep through the corridor
t_road = np.linspace(0, 1, 80)
# parameterise as slight S-curve starting at 203° and ending at 23° (reversed)
# reversed: road goes FROM Bethlehem (0,0) TOWARD Jerusalem
az_road_deg = 23.0 + (203.0 - 23.0) * (1 - t_road)  # 23 → 203 in reversed direction
# actually we want the road from Bethlehem toward Jerusalem
# road goes north (bearing ~23°) but deviates to accommodate terrain
# we model the portion nearest Bethlehem sweeping 190°-210° (looking south from Jerusalem)
# so from Bethlehem NORTHWARD, the road sweeps roughly 0°-30° azimuth
# We just want a schematic gentle curve between the two points
n_pts = 60
s = np.linspace(0, 1, n_pts)
# road x,y from (0,0) to (jx,jy) with a slight curve
ctrl_x, ctrl_y = jx * 0.4 + 0.18, jy * 0.4  # control point for quadratic bezier
rx = (1-s)**2 * 0 + 2*(1-s)*s * ctrl_x + s**2 * jx
ry = (1-s)**2 * 0 + 2*(1-s)*s * ctrl_y + s**2 * jy
ax2.plot(rx, ry, '-', color=C_ROAD, lw=2.2, alpha=0.85, zorder=4)

# North arrow (top of diagram)
ax2.annotate('', xy=(0, 1.35), xytext=(0, 1.10),
             arrowprops=dict(arrowstyle='->', color='black', lw=1.5,
                             mutation_scale=12))
ax2.text(0, 1.42, 'N', ha='center', va='bottom',
         fontsize=10, fontweight='bold')

# corridor annotation
ax2.text(np.sin(np.radians(200)) * 0.80,
         np.cos(np.radians(200)) * 0.80,
         'Azimuth\ncorridor\n190°–210°',
         ha='center', va='center', fontsize=8.0, color=C_CORRIDOR,
         bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                   edgecolor=C_CORRIDOR, alpha=0.9), zorder=8)

ax2.set_title('(b)', fontsize=12, fontweight='bold', loc='left')

# ── save ─────────────────────────────────────────────────────────────────────
plt.tight_layout(pad=1.5)
outpath = OUTDIR + 'figure_paper1_sky_geometry.png'
fig.savefig(outpath, dpi=600, bbox_inches='tight')
plt.close()
print(f"  Saved figure_paper1_sky_geometry.png  (600 DPI)")
print("Done.")
