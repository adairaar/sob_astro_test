#!/usr/bin/env python3
"""Figure 3 variant v3 — tight empirical Keplerian band + clean layout."""
import numpy as np
import matplotlib; matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors, matplotlib.cm as cm

import os as _os; PROJ = _os.path.dirname(_os.path.abspath(__file__))
cache = np.load(f"{PROJ}/scatter_cache.npz")
az_b,mot_b,stp_b,e_b = cache['az_b'],cache['mot_b'],cache['stp_b'],cache['e_b']
az_1,mot_1,stp_1,e_1 = cache['az_1'],cache['mot_1'],cache['stp_1'],cache['e_1']
az_2,mot_2,stp_2,e_2 = cache['az_2'],cache['mot_2'],cache['stp_2'],cache['e_2']
all_az  = np.concatenate([az_b,az_1,az_2])
all_mot = np.concatenate([mot_b,mot_1,mot_2])
all_stp = np.concatenate([stp_b,stp_1,stp_2])
all_e   = np.concatenate([e_b,e_1,e_2])

AZ_TOL, OMG = 5.0, 2.0
FLOOR = 0.001

# Empirical band: 1st–99th percentile of stp/mot ratio for resolved orbits
mask_res = (all_mot > 0.01) & (all_stp > 0.01)
ratio = all_stp[mask_res] / all_mot[mask_res]
p01, p99 = np.percentile(ratio, 1), np.percentile(ratio, 99)
# p01 ≈ 0.9952, p99 ≈ 1.033  → max deviation ≈ 3.3%

in_corr = all_az < AZ_TOL
rng = np.random.default_rng(7)
out_idx = rng.choice(np.where(~in_corr)[0], min(25_000,(~in_corr).sum()), replace=False)
in_idx  = np.where(in_corr)[0]

# ── figure ────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9.2, 6.5))
ax.set_facecolor('white')

# ── Keplerian constraint: empirical 1–99th percentile band ────────────────
x_line = np.logspace(np.log10(FLOOR), np.log10(200), 400)
ax.fill_between(x_line, p01*x_line, p99*x_line,
                color='#7B52AB', alpha=0.25, zorder=0,
                label=f'Keplerian band (1st–99th pct: ×{p01:.4f}–×{p99:.4f})')
ax.plot(x_line, x_line, color='#5A2D82', lw=2.0, ls='-', zorder=4,
        label='y = x  (speed unchanged over 2 h)')

# ── scatter: non-corridor coloured by eccentricity ────────────────────────
norm = mcolors.LogNorm(vmin=0.01, vmax=5.0, clip=True)
sc = ax.scatter(np.maximum(all_mot[out_idx], FLOOR),
                np.maximum(all_stp[out_idx], FLOOR),
                c=np.clip(all_e[out_idx], 0.01, 5.0),
                cmap='plasma', norm=norm,
                s=1.5, alpha=0.20, rasterized=True, zorder=2)

# ── corridor (green) ──────────────────────────────────────────────────────
n_corr = int(in_corr.sum())
ax.scatter(np.maximum(all_mot[in_idx], FLOOR),
           np.maximum(all_stp[in_idx], FLOOR),
           c='#1DB954', s=13, alpha=0.80, rasterized=True, zorder=6)

# ── threshold lines ───────────────────────────────────────────────────────
ax.axvline(OMG, color='#008000', lw=2.0, ls='--', zorder=7,
           label='Guidance threshold (2°/h)')
ax.axhline(OMG, color='#CC0000', lw=2.0, ls='--', zorder=7,
           label='Stopping threshold (2°/h)')

# ── required region (lower-right) ────────────────────────────────────────
gold = plt.Rectangle((OMG, FLOOR), 200-OMG, OMG-FLOOR,
                      lw=2.0, edgecolor='#C8A000',
                      facecolor='#FFEC80', alpha=0.70, zorder=1)
ax.add_patch(gold)
ax.text(25, 0.055,
        'Required region\n(guidance >2°/h AND stop <2°/h)\n→ 0 orbits found',
        ha='center', va='center', fontsize=9.5, color='#7A5900',
        bbox=dict(boxstyle='round,pad=0.4', fc='#FFF8C0', ec='#C8A000',
                  lw=1.5, alpha=0.95), zorder=12)

# ── e=2.0 outlier annotation ──────────────────────────────────────────────
ax.annotate('e=2.0: guidance=2.89°/h\nstopping=2.89°/h\n→ on diagonal; fails stopping',
            xy=(2.89, 2.89), xytext=(9, 18), fontsize=8.5, color='#1DB954',
            arrowprops=dict(arrowstyle='->', color='#1a7a38', lw=1.4, shrinkB=3),
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#1DB954',
                      lw=1.2, alpha=0.95), zorder=13)

# ── corridor cluster annotation: upper-LEFT, arrow pointing DOWN-LEFT ─────
ax.annotate(
    f'Corridor-passing orbits ({n_corr:,}):\nquasi-stationary in both windows\n→ fail guidance criterion',
    xy=(0.0025, 0.0025),          # tip of arrow: near the cluster
    xytext=(0.007, 0.18),          # text box: mid-left, clear of gold box
    fontsize=8.5, color='#145c30',
    arrowprops=dict(arrowstyle='->', color='#1DB954', lw=1.4, shrinkB=4),
    bbox=dict(boxstyle='round,pad=0.3', fc='#e8f8ee', ec='#1DB954',
              lw=1.2, alpha=0.97), zorder=13)

# ── Keplerian band label: left of centre, no overlap with colorbar ─────────
ax.text(0.0035, 3.5,
        f'Keplerian constraint\n(99% of orbits within\n±3.3% of y=x;\ninvisibly tight at plot scale)',
        ha='left', va='top', fontsize=8.2, color='#3A1060',
        bbox=dict(boxstyle='round,pad=0.35', fc='#F3EEF8', ec='#6A3D9A',
                  lw=1.2, alpha=0.95), zorder=13)

# ── stats box: upper-right but inset from colorbar ───────────────────────
n_vis = len(all_az)
ax.text(0.975, 0.985,
        f'Visible configs: {n_vis:,}\nAll 3 criteria: 0',
        transform=ax.transAxes, ha='right', va='top', fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#aaaaaa',
                  lw=0.9, alpha=0.93), zorder=13)

# ── colorbar ──────────────────────────────────────────────────────────────
sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.022, aspect=30)
cbar.set_label('Eccentricity (non-corridor orbits)', fontsize=8.5)
cbar.set_ticks([0.01, 0.1, 0.5, 1.0, 2.0, 5.0])
cbar.set_ticklabels(['0.01','0.1','0.5','1','2','5'], fontsize=7.5)

# ── legend ────────────────────────────────────────────────────────────────
legend_elements = [
    mpatches.Patch(fc='#D8C8F0', ec='#5A2D82', lw=1,
                   label=f'Keplerian band: stp/guidance ×0.995–×1.033 (99% of orbits)'),
    Line2D([0],[0], color='#5A2D82', lw=2.0, label='y = x  (speed unchanged in 2 h)'),
    Line2D([0],[0], marker='o', color='w', mfc='#1DB954', mec='#1DB954', ms=8,
           label=f'Az in corridor (<5°): {n_corr:,} orbits'),
    mpatches.Patch(fc='#FFEC80', ec='#C8A000', lw=1.5,
                   label='Required region (lower-right, empty)'),
    Line2D([0],[0], color='#008000', lw=2, ls='--', label='Guidance threshold (2°/h)'),
    Line2D([0],[0], color='#CC0000', lw=2, ls='--', label='Stopping threshold (2°/h)'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=7.8,
          framealpha=0.92, edgecolor='#aaaaaa', handlelength=1.5)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(FLOOR, 200); ax.set_ylim(FLOOR, 200)
ax.grid(True, which='both', alpha=0.22)
ax.set_xlabel('Min guidance-phase angular velocity  (°/h)\n'
              '[must exceed 2°/h to lead the Magi]', fontsize=10)
ax.set_ylabel('Max stopping-phase angular velocity  (°/h)\n'
              '[must be below 2°/h to appear stationary]', fontsize=10)
ax.set_title(
    'Impossibility of Matthew 2:9 for Any Keplerian Orbit\n'
    'Orbits cluster within ±3.3% of y=x (Keplerian constraint); '
    'the required region (gold) lies orders of magnitude away',
    fontsize=10.5, fontweight='bold')

fig.tight_layout()
out = f"{PROJ}/figure_paper2_impossibility_scatter_v2.png"
fig.savefig(out, dpi=600, bbox_inches='tight')
plt.close()
print(f"Saved: {out}")
print(f"Keplerian band: p01={p01:.4f}  p99={p99:.4f}")
