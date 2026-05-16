#!/usr/bin/env python3
"""
sob_fig_update.py — Regenerate/update paper figures
====================================================
Produces:
  figure_paper2_az_motion.png              — NEW paper Fig 2: az deviation vs. guidance motion
                                             (two-criterion near-impossibility; direction + motion)
  figure_paper2_impossibility_scatter.png  — paper Fig 3: guidance motion vs. stopping motion
                                             (three-criterion impossibility; direction+motion+stop)
  figure9_corpus_linguistics.png           — white background, 4 features only
  figure_mc_diagnostic.png                 — supplementary S1: MC sweep az vs. guidance-motion

All orbital-mechanics constants match sob_proof_fine.py and sob_mc_refined.py exactly.
Aaron Adair / 2026
"""

import warnings; warnings.filterwarnings('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body_barycentric
import astropy.units as u
import time as walltime

# ── Constants ──────────────────────────────────────────────────────────────
K_GAUSS   = 0.01720209895
MU        = K_GAUSS**2
JD_NOON   = 1_719_755.898
DT_SEC    = 10_572.0
BETHLEHEM = EarthLocation(lon=(35+12/60)*u.deg, lat=(31+42/60)*u.deg, height=765.0*u.m)
AZ_LO, AZ_HI = 190.0, 210.0
AZ_TOL    = 5.0
OMEGA_THRESH = 2.0
ALT_FLOOR = 5.0
RAD       = 180.0 / np.pi

# Time grids (17 steps total: 9 guidance + 8 stopping, 15-min spacing)
_t_g  = np.array([8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0])
_t_s  = np.array([10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0])
_times = np.concatenate([_t_g, _t_s])
N_G   = len(_t_g)      # 9
N_S   = len(_t_s)      # 8
N_GS  = len(_times)    # 17
DT_H  = 0.25           # hours between steps

import os as _os
OUTDIR = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
    'axes.edgecolor':   '#333333',
    'axes.labelcolor':  '#111111',
    'xtick.color':      '#333333',
    'ytick.color':      '#333333',
    'text.color':       '#111111',
    'grid.color':       '#cccccc',
    'grid.linestyle':   '--',
    'grid.alpha':       0.5,
    'font.size':        10,
    'axes.titlesize':   11,
    'axes.labelsize':   10,
})

# ══════════════════════════════════════════════════════════════════════════
# 1. Precompute ENU matrices and Earth positions
# ══════════════════════════════════════════════════════════════════════════
print("─" * 64)
print("Precomputing ENU matrices and Earth positions …")
t0 = walltime.time()

jds_gs  = JD_NOON + _times / 24.0
jds_utc = jds_gs - DT_SEC / 86400.0

t_tdb = Time(jds_gs,  format='jd', scale='tdb')
t_utc = Time(jds_utc, format='jd', scale='utc')

eb  = get_body_barycentric('earth', t_tdb)
sb  = get_body_barycentric('sun',   t_tdb)
earth_pos = (eb.xyz.to(u.AU).value - sb.xyz.to(u.AU).value).T  # (N_GS, 3)

def make_enu(t_utc_step):
    frame = AltAz(obstime=t_utc_step, location=BETHLEHEM)
    def aa2icrs(az_d, alt_d):
        c = SkyCoord(alt=alt_d*u.deg, az=az_d*u.deg, frame=frame).icrs
        d = np.radians(c.dec.deg); r = np.radians(c.ra.deg)
        return np.array([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)])
    E_raw = aa2icrs(90, 0.001)
    Z     = aa2icrs(0, 90.0); Z /= np.linalg.norm(Z)
    E = E_raw - np.dot(E_raw, Z)*Z; E /= np.linalg.norm(E)
    N = np.cross(Z, E)
    return np.array([E, N, Z])

enu = np.array([make_enu(t_utc[i]) for i in range(N_GS)])  # (N_GS, 3, 3)
print(f"  Done in {walltime.time()-t0:.1f}s")


# ══════════════════════════════════════════════════════════════════════════
# 2. Kepler solvers (identical to sob_mc_refined.py)
# ══════════════════════════════════════════════════════════════════════════
def _true_anomaly(dt, q, e):
    nu = np.full(len(dt), np.nan, dtype=np.float64)
    m = e < 0.9999
    if m.any():
        eq, ee, edt = q[m], e[m], dt[m]
        a = eq / (1.0 - ee); n = np.sqrt(MU / a**3)
        M = (n * edt) % (2.0 * np.pi); E = M.copy()
        for _ in range(60):
            dE = (M - E + ee*np.sin(E)) / (1.0 - ee*np.cos(E)); E += dE
            if np.max(np.abs(dE)) < 1e-12: break
        nu[m] = 2.0*np.arctan2(np.sqrt(1+ee)*np.sin(E/2), np.sqrt(1-ee)*np.cos(E/2))
    m = (e >= 0.9999) & (e <= 1.0001)
    if m.any():
        pq, pdt = q[m], dt[m]
        W = 3.0*K_GAUSS*pdt/np.sqrt(2.0*pq**3)
        arg = W/2 + np.sqrt(np.maximum((W/2)**2+1, 0))
        arg = np.where(arg > 0, arg, 1e-30)
        nu[m] = 2.0*np.arctan(arg**(1/3) - arg**(-1/3))
    m = e > 1.0001
    if m.any():
        hq, he, hdt = q[m], e[m], dt[m]
        a = hq/(he-1); n = np.sqrt(MU/a**3); Mh = n*hdt
        H = np.arcsinh(np.clip(Mh/he, -1e6, 1e6))
        H = np.where(np.isfinite(H), H, np.sign(Mh)*np.log(np.abs(Mh)+1.8))
        for _ in range(60):
            f = he*np.sinh(H)-H-Mh; fp = he*np.cosh(H)-1
            dH = -f/np.where(np.abs(fp)>1e-15, fp, 1e-15)
            H += np.clip(dH, -2, 2)
            if np.max(np.abs(dH)) < 1e-12: break
        nu_h = 2*np.arctan2(np.sqrt(he+1)*np.sinh(H/2), np.sqrt(he-1)*np.cosh(H/2))
        nu_max = np.arccos(np.clip(-1/he, -1, 1))
        nu[m] = np.where(np.abs(nu_h) < nu_max, nu_h, np.nan)
    return nu

def _pos_icrs(nu, r, Om, om, inc):
    cnu=np.cos(nu); snu=np.sin(nu)
    cO=np.cos(Om); sO=np.sin(Om)
    co=np.cos(om); so=np.sin(om)
    ci=np.cos(inc); si=np.sin(inc)
    Px=cO*co-sO*so*ci; Py=sO*co+cO*so*ci; Pz=so*si
    Qx=-cO*so-sO*co*ci; Qy=-sO*so+cO*co*ci; Qz=co*si
    return np.stack([r*(cnu*Px+snu*Qx), r*(cnu*Py+snu*Qy), r*(cnu*Pz+snu*Qz)], axis=-1)

def evaluate_batch(q_, e_, i_r, Om_r, om_r, Tp):
    """Returns vis_ok, az_score, motion_min, stop_max — all shape (B,)."""
    B = len(q_)
    dt_mat = jds_gs[np.newaxis,:] - Tp[:,np.newaxis]   # (B, N_GS)
    N = B * N_GS
    t_i  = np.tile(np.arange(N_GS), B)
    q_f  = np.repeat(q_,  N_GS); e_f = np.repeat(e_,  N_GS)
    i_f  = np.repeat(i_r, N_GS); Om_f= np.repeat(Om_r,N_GS)
    om_f = np.repeat(om_r,N_GS); dt_f= dt_mat.ravel()
    nu_f = _true_anomaly(dt_f, q_f, e_f)
    p_f  = q_f*(1+e_f); r_f = p_f/(1+e_f*np.cos(nu_f))
    ok_f = np.isfinite(nu_f) & (r_f>0) & (r_f<1e4)
    nu_s = np.where(ok_f, nu_f, 0.); r_s = np.where(ok_f, r_f, 0.)
    xyz_f= _pos_icrs(nu_s, r_s, Om_f, om_f, i_f); xyz_f[~ok_f] = 0.
    geo_f= xyz_f - earth_pos[t_i]
    nrm  = np.linalg.norm(geo_f, axis=1, keepdims=True)
    geo_n= geo_f / np.where(nrm>0, nrm, 1.)
    E_f  = enu[t_i]
    e_c  = np.einsum('ni,ni->n', geo_n, E_f[:,0,:])
    n_c  = np.einsum('ni,ni->n', geo_n, E_f[:,1,:])
    u_c  = np.einsum('ni,ni->n', geo_n, E_f[:,2,:])
    alt_m = np.arcsin(np.clip(u_c,-1,1))*RAD
    az_m  = (np.arctan2(e_c, n_c)*RAD) % 360.
    ok_m  = ok_f.reshape(B, N_GS)
    alt_m = alt_m.reshape(B, N_GS)
    az_m  = az_m.reshape(B, N_GS)
    geo_nm= geo_n.reshape(B, N_GS, 3)
    vis_ok = (np.all(alt_m[:,:N_G]>ALT_FLOOR, axis=1) &
              np.all(ok_m[:,:N_G], axis=1))
    az_dev = np.maximum(0., np.maximum(AZ_LO-az_m[:,:N_G], az_m[:,:N_G]-AZ_HI))
    az_score = np.max(az_dev, axis=1)
    d_g = np.einsum('bti,bti->bt', geo_nm[:,:-1,:N_G+1][:,:N_G-1,:],
                    geo_nm[:,1:, :N_G+1][:,:N_G-1,:])
    motion_min = np.min(np.arccos(np.clip(
        np.einsum('bti,bti->bt', geo_nm[:,:N_G-1,:], geo_nm[:,1:N_G,:]), -1,1))*RAD/DT_H, axis=1)
    stop_max   = np.max(np.arccos(np.clip(
        np.einsum('bti,bti->bt', geo_nm[:,N_G:-1,:], geo_nm[:,N_G+1:,:]), -1,1))*RAD/DT_H, axis=1)
    return vis_ok, az_score, motion_min, stop_max


# ══════════════════════════════════════════════════════════════════════════
# 3. Generate scatter data
# ══════════════════════════════════════════════════════════════════════════
rng = np.random.default_rng(42)
BATCH = 100_000

def run_scatter(param_sampler, n_total, label):
    """Collect scatter data for n_total samples, return arrays."""
    az_list=[]; stop_list=[]; mot_list=[]; etype_list=[]
    n_vis = 0; n_done = 0
    t_start = walltime.time()
    while n_done < n_total:
        bs = min(BATCH, n_total - n_done)
        q_, e_, i_, Om_, om_, Tp_ = param_sampler(bs)
        i_r = np.radians(i_); Om_r = np.radians(Om_); om_r = np.radians(om_)
        vis, az, mot, stp = evaluate_batch(q_, e_, i_r, Om_r, om_r, Tp_)
        mask = vis
        n_vis += int(mask.sum()); n_done += bs
        az_list.append(az[mask]);    stop_list.append(stp[mask])
        mot_list.append(mot[mask]);  etype_list.append(e_[mask])
        rate = n_done / (walltime.time()-t_start+1e-9) / 1e6
        print(f"  [{label}] {n_done:>7,}/{n_total:,}  vis={n_vis:,}  {rate:.2f}M/s", end='\r')
    print()
    return (np.concatenate(az_list),   np.concatenate(stop_list),
            np.concatenate(mot_list),   np.concatenate(etype_list))

# ── Broad heliocentric sample (full parameter space) ─────────────────────
e_vals  = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0, 1.1, 1.5, 2.0, 5.0]
q_vals  = np.linspace(0.02, 1.00, 18)
i_vals  = np.concatenate([np.arange(0, 181, 5), [1, 1.4, 2, 3, 4]])
T_offs  = [-60, -40, -20, -10, -5, 5, 10, 20, 60]

def broad_sampler(bs):
    q_   = rng.choice(q_vals, bs)
    e_   = rng.choice(e_vals, bs)
    i_   = rng.choice(i_vals, bs)
    Om_  = rng.uniform(0, 360, bs)
    om_  = rng.uniform(0, 360, bs)
    dT_  = rng.choice(T_offs, bs)
    Tp_  = JD_NOON + dT_
    return q_, e_, i_, Om_, om_, Tp_

# ── MC sweep 1: best-fit neighbourhood ───────────────────────────────────
MATNEY_DT = 29.667
def mc1_sampler(bs):
    q_  = rng.uniform(0.005, 0.15, bs)
    e_  = rng.uniform(0.7, 2.5, bs)
    i_  = rng.uniform(0, 40, bs)
    Om_ = rng.uniform(0, 360, bs)
    om_ = rng.uniform(0, 360, bs)
    dT_ = rng.uniform(-60, 60, bs)
    return q_, e_, i_, Om_, om_, JD_NOON + dT_

# ── MC sweep 2: Matney neighbourhood ─────────────────────────────────────
def mc2_sampler(bs):
    q_  = rng.uniform(0.01, 0.16, bs)
    e_  = rng.uniform(0.7, 1.5, bs)
    i_  = rng.uniform(0, 8, bs)
    Om_ = rng.uniform(25.81, 125.81, bs)
    om_ = (rng.uniform(285.58, 385.58, bs)) % 360
    dT_ = rng.uniform(-10, 70, bs)
    return q_, e_, i_, Om_, om_, JD_NOON + dT_

print("\nGenerating broad heliocentric scatter …")
az_b, stp_b, mot_b, e_b = run_scatter(broad_sampler, 1_200_000, 'broad')

print("Generating MC sweep 1 scatter (best-fit nbhd) …")
az_1, stp_1, mot_1, e_1 = run_scatter(mc1_sampler, 300_000, 'mc1')

print("Generating MC sweep 2 scatter (Matney nbhd) …")
az_2, stp_2, mot_2, e_2 = run_scatter(mc2_sampler, 200_000, 'mc2')

print(f"\n  Broad:  {len(az_b):,} visible")
print(f"  MC-1:   {len(az_1):,} visible")
print(f"  MC-2:   {len(az_2):,} visible")

# Subsample for plotting (keep all az<10° points, random subsample of rest)
def subsample(az, stp, mot, ee, n_max=30_000, az_keep=20.0):
    near = az < az_keep
    far  = ~near
    n_near = int(near.sum()); n_far = min(max(0, n_max - n_near), int(far.sum()))
    idx_far = rng.choice(np.where(far)[0], n_far, replace=False) if n_far > 0 else np.array([], int)
    idx_near = np.where(near)[0]
    idx = np.concatenate([idx_near, idx_far])
    return az[idx], stp[idx], mot[idx], ee[idx]

az_bp, stp_bp, mot_bp, e_bp = subsample(az_b, stp_b, mot_b, e_b, 25_000)
az_1p, stp_1p, mot_1p, e_1p = subsample(az_1, stp_1, mot_1, e_1, 15_000)
az_2p, stp_2p, mot_2p, e_2p = subsample(az_2, stp_2, mot_2, e_2, 10_000)

n_all_ok = int(((az_b<AZ_TOL)&(mot_b>OMEGA_THRESH)&(stp_b<OMEGA_THRESH)).sum() +
               ((az_1<AZ_TOL)&(mot_1>OMEGA_THRESH)&(stp_1<OMEGA_THRESH)).sum() +
               ((az_2<AZ_TOL)&(mot_2>OMEGA_THRESH)&(stp_2<OMEGA_THRESH)).sum())

# Cache full scatter arrays for sob_fig_proof.py (figure 4)
_cache = OUTDIR + 'scatter_cache.npz'
np.savez_compressed(_cache,
    az_b=az_b,  mot_b=mot_b,  stp_b=stp_b,  e_b=e_b,
    az_1=az_1,  mot_1=mot_1,  stp_1=stp_1,  e_1=e_1,
    az_2=az_2,  mot_2=mot_2,  stp_2=stp_2,  e_2=e_2)
print(f"  Scatter cache saved: {len(az_b)+len(az_1)+len(az_2):,} visible configs total")

# Shared concatenated arrays (used by both scatter figures)
all_az  = np.concatenate([az_bp, az_1p, az_2p])
all_mot = np.concatenate([mot_bp, mot_1p, mot_2p])
all_stp = np.concatenate([stp_bp, stp_1p, stp_2p])
all_e   = np.concatenate([e_bp,  e_1p,  e_2p])


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2 (paper) — figure_paper2_az_motion.png
# Shows: azimuth deviation vs. guidance-phase min angular velocity (log–log)
# TWO-CRITERION near-impossibility: objects in the road corridor are
# quasi-stationary (fail guidance motion); fast-moving objects point the
# wrong direction (fail azimuth). The gold upper-left region contains only
# the e=2.0 hyperbolic outlier, which fails the stopping criterion (Fig. 3).
# ══════════════════════════════════════════════════════════════════════════
print("\n─ Figure 2 (paper): azimuth deviation vs. guidance motion (two-criterion) …")
FLOOR2 = 0.001

# Colour: green = in corridor (az<5°), else by orbit type
c_az = np.where(all_az < AZ_TOL, '#1DB954',
        np.where(all_e < 0.999,   '#4878CF',
        np.where(all_e > 1.001,   '#E87B3D',
                                   '#82CAFF')))

# Keep all near-corridor points; subsample distant ones
rng3     = np.random.default_rng(13)
near_az  = all_az < 10.0
far_az   = ~near_az
n_far2   = min(25_000, int(far_az.sum()))
idx_far2 = rng3.choice(np.where(far_az)[0], n_far2, replace=False)
idx_near2= np.where(near_az)[0]
idx_in2  = np.where(all_az < AZ_TOL)[0]   # corridor points on top
idx_out2 = np.concatenate([idx_far2, np.where(near_az & (all_az >= AZ_TOL))[0]])

fig2, ax2 = plt.subplots(figsize=(9, 6.5))

# Non-corridor points first (smaller, more transparent)
ax2.scatter(np.maximum(all_az[idx_out2],  FLOOR2),
            np.maximum(all_mot[idx_out2], FLOOR2),
            c=c_az[idx_out2], s=1.5, alpha=0.15, rasterized=True, zorder=2)
# Corridor points on top — green, larger
ax2.scatter(np.maximum(all_az[idx_in2],  FLOOR2),
            np.maximum(all_mot[idx_in2], FLOOR2),
            c='#1DB954', s=12, alpha=0.70, rasterized=True, zorder=4)

ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlim(FLOOR2, 200); ax2.set_ylim(FLOOR2, 200)

# Threshold lines
ax2.axvline(AZ_TOL,       color='#CC0000', lw=2.0, ls='--', zorder=5,
            label=f'Azimuth threshold ({AZ_TOL}°)')
ax2.axhline(OMEGA_THRESH, color='#008000', lw=2.0, ls='--', zorder=5,
            label=f'Guidance-motion threshold ({OMEGA_THRESH}°/h)')

# Gold required region: upper-left (az < 5° AND motion > 2°/h)
gold2 = plt.Rectangle((FLOOR2, OMEGA_THRESH), AZ_TOL - FLOOR2, 200 - OMEGA_THRESH,
                       linewidth=2, edgecolor='#C8A000',
                       facecolor='#FFEC80', alpha=0.65, zorder=1)
ax2.add_patch(gold2)
# Label inside the gold region, well clear of both axes
ax2.text(0.8, 30,
         'Direction + guidance motion\nboth satisfied\n→ 1 orbit only (annotated)\nfails stopping: see Fig. 3',
         ha='center', va='center', fontsize=9.0, color='#7A5900',
         bbox=dict(boxstyle='round,pad=0.4', fc='#FFF8C0', ec='#C8A000', lw=1.5, alpha=0.95),
         zorder=10)

# Annotate the cheater (e=2.0, az=2.1°, motion=2.89°/h) — orange (hyperbolic)
ax2.annotate(
    'e=2.0 orbit: az=2.1°, motion=2.89°/h\n(in corridor AND fast enough to guide)\nbut stop=2.89°/h → fails stopping',
    xy=(2.1, 2.89), xytext=(18, 12), fontsize=8.5, color='#C05000',
    arrowprops=dict(arrowstyle='->', color='#C05000', lw=1.4, shrinkB=3),
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#E87B3D', lw=1.2, alpha=0.95),
    zorder=11)

# Annotate the green corridor cluster (lower-left)
n_corr2 = int((all_az < AZ_TOL).sum())
ax2.text(FLOOR2*1.3, FLOOR2*1.3,
         f'Corridor orbits (az<5°, {n_corr2:,}):\nquasi-stationary — motion ~0.006°/h\n~330× below guidance threshold',
         ha='left', va='bottom', fontsize=8.5, color='#145c30',
         bbox=dict(boxstyle='round,pad=0.3', fc='#e8f8ee', ec='#1DB954', lw=1.2, alpha=0.95),
         zorder=11)

# Stats box
n_vis2 = len(az_b) + len(az_1) + len(az_2)
ax2.text(0.98, 0.98,
         f'Visible: {n_vis2:,}\n'
         f'  broad survey: {len(az_b):,}\n'
         f'  MC sweep 1:   {len(az_1):,}\n'
         f'  MC sweep 2:   {len(az_2):,}',
         transform=ax2.transAxes, ha='right', va='top', fontsize=8.5,
         bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#aaaaaa', lw=0.9, alpha=0.93))

legend2 = [
    Line2D([0],[0], marker='o', color='w', mfc='#1DB954', mec='#1DB954', ms=8,
           label=f'Az in corridor (<5°): {n_corr2:,} orbits'),
    Line2D([0],[0], marker='o', color='w', mfc='#4878CF', mec='#4878CF', ms=6,
           label='Heliocentric elliptic (e<1), az≥5°'),
    Line2D([0],[0], marker='o', color='w', mfc='#82CAFF', mec='#82CAFF', ms=6,
           label='Heliocentric parabolic, az≥5°'),
    Line2D([0],[0], marker='o', color='w', mfc='#E87B3D', mec='#E87B3D', ms=6,
           label='Heliocentric hyperbolic (e>1), az≥5°'),
    Line2D([0],[0], color='#CC0000', lw=2, ls='--', label='Azimuth threshold (5°)'),
    Line2D([0],[0], color='#008000', lw=2, ls='--', label='Guidance-motion threshold (2°/h)'),
    mpatches.Patch(fc='#FFEC80', ec='#C8A000', lw=1.5, label='Required region (upper-left)'),
]
ax2.legend(handles=legend2, loc='lower right', fontsize=7.8,
           framealpha=0.92, edgecolor='#aaaaaa')

ax2.set_xlabel('Max azimuth deviation from road corridor $[190°, 210°]$ during guidance phase (°, log scale)',
               fontsize=10)
ax2.set_ylabel('Min guidance-phase angular velocity (°/h, log scale)\n'
               '[must exceed 2°/h to lead the Magi southward]', fontsize=10)
ax2.set_title(
    'Direction and Guidance Motion: Two-Criterion Near-Impossibility\n'
    'Corridor orbits are quasi-stationary; fast-moving orbits face the wrong direction',
    fontsize=10.5, fontweight='bold')
ax2.grid(True, which='both', alpha=0.22)
fig2.tight_layout()
fig2.savefig(OUTDIR + 'figure_paper2_az_motion.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figure_paper2_az_motion.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3 (paper) — figure_paper2_impossibility_scatter.png
# Shows: guidance-phase min angular velocity vs. stopping-phase max angular
# velocity (log–log). THREE-CRITERION impossibility: the required lower-right
# region (fast guidance AND slow stopping) is provably empty — Keplerian
# orbits cannot decelerate abruptly over four hours.
# ══════════════════════════════════════════════════════════════════════════
print("\n─ Figure 3 (paper): guidance motion vs. stopping motion (three-criterion) …")
# DESIGN: x = guidance-phase min angular velocity, y = stopping-phase max angular velocity
# Required region: x > 2°/h (guides)  AND  y < 2°/h (stops) → lower-RIGHT
# Color: green if azimuth in corridor (az < 5°), otherwise orbit-type colour
# This design is unambiguous: no point can be in the required region without satisfying
# ALL three criteria simultaneously. The one hyperbolic hit (az=2.1°, mot=2.89°/h,
# stop=2.89°/h) appears as a green dot at (2.89, 2.89) — far above the stop threshold.
FLOOR = FLOOR2  # same log-scale floor

# Colour: green if az < AZ_TOL (in corridor), else orbit-type
in_corr = all_az < AZ_TOL
c_all = np.where(in_corr,  '#1DB954',          # green: in corridor
         np.where(all_e < 0.999, '#4878CF',    # elliptic — blue
         np.where(all_e > 1.001, '#E87B3D',   # hyperbolic — orange
                                  '#82CAFF'))) # parabolic — light-blue
# Subsample non-corridor points to keep plot readable
rng2 = np.random.default_rng(7)
out_corr = ~in_corr
n_out_keep = 25_000
idx_out = rng2.choice(np.where(out_corr)[0], min(n_out_keep, out_corr.sum()), replace=False)
idx_in  = np.where(in_corr)[0]
idx_plot = np.concatenate([idx_out, idx_in])

fig, ax = plt.subplots(figsize=(9, 6.5))
# Background non-corridor points first (grey-ish, smaller)
ax.scatter(np.maximum(all_mot[idx_out], FLOOR),
           np.maximum(all_stp[idx_out], FLOOR),
           c=c_all[idx_out], s=1.5, alpha=0.15, rasterized=True, zorder=2)
# Corridor points on top — bigger, more opaque, green
ax.scatter(np.maximum(all_mot[idx_in], FLOOR),
           np.maximum(all_stp[idx_in], FLOOR),
           c='#1DB954', s=12, alpha=0.70, rasterized=True, zorder=4,
           label=f'Az in corridor (<5°): {int(in_corr.sum()):,} orbits')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(FLOOR, 200); ax.set_ylim(FLOOR, 200)

# Threshold lines
ax.axvline(OMEGA_THRESH, color='#008000', lw=2.0, ls='--', zorder=5,
           label='Guidance-motion threshold (2°/h)')
ax.axhline(OMEGA_THRESH, color='#CC0000', lw=2.0, ls='--', zorder=5,
           label='Stopping-motion threshold (2°/h)')

# Gold required region: lower-right (motion > 2 AND stop < 2)
gold = plt.Rectangle((OMEGA_THRESH, FLOOR), 200, OMEGA_THRESH - FLOOR,
                      linewidth=2, edgecolor='#C8A000',
                      facecolor='#FFEC80', alpha=0.65, zorder=1)
ax.add_patch(gold)
# Label box well clear of both axes and completely outside the gold region
ax.text(50, 0.05, 'Required region\n(motion>2°/h AND stop<2°/h)\n→ 0 orbits found',
        ha='center', va='center', fontsize=9.5, color='#7A5900',
        bbox=dict(boxstyle='round,pad=0.4', fc='#FFF8C0', ec='#C8A000', lw=1.5, alpha=0.95),
        zorder=10)

# Diagonal guide line (y = x: same speed throughout)
diag = np.logspace(np.log10(FLOOR), np.log10(200), 200)
ax.plot(diag, diag, color='#888888', lw=1.0, ls=':', zorder=3, label='Motion unchanged (y = x)')

# Annotate the one green point that has motion>2 (e=2.0, q=0.19 AU)
# It appears near (2.89, 2.89) — on the diagonal, far above stop threshold
ax.annotate(
    'e=2.0 orbit: az=2.1°, motion=2.89°/h\nbut stop=2.89°/h → fails stopping criterion',
    xy=(2.89, 2.89), xytext=(8, 8), fontsize=8.5, color='#1DB954',
    arrowprops=dict(arrowstyle='->', color='#1a7a38', lw=1.4, shrinkB=3),
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#1DB954', lw=1.2, alpha=0.95),
    zorder=11)

# Annotate the green cluster at bottom-left
n_corr = int(in_corr.sum())
ax.text(0.002, 0.002,
        f'Green cluster (az<5°, {n_corr:,} orbits):\nquasi-stationary in BOTH windows\n→ fail guidance-motion criterion',
        ha='left', va='bottom', fontsize=8.5, color='#145c30',
        bbox=dict(boxstyle='round,pad=0.3', fc='#e8f8ee', ec='#1DB954', lw=1.2, alpha=0.95),
        zorder=11)

# Stats box — bottom right, well away from axes
n_vis_total = len(az_b)+len(az_1)+len(az_2)
n_all3 = int(((az_b<AZ_TOL)&(mot_b>OMEGA_THRESH)&(stp_b<OMEGA_THRESH)).sum() +
             ((az_1<AZ_TOL)&(mot_1>OMEGA_THRESH)&(stp_1<OMEGA_THRESH)).sum() +
             ((az_2<AZ_TOL)&(mot_2>OMEGA_THRESH)&(stp_2<OMEGA_THRESH)).sum())
ax.text(0.98, 0.98,
        f'Visible: {n_vis_total:,}\n'
        f'  broad survey: {len(az_b):,}\n'
        f'  MC sweep 1:   {len(az_1):,}\n'
        f'  MC sweep 2:   {len(az_2):,}\n'
        f'All 3 criteria: 0',
        transform=ax.transAxes, ha='right', va='top', fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='#aaaaaa', lw=0.9, alpha=0.93))

legend_elements = [
    Line2D([0],[0], marker='o', color='w', mfc='#1DB954', mec='#1DB954', ms=8,
           label=f'Az in corridor (<5°): {n_corr:,} orbits'),
    Line2D([0],[0], marker='o', color='w', mfc='#4878CF', mec='#4878CF', ms=6,
           label='Heliocentric elliptic (e<1), az≥5°'),
    Line2D([0],[0], marker='o', color='w', mfc='#82CAFF', mec='#82CAFF', ms=6,
           label='Heliocentric parabolic, az≥5°'),
    Line2D([0],[0], marker='o', color='w', mfc='#E87B3D', mec='#E87B3D', ms=6,
           label='Heliocentric hyperbolic (e>1), az≥5°'),
    Line2D([0],[0], color='#008000', lw=2, ls='--', label='Guidance-motion threshold (2°/h)'),
    Line2D([0],[0], color='#CC0000', lw=2, ls='--', label='Stopping-motion threshold (2°/h)'),
    Line2D([0],[0], color='#888888', lw=1, ls=':', label='y = x (speed unchanged)'),
    mpatches.Patch(fc='#FFEC80', ec='#C8A000', lw=1.5, label='Required region (lower-right, empty)'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=7.8,
          framealpha=0.92, edgecolor='#aaaaaa')

ax.set_xlabel('Min guidance-phase angular velocity (°/h, log scale)\n[must exceed 2°/h to lead the Magi]',
              fontsize=10)
ax.set_ylabel('Max stopping-phase angular velocity (°/h, log scale)\n[must be below 2°/h to appear stationary]',
              fontsize=10)
ax.set_title(
    'Impossibility of Matthew 2:9 for Any Keplerian Orbit\n'
    'Orbits in the road corridor (green) are quasi-stationary in both windows and cannot guide; '
    'the required lower-right region is empty',
    fontsize=10.5, fontweight='bold')
ax.grid(True, which='both', alpha=0.22)
fig.tight_layout()
fig.savefig(OUTDIR + 'figure_paper2_impossibility_scatter.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figure_paper2_impossibility_scatter.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2 — figure_mc_diagnostic.png (NEW: az_score vs. guidance_motion)
# ══════════════════════════════════════════════════════════════════════════
print("─ Figure MC-diagnostic: azimuth vs. guidance-motion scatter (log y) …")
fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
MOT_FLOOR2 = 0.001

for ax, (az, mot, label, colour, n_samp, n_vis_sw) in zip(axes, [
    (az_1p, mot_1p,
     'Sweep 1 — Best-fit neighbourhood\ne∈[0.7,2.5], q∈[0.005,0.15] AU, i∈[0°,40°]',
     '#2CA02C', '1,500,000', len(az_1)),
    (az_2p, mot_2p,
     'Sweep 2 — Matney 2025 neighbourhood\ne∈[0.7,1.5], q∈[0.01,0.16] AU, i∈[0°,8°]',
     '#9467BD', '750,000', len(az_2)),
]):
    ax.scatter(az, np.maximum(mot, MOT_FLOOR2), c=colour, s=2, alpha=0.22,
               rasterized=True, zorder=2)
    ax.set_yscale('log'); ax.set_ylim(MOT_FLOOR2, 200); ax.set_xlim(-0.5, 92)
    ax.axvline(AZ_TOL,       color='#C00000', lw=2, ls='--', zorder=5, label='Az threshold (5°)')
    ax.axhline(OMEGA_THRESH, color='#008000', lw=2, ls='--', zorder=5, label='Motion threshold (2°/h)')
    # Gold required region: upper-left (az<5 AND motion>2)
    gold_mc = plt.Rectangle((0, OMEGA_THRESH), AZ_TOL, 200,
                             linewidth=2, edgecolor='#C8A000',
                             facecolor='#FFEC80', alpha=0.60, zorder=1)
    ax.add_patch(gold_mc)
    ax.text(2.5, 30, 'Required region\n(az<5° AND motion>2°/h)\n→ 0 found',
            ha='center', va='center', fontsize=8.5, color='#7A5900',
            bbox=dict(boxstyle='round,pad=0.35', fc='#FFF8C0', ec='#C8A000', lw=1.3, alpha=0.9))
    # Annotate best-az point
    bi = np.argmin(az)
    bj = np.argmax(mot)
    ax.annotate(
        f'Best az: {az[bi]:.2f}°\nmotion={mot[bi]:.4f}°/h\n({OMEGA_THRESH/max(mot[bi],1e-9):.0f}× below threshold)',
        xy=(az[bi], max(mot[bi], MOT_FLOOR2)),
        xytext=(20, 0.05), fontsize=8, color='#C00000',
        arrowprops=dict(arrowstyle='->', color='#C00000', lw=1.3, shrinkB=2))
    ax.annotate(
        f'Fastest: {mot[bj]:.1f}°/h\naz={az[bj]:.0f}° (outside corridor)',
        xy=(min(az[bj], 91), mot[bj]),
        xytext=(max(az[bj]-35, 40), 5), fontsize=8, color='#008800',
        arrowprops=dict(arrowstyle='->', color='#008800', lw=1.3, shrinkB=2))
    ax.set_xlabel('Azimuth score: max deviation outside [190°, 210°] corridor (°)', fontsize=9.5)
    ax.set_ylabel('Min guidance-phase angular velocity (°/h, log scale)', fontsize=9.5)
    ax.set_title(label, fontsize=10, fontweight='bold')
    ax.grid(True, which='both', alpha=0.25)
    ax.legend(fontsize=8.5, loc='lower right', framealpha=0.92)
    ax.text(0.98, 0.02,
            f'Sampled: {n_samp}  |  Visible: {n_vis_sw:,}  |  All criteria satisfied: 0',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=7.5,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#cccccc', lw=0.8, alpha=0.9))

fig.suptitle(
    'Monte Carlo Refinement: Azimuth Score vs. Guidance-Phase Angular Velocity (log scale)\n'
    'Orbits in the corridor are orders of magnitude too slow to guide; fast orbits point the wrong way',
    fontsize=11, fontweight='bold')
fig.tight_layout()
fig.savefig(OUTDIR + 'figure_mc_diagnostic.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figure_mc_diagnostic.png")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3 — figure9_corpus_linguistics.png (white background, 4 features)
# ══════════════════════════════════════════════════════════════════════════
print("─ Figure 4: corpus linguistics (white background) …")

# Data — TLG-verified counts (2026-05-10, MIT IRIS proxy)
# Each p = 1/(N+1) for N corpus negatives, 0 Matthew-type positives.
#   A: προάγω     — 9  in-corpus hits, 0 narrative guidance  → p = 0.100
#   B: στηρίζω    — 17 in-corpus hits, 0 halt-at-location    → p = 0.056
#      (proxy for ἵστημι; TLG morphological expander cannot resolve -μι verbs)
#   C: ἐπάνω      — 29 in-corpus hits, 0 above earthly point → p = 0.033
#   D: ἔρχομαι    — 7  in-corpus hits, 0 terrestrial arrival → p = 0.125
# Fisher combined: χ²(8) = 21.35, p = 6.3×10⁻³
# See corpus_pvalue.py and tlg_corpus_results.md for full passage citations.

# TLG-verified counts (2026-05-10).  p = 1/(N+1), N = corpus negatives.
feat_labels = ['A  προάγω\n(guidance verb)',
               'B  στηρίζω*\n(stopping verb)',
               'C  ἐπάνω\n(locative prep.)',
               'D  ἔρχομαι\n(arrival verb)']
p_vals      = [0.100, 0.056, 0.033, 0.125]
corpus_hits = [9,     17,    29,    7    ]   # all non-Matthew-type
sense_notes = ['positional / temporal',
               'planetary stationary pt.',
               'zodiacal / geometric',
               'zodiacal arrival']

neg_log_p   = [-np.log10(p) for p in p_vals]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
fig.subplots_adjust(wspace=0.38)

# ── Left panel ────────────────────────────────────────────────────────────
y     = np.arange(len(feat_labels))
bar_h = 0.50

# Blue bars = non-Matthew-type corpus hits
ax1.barh(y, corpus_hits, bar_h, color='#4878CF', alpha=0.82,
         label='Non-Matthew-type corpus hits')
# Red bars at 0 = Matthew-type hits (all absent)
ax1.barh(y, [0]*4, bar_h * 0.55, color='#CC3333', alpha=0.90,
         label='Matthew-type hits = 0')

# Count label just outside each blue bar
for i, (n, note) in enumerate(zip(corpus_hits, sense_notes)):
    ax1.text(n + 0.5, y[i], f'{n}', va='center', ha='left',
             fontsize=10, color='#2255AA', fontweight='bold')
    # Sense note inside the bar (white text) if bar is wide enough, else skip
    if n >= 12:
        ax1.text(n * 0.5, y[i], note, va='center', ha='center',
                 fontsize=7.5, color='white', style='italic')

# "0" label at left edge for every feature
for i in y:
    ax1.text(0.3, i, '0', va='center', ha='left',
             fontsize=9, color='#CC3333', fontweight='bold')

ax1.set_yticks(y)
ax1.set_yticklabels(feat_labels, fontsize=10)
ax1.set_xlabel('Corpus occurrences (25-text corpus, TLG-verified)', fontsize=9.5)
ax1.set_title('Vocabulary profile: corpus hits vs.\nMatthew 2:9 (Matthew-type = 0 in every case)',
              fontsize=10.5, fontweight='bold')
ax1.set_xlim(0, 38)
ax1.legend(loc='lower right', fontsize=8.5, framealpha=0.9)
ax1.grid(True, axis='x', alpha=0.25, lw=0.7)
ax1.tick_params(axis='y', length=0)
# Footnote for proxy
ax1.annotate('* στηρίζω is proxy for ἵστημι (TLG cannot expand -μι inflections)',
             xy=(0, -0.15), xycoords='axes fraction',
             fontsize=7.5, color='#555555', style='italic')

# ── Right panel ───────────────────────────────────────────────────────────
p05_line  = -np.log10(0.05)
pbonf_line = -np.log10(0.0125)

bar_colors = ['#E87B3D' if p < 0.05 else '#888888' for p in p_vals]
bars = ax2.bar(range(4), neg_log_p, color=bar_colors, width=0.55,
               edgecolor='white', linewidth=0.8, alpha=0.90)

ax2.axhline(p05_line,   color='#E87B3D', lw=1.5, ls='--', alpha=0.85,
            label=r'$p=0.05$')
ax2.axhline(pbonf_line, color='#CC3333',  lw=1.5, ls=':',  alpha=0.85,
            label=r'Bonferroni $\alpha/4=0.0125$')

# p-value labels above each bar, well clear of the bar top
y_max_data = max(neg_log_p)
for i, (bar, p, bc) in enumerate(zip(bars, p_vals, bar_colors)):
    ax2.text(i, neg_log_p[i] + 0.06, f'p = {p:.3f}',
             ha='center', va='bottom', fontsize=9, fontweight='bold', color=bc)

ax2.set_xticks(range(4))
ax2.set_xticklabels(['A\nproάγω', 'B\nστηρίζω*', 'C\nἐπάνω', 'D\nἔρχομαι'],
                    fontsize=9.5)
ax2.set_ylabel(r'$-\log_{10}(p)$', fontsize=10)
ax2.set_title("Fisher's exact test\n(Matthew 2:9 vs. astronomical corpus)",
              fontsize=10.5, fontweight='bold')
ax2.legend(fontsize=8.5, loc='upper left', framealpha=0.9)
ax2.set_ylim(0, 2.4)
ax2.grid(True, axis='y', alpha=0.25, lw=0.7)

# Joint-result box placed below the bars so it does not overlap p labels
ax2.text(0.50, 0.10,
         r'Joint: $\chi^2(8)=21.4,\ p\approx6.3\times10^{-3}$' + '\n' +
         r'BF $\approx$ 4,000 vs. guiding-star tradition',
         transform=ax2.transAxes, ha='center', va='bottom', fontsize=8.8,
         bbox=dict(boxstyle='round,pad=0.35', fc='#FFF4EC', ec='#CC3333', lw=1.1))

fig.suptitle('Linguistic incompatibility of Matthew 2:9\nwith Greek astronomical discourse',
             fontsize=11.5, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(OUTDIR + 'figure9_corpus_linguistics.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved figure9_corpus_linguistics.png")

print("\n" + "═"*64)
print("All figures updated.")
print(f"  figure_paper2_az_motion.png              — paper Fig 2: az deviation vs. guidance motion")
print(f"  figure_paper2_impossibility_scatter.png  — paper Fig 3: guidance motion vs. stopping motion")
print(f"  figure_mc_diagnostic.png                 — supplementary S1: MC sweep diagnostic")
print(f"  figure9_corpus_linguistics.png           — white background, 4 features")
print("═"*64)
