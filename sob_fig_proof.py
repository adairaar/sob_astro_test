#!/usr/bin/env python3
"""
sob_fig_proof.py — Regenerate figure6_impossibility_proof.png
=============================================================
Loads the full scatter dataset cached by sob_fig_update.py
(scatter_cache.npz, 717k visible configurations) and produces
a fully updated four-panel proof figure with dense scatter data
in panel A.

Run AFTER sob_fig_update.py so the cache exists.

Panels:
  A: Guidance azimuth score vs. stopping motion — full 717k dataset
  B: Full-day azimuth tracks (best-fit orbit and Matney 2025 comet)
  C: Guidance score distributions by eccentricity — full dataset
  D: Required ΔV vs. available gravitational forces (analytical)

Aaron Adair / 2026
"""

import warnings; warnings.filterwarnings('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body_barycentric
import astropy.units as u
import time as walltime

import os as _os
OUTDIR = _os.path.dirname(_os.path.abspath(__file__)) + _os.sep

# ── Constants (must match sob_fig_update.py exactly) ──────────────────────
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

# Matney (2025) comet orbital elements
M_q, M_e, M_i   = 0.0407, 1.0, 1.39      # AU, -, deg
M_Om, M_om, M_T = 93.01, 346.25, 1_719_785.565  # deg, deg, JD

# Physical constants for panel D
MU_KM3_S2 = 1.32712440018e11
GM_EARTH   = 3.986004418e5
AU_KM      = 1.495978707e8

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'axes.edgecolor': '#333333', 'axes.labelcolor': '#111111',
    'xtick.color': '#333333', 'ytick.color': '#333333',
    'text.color': '#111111', 'grid.color': '#cccccc',
    'grid.linestyle': '--', 'grid.alpha': 0.5,
    'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10,
})

# ══════════════════════════════════════════════════════════════════════════
# 1. Load cached scatter data
# ══════════════════════════════════════════════════════════════════════════
print("─" * 64)
print("Loading scatter cache …")
t0 = walltime.time()
cache = np.load(OUTDIR + 'scatter_cache.npz')
az_b,  mot_b,  stp_b,  e_b  = cache['az_b'],  cache['mot_b'],  cache['stp_b'],  cache['e_b']
az_1,  mot_1,  stp_1,  e_1  = cache['az_1'],  cache['mot_1'],  cache['stp_1'],  cache['e_1']
az_2,  mot_2,  stp_2,  e_2  = cache['az_2'],  cache['mot_2'],  cache['stp_2'],  cache['e_2']
az_all  = np.concatenate([az_b,  az_1,  az_2])
mot_all = np.concatenate([mot_b, mot_1, mot_2])
stp_all = np.concatenate([stp_b, stp_1, stp_2])
e_all   = np.concatenate([e_b,   e_1,   e_2])
n_total = len(az_all)
print(f"  Loaded {n_total:,} visible configurations in {walltime.time()-t0:.1f}s")

# ══════════════════════════════════════════════════════════════════════════
# 2. ENU matrices — guidance+stopping window (for orbit-finding scan)
#    and full-day window (for panel B azimuth tracks)
# ══════════════════════════════════════════════════════════════════════════
print("Precomputing ENU matrices …")
t1 = walltime.time()

_t_g = np.array([8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0])
_t_s = np.array([10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0])
_times = np.concatenate([_t_g, _t_s])
N_G = len(_t_g); N_S = len(_t_s); N_GS = len(_times); DT_H = 0.25

jds_gs  = JD_NOON + _times / 24.0
jds_utc = jds_gs - DT_SEC / 86400.0
t_tdb   = Time(jds_gs,  format='jd', scale='tdb')
t_utc   = Time(jds_utc, format='jd', scale='utc')
eb      = get_body_barycentric('earth', t_tdb)
sb      = get_body_barycentric('sun',   t_tdb)
earth_pos = (eb.xyz.to(u.AU).value - sb.xyz.to(u.AU).value).T  # (17, 3)

def make_enu(t_utc_step):
    frame = AltAz(obstime=t_utc_step, location=BETHLEHEM)
    def aa2icrs(az_d, alt_d):
        c = SkyCoord(alt=alt_d*u.deg, az=az_d*u.deg, frame=frame).icrs
        d = np.radians(c.dec.deg); r = np.radians(c.ra.deg)
        return np.array([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)])
    E_raw = aa2icrs(90, 0.001)
    Z = aa2icrs(0, 90.0); Z /= np.linalg.norm(Z)
    E = E_raw - np.dot(E_raw, Z)*Z; E /= np.linalg.norm(E)
    N = np.cross(Z, E)
    return np.array([E, N, Z])

enu = np.array([make_enu(t_utc[i]) for i in range(N_GS)])

# Full-day grid: 04:00–14:00, 15-min steps (41 points)
_t_day  = np.arange(4.0, 14.25, 0.25)
jds_day = JD_NOON + _t_day / 24.0
jds_day_utc = jds_day - DT_SEC / 86400.0
t_day_tdb = Time(jds_day, format='jd', scale='tdb')
t_day_utc = Time(jds_day_utc, format='jd', scale='utc')
eb_d = get_body_barycentric('earth', t_day_tdb)
sb_d = get_body_barycentric('sun',   t_day_tdb)
epos_day = (eb_d.xyz.to(u.AU).value - sb_d.xyz.to(u.AU).value).T
enu_day  = np.array([make_enu(t_day_utc[i]) for i in range(len(_t_day))])
print(f"  ENU ready in {walltime.time()-t1:.1f}s")

# ══════════════════════════════════════════════════════════════════════════
# 3. Kepler solvers (identical to sob_fig_update.py)
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
    cO=np.cos(Om); sO=np.sin(Om); co=np.cos(om); so=np.sin(om)
    ci=np.cos(inc); si=np.sin(inc)
    Px=cO*co-sO*so*ci; Py=sO*co+cO*so*ci; Pz=so*si
    Qx=-cO*so-sO*co*ci; Qy=-sO*so+cO*co*ci; Qz=co*si
    return np.stack([r*(cnu*Px+snu*Qx), r*(cnu*Py+snu*Qy), r*(cnu*Pz+snu*Qz)], axis=-1)

def evaluate_batch(q_, e_, i_r, Om_r, om_r, Tp):
    B = len(q_)
    dt_mat = jds_gs[np.newaxis,:] - Tp[:,np.newaxis]
    N = B * N_GS
    t_i = np.tile(np.arange(N_GS), B)
    q_f = np.repeat(q_, N_GS); e_f = np.repeat(e_, N_GS)
    i_f = np.repeat(i_r, N_GS); Om_f = np.repeat(Om_r, N_GS)
    om_f = np.repeat(om_r, N_GS); dt_f = dt_mat.ravel()
    nu_f = _true_anomaly(dt_f, q_f, e_f)
    p_f = q_f*(1+e_f); r_f = p_f/(1+e_f*np.cos(nu_f))
    ok_f = np.isfinite(nu_f) & (r_f>0) & (r_f<1e4)
    nu_s = np.where(ok_f, nu_f, 0.); r_s = np.where(ok_f, r_f, 0.)
    xyz_f = _pos_icrs(nu_s, r_s, Om_f, om_f, i_f); xyz_f[~ok_f] = 0.
    geo_f = xyz_f - earth_pos[t_i]
    nrm = np.linalg.norm(geo_f, axis=1, keepdims=True)
    geo_n = geo_f / np.where(nrm>0, nrm, 1.)
    E_f = enu[t_i]
    e_c = np.einsum('ni,ni->n', geo_n, E_f[:,0,:])
    n_c = np.einsum('ni,ni->n', geo_n, E_f[:,1,:])
    u_c = np.einsum('ni,ni->n', geo_n, E_f[:,2,:])
    alt_m = np.arcsin(np.clip(u_c,-1,1))*RAD
    az_m  = (np.arctan2(e_c, n_c)*RAD) % 360.
    ok_m  = ok_f.reshape(B, N_GS)
    alt_m = alt_m.reshape(B, N_GS)
    az_m  = az_m.reshape(B, N_GS)
    geo_nm = geo_n.reshape(B, N_GS, 3)
    vis_ok = (np.all(alt_m[:,:N_G]>ALT_FLOOR, axis=1) & np.all(ok_m[:,:N_G], axis=1))
    az_dev = np.maximum(0., np.maximum(AZ_LO-az_m[:,:N_G], az_m[:,:N_G]-AZ_HI))
    az_score = np.max(az_dev, axis=1)
    motion_min = np.min(np.arccos(np.clip(
        np.einsum('bti,bti->bt', geo_nm[:,:N_G-1,:], geo_nm[:,1:N_G,:]), -1,1))*RAD/DT_H, axis=1)
    stop_max = np.max(np.arccos(np.clip(
        np.einsum('bti,bti->bt', geo_nm[:,N_G:-1,:], geo_nm[:,N_G+1:,:]), -1,1))*RAD/DT_H, axis=1)
    return vis_ok, az_score, motion_min, stop_max

# ══════════════════════════════════════════════════════════════════════════
# 4. Find best-fit guidance orbit
#    Broad search: multiple (e, q, i) combinations × all T_offs from fine
#    grid survey × Ω/ω at 5° steps.  Goal: reproduce the orbit that the
#    full 5°-grid survey found with az_score=0.015° (e=1.1, q=0.03, i=10°).
# ══════════════════════════════════════════════════════════════════════════
print("Finding best-fit guidance orbit …")
t2 = walltime.time()
Om_s = np.arange(0, 360, 5); om_s = np.arange(0, 360, 5)
Om_g = np.repeat(Om_s, len(om_s)).astype(float)
om_g = np.tile(om_s, len(Om_s)).astype(float)
B_s  = len(Om_g)

T_OFFS_SCAN = [-60, -40, -20, -10, -5, 5, 10, 20, 60]   # same as fine-grid survey
E_SCAN = [0.9, 0.95, 0.99, 1.0, 1.1, 1.5]
Q_SCAN = [0.02, 0.03, 0.05, 0.07, 0.10]
I_SCAN = [1.4, 5.0, 10.0, 20.0]

best_az_score = 999.; best_Om = 0.; best_om = 0.
best_e_found = 1.1; best_q_found = 0.03; best_i_found = 10.0; best_Tp_found = JD_NOON - 5.

for e_try in E_SCAN:
    for q_try in Q_SCAN:
        for i_try in I_SCAN:
            for T_off_try in T_OFFS_SCAN:
                Tp_try = JD_NOON + T_off_try
                vis_s, az_s, mot_s, stp_s = evaluate_batch(
                    np.full(B_s, q_try), np.full(B_s, e_try),
                    np.full(B_s, np.radians(i_try)),
                    np.radians(Om_g), np.radians(om_g),
                    np.full(B_s, Tp_try))
                if vis_s.any():
                    scores = np.where(vis_s, az_s, 999.)
                    sc = float(scores.min())
                    if sc < best_az_score:
                        best_az_score = sc
                        bi = int(scores.argmin())
                        best_Om, best_om = float(Om_g[bi]), float(om_g[bi])
                        best_e_found, best_q_found = e_try, q_try
                        best_i_found, best_Tp_found = i_try, Tp_try

print(f"  Done in {walltime.time()-t2:.1f}s")
print(f"  Best: az_score={best_az_score:.3f}°, e={best_e_found}, q={best_q_found} AU, "
      f"i={best_i_found}°, Ω={best_Om:.0f}°, ω={best_om:.0f}°, "
      f"T_off={best_Tp_found - JD_NOON:+.0f}d")

# ── Full-day azimuth tracks ────────────────────────────────────────────────
def _track(q_, e_, i_d, Om_d, om_d, Tp):
    N   = len(_t_day)
    dt  = jds_day - Tp
    nu  = _true_anomaly(dt, np.full(N, q_), np.full(N, e_))
    r   = q_ * (1 + e_) / (1 + e_ * np.cos(nu))
    ok  = np.isfinite(nu) & (r > 0) & (r < 1e4)
    nu_s = np.where(ok, nu, 0.); r_s = np.where(ok, r, 0.)
    xyz  = _pos_icrs(nu_s, r_s,
                     np.full(N, np.radians(Om_d)), np.full(N, np.radians(om_d)),
                     np.full(N, np.radians(i_d)))
    xyz[~ok] = 0.
    geo  = xyz - epos_day
    nrm  = np.linalg.norm(geo, axis=1, keepdims=True)
    geo_n = geo / np.where(nrm > 0, nrm, 1.)
    e_c  = np.einsum('ni,ni->n', geo_n, enu_day[:, 0, :])
    n_c  = np.einsum('ni,ni->n', geo_n, enu_day[:, 1, :])
    u_c  = np.einsum('ni,ni->n', geo_n, enu_day[:, 2, :])
    alt  = np.degrees(np.arcsin(np.clip(u_c, -1, 1)))
    az   = np.degrees(np.arctan2(e_c, n_c)) % 360.
    return az, alt

az_best, alt_best = _track(best_q_found, best_e_found, best_i_found, best_Om, best_om, best_Tp_found)
az_mat,  alt_mat  = _track(M_q,  M_e,  M_i,  M_Om,  M_om,  M_T)

# ══════════════════════════════════════════════════════════════════════════
# 5. Panel C eccentricity bins
# ══════════════════════════════════════════════════════════════════════════
e_bins   = [(0.00, 0.50, 'e≤0.5'),
            (0.50, 0.90, 'e=0.5–0.9'),
            (0.90, 1.00, 'e=0.9–1'),
            (1.00, 1.01, 'e≈1 (parabolic)'),
            (1.01, 2.00, 'e=1–2'),
            (2.00, 6.00, 'e>2')]
e_colors = plt.cm.viridis(np.linspace(0, 1, len(e_bins)))

# ══════════════════════════════════════════════════════════════════════════
# 6. Build figure
# ══════════════════════════════════════════════════════════════════════════
print("Drawing figure6_impossibility_proof.png …")
fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle(
    "Proof of Impossibility: No Keplerian Orbit Reproduces Matthew 2:9\n"
    f"(Guidance: Az in [190°,210°] corridor + motion >{OMEGA_THRESH}°/h; "
    f"Stopping: motion <{OMEGA_THRESH}°/h)",
    fontsize=12, fontweight='bold')

# ── Panel A: guidance azimuth score vs. stopping motion ──────────────────
ax = axes[0, 0]
# Subsample to ≤150k for scatter speed; keep all guiding (motion>2) orbits
rng_a = np.random.default_rng(11)
motion_ok = mot_all > OMEGA_THRESH
n_guide = int(motion_ok.sum()); n_rest = min(150_000 - n_guide, (~motion_ok).sum())
idx_rest = rng_a.choice(np.where(~motion_ok)[0], n_rest, replace=False)
idx_A = np.concatenate([np.where(motion_ok)[0], idx_rest])

ax.scatter(az_all[idx_A[~motion_ok[idx_A]]],
           stp_all[idx_A[~motion_ok[idx_A]]],
           s=1.2, alpha=0.10, c='#888888', rasterized=True,
           label=f'Motion ≤{OMEGA_THRESH}°/h (not guiding)')
ax.scatter(az_all[motion_ok], stp_all[motion_ok],
           s=5, alpha=0.50, c='#4878CF', rasterized=True,
           label=f'Motion >{OMEGA_THRESH}°/h (guiding orbits)')
ax.axvline(AZ_TOL,       color='red',   lw=2, ls='--',
           label=f'Azimuth tolerance ({AZ_TOL}°)')
ax.axhline(OMEGA_THRESH, color='green', lw=2, ls='--',
           label=f'Stopping criterion ({OMEGA_THRESH}°/h)')
ax.fill_between([0, AZ_TOL], [0, 0], [OMEGA_THRESH, OMEGA_THRESH],
                alpha=0.25, color='gold', label='Must satisfy ALL THREE → 0 orbits')
ax.set_xlabel('Guidance azimuth score: max degrees outside [190°, 210°] corridor', fontsize=10)
ax.set_ylabel('Stopping score: max apparent velocity (°/h)', fontsize=10)
ax.set_title('(A) Guidance vs. Stopping Scores\n(blue = orbits with motion >{:.0f}°/h to guide)'.format(OMEGA_THRESH),
             fontsize=10, fontweight='bold')
ax.set_xlim(0, 90); ax.set_ylim(0, 30)
ax.legend(fontsize=8, markerscale=4)
n_all_A = int(((az_all < AZ_TOL) & (mot_all > OMEGA_THRESH) & (stp_all < OMEGA_THRESH)).sum())
ax.text(0.02, 0.75, f"N visible: {n_total:,}\nSatisfy ALL: {n_all_A}",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
ax.grid(True, alpha=0.2)

# ── Panel B: full-day azimuth tracks ─────────────────────────────────────
ax = axes[0, 1]
vis_best = alt_best > 5; vis_mat = alt_mat > 5
ax.plot(_t_day[vis_best], az_best[vis_best], 'b-', lw=2,
        label=f'Best-fit guidance orbit (e={best_e_found}, q={best_q_found} AU, '
              f'i={best_i_found}°, az_score={best_az_score:.2f}°)')
ax.plot(_t_day[vis_mat],  az_mat[vis_mat],   'r--', lw=2,
        label='Matney (2025) 5 BCE comet')
ax.axhspan(AZ_LO, AZ_HI, alpha=0.12, color='green',
           label=f'Road corridor [{AZ_LO:.0f}°,{AZ_HI:.0f}°]')
ax.axvspan(8, 10,  alpha=0.08, color='green')
ax.axvspan(10, 12, alpha=0.08, color='gold', label='Stopping window')
ax.set_xlabel('Local solar time (h)', fontsize=10)
ax.set_ylabel('Azimuth (°)', fontsize=10)
ax.set_title('(B) Full-Day Az Track: Best Orbit vs. Matney\n'
             '(Neither stays in road corridor during stopping window)',
             fontsize=10, fontweight='bold')
ax.set_xlim(4, 14); ax.set_ylim(150, 270)
ax.set_xticks(range(4, 15, 2))
ax.legend(fontsize=8, loc='upper left'); ax.grid(True, alpha=0.2)

# ── Panel C: per-eccentricity distributions ───────────────────────────────
ax = axes[1, 0]
for (e_lo, e_hi, e_lab), col in zip(e_bins, e_colors):
    mask = (e_all >= e_lo) & (e_all < e_hi)
    data = az_all[mask]
    if len(data) > 0:
        d = data[data < 90]
        ax.hist(d, bins=np.linspace(0, 90, 31), color=col,
                alpha=0.5, label=f'{e_lab} (N={len(data):,})', density=True)
ax.axvline(AZ_TOL, color='red', lw=2.5, ls='--',
           label=f'{AZ_TOL}° outside-corridor tolerance')
ax.set_xlabel('Guidance score: max degrees outside road corridor [190°,210°]', fontsize=10)
ax.set_ylabel('Density', fontsize=10)
ax.set_title('(C) Guidance Score Distribution by Eccentricity\n'
             '(No eccentricity keeps orbit inside road corridor)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=7, loc='upper right', ncol=2); ax.grid(True, alpha=0.2)

# ── Panel D: required ΔV vs. available forces (analytical) ───────────────
ax = axes[1, 1]
r_AU = np.logspace(-3, 0, 200); r_km = r_AU * AU_KM
# Kinematic requirement from table (alt=45°): ω_guide=15.55°/h, ω_stop=14.76°/h
omg_guide_d, omg_stop_d = 15.55, 14.76
v_g   = omg_guide_d * np.pi/180 * r_km / 3600
v_s   = omg_stop_d  * np.pi/180 * r_km / 3600
dv    = v_g - v_s
dt_s  = 30 * 60
a_sol = 2 * MU_KM3_S2 * r_km / AU_KM**3
a_ear = GM_EARTH / r_km**2
dv_sol = a_sol * dt_s; dv_ear = a_ear * dt_s
ax.loglog(r_AU, dv,     'b-',     lw=2.5, label='Required Δv (guide→stop)')
ax.loglog(r_AU, dv_sol, color='orange', lw=2, ls='--', label='Solar tidal Δv in 30 min')
ax.loglog(r_AU, dv_ear, 'green',  lw=2, ls='-.', label="Earth gravity Δv in 30 min")
ax.axvline(0.0026, color='red', lw=1.5, ls=':',
           label="Matney's distance (0.0026 AU)")
mask_ok = dv_ear > dv
if mask_ok.any():
    ax.fill_between(r_AU, dv, dv_ear, where=mask_ok,
                    alpha=0.2, color='green', label='Earth gravity sufficient')
ax.set_xlabel('Geocentric distance r (AU)', fontsize=10)
ax.set_ylabel('Velocity change (km/s)', fontsize=10)
ax.set_title('(D) Required Δv for Transition vs. Available Forces\n'
             "(Earth gravity insufficient at Matney's proposed distance)",
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.2)

plt.tight_layout()
outpath = OUTDIR + 'figure6_impossibility_proof.png'
fig.savefig(outpath, dpi=600, bbox_inches='tight')
plt.close()
print(f"  Saved figure6_impossibility_proof.png")
print(f"  Panel A: {n_total:,} visible configs (was 50,347)")
print("═" * 64)
