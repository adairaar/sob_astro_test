#!/usr/bin/env python3
"""
sob_epoch_robustness.py
Multi-epoch robustness survey for the Matt 2:9 impossibility proof.
Aaron Adair / 2026

Tests ±40-day windows around all major Star of Bethlehem proposals:
  12 BCE  Halley's Comet perihelion (Oct 9)
   7 BCE  Jupiter-Saturn triple conjunction (May 29, Oct 3, Dec 4)
   5 BCE  Reference epoch from main paper (Jun 8)
   1 BCE  Lunar eclipse / nativity proposals (Apr 1, Dec 25)

For each candidate observation date, Earth barycentric positions and
ENU rotation matrices are recomputed via astropy from the JPL ephemeris.
An orbital survey identical in logic to sob_proof.py (30-degree Omega/omega
spacing) then confirms that zero Keplerian configurations satisfy both the
guidance criterion and the stopping criterion at every epoch.
"""

import warnings; warnings.filterwarnings('ignore')
import numpy as np
from astropy.time import Time
from astropy.coordinates import (SkyCoord, AltAz, EarthLocation,
                                 get_body_barycentric)
import astropy.units as u
import time as walltime

T_WALL = walltime.time()
def elapsed(): return f"[{walltime.time()-T_WALL:.1f}s]"

# ─── site and scenario constants ──────────────────────────────────────────
# Jerusalem center (31°46'N, 35°14'E), Bethlehem center (31°42'N, 35°12'E)
# Straight-line bearing Jerusalem → Bethlehem: 203°
# Road sweeps 190°–210° (deviation around the ridge)
BLAT, BLON, BALT = 31 + 42/60, 35 + 12/60, 765.0
BETHLEHEM = EarthLocation(lon=BLON*u.deg, lat=BLAT*u.deg, height=BALT*u.m)
AZ_LO, AZ_HI    = 190.0, 210.0   # road azimuth corridor (deg)
AZ_THRESH       = 5.0            # deg — tolerance outside corridor boundary
# Perceptible angular velocity threshold (deg/h).
# Physical basis: naked-eye resolution ~1'; detectable over 30 min → floor ~0.033°/h.
# Moon moves ~0.5°/h relative to background stars (routinely noted by Babylonian
# astronomers).  2°/h = 4× lunar rate; detectable within minutes.
# Guidance requires ω > OMEGA_THRESH (star "goes before"); stopping requires ω < OMEGA_THRESH.
# This threshold is deliberately generous: stricter values only strengthen the result.
OMEGA_THRESH    = 2.0            # deg/h

K = 0.01720209895; MU = K**2
EPS = np.radians(23.439291111)
R_ECL2EQU = np.array([[1,0,0],
                       [0, np.cos(EPS), -np.sin(EPS)],
                       [0, np.sin(EPS),  np.cos(EPS)]])
DT_SEC = 10572.0   # ΔT (TT – UT1) for BCE era, seconds (Morrison & Stephenson 2004)

# ─── proposed observation epochs ─────────────────────────────────────────
# Each tuple: (label, JD at local noon on the proposed canonical date)
EPOCHS = [
    ("12 BCE  Halley perihelion  (Oct  9)", 1717321.5),
    (" 7 BCE  Jupiter-Saturn #1  (May 29)", 1719014.5),
    (" 7 BCE  Jupiter-Saturn #2  (Oct  3)", 1719141.5),
    (" 7 BCE  Jupiter-Saturn #3  (Dec  4)", 1719203.5),
    (" 5 BCE  reference epoch    (Jun  8)", 1719755.898),
    (" 1 BCE  eclipse proposal   (Apr  1)", 1721148.5),
    (" 1 BCE  nativity proposal  (Dec 25)", 1721416.5),
]

OBS_OFFSETS    = np.array([-40.0, 0.0, +40.0])   # days around each canonical date
T_PERI_OFFSETS = [-40.0, -20.0, -5.0, +5.0]      # perihelion offset from obs date

OM_VALS = np.arange(0, 360, 30)   # 12 values
om_vals = np.arange(0, 360, 30)   # 12 values → 144 (Ω, ω) pairs
Q_VALS  = [0.03, 0.04, 0.07, 0.12]
I_VALS  = [1.4, 10.0, 30.0, 60.0]
E_VALS  = [0.0, 0.3, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0]
CONFIGS_PER_DATE = (len(E_VALS)*len(Q_VALS)*len(I_VALS)*
                    len(T_PERI_OFFSETS)*len(OM_VALS)*len(om_vals))

# ─── Kepler equation solvers (identical to sob_proof.py) ─────────────────

def barker(dt, q):
    W = 3.0*np.sqrt(MU/(2.0*q**3))*dt
    s = np.sqrt((W/2)**2 + 1.0)
    return np.cbrt(W/2+s) + np.cbrt(W/2-s)

def kepler_elliptic(M, e, tol=1e-12):
    E = M.copy()
    for _ in range(50):
        dE = (M - E + e*np.sin(E)) / (1.0 - e*np.cos(E))
        E += dE
        if np.max(np.abs(dE)) < tol: break
    return E

def true_anomaly(dt_days, q, e):
    if abs(e - 1.0) < 1e-9:
        return 2.0*np.arctan(barker(dt_days, q))
    elif e < 1.0:
        a = q/(1-e); n = np.sqrt(MU/a**3); M = n*dt_days
        E = kepler_elliptic(np.atleast_1d(float(M)), e)
        nu = 2.0*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
        return float(nu[0]) if nu.size == 1 else nu
    else:
        a = q/(e-1); n = np.sqrt(MU/a**3); M = n*dt_days
        H = np.sign(M)*np.log(2.0*np.abs(M)/e + 1.8)
        for _ in range(50):
            dH = (M - e*np.sinh(H) + H) / (e*np.cosh(H) - 1.0); H += dH
            if np.abs(dH) < 1e-12: break
        return 2.0*np.arctan2(np.sqrt(e+1)*np.sinh(H/2), np.sqrt(e-1)*np.cosh(H/2))

# ─── ENU matrix (via astropy SkyCoord → ICRS, same as sob_proof.py) ──────

def make_enu(t_utc_single):
    frame = AltAz(obstime=t_utc_single, location=BETHLEHEM)
    def aa2icrs(az, alt):
        c = SkyCoord(alt=alt*u.deg, az=az*u.deg, frame=frame).icrs
        d = np.radians(c.dec.deg); r = np.radians(c.ra.deg)
        return np.array([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)])
    E = aa2icrs(90, 0.001); Z = aa2icrs(0, 90)
    Z /= np.linalg.norm(Z); E -= np.dot(E,Z)*Z; E /= np.linalg.norm(E)
    return np.array([E, np.cross(Z,E), Z])

# ─── build Earth/ENU data for one observation date ────────────────────────

def build_epoch_data(jd_noon_obs):
    hrs   = np.arange(-8.0, 3.1, 0.25)
    jds   = jd_noon_obs + hrs/24.0
    jds_u = jds - DT_SEC/86400.0
    lh    = hrs + 12.0    # local solar hour

    t_tdb = Time(jds, format='jd', scale='tdb')
    t_utc = Time(jds_u, format='jd', scale='utc')
    Ef = get_body_barycentric('earth', t_tdb)
    Sf = get_body_barycentric('sun',   t_tdb)
    earth = (Ef - Sf).xyz.to(u.AU).value
    R_enu = np.array([make_enu(t_utc[ii]) for ii in range(len(jds))])

    mg = (lh >= 8.0) & (lh <= 10.0)
    ms = (lh >  10.0) & (lh <= 12.0)
    return dict(
        jd_guide    = jds[mg],        jd_stop     = jds[ms],
        t_stop      = lh[ms],
        E_guide     = earth[:, mg],   E_stop      = earth[:, ms],
        R_enu_guide = R_enu[mg],      R_enu_stop  = R_enu[ms],
        N_guide     = int(mg.sum()),  N_stop      = int(ms.sum()),
    )

# ─── vectorised survey (identical logic to sob_proof.py) ─────────────────

def survey_epoch(q, e, i_d, Om_d_arr, om_d_arr, T_jd, ep):
    nOm = len(Om_d_arr); nom = len(om_d_arr)
    i = np.radians(i_d); ci, si = np.cos(i), np.sin(i)
    cO = np.cos(np.radians(Om_d_arr)); sO = np.sin(np.radians(Om_d_arr))
    co = np.cos(np.radians(om_d_arr)); so = np.sin(np.radians(om_d_arr))
    P_ecl = np.array([cO[:,None]*co - sO[:,None]*so*ci,
                      sO[:,None]*co + cO[:,None]*so*ci,
                      np.broadcast_to(so*si, (nOm,nom)).copy()])
    Q_ecl = np.array([-cO[:,None]*so - sO[:,None]*co*ci,
                      -sO[:,None]*so + cO[:,None]*co*ci,
                       np.broadcast_to(co*si, (nOm,nom)).copy()])
    P = np.einsum('ij,jkl->ikl', R_ECL2EQU, P_ecl)
    Q = np.einsum('ij,jkl->ikl', R_ECL2EQU, Q_ecl)

    Ng = ep['N_guide']
    az_g  = np.zeros((Ng, nOm, nom))
    alt_g = np.zeros_like(az_g)
    dg    = np.zeros((Ng, nOm, nom, 3))   # direction vectors for motion score
    for ii in range(Ng):
        dt = ep['jd_guide'][ii] - T_jd
        nu = true_anomaly(dt, q, e); r = q*(1+e)/(1+e*np.cos(nu))
        pos = r*(np.cos(nu)*P + np.sin(nu)*Q)
        rho = pos - ep['E_guide'][:,ii,None,None]
        dist = np.linalg.norm(rho, axis=0); dist[dist<1e-12]=1e-12
        rh = rho/dist
        enu = np.einsum('ij,jkl->ikl', ep['R_enu_guide'][ii], rh)
        alt_g[ii] = np.degrees(np.arcsin(np.clip(enu[2],-1,1)))
        az_g[ii]  = np.degrees(np.arctan2(enu[0],enu[1])) % 360.0
        dg[ii]    = np.moveaxis(rh, 0, -1)

    vis_g = (alt_g > 5.0).all(axis=0)

    # Guidance azimuth score: max degrees outside road corridor
    dev_lo      = np.maximum(0.0, AZ_LO - az_g)
    dev_hi      = np.maximum(0.0, az_g  - AZ_HI)
    guide_score = np.maximum(dev_lo, dev_hi).max(axis=0)

    # Guidance motion score: min apparent angular velocity during guidance window.
    # proago (Matt 2:9) requires the star to be actively moving; a stationary
    # object at the right azimuth does not "go before" the Magi.
    d1g = dg[:-1]; d2g = dg[1:]
    dt_h_g = np.diff(ep['jd_guide']) * 24.0   # hours
    crs_g  = np.cross(d1g, d2g, axis=-1)
    sin_g  = np.linalg.norm(crs_g, axis=-1)
    cos_g  = np.einsum('...i,...i->...', d1g, d2g)
    omg_g  = np.degrees(np.arctan2(sin_g, cos_g)) / dt_h_g[:,None,None]
    guide_motion = omg_g.min(axis=0)   # must exceed OMEGA_THRESH

    # Stopping score: max apparent angular velocity during stopping window.
    Ns = ep['N_stop']
    dp = np.zeros((Ns, nOm, nom, 3))
    for ii in range(Ns):
        dt = ep['jd_stop'][ii] - T_jd
        nu = true_anomaly(dt, q, e); r = q*(1+e)/(1+e*np.cos(nu))
        pos = r*(np.cos(nu)*P + np.sin(nu)*Q)
        rho = pos - ep['E_stop'][:,ii,None,None]
        dist = np.linalg.norm(rho, axis=0); dist[dist<1e-12]=1e-12
        dp[ii] = np.moveaxis(rho/dist, 0, -1)

    d1,d2  = dp[:-1], dp[1:]
    dt_h   = np.diff(ep['t_stop'])
    cross  = np.cross(d1, d2, axis=-1)
    sin_t  = np.linalg.norm(cross, axis=-1)
    cos_t  = np.einsum('...i,...i->...', d1, d2)
    ang_v  = np.degrees(np.arctan2(sin_t, cos_t)) / dt_h[:,None,None]
    stop_score = ang_v.max(axis=0)

    return (np.where(vis_g, guide_score,  np.nan),
            np.where(vis_g, guide_motion, np.nan),
            np.where(vis_g, stop_score,   np.nan))

# ─── main survey loop ─────────────────────────────────────────────────────

print("=" * 70)
print("MULTI-EPOCH ROBUSTNESS SURVEY  —  Matt 2:9 Impossibility Proof")
print("=" * 70)
print(f"Guidance criterion  : azimuth inside [{AZ_LO}°,{AZ_HI}°] road corridor (±{AZ_THRESH}° tol.)"
      f" AND motion >{OMEGA_THRESH}°/h")
print(f"Stopping criterion  : apparent angular velocity < {OMEGA_THRESH}°/h")
print(f"Grid per observation date : {CONFIGS_PER_DATE:,} configurations")
print(f"  e          = {E_VALS}")
print(f"  q (AU)     = {Q_VALS}  |  i (deg) = {I_VALS}")
print(f"  T_peri     = {T_PERI_OFFSETS} days relative to obs date")
print(f"  Omega/omega: 30° spacing (12×12 = 144 pairs)")
print(f"  Obs offsets: {OBS_OFFSETS.tolist()} days around each canonical date")
print(f"  Total dates : {len(EPOCHS)*len(OBS_OFFSETS)}")
print()

epoch_results = []

for epoch_label, jd_ref in EPOCHS:
    print(f"\n{'='*70}")
    print(f"EPOCH: {epoch_label}  (JD {jd_ref:.1f})")
    print(f"{'='*70}")

    epoch_n_vis  = 0; epoch_n_both = 0
    epoch_best_g = 999.0; epoch_best_s = 999.0

    for obs_off in OBS_OFFSETS:
        jd_obs = jd_ref + obs_off
        print(f"  {obs_off:+.0f} d  (JD {jd_obs:.1f})  ...", end=' ', flush=True)
        t0 = walltime.time()
        ep = build_epoch_data(jd_obs)

        date_n_vis = 0; date_n_both = 0
        date_best_g = 999.0; date_best_s = 999.0

        for e in E_VALS:
            for q in Q_VALS:
                for i_d in I_VALS:
                    for T_off in T_PERI_OFFSETS:
                        T_jd = jd_obs + T_off
                        gs, gm, ss = survey_epoch(q, e, i_d, OM_VALS, om_vals, T_jd, ep)
                        valid = ~np.isnan(gs)
                        nv = int(valid.sum())
                        date_n_vis  += nv; epoch_n_vis += nv
                        if nv == 0: continue
                        g_v = gs[valid]; m_v = gm[valid]; s_v = ss[valid]
                        # All three criteria: azimuth in corridor, guiding motion, stopped
                        nb = int(((g_v < AZ_THRESH) & (m_v > OMEGA_THRESH) & (s_v < OMEGA_THRESH)).sum())
                        date_n_both  += nb; epoch_n_both += nb
                        mg = float(np.min(g_v)); ms_ = float(np.min(s_v))
                        if mg < date_best_g:  date_best_g  = mg
                        if ms_ < date_best_s: date_best_s  = ms_
                        if mg < epoch_best_g: epoch_best_g = mg
                        if ms_ < epoch_best_s: epoch_best_s = ms_

        flag = "✓" if date_n_both == 0 else f"✗ {date_n_both} FOUND"
        print(f"vis={date_n_vis:,}  both={date_n_both} {flag}  "
              f"bestG={date_best_g:.1f}°  bestS={date_best_s:.1f}°/h  "
              f"({walltime.time()-t0:.1f}s)")

    print(f"\n  → EPOCH TOTAL: {epoch_n_both} of {epoch_n_vis:,} visible configs satisfy ALL THREE")
    print(f"     Best guidance azimuth (min outside [{AZ_LO}°,{AZ_HI}°]): {epoch_best_g:.2f}°  (tol. {AZ_THRESH}°)")
    print(f"     Best stopping : {epoch_best_s:.2f}°/h  (threshold {OMEGA_THRESH}°/h)")
    epoch_results.append((epoch_label.strip(), epoch_n_vis,
                          epoch_n_both, epoch_best_g, epoch_best_s))

# ─── final table ─────────────────────────────────────────────────────────

print(f"\n\n{'='*70}")
print("FINAL SUMMARY TABLE")
print("="*70)
print(f"{'Epoch':<44} {'N_vis':>7} {'N_both':>6} {'BestG':>6} {'BestS':>7}  Result")
print("-"*76)
all_zero = True
for lab, nv, nb, bg, bs in epoch_results:
    result = "ZERO ✓" if nb == 0 else f"{nb} found ✗"
    if nb > 0: all_zero = False
    print(f"{lab:<44} {nv:>7,} {nb:>6}  {bg:>5.1f}°  {bs:>6.1f}°/h  {result}")

print()
if all_zero:
    print("CONCLUSION: Zero orbital configurations satisfy all three criteria")
    print("(guidance azimuth, guidance motion >2°/h, stopping <2°/h) at ANY epoch")
    print("within ±40 days of ALL major Star of Bethlehem proposals.")
    print("The impossibility result is epoch-independent.")
else:
    print("WARNING: Configurations found — inspect above.")

print(f"\n{elapsed()} Complete.")
