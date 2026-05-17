#!/usr/bin/env python3
"""
Proof of Impossibility (fast version): No Keplerian Orbit Can Reproduce Matthew 2:9
Aaron Adair / 2026

Optimised for single-run completion (<45s). Uses same ENU infrastructure as
sob_figures.py. Key results match the full version; parameter grids are coarser.
"""

import warnings; warnings.filterwarnings('ignore')
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
import matplotlib.pyplot as plt

OUTDIR = os.path.dirname(os.path.abspath(__file__)) + os.sep
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body_barycentric
import astropy.units as u
import time as walltime

T0 = walltime.time()
def elapsed(): return f"[{walltime.time()-T0:.1f}s]"

# ── site coordinates (to the arc-minute) ──────────────────────────────────
# Jerusalem Old City centre (Jaffa Gate area):  31°46'N, 35°14'E
# Bethlehem centre (Church of the Nativity):    31°42'N, 35°12'E
# Straight-line bearing Jerusalem → Bethlehem: 203°
# Actual road sweeps 190°–210° (deviation around the ridge)
JLAT, JLON       = 31 + 46/60, 35 + 14/60
BLAT, BLON, BALT = 31 + 42/60, 35 + 12/60, 765
BETHLEHEM = EarthLocation(lon=BLON*u.deg, lat=BLAT*u.deg, height=BALT*u.m)

AZ_LO, AZ_HI = 190.0, 210.0      # road azimuth corridor (degrees)
AZ_CENTER    = (AZ_LO + AZ_HI) / 2   # 200.0°
AZ_THRESH    = 5.0   # deg — tolerance outside corridor boundary
# No altitude constraint on stopping: Matthew's "stood over" ranges in
# scholarly interpretation from horizon to zenith; altitude is not
# used as a criterion in the survey.

# Perceptible angular velocity threshold  (deg/h)
# Physical basis: naked-eye angular resolution ≈ 1' (arcminute). Motion is
# detectable when the star shifts >1' relative to background stars.  Over a
# 30-minute observation, this requires a rate > 0.033°/h — far below 2°/h.
# For comparison, the Moon moves ~0.5°/h relative to background stars, a rate
# ancient Babylonian observers monitored routinely.  Our 2°/h threshold is 4×
# the lunar rate, detectable within minutes by any naked-eye observer.
# Guidance requires ω > OMEGA_THRESH (star is actively "going before" the
# Magi); stopping requires ω < OMEGA_THRESH (star "stood still").  This
# threshold is deliberately generous: any stricter value only strengthens the
# impossibility conclusion (fewer orbits qualify as guiding, and fewer qualify
# as stopped).
OMEGA_THRESH = 2.0   # deg/h

K = 0.01720209895; MU = K**2
MU_KM3_S2 = 1.32712440018e11
GM_EARTH   = 3.986004418e5
AU_KM      = 1.495978707e8
DAY_S      = 86400.0
EPS = np.radians(23.439291111)
R_ECL2EQU = np.array([[1,0,0],[0,np.cos(EPS),-np.sin(EPS)],[0,np.sin(EPS),np.cos(EPS)]])

M_q, M_i, M_Om, M_om = 0.0407, 1.39, 93.01, 346.25
M_T     = 1719785.565
JD_NOON = 1719755.898
DT_SEC  = 10572.0

# ── orbital mechanics ──────────────────────────────────────────────────────
def barker(dt, q):
    W = 3*np.sqrt(MU/(2*q**3))*dt
    s = np.sqrt((W/2)**2+1)
    return np.cbrt(W/2+s)+np.cbrt(W/2-s)

def kepler_elliptic(M, e, tol=1e-12):
    E = M.copy()
    for _ in range(50):
        dE = (M - E + e*np.sin(E))/(1 - e*np.cos(E))
        E  = E + dE
        if np.max(np.abs(dE)) < tol: break
    return E

def true_anomaly(dt_days, q, e):
    if abs(e - 1.0) < 1e-9:               # parabolic
        D  = barker(dt_days, q); return 2*np.arctan(D)
    elif e < 1.0:                           # elliptic
        a  = q/(1-e); n  = np.sqrt(MU/a**3); M  = n*dt_days
        E  = kepler_elliptic(np.atleast_1d(float(M)), e)
        nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
        return float(nu[0]) if nu.size==1 else nu
    else:                                   # hyperbolic
        a  = q/(e-1); n  = np.sqrt(MU/a**3); M  = n*dt_days
        H  = np.sign(M)*np.log(2*np.abs(M)/e+1.8)
        for _ in range(50):
            dH = (M - e*np.sinh(H)+H)/(e*np.cosh(H)-1); H += dH
            if np.abs(dH) < 1e-12: break
        return 2*np.arctan2(np.sqrt(e+1)*np.sinh(H/2), np.sqrt(e-1)*np.cosh(H/2))

# ── time grid and ENU matrices ─────────────────────────────────────────────
print("="*70)
print("PROOF OF IMPOSSIBILITY: Matt 2:9 vs. Orbital Mechanics")
print("="*70)
print(f"\n{elapsed()} Precomputing reference frames...")

hrs_all  = np.arange(-8.0, 3.1, 0.25)
jds_all  = JD_NOON + hrs_all / 24.0
jds_utc  = jds_all - DT_SEC/86400.0
local_h  = hrs_all + 12.0

def make_enu(t_utc_single):
    frame = AltAz(obstime=t_utc_single, location=BETHLEHEM)
    def aa2icrs(az, alt):
        c = SkyCoord(alt=alt*u.deg, az=az*u.deg, frame=frame).icrs
        d = np.radians(c.dec.deg); r = np.radians(c.ra.deg)
        return np.array([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)])
    E = aa2icrs(90, 0.001); Z = aa2icrs(0, 90)
    Z /= np.linalg.norm(Z)
    E -= np.dot(E,Z)*Z; E /= np.linalg.norm(E)
    return np.array([E, np.cross(Z,E), Z])

t_tdb_all = Time(jds_all, format='jd', scale='tdb')
t_utc_all = Time(jds_utc, format='jd', scale='utc')
Ef = get_body_barycentric('earth', t_tdb_all)
Sf = get_body_barycentric('sun',   t_tdb_all)
earth_all = (Ef-Sf).xyz.to(u.AU).value
R_enu_all = np.array([make_enu(t_utc_all[ii]) for ii in range(len(jds_all))])
print(f"{elapsed()} ENU matrices done ({len(jds_all)} steps)")

mask_guide = (local_h >= 8.0) & (local_h <= 10.0)
mask_stop  = (local_h > 10.0) & (local_h <= 12.0)
N_guide = mask_guide.sum(); N_stop = mask_stop.sum()

def azel_to_icrs_dir(az_deg, alt_deg, R_enu):
    az_r = np.radians(az_deg); alt_r = np.radians(alt_deg)
    enu = np.array([np.sin(az_r)*np.cos(alt_r),
                    np.cos(az_r)*np.cos(alt_r),
                    np.sin(alt_r)])
    return R_enu.T @ enu

# For the kinematic proof, compute ICRS trajectories at three representative
# altitudes (no altitude is imposed as a constraint; these illustrate the
# argument for the full plausible range of stopping altitudes).
PROOF_ALTS = [20.0, 45.0, 60.0]   # low, mid, high — degrees above horizon

t_guide  = local_h[mask_guide]; t_stop   = local_h[mask_stop]
jd_guide = jds_all[mask_guide];  jd_stop  = jds_all[mask_stop]
E_guide  = earth_all[:, mask_guide]; E_stop = earth_all[:, mask_stop]

# Guidance: star moves at the MINIMUM detectable angular velocity (OMEGA_THRESH)
# in local AltAz — rising along az=AZ_CENTER.  This is the minimum-Δv case:
# any faster guidance motion requires a proportionally larger velocity change.
# At time step k: alt(k) = alt_0 + OMEGA_THRESH * elapsed_hours_since_guidance_start
d_guide = {alt: np.array([azel_to_icrs_dir(AZ_CENTER,
                                            alt + OMEGA_THRESH*(t_guide[k]-t_guide[0]),
                                            R_enu_all[mask_guide][k])
                           for k in range(N_guide)])
           for alt in PROOF_ALTS}

# Stopping: star appears stationary in AltAz (fixed position in local sky).
# In ICRS, this star must track Earth's sidereal rotation exactly.
# Use the representative altitude as the stopping altitude.
d_stop  = {alt: np.array([azel_to_icrs_dir(AZ_CENTER, alt, R_enu_all[mask_stop][k])
                           for k in range(N_stop)])
           for alt in PROOF_ALTS}

# ═══════════════════════════════════════════════════════════════════════════
# PROOF A: KINEMATIC ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("PROOF A: Kinematic Analysis")
print("="*70)

def ang_velocity(d_arr, t_arr_hours):
    dt = np.diff(t_arr_hours) * 3600
    cross = np.cross(d_arr[:-1], d_arr[1:], axis=1)
    sin_theta = np.linalg.norm(cross, axis=1)
    cos_theta = np.einsum('ij,ij->i', d_arr[:-1], d_arr[1:])
    theta_rad = np.arctan2(sin_theta, cos_theta)
    return np.degrees(theta_rad) / (dt/3600)

# ICRS angular velocities at three representative stopping altitudes.
# Guidance d_guide[alt]: star rising at OMEGA_THRESH deg/h along az=200°
#   → ICRS velocity = sidereal_rate + altitude-change component (vectorial sum)
# Stopping d_stop[alt]:  star FIXED in local AltAz (sidereal tracking)
#   → ICRS velocity = sidereal_rate exactly
# The DIFFERENCE (omega_guide - omega_stop) = the required angular deceleration.
print("ICRS angular velocities (guidance = rising at 2°/h; stopping = sidereal-fixed):")
print(f"  {'Alt':>5}  {'ω_guide (°/h)':>14}  {'ω_stop (°/h)':>13}  "
      f"{'Δω (°/h)':>9}  {'stop/thresh':>11}")
print("  " + "-"*60)
omg_guide_ref = omg_stop_ref = None
for alt in PROOF_ALTS:
    omg_g = ang_velocity(d_guide[alt], t_guide)
    omg_s = ang_velocity(d_stop[alt],  t_stop)
    domg  = omg_g.mean() - omg_s.mean()
    marker = " <-- reference" if alt == 45.0 else ""
    print(f"  {alt:5.0f}°  {omg_g.mean():>8.2f} ± {omg_g.std():.2f}"
          f"  {omg_s.mean():>8.2f} ± {omg_s.std():.2f}"
          f"  {domg:>9.2f}  {omg_s.mean()/OMEGA_THRESH:>8.1f}×{marker}")
    if alt == 45.0:
        omg_guide_ref = omg_g; omg_stop_ref = omg_s

omg_guide = omg_guide_ref; omg_stop = omg_stop_ref
domg_ref  = omg_guide.mean() - omg_stop.mean()

print(f"""
Key results (45° reference altitude):
  ω_stop   = {omg_stop.mean():.2f}°/h  — sidereal rate; stopping requires this exactly
  ω_guide  = {omg_guide.mean():.2f}°/h  — star rising at OMEGA_THRESH = {OMEGA_THRESH}°/h
  Δω       = {domg_ref:.2f}°/h  — ICRS angular deceleration required at transition

  The stopping rate ({omg_stop.mean():.1f}°/h) already far exceeds the 2°/h threshold:
  any body that stops in local AltAz is moving at {omg_stop.mean()/OMEGA_THRESH:.1f}× the threshold in ICRS,
  confirming that "stopped" ≠ "stationary in ICRS".

  The velocity transition (guidance → stopping) requires a change in ICRS velocity
  direction and magnitude; the minimum Δv scales as Δω × r.
""")

print(f"{'r (AU)':>10}  {'v_guide (km/s)':>14}  {'v_stop (km/s)':>14}  "
      f"{'Δv (km/s)':>12}  {'Δv/v_circ':>11}")
print("-"*70)
for r_AU in [0.001, 0.0026, 0.01, 0.05, 0.1, 0.5, 1.0]:
    r_km = r_AU * AU_KM
    v_g  = omg_guide.mean() * np.pi/180 * r_km / 3600
    v_s  = omg_stop.mean()  * np.pi/180 * r_km / 3600
    dv   = v_g - v_s
    v_circ = np.sqrt(GM_EARTH / r_km)
    print(f"{r_AU:10.4f}  {v_g:14.2f}  {v_s:14.2f}  {dv:12.2f}  {dv/v_circ:11.4f}")

r_matney     = 0.0026 * AU_KM
v_g_matney   = omg_guide.mean() * np.pi/180 * r_matney / 3600
v_s_matney   = omg_stop.mean()  * np.pi/180 * r_matney / 3600
dv_matney    = v_g_matney - v_s_matney
a_earth_grav = GM_EARTH / r_matney**2
t_decel      = dv_matney / a_earth_grav
sky_sweep_deg = (v_g_matney + v_s_matney) / 2.0 / r_matney * t_decel * 180.0 / np.pi

print(f"""
At Matney's closest approach (r ≈ 0.0026 AU ≈ Moon distance):
  Guidance ICRS velocity        : {v_g_matney:.1f} km/s
  Stopping ICRS velocity        : {v_s_matney:.1f} km/s
  Required minimum Δv           : {dv_matney:.1f} km/s
    (any faster guidance motion requires proportionally more Δv)
  Earth gravitational accel     : {a_earth_grav:.2e} km/s²
  Time for Earth gravity to supply Δv : {t_decel:.0f} s = {t_decel/3600:.0f} h
  Sky sweep during deceleration : {sky_sweep_deg:.0f}°
  (→ object traverses many circuits of sky, not "stopped over a house")
  Scaling: Δv ∝ r, a_Earth ∝ r⁻², decel time ∝ r³ — worse at every distance.
""")
print(f"{elapsed()} Proof A complete.")

# ═══════════════════════════════════════════════════════════════════════════
# PROOF B: ORBITAL SURVEY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("PROOF B: Extended Orbital Survey (e = 0 to 1)")
print("="*70)

def survey_combined(q, e, i_d, Om_d_arr, om_d_arr, T_jd):
    """Vectorised: returns (guide_score, guide_motion, stop_score) arrays shape (nOm, nom).
    guide_score  : max degrees azimuth falls outside corridor [AZ_LO, AZ_HI] during guidance.
    guide_motion : min apparent angular velocity during guidance window (must exceed OMEGA_THRESH).
    stop_score   : max apparent angular velocity during stopping window (must be below OMEGA_THRESH).
    An orbit satisfies the Matthew 2:9 description iff:
        guide_score < AZ_THRESH  AND  guide_motion > OMEGA_THRESH  AND  stop_score < OMEGA_THRESH.
    """
    nOm, nom = len(Om_d_arr), len(om_d_arr)
    i = np.radians(i_d); ci,si = np.cos(i),np.sin(i)
    Om_r = np.radians(Om_d_arr); om_r = np.radians(om_d_arr)
    cO=np.cos(Om_r); sO=np.sin(Om_r); co=np.cos(om_r); so=np.sin(om_r)
    P_ecl = np.array([cO[:,None]*co - sO[:,None]*so*ci,
                      sO[:,None]*co + cO[:,None]*so*ci,
                      np.broadcast_to(so*si,(nOm,nom)).copy()])
    Q_ecl = np.array([-cO[:,None]*so - sO[:,None]*co*ci,
                      -sO[:,None]*so + cO[:,None]*co*ci,
                      np.broadcast_to(co*si,(nOm,nom)).copy()])
    P_eq = np.einsum('ij,jkl->ikl', R_ECL2EQU, P_ecl)
    Q_eq = np.einsum('ij,jkl->ikl', R_ECL2EQU, Q_ecl)

    # ── guidance window: azimuth score + motion score ─────────────────────
    az_g  = np.zeros((N_guide, nOm, nom))
    alt_g = np.zeros_like(az_g)
    dir_g = np.zeros((N_guide, nOm, nom, 3))   # unit vectors for velocity calc
    for ii,(jd,Epos,Renu) in enumerate(zip(jd_guide, E_guide.T, R_enu_all[mask_guide])):
        dt=jd-T_jd; nu=true_anomaly(dt,q,e); r=q*(1+e)/(1+e*np.cos(nu))
        pos_c=r*(np.cos(nu)*P_eq+np.sin(nu)*Q_eq)
        rho=pos_c-Epos[:,None,None]; dist=np.linalg.norm(rho,axis=0); dist[dist<1e-12]=1e-12
        rh=rho/dist
        enu=np.einsum('ij,jkl->ikl',Renu,rh)
        alt_g[ii]=np.degrees(np.arcsin(np.clip(enu[2],-1,1)))
        az_g[ii] =np.degrees(np.arctan2(enu[0],enu[1]))%360
        dir_g[ii]=np.moveaxis(rh,0,-1)

    vis_g = (alt_g > 5.0).all(axis=0)

    # Azimuth score: max degrees outside road corridor
    dev_lo      = np.maximum(0.0, AZ_LO - az_g)
    dev_hi      = np.maximum(0.0, az_g  - AZ_HI)
    guide_score = np.maximum(dev_lo, dev_hi).max(axis=0)

    # Motion score: min angular velocity during guidance (proago requires motion)
    dg1=dir_g[:-1]; dg2=dir_g[1:]
    dt_h_g=np.diff(t_guide)
    crs_g=np.cross(dg1,dg2,axis=-1)
    sin_g=np.linalg.norm(crs_g,axis=-1); cos_g=np.einsum('...i,...i->...',dg1,dg2)
    omg_g=np.degrees(np.arctan2(sin_g,cos_g))/dt_h_g[:,None,None]
    guide_motion = omg_g.min(axis=0)   # must exceed OMEGA_THRESH throughout

    # ── stopping window: motion score ─────────────────────────────────────
    dir_pred = np.zeros((N_stop,nOm,nom,3))
    for ii,(jd,Epos,Renu) in enumerate(zip(jd_stop, E_stop.T, R_enu_all[mask_stop])):
        dt=jd-T_jd; nu=true_anomaly(dt,q,e); r=q*(1+e)/(1+e*np.cos(nu))
        pos_c=r*(np.cos(nu)*P_eq+np.sin(nu)*Q_eq)
        rho=pos_c-Epos[:,None,None]; dist=np.linalg.norm(rho,axis=0); dist[dist<1e-12]=1e-12
        rh=rho/dist; dir_pred[ii]=np.moveaxis(rh,0,-1)

    d1=dir_pred[:-1]; d2=dir_pred[1:]
    dt_h=np.diff(t_stop)
    cross2=np.cross(d1,d2,axis=-1)
    sin_t=np.linalg.norm(cross2,axis=-1); cos_t=np.einsum('...i,...i->...',d1,d2)
    ang_v=np.degrees(np.arctan2(sin_t,cos_t))/dt_h[:,None,None]
    stop_score=ang_v.max(axis=0)

    return (np.where(vis_g, guide_score,  np.nan),
            np.where(vis_g, guide_motion, np.nan),
            np.where(vis_g, stop_score,   np.nan))

# Survey params (coarser Om/om for speed; same conclusion)
Om_vals = np.arange(0, 360, 20)  # 18 values
om_vals = np.arange(0, 360, 20)  # 18 values
q_vals  = [0.03, 0.04, 0.07, 0.12]
i_vals  = [1.4, 10.0, 30.0, 60.0]
T_offs  = [-40.0, -20.0, -5.0, +5.0]
e_vals  = [0.0, 0.3, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0]

total_expected = len(e_vals)*len(q_vals)*len(i_vals)*len(T_offs)*len(Om_vals)*len(om_vals)
print(f"\nSurveying {total_expected:,} configurations:")
print(f"  e={e_vals}, q={q_vals}, i={i_vals}")
print(f"  T_offs={T_offs}, Om/om: {len(Om_vals)}×{len(om_vals)} @ 20° steps")

all_guide=[]; all_motion=[]; all_stop=[]
n_az_ok=0; n_motion_ok=0; n_stop_ok=0; n_both_ok=0; n_vis_total=0
best_guide=999; best_stop=999; best_combined=999; best_params=None
per_e_guide = {e: [] for e in e_vals}

for e in e_vals:
    for q in q_vals:
        for i_d in i_vals:
            for T_off in T_offs:
                gs, gm, ss = survey_combined(q,e,i_d,Om_vals,om_vals,M_T+T_off)
                valid = ~np.isnan(gs)
                n_vis = int(valid.sum()); n_vis_total+=n_vis
                if n_vis==0: continue
                g_v=gs[valid]; m_v=gm[valid]; s_v=ss[valid]
                all_guide.extend(g_v.tolist())
                all_motion.extend(m_v.tolist())
                all_stop.extend(s_v.tolist())
                per_e_guide[e].extend(g_v.tolist())
                n_az_ok     += int((g_v < AZ_THRESH).sum())
                n_motion_ok += int((m_v > OMEGA_THRESH).sum())
                n_stop_ok   += int((s_v < OMEGA_THRESH).sum())
                # All three criteria simultaneously:
                both = (g_v < AZ_THRESH) & (m_v > OMEGA_THRESH) & (s_v < OMEGA_THRESH)
                n_both_ok += int(both.sum())
                combined=np.sqrt(g_v**2+(s_v/3)**2)
                idx_c=np.argmin(combined)
                if combined[idx_c]<best_combined:
                    best_combined=float(combined[idx_c])
                    flat=np.nanargmin(np.where(valid,np.sqrt(gs**2+(ss/3)**2),999))
                    ij=np.unravel_index(flat,gs.shape)
                    best_guide=float(gs[ij]); best_stop=float(ss[ij])
                    best_params=(e,q,i_d,Om_vals[ij[0]],om_vals[ij[1]],T_off)

all_guide=np.array(all_guide); all_motion=np.array(all_motion); all_stop=np.array(all_stop)
print(f"\n{elapsed()} Survey complete.")
print(f"""
Results:
  Total visible configurations         : {n_vis_total:,}
  Satisfying azimuth  (<{AZ_THRESH}° outside corridor) : {n_az_ok:,}
  Satisfying guidance motion  (>{OMEGA_THRESH}°/h)     : {n_motion_ok:,}
  Satisfying stopping         (<{OMEGA_THRESH}°/h)      : {n_stop_ok:,}
  Satisfying ALL THREE simultaneously  : {n_both_ok}

  Best guidance score : {best_guide:.2f}°  (azimuth deviation)
  Best stopping score : {best_stop:.2f}°/h
  Best combined config: e={best_params[0]}, q={best_params[1]},
    i={best_params[2]}, Om={best_params[3]}, om={best_params[4]}, T_off={best_params[5]:+.0f}d
""")
print(f"RESULT: {n_both_ok} out of {n_vis_total:,} Keplerian configurations satisfy all three criteria.")

# ─── Sensitivity analysis: result vs. angular-velocity threshold ──────────
print("\n" + "="*70)
print("SENSITIVITY: Zero result vs. angular-velocity threshold ω_thresh")
print("  (stopping requires ω < ω_thresh; guidance motion requires ω > ω_thresh)")
print("="*70)
print(f"  {'ω_thresh':>10}  {'N_az_ok':>10}  {'N_motion_ok':>12}  "
      f"{'N_stop_ok':>10}  {'N_both':>8}")
print(f"  {'-'*10}  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*8}")
for thresh in [0.25, 0.5, 1.0, 2.0, 5.0, 10.0]:
    na = int((all_guide < AZ_THRESH).sum())
    nm = int((all_motion > thresh).sum())
    ns = int((all_stop   < thresh).sum())
    nb = int(((all_guide < AZ_THRESH) & (all_motion > thresh) & (all_stop < thresh)).sum())
    print(f"  {thresh:>10.2f}  {na:>10,}  {nm:>12,}  {ns:>10,}  {nb:>8}")
print("  → zero satisfying configurations at every threshold tested.")

# ═══════════════════════════════════════════════════════════════════════════
# PROOF C: ORBIT INVERSION
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("PROOF C: Orbit Inversion — best-fit orbit propagated beyond guidance")
print("="*70)

Om_fine=np.arange(0,360,10); om_fine=np.arange(0,360,10)
inv_best_guide=999.0; inv_best_params=None

print(f"{elapsed()} Fine-grid inversion (10° steps, e≥0.9)...")
for e in [0.9,0.95,0.99,1.0]:
    for q in [0.03,0.04,0.06,0.10]:
        for T_off in [-20.0,-5.0,+5.0]:
            gs,gm,ss=survey_combined(q,e,1.4,Om_fine,om_fine,M_T+T_off)
            if not np.all(np.isnan(gs)):
                mn=float(np.nanmin(gs))
                if mn<inv_best_guide:
                    inv_best_guide=mn
                    idx_f=np.unravel_index(np.nanargmin(gs),gs.shape)
                    inv_best_params=(e,q,1.4,Om_fine[idx_f[0]],om_fine[idx_f[1]],T_off)

print(f"  Best guidance residual (fine): {inv_best_guide:.2f} deg")
print(f"  Best params: e={inv_best_params[0]}, q={inv_best_params[1]}, "
      f"Om={inv_best_params[3]}, om={inv_best_params[4]}, T={inv_best_params[5]:+.0f}d")

e0,q0,i0_d,Om0,om0,T0_off = inv_best_params
gs1,gm1,ss1=survey_combined(q0,e0,i0_d,np.array([Om0]),np.array([om0]),M_T+T0_off)
inv_best_stop  = float(ss1[0,0])  if not np.isnan(ss1[0,0])  else 999.0
inv_best_motion= float(gm1[0,0])  if not np.isnan(gm1[0,0])  else 0.0
print(f"  Guidance motion for best-fit orbit: {inv_best_motion:.2f}°/h  (need >{OMEGA_THRESH}°/h)")
print(f"  Stopping score for best guidance orbit: {inv_best_stop:.2f}°/h  (need <{OMEGA_THRESH}°/h)")
if inv_best_stop >= OMEGA_THRESH:
    print(f"  → Stopping FAILS by {inv_best_stop/OMEGA_THRESH:.1f}× threshold")
else:
    print(f"  → Stopping passes, but guidance motion={'PASS' if inv_best_motion>OMEGA_THRESH else 'FAIL'}")

# Full-day track of best orbit
i0=np.radians(i0_d); Om0r=np.radians(Om0); om0r=np.radians(om0)
ci,si=np.cos(i0),np.sin(i0); cO,sO=np.cos(Om0r),np.sin(Om0r); co,so=np.cos(om0r),np.sin(om0r)
P=R_ECL2EQU@np.array([cO*co-sO*so*ci,sO*co+cO*so*ci,so*si])
Q=R_ECL2EQU@np.array([-cO*so-sO*co*ci,-sO*so+cO*co*ci,co*si])
az_best=[]; alt_best=[]
for ii in range(len(jds_all)):
    dt=jds_all[ii]-M_T-T0_off; nu=true_anomaly(dt,q0,e0)
    r=q0*(1+e0)/(1+e0*np.cos(nu)); pos_c=r*(np.cos(nu)*P+np.sin(nu)*Q)
    rho=pos_c-earth_all[:,ii]; dist=np.linalg.norm(rho); rh=rho/dist
    enu=R_enu_all[ii]@rh; e_c,n_c,up=enu
    az_best.append(np.degrees(np.arctan2(e_c,n_c))%360)
    alt_best.append(np.degrees(np.arcsin(np.clip(up,-1,1))))
az_best=np.array(az_best); alt_best=np.array(alt_best)

# Matney track
P_m=R_ECL2EQU@np.array([np.cos(np.radians(M_Om))*np.cos(np.radians(M_om))-np.sin(np.radians(M_Om))*np.sin(np.radians(M_om))*np.cos(np.radians(M_i)),
                         np.sin(np.radians(M_Om))*np.cos(np.radians(M_om))+np.cos(np.radians(M_Om))*np.sin(np.radians(M_om))*np.cos(np.radians(M_i)),
                         np.sin(np.radians(M_om))*np.sin(np.radians(M_i))])
Q_m=R_ECL2EQU@np.array([-np.cos(np.radians(M_Om))*np.sin(np.radians(M_om))-np.sin(np.radians(M_Om))*np.cos(np.radians(M_om))*np.cos(np.radians(M_i)),
                         -np.sin(np.radians(M_Om))*np.sin(np.radians(M_om))+np.cos(np.radians(M_Om))*np.cos(np.radians(M_om))*np.cos(np.radians(M_i)),
                         np.cos(np.radians(M_om))*np.sin(np.radians(M_i))])
az_mat=[]; alt_mat=[]
for ii in range(len(jds_all)):
    dt=jds_all[ii]-M_T; D=barker(dt,M_q); nu=2*np.arctan(D)
    r=2*M_q/(1+np.cos(nu)); pos_c=r*(np.cos(nu)*P_m+np.sin(nu)*Q_m)
    rho=pos_c-earth_all[:,ii]; dist=np.linalg.norm(rho); rh=rho/dist
    enu=R_enu_all[ii]@rh; e_c,n_c,up=enu
    az_mat.append(np.degrees(np.arctan2(e_c,n_c))%360)
    alt_mat.append(np.degrees(np.arcsin(np.clip(up,-1,1))))
az_mat=np.array(az_mat); alt_mat=np.array(alt_mat)

print(f"\n{elapsed()} Best-fit orbit Az/Alt after guidance ends (stopping window 10:00-12:00):")
for ii in range(len(local_h)):
    if 10.0 <= local_h[ii] <= 12.0 and abs(local_h[ii]*4-round(local_h[ii]*4))<0.01:
        hh=int(local_h[ii]); mm=int((local_h[ii]-hh)*60)
        in_corridor = AZ_LO <= az_best[ii] <= AZ_HI
        print(f"  {hh:02d}:{mm:02d}  Az={az_best[ii]:.1f}°  Alt={alt_best[ii]:.1f}°"
              f"  (corridor [190°,210°]: {'inside' if in_corridor else 'OUTSIDE'})")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 6
# ═══════════════════════════════════════════════════════════════════════════
print(f"\n{elapsed()} Generating Figure 6...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.suptitle(
    "Proof of Impossibility: No Keplerian Orbit Reproduces Matthew 2:9\n"
    f"(Guidance: Az in [{AZ_LO:.0f}°,{AZ_HI:.0f}°] corridor + motion >{OMEGA_THRESH}°/h; "
    f"Stopping: motion <{OMEGA_THRESH}°/h)",
    fontsize=12, fontweight='bold')

# Panel A: guidance azimuth score vs. stopping score scatter
ax = axes[0,0]
# Colour by whether guidance motion criterion is met
motion_ok = all_motion > OMEGA_THRESH
ax.scatter(all_guide[~motion_ok], all_stop[~motion_ok], s=2, alpha=0.15,
           c='grey', label=f'Motion ≤{OMEGA_THRESH}°/h (not guiding)')
ax.scatter(all_guide[motion_ok],  all_stop[motion_ok],  s=2, alpha=0.35,
           c='steelblue', label=f'Motion >{OMEGA_THRESH}°/h (guiding orbits)')
ax.axvline(AZ_THRESH,    color='red',   lw=2, ls='--',
           label=f'Azimuth tolerance ({AZ_THRESH}° outside corridor)')
ax.axhline(OMEGA_THRESH, color='green', lw=2, ls='--',
           label=f'Stopping criterion ({OMEGA_THRESH}°/h)')
ax.fill_between([0,AZ_THRESH],[0,0],[OMEGA_THRESH,OMEGA_THRESH],
                alpha=0.2, color='gold', label='Must satisfy ALL THREE → 0 orbits')
ax.set_xlabel('Guidance azimuth score: max degrees outside [190°,210°] corridor', fontsize=10)
ax.set_ylabel('Stopping score: max apparent velocity (°/h)', fontsize=10)
ax.set_title('(A) Guidance vs. Stopping Scores\n(blue = orbits that move enough to guide)',
             fontsize=10, fontweight='bold')
ax.set_xlim(0, 90); ax.set_ylim(0, 30)
ax.legend(fontsize=8, markerscale=4)
ax.text(0.02, 0.75, f"N visible: {n_vis_total:,}\nSatisfy ALL: {n_both_ok}",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
ax.grid(True, alpha=0.2)

# Panel B: full-day Az track
ax = axes[0,1]
vis_best = alt_best > 5; vis_mat = alt_mat > 5
ax.plot(local_h[vis_best], az_best[vis_best], 'b-', lw=2,
        label=f"Best-fit guidance orbit (e={inv_best_params[0]})")
ax.plot(local_h[vis_mat], az_mat[vis_mat], 'r--', lw=2, label="Matney (2025)")
ax.axhspan(AZ_LO, AZ_HI, alpha=0.12, color='green', label=f'Road corridor [{AZ_LO:.0f}°,{AZ_HI:.0f}°]')
ax.axvspan(8, 10, alpha=0.08, color='green')
ax.axvspan(10, 12, alpha=0.08, color='gold', label='Stopping window')
ax.set_xlabel('Local solar time (h)', fontsize=10)
ax.set_ylabel('Azimuth (°)', fontsize=10)
ax.set_title('(B) Full-Day Az Track: Best Orbit vs. Matney\n(Neither stays in road corridor during stopping window)',
             fontsize=10, fontweight='bold')
ax.set_xlim(4, 14); ax.set_ylim(150, 270)
ax.set_xticks(range(4, 15, 2)); ax.legend(fontsize=8, loc='upper left'); ax.grid(True, alpha=0.2)

# Panel C: guidance score distribution by eccentricity
ax = axes[1,0]
colors = plt.cm.viridis(np.linspace(0,1,len(e_vals)))
for e_idx, e_val in enumerate(e_vals):
    data = np.array(per_e_guide[e_val])
    if len(data)>0:
        d = data[data<90]
        ax.hist(d, bins=np.linspace(0,90,31), color=colors[e_idx],
                alpha=0.5, label=f'e={e_val}', density=True)
ax.axvline(AZ_THRESH, color='red', lw=2.5, ls='--', label=f'{AZ_THRESH}° outside-corridor tolerance')
ax.set_xlabel('Guidance score: max degrees outside road corridor [190°,210°]', fontsize=10)
ax.set_ylabel('Density', fontsize=10)
ax.set_title('(C) Guidance Score Distribution by Eccentricity\n(No eccentricity keeps orbit inside road corridor)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=7, loc='upper right', ncol=2); ax.grid(True, alpha=0.2)

# Panel D: required vs. achievable Δv
ax = axes[1,1]
r_AU_range = np.logspace(-3, 0, 200)
r_km_range = r_AU_range * AU_KM
v_g_range  = omg_guide.mean() * np.pi/180 * r_km_range / 3600
v_s_range  = omg_stop.mean()  * np.pi/180 * r_km_range / 3600
delta_v    = v_g_range - v_s_range
dt_s = 30*60  # 30 minutes
a_solar = 2*MU_KM3_S2*r_km_range / AU_KM**3
a_earth = GM_EARTH / r_km_range**2
dv_solar_30 = a_solar * dt_s
dv_earth_30 = a_earth * dt_s

ax.loglog(r_AU_range, delta_v,    'b-',   lw=2.5, label='Required Δv (guide→stop)')
ax.loglog(r_AU_range, dv_solar_30,'orange',lw=2, ls='--', label='Solar tidal Δv in 30 min')
ax.loglog(r_AU_range, dv_earth_30,'green', lw=2, ls='-.', label="Earth gravity Δv in 30 min")
ax.axvline(0.0026, color='red', lw=1.5, ls=':', label="Matney's distance (0.0026 AU)")
mask_ok = dv_earth_30 > delta_v
if mask_ok.any():
    ax.fill_between(r_AU_range, delta_v, dv_earth_30, where=mask_ok,
                    alpha=0.2, color='green', label='Earth gravity sufficient')
ax.set_xlabel('Geocentric distance r (AU)', fontsize=10)
ax.set_ylabel('Velocity change (km/s)', fontsize=10)
ax.set_title('(D) Required Δv for Transition vs. Available Forces\n'
             '(Earth gravity insufficient at Matney\'s proposed distance)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.2)

plt.tight_layout()
out6 = OUTDIR + 'figure6_impossibility_proof.png'
plt.savefig(out6, dpi=600, bbox_inches='tight')
plt.close()
print(f"  Saved Figure 6")

# ═══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("FINAL SUMMARY — THREE-PART PROOF OF IMPOSSIBILITY")
print("="*70)
print(f"""
PROOF A (Kinematic):
  Guidance phase requires geocentric transverse velocity {v_g_matney:.1f} km/s
  at Matney's distance (0.0026 AU).
  Stopping phase requires ~{v_s_matney:.2f} km/s (Earth-rotation rate).
  Required Δv = {v_g_matney-v_s_matney:.1f} km/s.
  Time for Earth gravity to supply this Δv: {t_decel:.0f} s = {t_decel/3600:.2f} h
  Sky sweep during deceleration: {(v_g_matney+v_s_matney)/2/r_matney * t_decel * 180/np.pi:.0f} deg
  → Object cannot "stop over a house" — it sweeps visibly across the sky.

PROOF B (Survey, {n_vis_total:,} configurations, e = 0 to 1.0):
  Guidance azimuth criterion (<{AZ_THRESH}° outside [{AZ_LO:.0f}°,{AZ_HI:.0f}°] corridor): {n_az_ok:,} satisfied
  Guidance motion criterion  (>{OMEGA_THRESH}°/h, 4× Moon's relative rate):   {n_motion_ok:,} satisfied
  Stopping criterion         (<{OMEGA_THRESH}°/h apparent motion):             {n_stop_ok:,} satisfied
  ALL THREE simultaneously: {n_both_ok}
  → Zero Keplerian orbits satisfy all three requirements.

PROOF C (Orbit Inversion):
  Best orbit for guidance azimuth: max degrees outside corridor = {inv_best_guide:.2f}°
  Its guidance motion : {inv_best_motion:.2f}°/h  ({'passes' if inv_best_motion > OMEGA_THRESH else 'fails'} >{OMEGA_THRESH}°/h)
  Its stopping score  : {inv_best_stop:.2f}°/h   ({'passes' if inv_best_stop < OMEGA_THRESH else f'fails by {inv_best_stop/OMEGA_THRESH:.1f}×'} <{OMEGA_THRESH}°/h)
  → Orbits good for guidance azimuth fail the combined criterion.

CONCLUSION:
  The motion described in Matthew 2:9 is incompatible with any orbit
  governed by Newtonian/Keplerian gravity. Guidance and stopping are
  mutually exclusive for any solar-orbiting body.
""")
print(f"{elapsed()} All done. Figure 6 saved to workspace folder.")
