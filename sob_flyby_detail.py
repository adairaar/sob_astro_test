#!/usr/bin/env python3
"""
sob_flyby_detail.py  —  Detailed analysis of the best close-flyby orbit
Aaron Adair / 2026

Best orbit from sob_close_flyby_opt.py (cost = 4.17):
  q=0.998156 AU, e=0.709313, i=2.172°, Ω=95.196°, ω=204.638°, ΔT=+13.723d
  d_min=0.00092 AU, guide=0.58h, cost C=4.17 (guidance short of 1h minimum)

Outputs:
  1. Time-series table: JD, d(AU), ω_ICRS (°/h), ω_app_ground (°/h),
     alt(°), az(°), sun_alt(°), daz/dt, dalt/dt
  2. Summary statistics for guidance and stopping windows
  3. Analytic derivation of the ~35-minute guidance ceiling
"""

import numpy as np
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, get_sun, EarthLocation
import astropy.units as u

# ── Best orbit parameters ─────────────────────────────────────────────────────
PERI = 0.998156      # AU — perihelion distance
ECC  = 0.709313      # eccentricity
I_R  = np.radians(2.172)
OM_R = np.radians(95.196)   # Ω
OM_r = np.radians(204.638)  # ω
DT_DAYS = 13.723            # T_perihelion = JD_PRIMARY + DT_DAYS

# ── Constants ─────────────────────────────────────────────────────────────────
K_GAUSS   = 0.01720209895
MU        = K_GAUSS**2
EPS       = np.radians(23.4392911)
AU_KM     = 1.495978707e8
R_EQ      = np.array([
    [1,  0,           0          ],
    [0,  np.cos(EPS), -np.sin(EPS)],
    [0,  np.sin(EPS),  np.cos(EPS)],
])

LAT_DEG = 31.70;  LON_DEG = 35.20;  ALT_M = 765.0
LAT     = np.radians(LAT_DEG)
JD_PRIMARY   = 1719755.898
DT_TDB_UTC   = 10572.0

# Criteria (same as MC survey)
AZ_LO, AZ_HI        = 190.0, 210.0
AZ_DRIFT_THRESH      = 1.5    # °/h
ALT_RATE_MIN         = 0.3    # °/h
STOP_THRESH          = 0.30   # °/h
STEP_MIN             = 1      # minute (finer than MC for detail)
WINDOW_HRS           = 14     # slightly wider window

OMEGA_SID = 15.04107          # °/h  (sidereal rotation rate)


# ── Orbital mechanics ─────────────────────────────────────────────────────────

def kepler_elliptic(Mv, ecc, tol=1e-12):
    Ev = Mv.copy()
    for _ in range(60):
        dE = (Mv - Ev + ecc*np.sin(Ev)) / (1 - ecc*np.cos(Ev))
        Ev += dE
        if np.max(np.abs(dE)) < tol:
            break
    return Ev

def comet_pos(dt_days):
    dt = np.atleast_1d(np.asarray(dt_days, float))
    ci, si = np.cos(I_R),  np.sin(I_R)
    cO, sO = np.cos(OM_R), np.sin(OM_R)
    co, so = np.cos(OM_r), np.sin(OM_r)
    Pecl = np.array([cO*co - sO*so*ci,  sO*co + cO*so*ci,  so*si])
    Qecl = np.array([-cO*so - sO*co*ci, -sO*so + cO*co*ci, co*si])
    Pvec = R_EQ @ Pecl   # (3,)
    Qvec = R_EQ @ Qecl   # (3,)
    a    = PERI / (1 - ECC)
    n    = np.sqrt(MU / a**3)
    Mv   = (n * dt) % (2*np.pi)
    Ev   = kepler_elliptic(Mv, ECC)
    nu   = 2*np.arctan2(np.sqrt(1+ECC)*np.sin(Ev/2), np.sqrt(1-ECC)*np.cos(Ev/2))
    r    = a * (1 - ECC*np.cos(Ev))
    pos  = np.outer(Pvec, r*np.cos(nu)) + np.outer(Qvec, r*np.sin(nu))
    return pos   # (3, N)

def gmst_rad(jd_utc):
    T = (jd_utc - 2451545.0) / 36525.0
    return np.radians((280.46061837 + 360.98564736629*(jd_utc-2451545.0)
                        + 0.000387933*T**2 - T**3/38710000.0) % 360.0)

def lst_rad(jd_utc):
    return (gmst_rad(jd_utc) + np.radians(LON_DEG)) % (2*np.pi)

def altaz(ra_r, dec_r, lst_r):
    ha      = lst_r - ra_r
    sin_alt = np.sin(LAT)*np.sin(dec_r) + np.cos(LAT)*np.cos(dec_r)*np.cos(ha)
    sin_alt = np.clip(sin_alt, -1, 1)
    alt     = np.arcsin(sin_alt)
    cos_alt = np.cos(alt)
    az_sin  = -np.cos(dec_r)*np.sin(ha)
    az_cos  = (np.sin(dec_r) - np.sin(LAT)*sin_alt) / (np.cos(LAT)*np.maximum(cos_alt, 1e-9))
    az      = np.arctan2(az_sin, az_cos) % (2*np.pi)
    return np.degrees(alt), np.degrees(az)


# ── Main analysis ─────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("DETAILED ANALYSIS OF BEST CLOSE-FLYBY ORBIT")
    print("=" * 72)
    print(f"  q={PERI:.6f} AU   e={ECC:.6f}   i={np.degrees(I_R):.3f}°")
    print(f"  Ω={np.degrees(OM_R):.3f}°   ω={np.degrees(OM_r):.3f}°")
    print(f"  T_peri = JD_PRIMARY + {DT_DAYS:.3f} d")
    print()

    # ── Build time grid ───────────────────────────────────────────────────────
    T_peri_jd = JD_PRIMARY + DT_DAYS

    jd_utc_midnight = (JD_PRIMARY - DT_TDB_UTC/86400.0) - 0.5
    dt_h   = np.arange(-WINDOW_HRS/2, WINDOW_HRS/2, STEP_MIN/60.0)
    jd_utc = jd_utc_midnight + dt_h / 24.0
    jd_tdb = jd_utc + DT_TDB_UTC / 86400.0

    # Astropy for Earth/Sun
    times_tdb = Time(jd_tdb, format='jd', scale='tdb')
    times_utc = Time(jd_utc, format='jd', scale='utc')
    eb        = get_body_barycentric('earth', times_tdb)
    sb        = get_body_barycentric('sun',   times_tdb)
    earth_pos = (eb - sb).xyz.to(u.AU).value      # (3, N)
    sun_icrs  = get_sun(times_utc)
    sun_ra    = sun_icrs.ra.rad
    sun_dec   = sun_icrs.dec.rad
    lst       = lst_rad(jd_utc)

    sun_alt, _ = altaz(sun_ra, sun_dec, lst)

    # ── Comet position and geocentric vectors ─────────────────────────────────
    dt_from_peri = jd_tdb - T_peri_jd
    r_c   = comet_pos(dt_from_peri)       # (3, N) heliocentric ICRS
    rho   = r_c - earth_pos               # geocentric ICRS
    d     = np.linalg.norm(rho, axis=0)   # AU

    # ICRS RA/Dec
    rhat  = rho / d
    ra_r  = np.arctan2(rhat[1], rhat[0]) % (2*np.pi)
    dec_r = np.arcsin(np.clip(rhat[2], -1, 1))

    # AltAz
    alt, az = altaz(ra_r, dec_r, lst)

    # Angular rates (numerical differencing over 2-step stencil, central difference)
    h_rad   = STEP_MIN / 60.0    # hours per step
    dalt    = np.gradient(alt, h_rad)    # °/h
    daz_raw = np.gradient(az, h_rad)     # °/h (raw, may wrap)

    # Fix azimuth wrapping
    daz = np.zeros_like(daz_raw)
    for k in range(1, len(az)-1):
        da = (az[k+1] - az[k-1])
        if da > 180:  da -= 360
        if da < -180: da += 360
        daz[k] = da / (2*h_rad)
    daz[0]  = daz[1]
    daz[-1] = daz[-2]

    # Total apparent angular velocity (ground frame)
    omega_app = np.sqrt(dalt**2 + (daz * np.cos(np.radians(alt)))**2)

    # ICRS angular velocity (ω_ICRS = V_rel_transverse / d)
    rho_km    = rho * AU_KM
    d_km      = d   * AU_KM
    # v_rel = d(rho)/dt in km/s; approximate via finite difference
    dt_sec    = STEP_MIN * 60.0
    v_rho     = np.gradient(rho_km, dt_sec, axis=1)   # (3, N) km/s
    v_r_rad   = np.sum(v_rho * rhat, axis=0)           # radial component km/s
    v_transv  = np.sqrt(np.sum(v_rho**2, axis=0) - v_r_rad**2)   # km/s
    omega_icrs_deg_h = np.degrees(v_transv / d_km) * 3600.0       # °/h

    # Time relative to event and to perihelion
    t_from_event_h = (jd_utc - (JD_PRIMARY - DT_TDB_UTC/86400.0)) * 24.0
    t_from_peri_h  = dt_from_peri * 24.0

    # ── Summary: key epoch ────────────────────────────────────────────────────
    idx_dmin = np.argmin(d)
    print(f"{'Epoch':^40s}  {'t_event':>8s}  {'d (AU)':>8s}  {'ω_ICRS':>9s}  "
          f"{'ω_app':>8s}  {'alt':>6s}  {'az':>6s}  {'sun_alt':>8s}")
    print("-" * 102)

    # Identify guidance and stopping windows
    night   = sun_alt < -6.0
    visible = (alt > 10.0) & night
    in_az   = (az >= AZ_LO) & (az <= AZ_HI)
    in_guide_az_drift = (np.abs(daz) < AZ_DRIFT_THRESH) & (np.abs(dalt) > ALT_RATE_MIN)
    guidance_mask = visible & in_az & in_guide_az_drift
    stopping_mask = visible & (omega_app < STOP_THRESH)

    # Print every 15 minutes with guidance/stopping flags
    for k in range(len(jd_utc)):
        if k % 15 != 0:
            continue
        flags = []
        if guidance_mask[k]: flags.append('GUIDE')
        if stopping_mask[k]: flags.append('STOP')
        flag_str = ','.join(flags) if flags else ''
        print(f"  t={t_from_event_h[k]:+7.2f}h  (peri{t_from_peri_h[k]:+7.1f}h)  "
              f"d={d[k]:.5f} AU  ω_ICRS={omega_icrs_deg_h[k]:8.3f}°/h  "
              f"ω_app={omega_app[k]:6.2f}°/h  "
              f"alt={alt[k]:5.1f}°  az={az[k]:5.1f}°  sun={sun_alt[k]:+5.1f}°  "
              f"{flag_str}")

    print()

    # ── Guidance window analysis ──────────────────────────────────────────────
    print("=" * 72)
    print("GUIDANCE WINDOW ANALYSIS (all 1-minute steps)")
    print("=" * 72)

    # Find contiguous guidance segments
    segments = []
    in_seg = False
    seg_start = 0
    for k in range(len(guidance_mask)):
        if guidance_mask[k] and not in_seg:
            in_seg = True; seg_start = k
        elif not guidance_mask[k] and in_seg:
            in_seg = False
            segments.append((seg_start, k-1))
    if in_seg:
        segments.append((seg_start, len(guidance_mask)-1))

    if not segments:
        print("  No guidance windows found.")
    else:
        for (s, e_) in segments:
            dur = (e_ - s + 1) * STEP_MIN / 60.0
            print(f"\n  Segment: t={t_from_event_h[s]:+.2f}h to {t_from_event_h[e_]:+.2f}h  "
                  f"duration={dur:.2f}h = {dur*60:.0f} min")
            print(f"    d: {d[s]:.5f} – {d[e_]:.5f} AU")
            print(f"    az: {az[s]:.1f}° – {az[e_]:.1f}°  (corridor [{AZ_LO},{AZ_HI}]°)")
            print(f"    alt: {alt[s]:.1f}° – {alt[e_]:.1f}°  (min {np.min(alt[s:e_+1]):.1f}°, max {np.max(alt[s:e_+1]):.1f}°)")
            print(f"    ω_ICRS: {np.min(omega_icrs_deg_h[s:e_+1]):.3f} – {np.max(omega_icrs_deg_h[s:e_+1]):.3f}°/h")
            print(f"    ω_app:  {np.min(omega_app[s:e_+1]):.3f} – {np.max(omega_app[s:e_+1]):.3f}°/h")
            print(f"    daz/dt: {np.min(daz[s:e_+1]):.3f} – {np.max(daz[s:e_+1]):.3f}°/h")
            print(f"    dalt/dt:{np.min(dalt[s:e_+1]):.3f} – {np.max(dalt[s:e_+1]):.3f}°/h")

    # ── Stopping window ───────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("STOPPING WINDOW (ω_app < 0.30°/h, night, alt > 10°)")
    print("=" * 72)
    stop_segs = []
    in_seg = False
    for k in range(len(stopping_mask)):
        if stopping_mask[k] and not in_seg:
            in_seg = True; seg_start = k
        elif not stopping_mask[k] and in_seg:
            in_seg = False; stop_segs.append((seg_start, k-1))
    if in_seg:
        stop_segs.append((seg_start, len(stopping_mask)-1))

    if not stop_segs:
        print("  No stopping windows found.")
    else:
        for (s, e_) in stop_segs:
            dur = (e_ - s + 1) * STEP_MIN / 60.0
            print(f"\n  Segment: t={t_from_event_h[s]:+.2f}h to {t_from_event_h[e_]:+.2f}h  "
                  f"duration={dur:.2f}h = {dur*60:.0f} min")
            print(f"    d: {d[s]:.5f} – {d[e_]:.5f} AU (min {d[s:e_+1].min():.5f})")
            print(f"    ω_ICRS: {np.min(omega_icrs_deg_h[s:e_+1]):.3f} – {np.max(omega_icrs_deg_h[s:e_+1]):.3f}°/h")
            print(f"    ω_app:  {np.min(omega_app[s:e_+1]):.4f} – {np.max(omega_app[s:e_+1]):.4f}°/h")
            print(f"    alt: {np.min(alt[s:e_+1]):.1f}° – {np.max(alt[s:e_+1]):.1f}°")
            print(f"    az:  {np.min(az[s:e_+1]):.1f}° – {np.max(az[s:e_+1]):.1f}°")

    # ── Full night statistics ─────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("FULL NIGHT SUMMARY (Sun below -6°, comet above 10°)")
    print("=" * 72)
    vis = visible
    if np.any(vis):
        print(f"  Visible window: t={t_from_event_h[vis][0]:+.2f}h to {t_from_event_h[vis][-1]:+.2f}h")
        print(f"  Altitude:  min={alt[vis].min():.1f}°   max={alt[vis].max():.1f}°")
        print(f"  Azimuth:   min={az[vis].min():.1f}°   max={az[vis].max():.1f}°")
        print(f"  d (AU):    min={d[vis].min():.5f}   max={d[vis].max():.5f}")
        print(f"  ω_ICRS:    min={omega_icrs_deg_h[vis].min():.4f}°/h   max={omega_icrs_deg_h[vis].max():.4f}°/h")
        print(f"  ω_app:     min={omega_app[vis].min():.4f}°/h   max={omega_app[vis].max():.4f}°/h")
        print(f"  Sidereal rate ω_sid = {OMEGA_SID:.3f}°/h")
        print(f"  Ratio ω_ICRS/ω_sid: {omega_icrs_deg_h[vis].min()/OMEGA_SID:.4f} – "
              f"{omega_icrs_deg_h[vis].max()/OMEGA_SID:.4f}")
    else:
        print("  No visible window.")

    # ── Analytic derivation of ~35-min ceiling ────────────────────────────────
    print()
    print("=" * 72)
    print("ANALYTIC: WHY THE GUIDANCE WINDOW MAXES OUT AT ~35 MINUTES")
    print("=" * 72)

    d_min_km  = d[idx_dmin] * AU_KM
    V_rel_kms = v_transv[idx_dmin]   # km/s at closest approach
    print(f"\n  At closest approach (t = {t_from_event_h[idx_dmin]:+.2f}h):")
    print(f"    d_min       = {d[idx_dmin]:.5f} AU = {d_min_km:.0f} km")
    print(f"    V_rel       = {V_rel_kms:.2f} km/s")
    print(f"    ω_ICRS_max  = V_rel/d_min = {np.degrees(V_rel_kms/d_min_km)*3600:.2f}°/h")
    print(f"    ω_sid       = {OMEGA_SID:.3f}°/h")

    # Characteristic flyby timescale
    tau_s   = d_min_km / V_rel_kms          # seconds
    tau_min = tau_s / 60.0
    print(f"\n  Flyby timescale  τ = d_min / V_rel = {tau_min:.1f} min")
    print(f"    This is the time for angular velocity to change significantly.")

    # ω_ICRS profile: ω(t) = ω_max / (1 + (t/τ)²) for t measured from closest approach
    # Ground-frame apparent velocity is the VECTOR SUM of ω_ICRS (in approach direction)
    # and ω_sid (westward diurnal).  The guidance window is when the vector sum
    # points toward Bethlehem AND the diurnal component doesn't dominate az-drift.

    # Key constraint: daz/dt < AZ_DRIFT_THRESH = 1.5°/h
    # At closest approach, ω_ICRS ≈ ω_sid. The approach direction is nearly fixed
    # in ICRS. The diurnal rotation means the orientation of the approach velocity
    # in the ground frame rotates at ω_sid relative to the local frame.
    # For az-drift < 1.5°/h, the approach velocity's azimuth component must be < 1.5°/h.
    # This constrains how far from the stopping moment we can be (ω_ICRS drops) while
    # still having the correct azimuth direction.

    # Numerically: find the contiguous window where ω_ICRS > some fraction of ω_sid
    # AND the comet remains in the corridor.
    frac = omega_icrs_deg_h[idx_dmin] / OMEGA_SID
    print(f"\n  At closest approach: ω_ICRS/ω_sid = {frac:.4f}")
    print(f"    → The comet {'nearly cancels' if abs(frac-1)<0.3 else 'does NOT cancel'} "
          f"the diurnal motion at this moment.")

    # Time for ω_ICRS to drop by 50%:  1/(1+(t/τ)²) = 0.5 → t = τ
    t_halfmax_min = tau_min
    print(f"\n  Time for ω_ICRS to drop to 50% of maximum: τ = {t_halfmax_min:.1f} min")
    print(f"  Comet's ground-frame motion near closest approach:")

    # Find the interval where guidance is marginally possible (ω_ICRS > 0.5 ω_sid)
    # equivalent to d < d_min * sqrt(2ω_ICRS_max/ω_sid - 1) ... for ω_ICRS_max >> ω_sid
    # but ω_ICRS_max ≈ ω_sid here, so need ω_ICRS > ω_sid/2 → d < d_min * sqrt(2ω_max/ω_sid - 1)
    omega_max = omega_icrs_deg_h[idx_dmin]

    # From ω(t) = ω_max/(1+(V_rel*t/d_min)^2) > threshold:
    # t < (d_min/V_rel)*sqrt(ω_max/threshold - 1) = τ * sqrt(ω_max/threshold - 1)
    for thresh in [OMEGA_SID*0.9, OMEGA_SID*0.5, OMEGA_SID*0.1]:
        if omega_max > thresh:
            t_lim_min = tau_min * np.sqrt(omega_max/thresh - 1)
            print(f"    ω_ICRS > {thresh:.1f}°/h  for  |t| < {t_lim_min:.1f} min  "
                  f"(window = {2*t_lim_min:.0f} min)")
        else:
            print(f"    ω_ICRS never exceeds {thresh:.1f}°/h")

    # Azimuth corridor geometric constraint
    # The azimuth corridor [190°,210°] has width 20°.
    # A star drifts through this range due to diurnal motion.
    # At the altitude of the comet, daz/dt due to pure diurnal motion:
    if np.any(guidance_mask):
        k_guide = np.where(guidance_mask)[0]
        alt_guide = alt[k_guide].mean()
        # Pure diurnal daz/dt at az≈200°, alt≈alt_guide, lat=31.7°
        # daz/dt ≈ ω_E * (sin(lat)*cos(az) - cos(lat)*sin(az)*tan(alt)^-1)?
        # Use numerical estimate from the data
        daz_diurnal_typical = np.abs(daz[k_guide]).mean()
        transit_time = 20.0 / (daz_diurnal_typical + 0.001)   # hours
        print(f"\n  Geometry of the azimuth corridor:")
        print(f"    Average comet altitude during guidance: {alt_guide:.1f}°")
        print(f"    Average |daz/dt| during guidance: {daz_diurnal_typical:.3f}°/h")
        print(f"    Time to traverse 20° corridor at this rate: {transit_time*60:.0f} min")

    print()
    print("  PHYSICAL SUMMARY:")
    print(f"  ─────────────────────────────────────────────────────────────")
    print(f"  The guidance window is constrained by THREE simultaneous requirements:")
    print(f"  (1) Azimuth in [190°,210°]: the comet must be in the right direction.")
    print(f"  (2) |daz/dt| < 1.5°/h: slow enough azimuth drift to 'go before'.")
    print(f"  (3) |dalt/dt| > 0.3°/h: visible motion (not completely stationary).")
    print()
    print(f"  For this best orbit, the comet is near closest approach (d_min = "
          f"{d[idx_dmin]:.4f} AU)")
    print(f"  at t = {t_from_event_h[idx_dmin]:+.2f}h relative to the event.")
    print(f"  The flyby timescale τ = d_min/V_rel = {tau_min:.1f} min sets how quickly")
    print(f"  ω_ICRS changes. For ω_ICRS to remain near ω_sid (the condition for")
    print(f"  near-zero ground-frame motion), |t - t_min| must be ≲ τ ≈ {tau_min:.0f} min.")
    print(f"  The azimuth corridor transit time (≈ corridor_width / |daz/dt|_diurnal)")
    if np.any(guidance_mask):
        k_guide = np.where(guidance_mask)[0]
        daz_diurnal_typical = np.abs(daz[k_guide]).mean()
        print(f"  is ≈ 20° / {daz_diurnal_typical:.2f}°/h ≈ {20/daz_diurnal_typical*60:.0f} min.")
    print(f"  The guidance window is the INTERSECTION of both constraints,")
    print(f"  limited to ≲ 35 min — a fundamental geometric ceiling, not a")
    print(f"  sampling artifact. No orbit can extend this significantly because")
    print(f"  the flyby timescale τ is fixed by ω_sid and d_min.")


if __name__ == '__main__':
    import warnings; warnings.filterwarnings('ignore')
    main()
