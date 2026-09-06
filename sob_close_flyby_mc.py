#!/usr/bin/env python3
"""
sob_close_flyby_mc.py — Monte Carlo survey of close-flyby heliocentric comets
Aaron Adair / 2026

Tests N_MC randomly sampled comet orbits for whether any can simultaneously:
  (A) Pass within D_FLYBY AU of Earth near the event date
  (B) Be visible at night from Bethlehem (Sun below -6°, comet alt > 10°)
  (C) Appear to GUIDE toward Bethlehem: azimuth maintained near 190–210°
      with low azimuth drift (|daz/dt| < AZ_DRIFT_THRESH) for ≥ GUIDE_DUR hours
  (D) Come to apparent rest (total ω_app < STOP_THRESH °/h) within
      STOP_GAP hours after the guidance window

Criteria are the strict ground-frame (literal) interpretation of Matt 2:9.
Guidance = object appears nearly fixed in the direction of Bethlehem.
Stopping = object appears nearly motionless in both azimuth and altitude.

Key insight: for a close flyby, ω_ICRS = V_perp × d_min / d(t)^2, which is
maximised at closest approach. Stopping requires ω_ICRS ≈ ω_sid there.
Guidance requires a DIFFERENT ω_ICRS (just enough to cancel azimuth drift)
at some earlier time. The survey tests whether both can be satisfied
simultaneously for any physically plausible orbit.
"""

import numpy as np
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation, AltAz, SkyCoord, get_sun, get_body_barycentric
)
import astropy.units as u

# ── Physical constants ────────────────────────────────────────────────────────
K_GAUSS   = 0.01720209895
MU        = K_GAUSS**2            # GM_sun [AU³/day²]
EPS       = np.radians(23.4392911)  # J2000 obliquity
AU_KM     = 1.495978707e8
R_EQ = np.array([                  # ecliptic → equatorial J2000
    [1,  0,           0          ],
    [0,  np.cos(EPS), -np.sin(EPS)],
    [0,  np.sin(EPS),  np.cos(EPS)],
])

# ── Precession: J2000 equatorial → equatorial of date ────────────────────────
# Comet and Earth positions are built in the J2000 frame, but sidereal time
# (gmst_rad below) is referred to the equinox OF DATE.  Combining the two
# directly leaves the ~28° of general precession accumulated since 5 BCE
# uncorrected, which enters as a spurious rotation of ~23° in azimuth.
# Positions are therefore precessed to the equinox of date before any hour
# angle is formed.
#
# Composition follows Meeus, Astronomical Algorithms ch. 21:
#     P = Rz(z) · Ry(−θ) · Rz(ζ)
# Note this is NOT the Rz(−z)Ry(θ)Rz(−ζ) form quoted in many references,
# which assumes rotation of the axes rather than of the vector and yields the
# inverse rotation.  Because ζ, z, θ are all negative for epochs before J2000,
# using the wrong form does not produce an obvious sign flip — it produces a
# doubled offset (~60° in RA at 5 BCE).  Validated in sob_frames.validate().

def _rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

def _ry(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])

def precession_matrix(jd_tt):
    """Mean J2000 equator/equinox → mean equator/equinox of date."""
    T = (jd_tt - 2451545.0) / 36525.0
    zeta  = (2306.2181*T + 0.30188*T**2 + 0.017998*T**3) / 3600.0 * np.pi/180.0
    z     = (2306.2181*T + 1.09468*T**2 + 0.018203*T**3) / 3600.0 * np.pi/180.0
    theta = (2004.3109*T - 0.42665*T**2 - 0.041833*T**3) / 3600.0 * np.pi/180.0
    return _rz(z) @ _ry(-theta) @ _rz(zeta)

# ── Site ─────────────────────────────────────────────────────────────────────
LAT_DEG = 31.70;  LON_DEG = 35.20;  ALT_M = 765.0
LAT     = np.radians(LAT_DEG)
BETHLEHEM = EarthLocation(lon=LON_DEG*u.deg, lat=LAT_DEG*u.deg, height=ALT_M*u.m)
DT_TDB_UTC = 10572.0              # TDB - UTC in 5 BCE [s]

# ── Event dates (JD TDB) ─────────────────────────────────────────────────────
# Primary = Matney's 5 BCE Jun 8 local noon. Perihelion ± 90 days sampled.
JD_PRIMARY = 1719755.898

# ── Survey parameters ─────────────────────────────────────────────────────────
import os
N_MC        = int(os.environ.get("SOB_N_MC", 10_000_000))  # Monte Carlo samples
D_REJECT    = 0.05        # AU — quick-reject threshold (generous)
D_FLYBY     = 0.020       # AU — close-flyby threshold (criterion A)
STEP_MIN    = 5           # minutes — timestep for detailed analysis
WINDOW_HRS  = 13          # hours of night window to search (centred on midnight)

# Criteria thresholds
AZ_LO           = 190.0   # °  — guidance azimuth corridor low
AZ_HI           = 210.0   # °  — guidance azimuth corridor high
ALT_MIN         = 10.0    # °  — minimum comet altitude
AZ_DRIFT_THRESH = 1.5     # °/h — max |daz/dt| during guidance
ALT_RATE_MIN    = 0.3     # °/h — min |dalt/dt| during guidance (must be moving)
GUIDE_DUR       = 1.0     # h   — minimum continuous guidance duration
STOP_THRESH     = 0.30    # °/h — max total angular velocity for "stopping"
STOP_GAP        = 4.0     # h   — max time from guidance end to stopping
SUN_ALT_MAX     = -6.0    # °   — astronomical twilight limit

# ── Orbital mechanics (vectorised) ───────────────────────────────────────────

def _barker_vec(dt, q):
    """Barker equation for parabolic orbit → tan(ν/2). dt in days."""
    W = 3.0 * np.sqrt(MU / (2.0 * q**3)) * dt
    s = np.sqrt((W / 2.0)**2 + 1.0)
    return np.cbrt(W / 2.0 + s) + np.cbrt(W / 2.0 - s)

def _kepler_elliptic(M, e, tol=1e-11, nmax=50):
    """Newton–Raphson for M = E − e sin E. Vectorised."""
    E = M.copy()
    for _ in range(nmax):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))
        E += dE
        if np.max(np.abs(dE)) < tol:
            break
    return E

def _kepler_hyperbolic(Mh, e, tol=1e-11, nmax=50):
    """Newton–Raphson for Mh = e sinh H − H. Vectorised."""
    H = np.arcsinh(Mh / e)
    for _ in range(nmax):
        dH = (Mh - e * np.sinh(H) + H) / (e * np.cosh(H) - 1.0)
        H += dH
        if np.max(np.abs(dH)) < tol:
            break
    return H

def comet_pos_scalar(dt_days, q, e, i_r, Om_r, om_r):
    """
    Heliocentric equatorial (ICRS J2000) position of comet [AU].
    dt_days: scalar or 1-D array (days from perihelion).
    Returns shape (3,) or (3, N).
    """
    dt = np.atleast_1d(np.asarray(dt_days, float))
    ci, si = np.cos(i_r),  np.sin(i_r)
    cO, sO = np.cos(Om_r), np.sin(Om_r)
    co, so = np.cos(om_r), np.sin(om_r)

    # Perifocal unit vectors in ecliptic frame → equatorial
    Pecl = np.array([cO*co - sO*so*ci,  sO*co + cO*so*ci,  so*si])
    Qecl = np.array([-cO*so - sO*co*ci, -sO*so + cO*co*ci, co*si])
    P = R_EQ @ Pecl   # shape (3,)
    Q = R_EQ @ Qecl

    if abs(e - 1.0) < 1e-7:
        tanv2 = _barker_vec(dt, q)
        nu = 2.0 * np.arctan(tanv2)
        r  = 2.0 * q / (1.0 + np.cos(nu))
    elif e < 1.0:
        a  = q / (1.0 - e)
        n  = np.sqrt(MU / a**3)
        M  = (n * dt) % (2 * np.pi)
        E  = _kepler_elliptic(M, e)
        nu = 2.0 * np.arctan2(np.sqrt(1+e) * np.sin(E/2), np.sqrt(1-e) * np.cos(E/2))
        r  = a * (1.0 - e * np.cos(E))
    else:  # hyperbolic
        a  = q / (1.0 - e)           # a < 0
        n  = np.sqrt(MU / (-a)**3)
        Mh = n * dt
        H  = _kepler_hyperbolic(Mh, e)
        nu = 2.0 * np.arctan2(np.sqrt(e+1) * np.sinh(H/2),
                               np.sqrt(e-1) * np.cosh(H/2))
        r  = a * (1.0 - e * np.cosh(H))   # > 0

    pos = np.outer(P, r * np.cos(nu)) + np.outer(Q, r * np.sin(nu))  # (3, N)
    return pos.squeeze() if dt.size == 1 else pos


def comet_pos_parabolic_batch(dt_arr, q_arr, i_arr, Om_arr, om_arr):
    """
    Parabolic approximation, fully vectorised over N orbits evaluated at
    one dt each.  dt_arr, q_arr, ... all shape (N,). Returns (3, N).
    """
    tanv2 = _barker_vec(dt_arr, q_arr)
    nu    = 2.0 * np.arctan(tanv2)
    r     = 2.0 * q_arr / (1.0 + np.cos(nu))

    ci, si = np.cos(i_arr),  np.sin(i_arr)
    cO, sO = np.cos(Om_arr), np.sin(Om_arr)
    co, so = np.cos(om_arr), np.sin(om_arr)

    Px = cO*co - sO*so*ci;  Py = sO*co + cO*so*ci;  Pz = so*si
    Qx = -cO*so - sO*co*ci; Qy = -sO*so + cO*co*ci; Qz = co*si

    cnu, snu = np.cos(nu), np.sin(nu)
    xe = r*(cnu*Px + snu*Qx)
    ye = r*(cnu*Py + snu*Qy)
    ze = r*(cnu*Pz + snu*Qz)

    # Ecliptic → equatorial
    x = xe
    y = np.cos(EPS)*ye - np.sin(EPS)*ze
    z = np.sin(EPS)*ye + np.cos(EPS)*ze
    return np.array([x, y, z])   # (3, N)


# ── AltAz helpers (analytic, no astropy dependency in inner loop) ─────────────

def altaz_analytic(ra_r, dec_r, lst_r, lat=LAT):
    """Return (alt°, az°) for arrays of RA, Dec, LST (radians)."""
    ha       = lst_r - ra_r
    sin_alt  = np.sin(lat)*np.sin(dec_r) + np.cos(lat)*np.cos(dec_r)*np.cos(ha)
    sin_alt  = np.clip(sin_alt, -1.0, 1.0)
    alt      = np.arcsin(sin_alt)
    cos_alt  = np.cos(alt)
    # Azimuth: 0 = N, 90 = E (standard astronomical convention)
    az_sin   = -np.cos(dec_r) * np.sin(ha)
    az_cos   = (np.sin(dec_r) - np.sin(lat)*sin_alt) / (np.cos(lat) * np.maximum(cos_alt, 1e-9))
    az       = np.arctan2(az_sin, az_cos) % (2*np.pi)
    return np.degrees(alt), np.degrees(az)


def gmst_rad(jd_utc):
    """Approximate GMST [radians] for array of JD UTC."""
    T = (jd_utc - 2451545.0) / 36525.0
    gmst_deg = (280.46061837 + 360.98564736629*(jd_utc - 2451545.0)
                + 0.000387933*T**2 - T**3/38710000.0) % 360.0
    return np.radians(gmst_deg)

def lst_rad(jd_utc, lon_deg=LON_DEG):
    return (gmst_rad(jd_utc) + np.radians(lon_deg)) % (2*np.pi)


# ── Precompute Earth and Sun positions ────────────────────────────────────────

def precompute_sky(jd_event_tdb, window_hrs=WINDOW_HRS, step_min=STEP_MIN):
    """
    Build arrays over a night window centred on local midnight nearest the event.
    Returns: jd_utc[N], earth_pos[3,N], sun_ra[N], sun_dec[N], lst[N]
    """
    # Local midnight: event noon - 0.5 day
    jd_utc_midnight = (jd_event_tdb - DT_TDB_UTC/86400.0) - 0.5
    dt_h = np.arange(-window_hrs/2, window_hrs/2, step_min/60.0)
    jd_utc = jd_utc_midnight + dt_h / 24.0
    jd_tdb = jd_utc + DT_TDB_UTC / 86400.0

    # Earth from VSOP87, not astropy's bundled ephemeris.  ERFA's epv00 is only
    # claimed valid 1900-2100 and is off by ~43" at 5 BCE.  Because this survey
    # tests near-stationary configurations of objects a few hundred thousand km
    # away, a 43" error in Earth's position displaces the observer by ~3e4 km
    # and swings the sightline by degrees -- enough to create or destroy an
    # apparent standstill.  Substituting the bundled ephemeris here produced a
    # spurious "solution" in the companion optimiser.
    import sob_frames as _sf
    earth_ecl = _sf.earth_helio_ecl_j2000(jd_tdb)
    earth_pos = _sf.ecl_to_equ(earth_ecl)          # (3, N) J2000 equatorial

    sun_vec  = _sf.ecl_to_equ(-earth_ecl)
    sun_n    = np.linalg.norm(sun_vec, axis=0)
    sun_ra   = np.arctan2(sun_vec[1], sun_vec[0]) % (2*np.pi)
    sun_dec  = np.arcsin(np.clip(sun_vec[2]/sun_n, -1, 1))

    lst = lst_rad(jd_utc)
    return jd_utc, earth_pos, sun_ra, sun_dec, lst


# ── Single-orbit detailed evaluation ─────────────────────────────────────────

def evaluate_orbit(q, e, i_r, Om_r, om_r, T_jd,
                   jd_utc, earth_pos, sun_ra, sun_dec, lst):
    """
    Evaluate one orbit against all Matt 2:9 criteria.
    Returns dict with d_min, pass flag, and failure reason.
    """
    step_h = (jd_utc[1] - jd_utc[0]) * 24.0

    jd_tdb = jd_utc + DT_TDB_UTC / 86400.0
    dt     = jd_tdb - T_jd

    try:
        r_c = comet_pos_scalar(dt, q, e, i_r, Om_r, om_r)   # (3, N)
    except Exception:
        return {'pass': False, 'reason': 'orbit_error', 'd_min': np.inf}

    rho  = r_c - earth_pos       # (3, N)
    d    = np.linalg.norm(rho, axis=0)
    d_min = float(np.min(d))

    if d_min > D_FLYBY:
        return {'pass': False, 'reason': f'd_min={d_min:.4f}>{D_FLYBY}', 'd_min': d_min}

    # RA, Dec of comet — precessed from J2000 to the equinox of date so that
    # the hour angle formed against `lst` (which is of-date) is frame-consistent.
    rhat = rho / np.maximum(d, 1e-12)
    P    = precession_matrix(float(np.mean(jd_tdb)))
    rhat = P @ rhat
    dec  = np.degrees(np.arcsin(np.clip(rhat[2], -1, 1)))
    ra   = np.degrees(np.arctan2(rhat[1], rhat[0])) % 360.0

    alt, az = altaz_analytic(np.radians(ra), np.radians(dec), lst)

    # Sun altitude (analytic) — precessed to the same frame
    sun_v = np.array([np.cos(sun_dec)*np.cos(sun_ra),
                      np.cos(sun_dec)*np.sin(sun_ra),
                      np.sin(sun_dec)])
    sun_v = P @ sun_v
    sun_ra_d  = np.arctan2(sun_v[1], sun_v[0]) % (2*np.pi)
    sun_dec_d = np.arcsin(np.clip(sun_v[2], -1, 1))
    sun_alt, _ = altaz_analytic(sun_ra_d, sun_dec_d, lst)

    # Apparent angular velocity in ground frame (finite differences)
    dt_h  = np.gradient(jd_utc) * 24.0          # hours per step (varies at edges)
    dalt  = np.gradient(alt) / dt_h              # °/h
    # Unwrap azimuth before differencing
    az_u  = np.degrees(np.unwrap(np.radians(az)))
    daz   = np.gradient(az_u) / dt_h             # °/h

    # Total apparent angular velocity (in AltAz metric)
    cos_alt = np.cos(np.radians(alt))
    omega_app = np.sqrt(dalt**2 + (daz * cos_alt)**2)   # °/h

    # ── Masks ────────────────────────────────────────────────────────────────
    night     = sun_alt < SUN_ALT_MAX
    high      = alt > ALT_MIN
    corridor  = (az >= AZ_LO) & (az <= AZ_HI)

    # Guidance mask: night, high, in corridor, azimuth nearly fixed,
    # and object actually moving (not already stopped)
    guide_mask = (night & high & corridor
                  & (np.abs(daz) < AZ_DRIFT_THRESH)
                  & (np.abs(dalt) > ALT_RATE_MIN))

    if not np.any(guide_mask):
        return {'pass': False, 'reason': 'no_guidance_window', 'd_min': d_min}

    # Require guidance duration ≥ GUIDE_DUR hours
    guide_dur = float(np.sum(guide_mask)) * step_h
    if guide_dur < GUIDE_DUR:
        return {'pass': False,
                'reason': f'guide_dur={guide_dur:.2f}h<{GUIDE_DUR}h', 'd_min': d_min}

    # End of last guidance moment
    guide_indices = np.where(guide_mask)[0]
    guide_end_idx = int(guide_indices[-1])

    # Look for stopping within STOP_GAP hours after guidance ends
    n_gap     = int(np.ceil(STOP_GAP / step_h))
    stop_end  = min(guide_end_idx + n_gap + 1, len(omega_app))
    stop_slice = slice(guide_end_idx, stop_end)

    # Must be nighttime during stop window
    stop_night = night[stop_slice]
    if not np.any(stop_night):
        return {'pass': False, 'reason': 'stop_in_daylight', 'd_min': d_min}

    omega_stop_region = omega_app[stop_slice][stop_night]
    min_stop = float(np.min(omega_stop_region))

    if min_stop > STOP_THRESH:
        return {'pass': False,
                'reason': f'min_ω={min_stop:.3f}>stop_thresh', 'd_min': d_min}

    # Passed all criteria
    gi = guide_indices[0]
    return {
        'pass':         True,
        'd_min':        d_min,
        'guide_az':     float(az[gi]),
        'guide_alt':    float(alt[gi]),
        'guide_dur_h':  guide_dur,
        'min_stop_omega': min_stop,
        'e': e, 'q': q, 'i_deg': np.degrees(i_r),
        'Om_deg': np.degrees(Om_r), 'om_deg': np.degrees(om_r),
    }


# ── Main MC survey ────────────────────────────────────────────────────────────

def run_survey():
    print("=" * 70)
    print("CLOSE-FLYBY MONTE CARLO SURVEY  —  Matt 2:9 strict ground-frame")
    print(f"N_MC = {N_MC:,}  |  d_flyby < {D_FLYBY} AU  |  step = {STEP_MIN} min")
    print("=" * 70)

    print("\n[1] Precomputing Earth/Sun over night window...")
    jd_utc, earth_pos, sun_ra, sun_dec, lst = precompute_sky(JD_PRIMARY)
    step_h = (jd_utc[1] - jd_utc[0]) * 24.0
    print(f"    {len(jd_utc)} steps × {step_h*60:.0f} min  "
          f"({jd_utc[0]:.3f} – {jd_utc[-1]:.3f} JD UTC)")

    earth_mid = earth_pos[:, len(jd_utc)//2]   # Earth pos at midnight

    # ── Sample orbital elements ──────────────────────────────────────────────
    print(f"\n[2] Drawing {N_MC:,} random orbits...")
    rng = np.random.default_rng(2026)

    # ── Perihelion distance q ────────────────────────────────────────────────
    # Concentrate near Earth-crossing orbits (q ≈ 1 AU) where close flybys
    # are geometrically most likely.  Stratified:
    #   60 %  in [0.80, 1.20] AU  — optimal Earth-crossing zone
    #   25 %  in [0.50, 0.80] ∪ [1.20, 1.40] AU — mid range
    #   15 %  in [0.30, 0.50] AU  — inner comets (completeness)
    N_q0 = int(0.60 * N_MC)
    N_q1 = int(0.25 * N_MC)
    N_q2 = N_MC - N_q0 - N_q1
    q_arr = np.concatenate([
        rng.uniform(0.80, 1.20, N_q0),
        np.where(rng.integers(0, 2, N_q1).astype(bool),
                 rng.uniform(0.50, 0.80, N_q1),
                 rng.uniform(1.20, 1.40, N_q1)),
        rng.uniform(0.30, 0.50, N_q2),
    ])

    # ── Eccentricity: focus on cometary (high-e) orbits ──────────────────────
    e_arr  = np.concatenate([
        rng.uniform(0.70, 1.00, N_MC // 3),                         # high-elliptic
        rng.uniform(1.00, 1.01, N_MC // 6),                         # near-parabolic
        rng.uniform(1.01, 5.00, N_MC - N_MC // 3 - N_MC // 6),     # hyperbolic
    ])

    # ── Inclination: stratified toward low-i (near-ecliptic) orbits ──────────
    # First-principles argument: for a close flyby to cancel diurnal drift and
    # produce ground-frame guidance → stopping, the comet's geocentric velocity
    # vector must lie close to the celestial equatorial plane.  The optimizer
    # confirmed the best orbit has i = 2.2°.  Isotropic sampling wastes ~97% of
    # draws on high-i orbits that are geometrically excluded.  Stratification:
    #   60 %  in [0°, 20°]   — prime candidate zone
    #   25 %  in [20°, 45°]  — moderate inclination
    #   15 %  in [45°, 180°] — full-sky tail (probabilistic statement + completeness)
    N_i0 = int(0.60 * N_MC)
    N_i1 = int(0.25 * N_MC)
    N_i2 = N_MC - N_i0 - N_i1
    # Convert degree bounds to cos(i) for uniform-in-angle draws
    i_low  = np.arccos(rng.uniform(np.cos(np.radians(20)),  1.0,                   N_i0))
    i_mid  = np.arccos(rng.uniform(np.cos(np.radians(45)),  np.cos(np.radians(20)), N_i1))
    i_high = np.arccos(rng.uniform(-1.0,                    np.cos(np.radians(45)), N_i2))
    i_arr  = np.concatenate([i_low, i_mid, i_high])

    Om_arr = rng.uniform(0.0, 2*np.pi, N_MC)
    om_arr = rng.uniform(0.0, 2*np.pi, N_MC)
    # Perihelion time within ±90 days of event
    T_arr  = JD_PRIMARY + rng.uniform(-90.0, 90.0, N_MC)

    # ── Vectorised quick reject ──────────────────────────────────────────────
    print("[3] Quick-reject phase (parabolic approximation)...")
    dt_mid = JD_PRIMARY - T_arr         # dt at event time (TDB)
    pos_approx = comet_pos_parabolic_batch(dt_mid, q_arr, i_arr, Om_arr, om_arr)
    dx = pos_approx[0] - earth_mid[0]
    dy = pos_approx[1] - earth_mid[1]
    dz = pos_approx[2] - earth_mid[2]
    d_approx = np.sqrt(dx**2 + dy**2 + dz**2)

    keep  = d_approx < D_REJECT
    n_keep = int(np.sum(keep))
    print(f"    Survived: {n_keep:,} / {N_MC:,} ({100*n_keep/N_MC:.4f}%)")

    if n_keep == 0:
        print("\n  No close-flyby candidates even after generous reject. Done.")
        return []

    # ── Detailed analysis ────────────────────────────────────────────────────
    print(f"\n[4] Detailed evaluation of {n_keep:,} candidates...")
    q_k  = q_arr[keep];  e_k  = e_arr[keep]
    i_k  = i_arr[keep];  Om_k = Om_arr[keep]
    om_k = om_arr[keep]; T_k  = T_arr[keep]

    n_flyby = 0
    n_pass  = 0
    failures = {}
    passes   = []

    for j in range(n_keep):
        res = evaluate_orbit(q_k[j], e_k[j], i_k[j], Om_k[j], om_k[j], T_k[j],
                             jd_utc, earth_pos, sun_ra, sun_dec, lst)
        if res['d_min'] <= D_FLYBY:
            n_flyby += 1
            reason = res.get('reason', 'PASS')
            failures[reason] = failures.get(reason, 0) + 1
            if res['pass']:
                n_pass += 1
                passes.append(res)

        if (j + 1) % 5000 == 0:
            pct = 100*(j+1)/n_keep
            print(f"    {j+1:6d}/{n_keep}  ({pct:.0f}%)  "
                  f"flyby={n_flyby}  pass={n_pass}", end='\r')

    print()

    # ── Results ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"MC samples drawn        : {N_MC:,}")
    print(f"Survived quick reject   : {n_keep:,}")
    print(f"True close-flyby orbits : {n_flyby:,}  (d_min < {D_FLYBY} AU)")
    print(f"Passed ALL criteria     : {n_pass}")

    if n_pass > 0:
        print("\n*** VIABLE CLOSE-FLYBY ORBITS FOUND ***")
        for r in passes:
            print(f"  e={r['e']:.3f}  q={r['q']:.4f}AU  i={r['i_deg']:.1f}°  "
                  f"Ω={r['Om_deg']:.1f}°  ω={r['om_deg']:.1f}°  "
                  f"d_min={r['d_min']:.4f}AU  "
                  f"az={r['guide_az']:.1f}°  guide={r['guide_dur_h']:.1f}h  "
                  f"ω_stop={r['min_stop_omega']:.3f}°/h")
    else:
        print("\nZero orbits satisfy guidance→stopping in nighttime corridor.")

    if failures:
        print(f"\nFailure breakdown (among {n_flyby} close-flyby orbits):")
        for reason, count in sorted(failures.items(), key=lambda x: -x[1]):
            label = reason if reason != 'PASS' else 'pass'
            print(f"  {count:6d}  {label}")

    return passes


if __name__ == '__main__':
    results = run_survey()
