#!/usr/bin/env python3
"""
sob_mc_refined.py
=================
Monte Carlo refinement around best-case orbital candidates found by the
exhaustive grid survey (sob_proof_fine.py).

BACKGROUND
----------
The grid survey (sob_proof_fine.py) sampled 529,079,040 heliocentric
configurations at 5° resolution in Ω and ω and found zero satisfying all
three criteria simultaneously.  The closest approach to feasibility was
the orbit with the best guidance score (0.015° inside the road corridor):

    e = 1.1, q = 0.03 AU, i = 10°

However, that orbit's angular velocity (0.23°/h) falls an order of
magnitude below the required 2°/h guidance motion.

A critic might ask: could fine-tuning around these parameters, or around
Matney's (2025) specifically proposed comet, reveal a viable orbit that
the 5° grid spacing missed?  This script answers that question by densely
sampling continuous parameter space around two candidate regions.

TWO SWEEPS
----------
SWEEP 1 — Best-fit guidance orbit neighbourhood
  Broad uniform random search centred on e=1.1, q=0.03 AU, i=10°,
  with Ω and ω drawn uniformly over all 360°.
  Purpose: exhaustively probe the vicinity of the most promising grid
  orbit, ruling out that any gap in the 5° grid hides a viable solution.

SWEEP 2 — Matney (2025) proposed 5 BCE comet neighbourhood
  Published 5 BCE osculating elements (Matney, JBAA 135, 2025):
    q = 0.0407 AU,  e = 1.0,   i = 1.43°
    ω₅ᴮᶜᴱ = 335.58°,  Ω₅ᴮᶜᴱ = 75.81°,  T_peri = JD 1719785.565
    → ΔT ≈ +29.7 days from reference epoch
  Purpose: confirm that even the one specifically proposed orbit, and its
  generous neighbourhood, contains no viable candidate.

CRITERIA (identical to sob_proof_fine.py)
-----------------------------------------
  Guidance window : 08:00–10:00 local solar time (JD_NOON = 1719755.898)
  Stopping window : 10:15–12:00 local solar time
  (A) Azimuth in road corridor [190°, 210°] ± 5° throughout guidance window
  (B) Angular velocity > 2°/h throughout guidance window
  (C) Angular velocity < 2°/h throughout stopping window
  Object must also be above 5° altitude throughout the guidance window.

USAGE
-----
  python sob_mc_refined.py [--n1 N] [--n2 N] [--batch B] [--seed S]

  --n1 N    Samples for sweep 1 (default: 10,000,000)
  --n2 N    Samples for sweep 2 (default:  5,000,000)
  --batch B Batch size for vectorised evaluation (default: 200,000)
  --seed S  Random seed for reproducibility (default: 42)

OUTPUT
------
  sob_mc_refined.log      — timestamped progress log
  sob_mc_results.json     — full statistics and best-found orbit parameters
"""

import argparse
import json
import os
import sys
import time as walltime

import numpy as np
from astropy.time import Time
from astropy.coordinates import (
    SkyCoord, AltAz, EarthLocation, get_body_barycentric
)
import astropy.units as u

# ── Physical constants ────────────────────────────────────────────────────
K_GAUSS = 0.01720209895      # Gaussian gravitational constant
MU      = K_GAUSS**2         # AU³ day⁻²
RAD     = 180.0 / np.pi      # radians → degrees
DEG     = np.pi / 180.0      # degrees → radians

# ── Observer ──────────────────────────────────────────────────────────────
BETHLEHEM = EarthLocation(
    lon=(35 + 12/60)*u.deg,
    lat=(31 + 42/60)*u.deg,
    height=765.0*u.m
)

# ── Reference epoch ───────────────────────────────────────────────────────
JD_NOON = 1_719_755.898   # local noon at Bethlehem, 5 BCE (TDB)
DT_SEC  = 10_572.0        # ΔT = UTC − TDB for 5 BCE (s)

# ── Criteria ──────────────────────────────────────────────────────────────
AZ_LO, AZ_HI = 190.0, 210.0
AZ_TOL       = 5.0    # ° — max deviation outside corridor
OMEGA_THRESH = 2.0    # °/h
ALT_FLOOR    = 5.0    # ° — minimum altitude for visibility

# ── Matney (2025) central parameters ─────────────────────────────────────
MATNEY = dict(
    q   = 0.0407,              # AU
    e   = 1.0,
    i   = 1.43,                # degrees (5 BCE epoch)
    om  = 335.58,              # ω degrees (5 BCE epoch)
    asc = 75.81,               # Ω degrees (5 BCE epoch)
    jd  = 1_719_785.565,       # T_perihelion (TDB)
)
MATNEY['dT'] = MATNEY['jd'] - JD_NOON   # ≈ +29.67 days

# ── Output files ──────────────────────────────────────────────────────────
_DIR     = os.path.dirname(os.path.abspath(__file__))
LOGFILE  = os.path.join(_DIR, 'sob_mc_refined.log')
OUTFILE  = os.path.join(_DIR, 'sob_mc_results.json')
T0       = walltime.time()

open(LOGFILE, 'w').close()   # reset log

def log(msg):
    line = f"[{walltime.time()-T0:8.1f}s] {msg}"
    print(line, flush=True)
    with open(LOGFILE, 'a') as f:
        f.write(line + '\n')


# ─────────────────────────────────────────────────────────────────────────
# Pre-computation: ENU matrices and Earth positions
# ─────────────────────────────────────────────────────────────────────────

def precompute_environment():
    """
    Build ENU rotation matrices and heliocentric Earth positions at each of
    the 17 guidance+stopping time steps (same grid as sob_proof_fine.py).
    Returns
    -------
    enu       : (N_GS, 3, 3)  — rows: [East, North, Up] unit vectors in ICRS
    earth_pos : (N_GS, 3)     — heliocentric Earth position in AU
    jds_gs    : (N_GS,)       — Julian dates (TDB) of each step
    N_G       : int            — number of guidance steps
    """
    log("Pre-computing ENU matrices and Earth positions …")

    hrs_all = np.arange(-8.0, 3.1, 0.25)       # hours from local noon
    local_h = hrs_all + 12.0                    # local solar time (hours)
    jds_tdb = JD_NOON + hrs_all / 24.0
    jds_utc = jds_tdb  - DT_SEC / 86400.0

    mask_g = (local_h >= 8.0)  & (local_h <= 10.0)   # guidance 08:00–10:00
    mask_s = (local_h >  10.0) & (local_h <= 12.0)   # stopping 10:15–12:00
    N_G = int(mask_g.sum())    # 9 steps
    N_S = int(mask_s.sum())    # 8 steps

    jds_gs  = np.concatenate([jds_tdb[mask_g], jds_tdb[mask_s]])   # (17,)
    utc_gs  = np.concatenate([jds_utc[mask_g], jds_utc[mask_s]])
    N_GS    = N_G + N_S

    t_tdb = Time(jds_gs, format='jd', scale='tdb')
    t_utc = Time(utc_gs, format='jd', scale='utc')

    # Heliocentric Earth positions (Earth – Sun, AU)
    Ef = get_body_barycentric('earth', t_tdb)
    Sf = get_body_barycentric('sun',   t_tdb)
    earth_pos = (Ef - Sf).xyz.to(u.AU).value.T    # (N_GS, 3)

    # ENU matrices
    def aa2icrs(t_utc_single, az_deg, alt_deg):
        frame = AltAz(obstime=t_utc_single, location=BETHLEHEM)
        c = SkyCoord(alt=alt_deg*u.deg, az=az_deg*u.deg, frame=frame).icrs
        d, r = np.radians(c.dec.deg), np.radians(c.ra.deg)
        return np.array([np.cos(d)*np.cos(r),
                         np.cos(d)*np.sin(r),
                         np.sin(d)])

    enu = np.zeros((N_GS, 3, 3))
    for k in range(N_GS):
        E = aa2icrs(t_utc[k], 90.0,  0.001)
        Z = aa2icrs(t_utc[k],  0.0, 90.0)
        Z /= np.linalg.norm(Z)
        E -= np.dot(E, Z) * Z;  E /= np.linalg.norm(E)
        enu[k] = [E, np.cross(Z, E), Z]

    log(f"  {N_G} guidance + {N_S} stopping steps.  "
        f"Earth x-range [{earth_pos[:,0].min():.3f}, {earth_pos[:,0].max():.3f}] AU")
    return enu, earth_pos, jds_gs, N_G


# ─────────────────────────────────────────────────────────────────────────
# Vectorised orbital mechanics
# ─────────────────────────────────────────────────────────────────────────

def _true_anomaly(dt, q, e):
    """
    Solve Kepler's equation for all three conic types in one vectorised call.

    Parameters
    ----------
    dt, q, e : 1-D arrays of length N — time from perihelion (days),
               perihelion distance (AU), eccentricity.

    Returns
    -------
    nu : 1-D array of length N — true anomaly in radians.
         NaN where the orbit is invalid (object past asymptote, or r < 0).
    """
    nu = np.full(len(dt), np.nan, dtype=np.float64)

    # ── Elliptic (e < 0.9999) ────────────────────────────────────────────
    m = e < 0.9999
    if m.any():
        eq, ee, edt = q[m], e[m], dt[m]
        a = eq / (1.0 - ee)
        n = np.sqrt(MU / a**3)
        M = (n * edt) % (2.0 * np.pi)
        E = M.copy()
        for _ in range(60):
            dE = (M - E + ee*np.sin(E)) / (1.0 - ee*np.cos(E))
            E += dE
            if np.max(np.abs(dE)) < 1e-12:
                break
        nu[m] = 2.0 * np.arctan2(
            np.sqrt(1.0 + ee) * np.sin(E / 2.0),
            np.sqrt(1.0 - ee) * np.cos(E / 2.0)
        )

    # ── Parabolic (0.9999 ≤ e ≤ 1.0001) ─────────────────────────────────
    m = (e >= 0.9999) & (e <= 1.0001)
    if m.any():
        pq, pdt = q[m], dt[m]
        W   = 3.0 * K_GAUSS * pdt / np.sqrt(2.0 * pq**3)
        arg = W / 2.0 + np.sqrt(np.maximum((W / 2.0)**2 + 1.0, 0.0))
        arg = np.where(arg > 0, arg, 1e-30)
        Y   = arg ** (1.0 / 3.0)
        nu[m] = 2.0 * np.arctan(Y - 1.0 / Y)

    # ── Hyperbolic (e > 1.0001) ──────────────────────────────────────────
    m = e > 1.0001
    if m.any():
        hq, he, hdt = q[m], e[m], dt[m]
        a   = hq / (he - 1.0)
        n   = np.sqrt(MU / a**3)
        Mh  = n * hdt
        # Initial guess for hyperbolic anomaly H
        H   = np.arcsinh(np.clip(Mh / he, -1e6, 1e6))
        H   = np.where(np.isfinite(H), H, np.sign(Mh) * np.log(np.abs(Mh) + 1.8))
        for _ in range(60):
            f  = he * np.sinh(H) - H - Mh
            fp = he * np.cosh(H) - 1.0
            dH = -f / np.where(np.abs(fp) > 1e-15, fp, 1e-15)
            dH = np.clip(dH, -2.0, 2.0)
            H  += dH
            if np.max(np.abs(dH)) < 1e-12:
                break
        nu_h    = 2.0 * np.arctan2(
            np.sqrt(he + 1.0) * np.sinh(H / 2.0),
            np.sqrt(he - 1.0) * np.cosh(H / 2.0)
        )
        nu_max  = np.arccos(np.clip(-1.0 / he, -1.0, 1.0))
        nu[m]   = np.where(np.abs(nu_h) < nu_max, nu_h, np.nan)

    return nu


def _pos_icrs(nu, r, Om, om, inc):
    """
    Transform perifocal coordinates to ICRS for a batch of N orbits.
    All inputs shape (N,); output shape (N, 3).
    """
    cnu = np.cos(nu);  snu = np.sin(nu)
    cO  = np.cos(Om);  sO  = np.sin(Om)
    co  = np.cos(om);  so  = np.sin(om)
    ci  = np.cos(inc); si  = np.sin(inc)
    # P-hat (perihelion direction)
    Px = cO*co - sO*so*ci;   Py = sO*co + cO*so*ci;   Pz = so*si
    # Q-hat (90° ahead of perihelion in orbital plane)
    Qx = -cO*so - sO*co*ci;  Qy = -sO*so + cO*co*ci;  Qz = co*si
    x = r * (cnu*Px + snu*Qx)
    y = r * (cnu*Py + snu*Qy)
    z = r * (cnu*Pz + snu*Qz)
    return np.stack([x, y, z], axis=-1)


# ─────────────────────────────────────────────────────────────────────────
# Batch evaluation
# ─────────────────────────────────────────────────────────────────────────

def evaluate_batch(q_, e_, i_r, Om_r, om_r, Tp,
                   enu, earth_pos, jds_gs, N_G, DT_H):
    """
    Evaluate B orbital configurations against all three criteria.

    Parameters
    ----------
    q_, e_           : (B,) — perihelion distance (AU), eccentricity
    i_r, Om_r, om_r  : (B,) — inclination, ascending node, arg. perihelion (rad)
    Tp               : (B,) — perihelion epoch (JD, TDB)
    enu              : (N_GS, 3, 3) — ENU rotation matrices
    earth_pos        : (N_GS, 3)   — heliocentric Earth positions (AU)
    jds_gs           : (N_GS,)     — Julian dates of evaluation steps
    N_G              : int          — number of guidance steps
    DT_H             : float        — time spacing between steps (hours)

    Returns (all shape (B,))
    -------
    vis_ok    : bool  — object above ALT_FLOOR throughout guidance window
    az_score  : float — max° outside [190°, 210°] during guidance (0 = in corridor)
    motion_min: float — min angular velocity °/h during guidance
    stop_max  : float — max angular velocity °/h during stopping
    all_ok    : bool  — all three criteria satisfied
    """
    B   = len(q_)
    N_GS = len(jds_gs)

    # dt[b, t] = time from perihelion (days) for sample b at step t
    dt_mat = jds_gs[np.newaxis, :] - Tp[:, np.newaxis]   # (B, N_GS)

    # Flatten to (B*N_GS,) for vectorised Kepler solution
    N    = B * N_GS
    t_i  = np.tile(np.arange(N_GS), B)   # time-step index for each flat entry
    q_f  = np.repeat(q_,  N_GS)
    e_f  = np.repeat(e_,  N_GS)
    i_f  = np.repeat(i_r, N_GS)
    Om_f = np.repeat(Om_r, N_GS)
    om_f = np.repeat(om_r, N_GS)
    dt_f = dt_mat.ravel()

    # True anomaly
    nu_f  = _true_anomaly(dt_f, q_f, e_f)

    # Heliocentric distance
    p_f   = q_f * (1.0 + e_f)
    r_f   = p_f / (1.0 + e_f * np.cos(nu_f))
    ok_f  = np.isfinite(nu_f) & (r_f > 0) & (r_f < 1e4)

    # Heliocentric ICRS position (use safe values to avoid NaN propagation)
    nu_s  = np.where(ok_f, nu_f, 0.0)
    r_s   = np.where(ok_f, r_f,  0.0)
    xyz_f = _pos_icrs(nu_s, r_s, Om_f, om_f, i_f)    # (N, 3)
    xyz_f[~ok_f] = 0.0

    # Geocentric unit vector
    geo_f  = xyz_f - earth_pos[t_i]                    # (N, 3)
    nrm    = np.linalg.norm(geo_f, axis=1, keepdims=True)
    geo_n  = geo_f / np.where(nrm > 0, nrm, 1.0)      # (N, 3)

    # AltAz via ENU rotation
    E_f = enu[t_i]                                     # (N, 3, 3)
    e_c = np.einsum('ni,ni->n', geo_n, E_f[:, 0, :])  # East component
    n_c = np.einsum('ni,ni->n', geo_n, E_f[:, 1, :])  # North component
    u_c = np.einsum('ni,ni->n', geo_n, E_f[:, 2, :])  # Up component

    alt_f = np.arcsin(np.clip(u_c, -1.0, 1.0)) * RAD
    az_f  = (np.arctan2(e_c, n_c) * RAD) % 360.0

    alt_m   = alt_f.reshape(B, N_GS)
    az_m    = az_f.reshape(B, N_GS)
    ok_m    = ok_f.reshape(B, N_GS)
    geo_nm  = geo_n.reshape(B, N_GS, 3)

    # ── Visibility: altitude > ALT_FLOOR for all guidance steps ──────────
    vis_ok = (np.all(alt_m[:, :N_G] > ALT_FLOOR, axis=1) &
              np.all(ok_m[:,  :N_G],              axis=1))

    # ── Azimuth score (guidance window) ──────────────────────────────────
    az_dev   = np.maximum(0.0,
                   np.maximum(AZ_LO - az_m[:, :N_G],
                               az_m[:, :N_G] - AZ_HI))
    az_score = np.max(az_dev, axis=1)    # 0 = inside corridor

    # ── Angular velocities via great-circle separations ──────────────────
    def _angvel(seg):
        """seg : (B, K, 3) consecutive unit vectors → (B, K-1) °/h"""
        d = np.einsum('bti,bti->bt', seg[:, :-1, :], seg[:, 1:, :])
        return np.arccos(np.clip(d, -1.0, 1.0)) * RAD / DT_H

    motion_min = np.min(_angvel(geo_nm[:, :N_G,  :]), axis=1)  # guidance
    stop_max   = np.max(_angvel(geo_nm[:,  N_G:, :]), axis=1)   # stopping

    # Null non-visible samples
    az_score   = np.where(vis_ok, az_score,   np.inf)
    motion_min = np.where(vis_ok, motion_min, np.nan)
    stop_max   = np.where(vis_ok, stop_max,   np.nan)

    all_ok = (vis_ok
              & (az_score   < AZ_TOL)
              & (motion_min > OMEGA_THRESH)
              & (stop_max   < OMEGA_THRESH))

    return vis_ok, az_score, motion_min, stop_max, all_ok


# ─────────────────────────────────────────────────────────────────────────
# Statistics tracker
# ─────────────────────────────────────────────────────────────────────────

class SweepStats:
    """Accumulates best-found results across batches in one MC sweep."""

    def __init__(self, name, desc):
        self.name = name
        self.desc = desc
        self.n_sampled = self.n_visible = self.n_all_ok = 0
        # Best azimuth score (lowest = most corridor-centred)
        self.best_az       = np.inf
        self.best_az_p     = None          # params at best az
        self.best_az_mot   = np.nan        # motion of best-az orbit
        self.best_az_stop  = np.nan        # stopping score of best-az orbit
        # Best guidance motion (highest = most dynamic)
        self.best_mot      = 0.0
        self.best_mot_p    = None
        self.best_mot_az   = np.inf
        # Best stopping score (lowest = most stationary)
        self.best_stop     = np.inf
        self.best_stop_p   = None
        # First all-criteria-satisfying orbit (if any)
        self.best_all_p    = None

    def update(self, params, vis_ok, az_score, motion_min, stop_max, all_ok):
        """params : (B, 6) — [q, e, i_rad, Om_rad, om_rad, Tp_jd]"""
        B = params.shape[0]
        self.n_sampled += B
        self.n_visible += int(np.sum(vis_ok))
        self.n_all_ok  += int(np.sum(all_ok))

        if not np.any(vis_ok):
            return

        # Best azimuth score (for visible orbits)
        az_vis = np.where(vis_ok, az_score, np.inf)
        idx    = np.argmin(az_vis)
        if az_vis[idx] < self.best_az:
            self.best_az      = float(az_vis[idx])
            self.best_az_p    = params[idx].copy()
            self.best_az_mot  = float(motion_min[idx]) if np.isfinite(motion_min[idx]) else np.nan
            self.best_az_stop = float(stop_max[idx])   if np.isfinite(stop_max[idx])   else np.nan

        # Best guidance motion
        mot_vis = np.where(vis_ok, motion_min, np.nan)
        idx2    = np.nanargmax(mot_vis) if np.any(np.isfinite(mot_vis)) else None
        if idx2 is not None and np.isfinite(mot_vis[idx2]) and mot_vis[idx2] > self.best_mot:
            self.best_mot    = float(mot_vis[idx2])
            self.best_mot_p  = params[idx2].copy()
            self.best_mot_az = float(az_score[idx2])

        # Best stopping score (for visible orbits)
        stp_vis = np.where(vis_ok, stop_max, np.inf)
        idx3    = np.argmin(stp_vis)
        if stp_vis[idx3] < self.best_stop:
            self.best_stop    = float(stp_vis[idx3])
            self.best_stop_p  = params[idx3].copy()

        # Any all-criteria orbit?
        if self.best_all_p is None and np.any(all_ok):
            self.best_all_p = params[np.where(all_ok)[0][0]].copy()

    def _fmt(self, p):
        if p is None:
            return None
        return dict(
            q_au    = float(p[0]),
            e       = float(p[1]),
            i_deg   = float(p[2] * RAD),
            Om_deg  = float(p[3] * RAD),
            om_deg  = float(p[4] * RAD),
            Tp_jd   = float(p[5]),
            dT_days = float(p[5] - JD_NOON),
        )

    def report(self):
        log(f"\n{'━'*70}")
        log(f"SWEEP: {self.name}")
        log(f"  {self.desc}")
        log(f"{'━'*70}")
        log(f"  Sampled              : {self.n_sampled:>13,}")
        log(f"  Visible              : {self.n_visible:>13,}")
        log(f"  All three criteria   : {self.n_all_ok:>13,}"
            + ("  ◄── VIABLE ORBIT FOUND" if self.n_all_ok > 0 else ""))
        log("")
        log(f"  Best azimuth score   : {self.best_az:.5f}°"
            f"  (threshold < {AZ_TOL}°;  0 = corridor centre)")
        if self.best_az_p is not None:
            p = self.best_az_p
            log(f"    orbit: q={p[0]:.5f} AU  e={p[1]:.4f}  i={p[2]*RAD:.3f}°"
                f"  Ω={p[3]*RAD:.2f}°  ω={p[4]*RAD:.2f}°  ΔT={p[5]-JD_NOON:+.2f} d")
            log(f"    at best-az orbit: guidance motion = {self.best_az_mot:.5f} °/h"
                f"  (need > {OMEGA_THRESH})  stopping = {self.best_az_stop:.3f} °/h")
        log("")
        log(f"  Best guidance motion : {self.best_mot:.5f} °/h"
            f"  (threshold > {OMEGA_THRESH} °/h)")
        if self.best_mot_p is not None:
            p = self.best_mot_p
            log(f"    orbit: q={p[0]:.5f} AU  e={p[1]:.4f}  i={p[2]*RAD:.3f}°"
                f"  Ω={p[3]*RAD:.2f}°  ω={p[4]*RAD:.2f}°  ΔT={p[5]-JD_NOON:+.2f} d")
            log(f"    azimuth score at this orbit: {self.best_mot_az:.4f}° "
                f"({'in' if self.best_mot_az < AZ_TOL else 'OUTSIDE'} corridor)")
        log("")
        log(f"  Best stopping score  : {self.best_stop:.5f} °/h"
            f"  (threshold < {OMEGA_THRESH} °/h)")
        log(f"{'━'*70}")

    def to_dict(self):
        return dict(
            name                = self.name,
            description         = self.desc,
            n_sampled           = int(self.n_sampled),
            n_visible           = int(self.n_visible),
            n_all_ok            = int(self.n_all_ok),
            best_az_score       = float(self.best_az),
            best_az_guidance_motion = float(self.best_az_mot),
            best_guidance_motion= float(self.best_mot),
            best_mot_az_score   = float(self.best_mot_az),
            best_stopping_score = float(self.best_stop),
            best_az_orbit       = self._fmt(self.best_az_p),
            best_mot_orbit      = self._fmt(self.best_mot_p),
            best_stop_orbit     = self._fmt(self.best_stop_p),
            best_all_orbit      = self._fmt(self.best_all_p),
        )


# ─────────────────────────────────────────────────────────────────────────
# MC sweep runner
# ─────────────────────────────────────────────────────────────────────────

def run_sweep(sweep_def, n_samples, batch_size,
              enu, earth_pos, jds_gs, N_G, DT_H, rng):
    """
    Draw n_samples random configurations from sweep_def ranges,
    evaluate in batches of batch_size, and return a SweepStats object.

    sweep_def keys (all ranges [lo, hi]):
        q_range, e_range, i_range (degrees), Om_range, om_range (degrees),
        dT_range (days from JD_NOON → Tp = JD_NOON + dT)
    """
    stats     = SweepStats(sweep_def['name'], sweep_def['desc'])
    n_batches = (n_samples + batch_size - 1) // batch_size
    t_sw      = walltime.time()

    for b in range(n_batches):
        bsz = min(batch_size, n_samples - b * batch_size)
        if bsz <= 0:
            break

        # Sample parameters uniformly within the defined ranges
        q_   = rng.uniform(*sweep_def['q_range'],   bsz)
        e_   = rng.uniform(*sweep_def['e_range'],   bsz)
        i_d  = np.clip(rng.uniform(*sweep_def['i_range'],  bsz), 0.0, 180.0)
        Om_d = rng.uniform(*sweep_def['Om_range'],  bsz) % 360.0
        om_d = rng.uniform(*sweep_def['om_range'],  bsz) % 360.0
        dT   = rng.uniform(*sweep_def['dT_range'],  bsz)
        Tp   = JD_NOON + dT

        # Pack into (B, 6) for SweepStats
        params = np.column_stack([q_, e_, i_d*DEG, Om_d*DEG, om_d*DEG, Tp])

        vis_ok, az_score, motion_min, stop_max, all_ok = evaluate_batch(
            q_, e_, i_d*DEG, Om_d*DEG, om_d*DEG, Tp,
            enu, earth_pos, jds_gs, N_G, DT_H
        )

        stats.update(params, vis_ok, az_score, motion_min, stop_max, all_ok)

        # Progress every ~5%
        if (b + 1) % max(1, n_batches // 20) == 0 or b == n_batches - 1:
            elapsed = walltime.time() - t_sw
            rate    = stats.n_sampled / elapsed if elapsed > 0 else 0
            pct     = 100.0 * (b + 1) / n_batches
            log(f"  [{sweep_def['name'][:10]}] {pct:5.1f}%  "
                f"sampled={stats.n_sampled:>9,}  "
                f"vis={stats.n_visible:>8,}  "
                f"all_ok={stats.n_all_ok}  "
                f"best_az={stats.best_az:.4f}°  "
                f"best_mot={stats.best_mot:.3f}°/h  "
                f"{rate/1e6:.2f}M/s")

    return stats


# ─────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Monte Carlo refinement around best-case SoB orbital candidates')
    parser.add_argument('--n1',    type=int, default=10_000_000,
                        help='Samples for sweep 1 — best-fit neighbourhood (default: 10M)')
    parser.add_argument('--n2',    type=int, default=5_000_000,
                        help='Samples for sweep 2 — Matney neighbourhood (default: 5M)')
    parser.add_argument('--batch', type=int, default=200_000,
                        help='Vectorised batch size (default: 200k)')
    parser.add_argument('--seed',  type=int, default=42,
                        help='Random seed (default: 42)')
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    log('═' * 70)
    log('sob_mc_refined.py — Monte Carlo refinement of best-case SoB orbits')
    log(f'  Sweep 1 (best-fit neighbourhood): {args.n1:,} samples')
    log(f'  Sweep 2 (Matney neighbourhood)  : {args.n2:,} samples')
    log(f'  Batch size: {args.batch:,}   Seed: {args.seed}')
    log('═' * 70)

    # Pre-compute environment (shared by both sweeps)
    enu, earth_pos, jds_gs, N_G = precompute_environment()
    N_GS = len(jds_gs)
    DT_H = 0.25   # hours between time steps

    # ── Sweep 1: best-fit guidance orbit neighbourhood ────────────────────
    # The grid survey's best guidance orbit had e=1.1, q=0.03 AU, i=10°.
    # We search broadly around these values, with Ω and ω over full 360°,
    # to rule out that any nearby continuous-space orbit is viable.
    sweep1 = dict(
        name     = 'best_fit_nbhd',
        desc     = ('Neighbourhood of best-fit grid orbit: '
                    'e=1.1, q=0.03 AU, i=10°, full Ω/ω, ΔT ∈ [−60,+60] d'),
        q_range  = (0.005, 0.15),    # AU: generous range around 0.03
        e_range  = (0.7,   2.50),    # elliptic through strongly hyperbolic
        i_range  = (0.0,   40.0),    # degrees: centred on 10°, broad
        Om_range = (0.0,  360.0),    # full range
        om_range = (0.0,  360.0),    # full range
        dT_range = (-60.0, 60.0),    # days: perihelion offset from JD_NOON
    )

    log(f"\n── Sweep 1: {sweep1['desc']}")
    s1 = run_sweep(sweep1, args.n1, args.batch,
                   enu, earth_pos, jds_gs, N_G, DT_H, rng)
    s1.report()

    # ── Sweep 2: Matney (2025) comet neighbourhood ────────────────────────
    # Matney's 5 BCE elements: q=0.0407 AU, e=1.0, i=1.43°,
    # ω=335.58°, Ω=75.81°, T_peri=JD 1719785.565 (ΔT≈+29.7 d).
    # We search within ±40° of Ω and ω, generous range in q, e, i,
    # and ΔT bracketing Matney's perihelion.
    sweep2 = dict(
        name     = 'matney_nbhd',
        desc     = (f'Neighbourhood of Matney 2025 comet: '
                    f'q={MATNEY["q"]} AU, e={MATNEY["e"]}, '
                    f'i={MATNEY["i"]}°, ΔT≈+{MATNEY["dT"]:.1f} d'),
        q_range  = (MATNEY['q'] * 0.25, MATNEY['q'] * 4.0),  # 0.010–0.163 AU
        e_range  = (0.70, 1.50),
        i_range  = (0.0,  8.0),                               # 0°–8°
        Om_range = (MATNEY['asc'] - 50, MATNEY['asc'] + 50),  # 25.81°–125.81°
        om_range = (MATNEY['om']  - 50, MATNEY['om']  + 50),  # 285.58°–385.58°
        dT_range = (MATNEY['dT'] - 40,  MATNEY['dT'] + 40),   # ≈ −10 to +70 d
    )

    log(f"\n── Sweep 2: {sweep2['desc']}")
    s2 = run_sweep(sweep2, args.n2, args.batch,
                   enu, earth_pos, jds_gs, N_G, DT_H, rng)
    s2.report()

    # ── Save results ──────────────────────────────────────────────────────
    output = dict(
        description    = 'Monte Carlo refinement around best-case SoB orbital candidates',
        reference_jd   = JD_NOON,
        omega_thresh   = OMEGA_THRESH,
        az_corridor    = [AZ_LO, AZ_HI],
        az_tolerance   = AZ_TOL,
        alt_floor      = ALT_FLOOR,
        matney_params  = MATNEY,
        sweeps         = [s1.to_dict(), s2.to_dict()],
    )
    with open(OUTFILE, 'w') as f:
        json.dump(output, f, indent=2)
    log(f"\nResults saved → {OUTFILE}")

    # ── Final verdict ─────────────────────────────────────────────────────
    log('\n' + '═' * 70)
    log('FINAL VERDICT')
    log('═' * 70)
    total_sampled = s1.n_sampled + s2.n_sampled
    total_visible = s1.n_visible + s2.n_visible
    total_all_ok  = s1.n_all_ok  + s2.n_all_ok
    log(f"  Total configurations sampled  : {total_sampled:>13,}")
    log(f"  Total visible                 : {total_visible:>13,}")
    log(f"  Satisfying all three criteria : {total_all_ok:>13,}")
    log("")
    for s in (s1, s2):
        log(f"  [{s.name}]")
        log(f"    Best az score:    {s.best_az:.5f}°  "
            f"(threshold < {AZ_TOL}°)")
        log(f"    Motion at best-az: {s.best_az_mot:.5f} °/h  "
            f"(need > {OMEGA_THRESH}; ratio = "
            + (f"{s.best_az_mot/OMEGA_THRESH:.3f}×)" if np.isfinite(s.best_az_mot) else "N/A)"))
        log(f"    Best guidance motion anywhere: {s.best_mot:.5f} °/h  "
            f"az_score of that orbit: {s.best_mot_az:.4f}°")
    log("")
    if total_all_ok == 0:
        log("  CONCLUSION: Zero configurations satisfying all three criteria")
        log("  found in either Monte Carlo sweep.")
        log("  The grid-survey null result is confirmed in continuous")
        log("  parameter space: no fine-tuning around the best-fit guidance")
        log("  orbit or Matney's proposed comet reveals a viable candidate.")
        log("  The impossibility is not an artefact of grid spacing.")
    else:
        log("  *** WARNING: viable orbit(s) found — see JSON for parameters ***")
    log('═' * 70)


if __name__ == '__main__':
    main()
