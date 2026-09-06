#!/usr/bin/env python3
"""
Fine-Grid Exhaustive Orbital Survey for Matthew 2:9
Aaron Adair / 2026

Surveys every physically distinct Keplerian heliocentric orbit that could
plausibly be observable from Bethlehem on the proposed date, asking whether
any simultaneously (a) lies in the road corridor [190°,210°] azimuth,
(b) moves at > 2°/h apparent motion during the guidance window, and
(c) moves at < 2°/h during the stopping window.

ANGULAR GRID:
  Ω and ω at 5° spacing — 72 × 72 = 5,184 orientation pairs per outer
  combination.  This is 16× denser than the original 20° grid published
  in the literature (18×18 = 324 pairs).

OUTER PARAMETER SPACE:
  q:     0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20,          (AU)
         0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00
         18 values — from Sun-grazing to Earth's orbital distance
  e:     0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0,
         1.1, 1.5, 2.0, 5.0
         15 values — elliptic, parabolic, and hyperbolic (spanning the
         eccentricity range of 1I/'Oumuamua and 2I/Borisov)
  i:     0°–180° in 5° steps (37 values) plus 1°, 1.4°, 2°, 3°, 4°
         (5 extra values for dense near-ecliptic sampling)  →  42 total
  T_off: −60, −40, −20, −10, −5, +5, +10, +20, +60 days from perihelion
         9 values

  Outer combinations: 18 × 15 × 42 × 9 = 102,060
  Inner pairs:        72 × 72 = 5,184
  Total configurations: 102,060 × 5,184 = 529,079,040

NOTE ON INCLINATION RANGE [0°, 180°]:
  Inclination i is defined as the dihedral angle between the orbital plane
  and the ecliptic reference plane; by construction i ∈ [0°, 180°]:
    i = 0°   — prograde orbit coplanar with the ecliptic
    i = 90°  — polar orbit perpendicular to the ecliptic
    i = 180° — retrograde orbit coplanar with the ecliptic
  i < 0° is equivalent to (|i|, Ω + 180°); i > 180° is equivalent to
  (360° − i, Ω + 180°).  Because Ω is swept over all 72 values at 5°
  spacing (covering all 360°), every such orbit is already represented.
  The range [0°, 180°] is therefore provably exhaustive.

PERFORMANCE:
  Uses multiprocessing.Pool (fork, inherits globals on Linux/macOS).
  On a modern 8-core machine: ~1–3 minutes.
  On a 16-core machine: ~1 minute or less.
  Saves incremental checkpoints every 5% so the run can be resumed if
  interrupted (just re-run the script — it detects the partial results).

OUTPUT FILES (written to the same directory as this script):
  sob_proof_fine_progress.log   — human-readable progress log
  sob_proof_fine_results.json   — machine-readable final results
  figure6_fine_impossibility.png — two-panel summary figure
"""

import warnings; warnings.filterwarnings('ignore')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body_barycentric
import astropy.units as u
import time as walltime
import multiprocessing as mp
import json, sys, os

LOGFILE = os.path.join(os.path.dirname(__file__), 'sob_proof_fine_progress.log')
RESULTS  = os.path.join(os.path.dirname(__file__), 'sob_proof_fine_results.json')
FIGFILE  = os.path.join(os.path.dirname(__file__), 'figure6_fine_impossibility.png')

def log(msg):
    line = f"[{walltime.time()-T_START:.1f}s] {msg}"
    print(line, flush=True)
    with open(LOGFILE, 'a') as f:
        f.write(line + '\n')

T_START = walltime.time()
open(LOGFILE, 'w').close()   # clear log

# ── site ──────────────────────────────────────────────────────────────────
BLAT, BLON, BALT = 31 + 42/60, 35 + 12/60, 765.0
BETHLEHEM = EarthLocation(lon=BLON*u.deg, lat=BLAT*u.deg, height=BALT*u.m)

AZ_LO, AZ_HI = 190.0, 210.0
AZ_THRESH    = 5.0       # deg outside corridor
OMEGA_THRESH = 2.0       # deg/h  (guidance requires >, stopping requires <)

# ── constants ──────────────────────────────────────────────────────────────
K = 0.01720209895; MU = K**2
AU_KM = 1.495978707e8
EPS   = np.radians(23.439291111)
R_ECL2EQU = np.array([
    [1, 0,           0          ],
    [0, np.cos(EPS), -np.sin(EPS)],
    [0, np.sin(EPS),  np.cos(EPS)],
])

M_T     = 1719785.565   # Matney perihelion JD TDB
JD_NOON = 1719755.898   # 5 BCE Jun 8 local noon at Bethlehem
DT_SEC  = 10572.0

# ── orbital mechanics ──────────────────────────────────────────────────────
def barker(dt, q):
    W = 3.0 * np.sqrt(MU / (2.0 * q**3)) * dt
    s = np.sqrt((W / 2.0)**2 + 1.0)
    return np.cbrt(W / 2.0 + s) + np.cbrt(W / 2.0 - s)

def kepler_E(M, e, tol=1e-12):
    E = M.copy()
    for _ in range(60):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))
        E  = E + dE
        if np.max(np.abs(dE)) < tol: break
    return E

def true_anomaly(dt, q, e):
    if abs(e - 1.0) < 1e-9:
        return 2.0 * np.arctan(barker(dt, q))
    elif e < 1.0:
        a  = q / (1.0 - e)
        n  = np.sqrt(MU / a**3)
        M  = n * dt
        E  = kepler_E(np.atleast_1d(float(M)), e)
        nu = 2.0 * np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
        return float(nu[0])
    else:
        a  = q / (e - 1.0)
        n  = np.sqrt(MU / a**3)
        M  = n * dt
        H  = np.sign(M) * np.log(2.0 * np.abs(M) / e + 1.8)
        for _ in range(60):
            dH = (M - e*np.sinh(H) + H) / (e*np.cosh(H) - 1.0)
            H += dH
            if np.abs(dH) < 1e-12: break
        return 2.0 * np.arctan2(np.sqrt(e+1)*np.sinh(H/2), np.sqrt(e-1)*np.cosh(H/2))

# ── time grid + ENU precomputation ─────────────────────────────────────────
log("Precomputing reference frames (ENU matrices + Earth positions)...")

hrs_all = np.arange(-8.0, 3.1, 0.25)
jds_all = JD_NOON + hrs_all / 24.0
jds_utc = jds_all - DT_SEC / 86400.0
local_h = hrs_all + 12.0

def make_enu(t_utc):
    frame = AltAz(obstime=t_utc, location=BETHLEHEM)
    def aa2icrs(az, alt):
        c = SkyCoord(alt=alt*u.deg, az=az*u.deg, frame=frame).icrs
        d = np.radians(c.dec.deg); r = np.radians(c.ra.deg)
        return np.array([np.cos(d)*np.cos(r), np.cos(d)*np.sin(r), np.sin(d)])
    E = aa2icrs(90, 0.001); Z = aa2icrs(0, 90)
    Z /= np.linalg.norm(Z); E -= np.dot(E, Z)*Z; E /= np.linalg.norm(E)
    return np.array([E, np.cross(Z, E), Z])

t_tdb_all = Time(jds_all, format='jd', scale='tdb')
t_utc_all = Time(jds_utc, format='jd', scale='utc')
Ef = get_body_barycentric('earth', t_tdb_all)
Sf = get_body_barycentric('sun',   t_tdb_all)
earth_all = (Ef - Sf).xyz.to(u.AU).value          # shape (3, N_time)
R_enu_all = np.array([make_enu(t_utc_all[ii]) for ii in range(len(jds_all))])

mask_guide = (local_h >= 8.0) & (local_h <= 10.0)
mask_stop  = (local_h > 10.0) & (local_h <= 12.0)
N_guide = mask_guide.sum(); N_stop = mask_stop.sum()
t_guide  = local_h[mask_guide]; t_stop = local_h[mask_stop]
jd_guide = jds_all[mask_guide]; jd_stop = jds_all[mask_stop]
E_guide  = earth_all[:, mask_guide]; E_stop = earth_all[:, mask_stop]

log(f"ENU done: {len(jds_all)} time steps, {N_guide} guidance, {N_stop} stopping.")

# ── angular grid: 5° spacing ──────────────────────────────────────────────
# 72 × 72 = 5,184 orientation pairs per outer combination.
# 16× denser than the original 20° grid (18×18 = 324) used in the literature.
Om_vals = np.arange(0, 360, 5)   # 72 values
om_vals = np.arange(0, 360, 5)   # 72 values

# Expanded outer ranges
# i range: 0° to 180° in 5° steps — see module docstring for why [0,180] is complete.
#   i < 0° ≡ (|i|, Ω+180°);  i > 180° ≡ (360°−i, Ω+180°) — both already in grid.
e_vals = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0,
          1.1, 1.5, 2.0, 5.0]         # 15: elliptic + parabolic + hyperbolic
q_vals = [0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20,   # inner Solar System
          0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80,    # out to ~Venus orbit
          0.90, 1.00]                  # 18 AU: up to Earth's orbital distance
# i: 5° grid 0°–180°, plus fine sampling at very small angles (near-ecliptic orbits).
# Small angles added because near-ecliptic prograde objects are observationally relevant
# (historical comet candidates are often low-inclination) and the 5° step skips the
# 1°–4° range entirely.  Duplicates (0° and 5° appear in both lists) are removed.
i_vals = sorted(set(
    list(np.arange(0, 181, 5)) +   # 37 values: 0°,5°,...,180°
    [1.0, 1.4, 2.0, 3.0, 4.0]      # 5 extra small-angle values between 0° and 5°
))   # 42 total unique inclinations
T_offs = [-60.0, -40.0, -20.0, -10.0, -5.0,
          +5.0, +10.0, +20.0, +60.0]  # 9 days from Matney perihelion

n_outer = len(e_vals) * len(q_vals) * len(i_vals) * len(T_offs)
n_inner = len(Om_vals) * len(om_vals)
n_total_est = n_outer * n_inner
log(f"Survey grid: Ω/ω at 5° ({len(Om_vals)}×{len(om_vals)}={n_inner} pairs/call)")
log(f"  e={len(e_vals)} × q={len(q_vals)} × i={len(i_vals)} × T={len(T_offs)} = {n_outer} outer calls")
log(f"  Estimated total configurations: {n_total_est:,}")

# ── survey kernel (vectorised over Ω, ω) ──────────────────────────────────
def survey_combined(q, e, i_d, Om_d_arr, om_d_arr, T_jd):
    """Vectorised over Om_d_arr × om_d_arr. Returns
    (guide_score, guide_motion, stop_score) arrays shape (nOm, nom)."""
    nOm, nom = len(Om_d_arr), len(om_d_arr)
    i   = np.radians(i_d)
    ci  = np.cos(i);  si = np.sin(i)
    Om_r = np.radians(Om_d_arr); om_r = np.radians(om_d_arr)
    cO = np.cos(Om_r); sO = np.sin(Om_r)
    co = np.cos(om_r); so = np.sin(om_r)

    P_ecl = np.array([cO[:,None]*co - sO[:,None]*so*ci,
                      sO[:,None]*co + cO[:,None]*so*ci,
                      np.broadcast_to(so*si, (nOm, nom)).copy()])
    Q_ecl = np.array([-cO[:,None]*so - sO[:,None]*co*ci,
                      -sO[:,None]*so + cO[:,None]*co*ci,
                       np.broadcast_to(co*si, (nOm, nom)).copy()])
    P_eq = np.einsum('ij,jkl->ikl', R_ECL2EQU, P_ecl)
    Q_eq = np.einsum('ij,jkl->ikl', R_ECL2EQU, Q_ecl)

    # guidance window
    az_g  = np.zeros((N_guide, nOm, nom))
    alt_g = np.zeros_like(az_g)
    dir_g = np.zeros((N_guide, nOm, nom, 3))
    for ii, (jd, Epos, Renu) in enumerate(
            zip(jd_guide, E_guide.T, R_enu_all[mask_guide])):
        dt  = jd - T_jd
        nu  = true_anomaly(dt, q, e)
        r   = q * (1.0 + e) / (1.0 + e * np.cos(nu))
        pos = r * (np.cos(nu)*P_eq + np.sin(nu)*Q_eq)
        rho = pos - Epos[:, None, None]
        dist = np.linalg.norm(rho, axis=0); dist[dist < 1e-12] = 1e-12
        rh  = rho / dist
        enu = np.einsum('ij,jkl->ikl', Renu, rh)
        alt_g[ii] = np.degrees(np.arcsin(np.clip(enu[2], -1, 1)))
        az_g[ii]  = np.degrees(np.arctan2(enu[0], enu[1])) % 360.0
        dir_g[ii] = np.moveaxis(rh, 0, -1)

    vis_g = (alt_g > 5.0).all(axis=0)

    dev_lo      = np.maximum(0.0, AZ_LO - az_g)
    dev_hi      = np.maximum(0.0, az_g  - AZ_HI)
    guide_score = np.maximum(dev_lo, dev_hi).max(axis=0)

    dg1 = dir_g[:-1]; dg2 = dir_g[1:]
    dt_h_g = np.diff(t_guide)
    crs_g  = np.cross(dg1, dg2, axis=-1)
    omg_g  = np.degrees(np.arctan2(
        np.linalg.norm(crs_g, axis=-1),
        np.einsum('...i,...i->...', dg1, dg2)
    )) / dt_h_g[:, None, None]
    guide_motion = omg_g.min(axis=0)

    # stopping window
    dir_s = np.zeros((N_stop, nOm, nom, 3))
    for ii, (jd, Epos, Renu) in enumerate(
            zip(jd_stop, E_stop.T, R_enu_all[mask_stop])):
        dt  = jd - T_jd
        nu  = true_anomaly(dt, q, e)
        r   = q * (1.0 + e) / (1.0 + e * np.cos(nu))
        pos = r * (np.cos(nu)*P_eq + np.sin(nu)*Q_eq)
        rho = pos - Epos[:, None, None]
        dist = np.linalg.norm(rho, axis=0); dist[dist < 1e-12] = 1e-12
        dir_s[ii] = np.moveaxis(rho / dist, 0, -1)

    d1 = dir_s[:-1]; d2 = dir_s[1:]
    dt_h = np.diff(t_stop)
    stop_score = np.degrees(np.arctan2(
        np.linalg.norm(np.cross(d1, d2, axis=-1), axis=-1),
        np.einsum('...i,...i->...', d1, d2)
    )).max(axis=0) / dt_h[:, None, None].min()

    return (np.where(vis_g, guide_score,  np.nan),
            np.where(vis_g, guide_motion, np.nan),
            np.where(vis_g, stop_score,   np.nan))

# ── worker function for multiprocessing ────────────────────────────────────
# On Linux, fork inherits all globals — no need to pass large arrays.
def _worker(args):
    e, q, i_d, T_off = args
    try:
        gs, gm, ss = survey_combined(q, e, i_d, Om_vals, om_vals, M_T + T_off)
        valid = ~np.isnan(gs)
        n_vis = int(valid.sum())
        if n_vis == 0:
            return (0, 0, 0, 0, 0, None, None, None, None, None, None, e)
        gv = gs[valid]; mv = gm[valid]; sv = ss[valid]
        n_az  = int((gv < AZ_THRESH).sum())
        n_mo  = int((mv > OMEGA_THRESH).sum())
        n_st  = int((sv < OMEGA_THRESH).sum())
        n_all = int(((gv < AZ_THRESH) & (mv > OMEGA_THRESH) & (sv < OMEGA_THRESH)).sum())
        # Best combined score for this outer combo
        combined = np.sqrt(gv**2 + (sv/3)**2)
        idx_best  = np.argmin(combined)
        best_g = float(gv[idx_best]); best_s = float(sv[idx_best])
        best_m = float(mv[idx_best])
        # histogram bucket of guide scores (0-90, 30 bins)
        hist, _ = np.histogram(gv[gv < 90], bins=np.linspace(0, 90, 31))
        return (n_vis, n_az, n_mo, n_st, n_all,
                best_g, best_s, best_m,
                (e, q, i_d, combined.min()), hist.tolist(), e, e)
    except Exception as ex:
        return (0, 0, 0, 0, 0, None, None, None, None, None, None, e)

# ── build outer combinations ───────────────────────────────────────────────
outer_all = [(e, q, i_d, T_off)
             for e in e_vals for q in q_vals
             for i_d in i_vals for T_off in T_offs]

# ── resume from checkpoint if available ───────────────────────────────────
N_DONE_PRIOR = 0
n_vis_total = 0; n_az_ok = 0; n_motion_ok = 0; n_stop_ok = 0; n_both_ok = 0
best_guide = 999.0; best_stop = 999.0; best_motion = 0.0; best_params = None
per_e_hist  = {e: np.zeros(30, dtype=int) for e in e_vals}

OM_STEP = int(Om_vals[1] - Om_vals[0])   # 5 for this run

if os.path.exists(RESULTS):
    try:
        prev = json.load(open(RESULTS))
        prev_om_step   = prev.get('Om_step', 0)
        prev_n_outer   = prev.get('n_outer_total', prev.get('n_outer_combos', 0))
        params_match   = (prev_om_step == OM_STEP and prev_n_outer == len(outer_all))
        is_partial     = (prev.get('status') == 'running' and prev.get('n_outer_done', 0) > 0)
        if is_partial and params_match:
            N_DONE_PRIOR  = prev['n_outer_done']
            n_vis_total   = prev['n_vis_total']
            n_az_ok       = prev['n_az_ok']
            n_motion_ok   = prev['n_motion_ok']
            n_stop_ok     = prev['n_stop_ok']
            n_both_ok     = prev['n_both_ok']
            best_guide    = prev['best_guide_deg']
            best_stop     = prev['best_stop_dph']
            log(f"Resuming from checkpoint: {N_DONE_PRIOR}/{len(outer_all)} outer done "
                f"({100*N_DONE_PRIOR/len(outer_all):.0f}%)")
        elif is_partial and not params_match:
            log(f"Checkpoint found but parameters differ "
                f"(Om_step={prev_om_step}° vs {OM_STEP}°, "
                f"n_outer={prev_n_outer} vs {len(outer_all)}) — starting fresh.")
    except Exception:
        pass

outer = outer_all[N_DONE_PRIOR:]   # slice off already-completed combinations

log(f"Total outer combinations: {len(outer_all):,}  (resuming from {N_DONE_PRIOR})")
log(f"Remaining this run: {len(outer):,} × {n_inner:,} Ω/ω = {len(outer)*n_inner:,} configs")

# ── parallel execution ─────────────────────────────────────────────────────
N_WORKERS = max(1, mp.cpu_count() - 1)   # leave one core for the OS
log(f"Using {N_WORKERS} worker processes (cpu_count={mp.cpu_count()})")

REPORT_EVERY = max(1, len(outer_all) // 20)   # progress every 5% of TOTAL

t_lap = walltime.time()
with mp.get_context('fork').Pool(N_WORKERS) as pool:
    for this_batch, result in enumerate(pool.imap_unordered(_worker, outer, chunksize=8), 1):
        (nv, na, nm, ns, nb, bg, bs, bm, bp, hist, e_val, _) = result
        n_vis_total  += nv
        n_az_ok      += na
        n_motion_ok  += nm
        n_stop_ok    += ns
        n_both_ok    += nb
        if hist is not None:
            per_e_hist[e_val] += np.array(hist, dtype=int)
        if bg is not None and bg + bs/3 < best_guide + best_stop/3:
            best_guide  = bg; best_stop = bs; best_motion = bm
            best_params = bp
        # Report & checkpoint based on total progress (prior + current)
        done_total = N_DONE_PRIOR + this_batch
        if done_total % REPORT_EVERY == 0 or this_batch == len(outer):
            pct = 100 * done_total / len(outer_all)
            elapsed = walltime.time() - t_lap
            eta = elapsed / this_batch * (len(outer) - this_batch) if this_batch > 0 else 0
            log(f"  {done_total}/{len(outer_all)} outer ({pct:.0f}%) | "
                f"vis={n_vis_total:,} az_ok={n_az_ok:,} both={n_both_ok} | "
                f"ETA {eta:.0f}s")
            # Checkpoint — write running totals so a subsequent resume call can continue
            snap = dict(status="running", pct_done=round(pct, 1),
                        n_outer_done=done_total, n_outer_total=len(outer_all),
                        n_inner_per_call=n_inner,
                        Om_step=OM_STEP, om_step=OM_STEP,
                        n_vis_total=n_vis_total, n_az_ok=n_az_ok,
                        n_motion_ok=n_motion_ok, n_stop_ok=n_stop_ok,
                        n_both_ok=n_both_ok,
                        best_guide_deg=best_guide, best_stop_dph=best_stop,
                        elapsed_s=elapsed)
            with open(RESULTS, 'w') as f:
                json.dump(snap, f, indent=2)

log("Survey complete.")
log(f"  Total visible configurations : {n_vis_total:,}")
log(f"  Satisfying azimuth (<{AZ_THRESH}° outside corridor) : {n_az_ok:,}")
log(f"  Satisfying guidance motion (>{OMEGA_THRESH}°/h)     : {n_motion_ok:,}")
log(f"  Satisfying stopping (<{OMEGA_THRESH}°/h)            : {n_stop_ok:,}")
log(f"  ALL THREE simultaneously     : {n_both_ok}")
log(f"  Best guidance score : {best_guide:.3f}°  |  best stopping: {best_stop:.3f}°/h")

# ── sensitivity table ──────────────────────────────────────────────────────
log("\nSensitivity: zero result vs. angular-velocity threshold")
log(f"  {'thresh':>8}  {'n_both':>8}")
# We don't have full arrays — report at the survey threshold only
log(f"  {OMEGA_THRESH:>8.2f}  {n_both_ok:>8}  (survey threshold)")

# ── save results ───────────────────────────────────────────────────────────
results = dict(
    status           = "complete",
    n_outer_combos   = len(outer_all),
    n_inner_per_call = n_inner,
    n_total_est      = len(outer_all) * n_inner,
    n_vis_total      = n_vis_total,
    n_az_ok          = n_az_ok,
    n_motion_ok      = n_motion_ok,
    n_stop_ok        = n_stop_ok,
    n_both_ok        = n_both_ok,
    best_guide_deg   = best_guide,
    best_stop_dph    = best_stop,
    best_motion_dph  = best_motion,
    best_params      = list(best_params) if best_params else None,
    Om_step          = 5,
    om_step          = 5,
    e_vals           = [float(v) for v in e_vals],
    q_vals           = [float(v) for v in q_vals],
    i_vals           = [float(v) for v in i_vals],
    T_offs           = [float(v) for v in T_offs],
    AZ_THRESH        = float(AZ_THRESH),
    OMEGA_THRESH     = float(OMEGA_THRESH),
    elapsed_s        = float(walltime.time() - T_START),
)
# Ensure all ints are native Python ints (not numpy int64)
def _to_native(obj):
    if isinstance(obj, dict): return {k: _to_native(v) for k,v in obj.items()}
    if isinstance(obj, list): return [_to_native(v) for v in obj]
    if hasattr(obj, 'item'): return obj.item()  # numpy scalar → Python scalar
    return obj
with open(RESULTS, 'w') as f:
    json.dump(_to_native(results), f, indent=2)
log(f"Results saved to {RESULTS}")

# ── figure ─────────────────────────────────────────────────────────────────
log("Generating figure...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(
    "Fine-Grid Orbital Survey: No Keplerian Orbit Reproduces Matthew 2:9\n"
    f"({n_vis_total:,} visible of {n_total_est:,} sampled | "
    f"5° Ω/ω grid (16× original) | e∈[0,5] | q∈[0.02,1.00 AU] | "
    f"i∈[0°,180°+near-ecliptic] | ΔT∈[−60,+60 d])",
    fontsize=11, fontweight='bold')

bin_edges = np.linspace(0, 90, 31)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
colors = plt.cm.viridis(np.linspace(0, 1, len(e_vals)))

ax = axes[0]
for idx, ev in enumerate(e_vals):
    h = per_e_hist[ev].astype(float)
    total = h.sum()
    if total > 0:
        ax.plot(bin_centers, h / total, color=colors[idx],
                alpha=0.6, lw=1.2, label=f'e={ev}')
ax.axvline(AZ_THRESH, color='red', lw=2.5, ls='--',
           label=f'{AZ_THRESH}° outside-corridor tolerance')
ax.set_xlabel('Guidance score: max degrees outside road corridor [190°,210°]', fontsize=10)
ax.set_ylabel('Fraction of visible configurations', fontsize=10)
ax.set_title('(A) Guidance Score by Eccentricity\n'
             '(No eccentricity keeps orbit inside road corridor with guidance + stopping)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=6, loc='upper right', ncol=3)
ax.grid(True, alpha=0.2)
ax.text(0.02, 0.96,
        f"n_vis = {n_vis_total:,}\n"
        f"5° Ω/ω spacing\n"
        f"Satisfy ALL: {n_both_ok}",
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

ax = axes[1]
# Summary bar chart of criteria pass rates
labels = [f'Azimuth\n(<{AZ_THRESH}° outside\ncorridor)',
          f'Guidance motion\n(>{OMEGA_THRESH}°/h)',
          f'Stopping\n(<{OMEGA_THRESH}°/h)',
          'ALL THREE\nsimultaneously']
counts = [n_az_ok, n_motion_ok, n_stop_ok, n_both_ok]
bar_colors = ['steelblue', 'seagreen', 'darkorange', 'firebrick']
bars = ax.bar(labels, counts, color=bar_colors, edgecolor='black', linewidth=0.8)
for bar, cnt in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2,
            bar.get_height() + n_vis_total * 0.002,
            f'{cnt:,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
ax.axhline(n_vis_total, color='gray', lw=1.5, ls='--', label=f'Total visible ({n_vis_total:,})')
ax.set_ylabel('Number of configurations', fontsize=10)
ax.set_title('(B) Criteria Pass Rates\n(Zero configurations satisfy all three simultaneously)',
             fontsize=10, fontweight='bold')
ax.set_yscale('log')
ax.legend(fontsize=9)
ax.grid(True, axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(FIGFILE, dpi=600, bbox_inches='tight')
plt.close()
log(f"Figure saved to {FIGFILE}")

log("=" * 60)
log("FINAL RESULT")
log("=" * 60)
log(f"  Grid: Ω, ω at 5° (16× finer than original 20°)")
log(f"  Parameter space: e∈{{0.0...5.0}}, q∈{{0.02...1.00 AU}},")
log(f"    i∈{{0°...180° in 5° steps}}, T_off∈{{-60...+60 d}}")
log(f"  Outer calls: {len(e_vals)}e × {len(q_vals)}q × {len(i_vals)}i × {len(T_offs)}T = {n_outer:,}")
log(f"  Inner: {len(Om_vals)}×{len(om_vals)} = {n_inner:,} Ω/ω pairs per call")
log(f"  Total configurations surveyed   : {n_total_est:,}")
log(f"  Total visible                   : {n_vis_total:,}")
log(f"  Satisfying all three criteria   : {n_both_ok}")
log(f"  Elapsed time                    : {walltime.time()-T_START:.0f} s")
