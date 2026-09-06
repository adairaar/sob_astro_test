#!/usr/bin/env python3
"""
sob_criteria_sweep.py — how far does the conclusion depend on the criteria?

Written in direct response to Reviewer 2's third round, which asks:

    "Rather than choosing a rigid screening criterion and seeing if orbits
     match, find out what are the tightest criteria that still allow such an
     orbit to exist or not. ... Maybe try stopping criteria of 1º/hr, 2º/hr,
     3º/hr, 4º/hr, and 5º/hr ... That way, it will become clear how dependent
     the orbit matching is on the detailed criterion chosen."

Two changes from sob_close_flyby_mc.py, both conceded to the referee.

1. THE STOPPING CRITERION IS NOW A FINITE ANGULAR SHIFT, NOT A RATE.
   R2 observes that the published thresholds for the eye's detection of
   absolute motion (~0.03 deg/s with a fixed reference; 0.5-1.0 deg/s in
   darkness) lie far above any rate at issue here, and concludes that "a
   finite angular shift over some finite time (as opposed to perceiving the
   actual motion) is probably a better metric".  That is correct, and it is
   adopted.  The object is held to have "stopped" if its total great-circle
   displacement over a window of SHIFT_WIN hours is below some threshold.

   The rate-based metric of the earlier version is retained alongside it so
   that the two can be compared, and so that nothing turns on the change.

2. THE AZIMUTH CORRIDOR IS SWEPT, NOT FIXED.
   R2 notes that Bethlehem lies east of the Jerusalem road, so that a party
   travelling to a house off the road -- he offers the Church of the Nativity
   as an example -- would have needed an azimuth nearer 170 deg, or further
   east still, for the final leg.  The corridor's lower edge is therefore
   swept from 190 deg down to 150 deg.

WHAT THIS SCRIPT REPORTS
    For every orbit surviving the close-flyby test that also admits a
    guidance window, the two stopping diagnostics are recorded.  The pass
    counts at every threshold then follow by post-processing, so the whole
    two-dimensional sweep costs one Monte Carlo run rather than twenty-four.

    The headline number the referee asked for is the MINIMUM achieved
    stopping metric over all candidates: the tightest criterion that any
    orbit satisfies.  If that minimum lies above 5 deg/h, the conclusion is
    insensitive to the whole range R2 proposes.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sob_close_flyby_mc as M                                 # noqa: E402

N_MC = int(os.environ.get("SOB_N_MC", 20_000_000))
SEED = int(os.environ.get("SOB_SEED", 2026))

# corridor lower edges to sweep; upper edge held at M.AZ_HI
AZ_LO_GRID = [float(x) for x in
              os.environ.get("AZ_LO_GRID", "190,180,170,160,150").split(",")]
# stopping thresholds, in deg/h for the rate metric and deg for the shift
# metric evaluated over SHIFT_WIN hours
STOP_GRID = [float(x) for x in
             os.environ.get("STOP_GRID", "0.3,1,2,3,4,5").split(",")]
# guidance azimuth-drift tolerances, deg/h. The published value is 1.5; the
# grid extends to 15, at which point a 2 h journey admits 30 deg of azimuth
# wander and the object can hardly be said to have "gone before" anyone.
DRIFT_GRID = [float(x) for x in
              os.environ.get("DRIFT_GRID", "1.5,3,5,10,15").split(",")]
SHIFT_WIN = float(os.environ.get("SHIFT_WIN", 1.0))    # hours


def great_circle(alt1, az1, alt2, az2):
    """Angular separation between two horizontal positions, in degrees."""
    a1, a2 = np.radians(alt1), np.radians(alt2)
    dz = np.radians(az1 - az2)
    c = np.sin(a1) * np.sin(a2) + np.cos(a1) * np.cos(a2) * np.cos(dz)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def track(q, e, i_r, Om_r, om_r, T_jd, sky):
    """Horizontal track of one orbit over the night window.

    Reproduces the geometry of sob_close_flyby_mc.evaluate_orbit exactly,
    including the precession of both comet and Sun to the equinox of date so
    that hour angles formed against the of-date sidereal time are frame
    consistent.  Returns (alt, az, night, d_min).
    """
    jd_utc, earth_pos, sun_ra, sun_dec, lst = sky
    jd_tdb = jd_utc + M.DT_TDB_UTC / 86400.0
    try:
        r_c = M.comet_pos_scalar(jd_tdb - T_jd, q, e, i_r, Om_r, om_r)
    except Exception:
        return None, None, None, np.inf

    rho = r_c - earth_pos
    d = np.linalg.norm(rho, axis=0)
    d_min = float(np.min(d))

    P = M.precession_matrix(float(np.mean(jd_tdb)))
    rhat = P @ (rho / np.maximum(d, 1e-12))
    dec = np.degrees(np.arcsin(np.clip(rhat[2], -1, 1)))
    ra = np.degrees(np.arctan2(rhat[1], rhat[0])) % 360.0
    alt, az = M.altaz_analytic(np.radians(ra), np.radians(dec), lst)

    sun_v = P @ np.array([np.cos(sun_dec) * np.cos(sun_ra),
                          np.cos(sun_dec) * np.sin(sun_ra),
                          np.sin(sun_dec)])
    s_ra = np.arctan2(sun_v[1], sun_v[0]) % (2 * np.pi)
    s_dec = np.arcsin(np.clip(sun_v[2], -1, 1))
    sun_alt, _ = M.altaz_analytic(s_ra, s_dec, lst)
    return alt, az, sun_alt < M.SUN_ALT_MAX, d_min


def sample_close_flybys(n_mc, sky, seed, batch=None):
    """Draw orbits on the same stratified priors as the published survey and
    return those achieving d_min < D_FLYBY.

    Sampling is batched: the element arrays are ~10 vectors of n_mc float64,
    so drawing 25 million at once exhausts memory. Batching leaves the priors
    and the result identical.
    """
    batch = batch or int(os.environ.get("SOB_BATCH", 4_000_000))
    if n_mc > batch:
        out = []
        done = 0
        k = 0
        while done < n_mc:
            n = min(batch, n_mc - done)
            out += sample_close_flybys(n, sky, seed + 1000 * k, batch=batch)
            done += n
            k += 1
            print(f"    ...{done:,}/{n_mc:,} sampled, {len(out)} flybys so far",
                  flush=True)
        return out

    jd_utc, earth_pos, _, _, _ = sky
    earth_mid = earth_pos[:, len(jd_utc) // 2]
    rng = np.random.default_rng(seed)

    n0, n1 = int(0.60 * n_mc), int(0.25 * n_mc)
    n2 = n_mc - n0 - n1
    q = np.concatenate([
        rng.uniform(0.80, 1.20, n0),
        np.where(rng.integers(0, 2, n1).astype(bool),
                 rng.uniform(0.50, 0.80, n1), rng.uniform(1.20, 1.40, n1)),
        rng.uniform(0.30, 0.50, n2)])
    e = np.concatenate([
        rng.uniform(0.70, 1.00, n_mc // 3),
        rng.uniform(1.00, 1.01, n_mc // 6),
        rng.uniform(1.01, 5.00, n_mc - n_mc // 3 - n_mc // 6)])
    m0, m1 = int(0.60 * n_mc), int(0.25 * n_mc)
    m2 = n_mc - m0 - m1
    i_r = np.concatenate([
        np.arccos(rng.uniform(np.cos(np.radians(20)), 1.0, m0)),
        np.arccos(rng.uniform(np.cos(np.radians(45)), np.cos(np.radians(20)), m1)),
        np.arccos(rng.uniform(-1.0, np.cos(np.radians(45)), m2))])
    Om = rng.uniform(0.0, 2 * np.pi, n_mc)
    om = rng.uniform(0.0, 2 * np.pi, n_mc)
    T = M.JD_PRIMARY + rng.uniform(-90.0, 90.0, n_mc)

    pos = M.comet_pos_parabolic_batch(M.JD_PRIMARY - T, q, i_r, Om, om)
    keep = np.linalg.norm(pos - earth_mid[:, None], axis=0) < M.D_REJECT
    del pos

    out = []
    idx = np.where(keep)[0]
    for j in idx:
        _, _, _, dm = track(q[j], e[j], i_r[j], Om[j], om[j], T[j], sky)
        if dm <= M.D_FLYBY:
            out.append((q[j], e[j], i_r[j], Om[j], om[j], T[j]))
    return out


def diagnose(q, e, i_r, Om_r, om_r, T_jd, sky):
    """Stopping diagnostics for one orbit over the (corridor, drift) grid.

    An earlier version of this sweep varied only the corridor and the
    stopping threshold, following the referee's request literally.  That was
    uninformative: no orbit reaches the stopping test at all, because the
    binding constraint is the GUIDANCE azimuth-drift tolerance, not the
    stopping threshold.  The drift tolerance is therefore swept as well.

    Sweeping it is warranted on the referee's own reasoning.  If what an
    observer can detect is a finite angular shift rather than a rate, then an
    azimuth drift of D deg/h sustained over a journey of GUIDE_DUR hours is
    "maintained" so long as D * GUIDE_DUR stays under the detectable shift.

    Returns {(az_lo, drift): (guidance_ok, omega_min, shift_min)}.
    """
    jd_utc = sky[0]
    step_h = (jd_utc[1] - jd_utc[0]) * 24.0

    alt, az, night, d_min = track(q, e, i_r, Om_r, om_r, T_jd, sky)
    if alt is None:
        return None, d_min

    dt_h = np.gradient(jd_utc) * 24.0
    dalt = np.gradient(alt) / dt_h
    az_u = np.degrees(np.unwrap(np.radians(az)))
    daz = np.gradient(az_u) / dt_h
    omega = np.sqrt(dalt**2 + (daz * np.cos(np.radians(alt)))**2)

    k = max(1, int(round(SHIFT_WIN / step_h)))
    shift = np.full_like(omega, np.inf)
    if len(alt) > k:
        shift[:-k] = great_circle(alt[:-k], az_u[:-k], alt[k:], az_u[k:])

    high = alt > M.ALT_MIN
    moving = np.abs(dalt) > M.ALT_RATE_MIN
    n_gap = int(np.ceil(M.STOP_GAP / step_h))

    out = {}
    for az_lo in AZ_LO_GRID:
        corridor = (az >= az_lo) & (az <= M.AZ_HI)
        base = night & high & corridor & moving
        for drift in DRIFT_GRID:
            guide = base & (np.abs(daz) < drift)
            if (not np.any(guide)) or float(np.sum(guide)) * step_h < M.GUIDE_DUR:
                out[(az_lo, drift)] = (False, np.inf, np.inf)
                continue
            gi = int(np.where(guide)[0][-1])
            sl = slice(gi, min(gi + n_gap + 1, len(omega)))
            nm = night[sl]
            if not np.any(nm):
                out[(az_lo, drift)] = (False, np.inf, np.inf)
                continue
            om_w = omega[sl][nm]
            sh_w = shift[sl][nm]
            sh_w = sh_w[np.isfinite(sh_w)]
            out[(az_lo, drift)] = (True, float(np.min(om_w)),
                                   float(np.min(sh_w)) if sh_w.size else np.inf)
    return out, d_min


def main():
    print("=" * 78)
    print("CRITERIA SENSITIVITY SWEEP — response to Reviewer 2, round 3")
    print("=" * 78)
    print(f"  N_MC              = {N_MC:,}")
    print(f"  close-flyby test  : d_min < {M.D_FLYBY} au")
    print(f"  corridor upper    : {M.AZ_HI:.0f} deg (fixed)")
    print(f"  corridor lower    : {AZ_LO_GRID} deg (swept)")
    print(f"  stopping grid     : {STOP_GRID}")
    print(f"  shift window      : {SHIFT_WIN:.2f} h")
    print(f"  guidance duration : >= {M.GUIDE_DUR} h, drift < {M.AZ_DRIFT_THRESH} deg/h")
    print()

    sky = M.precompute_sky(M.JD_PRIMARY)
    cands = sample_close_flybys(N_MC, sky, seed=SEED)
    print(f"  close flybys found: {len(cands):,}\n")
    if not len(cands):
        print("  no candidates; increase N_MC")
        return

    rows = {kk: {"n_guide": 0, "omega": [], "shift": []}
            for kk in ((a, d) for a in AZ_LO_GRID for d in DRIFT_GRID)}
    for c in cands:
        res, _ = diagnose(*c, sky)
        if res is None:
            continue
        for kk, (ok, om, sh) in res.items():
            if ok:
                rows[kk]["n_guide"] += 1
                rows[kk]["omega"].append(om)
                rows[kk]["shift"].append(sh)

    print("=" * 78)
    print("RESULT (1) — how many orbits reach the stopping test at all")
    print("=" * 78)
    hdr = "  corridor  " + "".join(f"{d:>8g}" for d in DRIFT_GRID)
    print("            guidance azimuth-drift tolerance, deg/h")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for a in AZ_LO_GRID:
        cells = "".join(f"{rows[(a, d)]['n_guide']:>8d}" for d in DRIFT_GRID)
        print(f"  {a:>5.0f}-{M.AZ_HI:<4.0f}{cells}")

    print("\n" + "=" * 78)
    print("RESULT (2) — orbits passing every criterion")
    print("=" * 78)
    for metric, key, unit in (("rate", "omega", "deg/h"),
                              (f"shift over {SHIFT_WIN:g} h", "shift", "deg")):
        print(f"\n  Stopping metric: {metric}")
        hdr = ("  corridor  drift " + "".join(f"{t:>8g}" for t in STOP_GRID))
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for a in AZ_LO_GRID:
            for d in DRIFT_GRID:
                v = np.array(rows[(a, d)][key])
                if not rows[(a, d)]["n_guide"]:
                    continue
                cells = "".join(f"{int(np.sum(v <= t)):>8d}" for t in STOP_GRID)
                print(f"  {a:>5.0f}-{M.AZ_HI:<4.0f}{d:>6g}{cells}")
        allv = np.concatenate([np.array(rows[kk][key]) for kk in rows
                               if len(rows[kk][key])]) \
            if any(len(rows[kk][key]) for kk in rows) else np.array([])
        best = float(np.min(allv)) if allv.size else np.inf
        print(f"\n    tightest value achieved by any orbit: {best:.3f} {unit}")
        if np.isfinite(best):
            print(f"    -> conclusion holds for every stopping criterion "
                  f"below {best:.2f} {unit}")

    print(f"""
==============================================================================
READING THIS
==============================================================================
  Each cell counts orbits satisfying the close-flyby test, the guidance test
  in that corridor, and the stopping test at that threshold.

  The number the referee asked for is the last line of each block: the
  tightest stopping criterion that any orbit in the sample achieves.  If it
  exceeds 5 deg/h (or 5 deg over {SHIFT_WIN:g} h), then no choice within the
  range he proposes changes the conclusion, and the disagreement over where
  exactly to set the threshold is immaterial.

  If instead some cell is non-zero, the table shows precisely which
  combination of corridor and threshold admits an orbit, which is the honest
  form of the result and is what should be reported.""")


if __name__ == "__main__":
    main()
