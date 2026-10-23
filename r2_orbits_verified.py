#!/usr/bin/env python3
"""
r2_orbits_verified.py — reviewer-proposed orbits, recomputed with verified frames.

Supersedes r2_reproduce.py, which contained a low-precision Sun series that
returned mean-equinox-of-date longitude into a J2000 pipeline, misplacing Earth
by ~28 deg.  Every number this script prints replaces the corresponding figure
from that script.

Uses sob_frames.py, validated against the reviewer's own independent Stellarium
output (see sob_frames.validate()).

Two distinct rate quantities are reported, because the manuscript uses both and
they are NOT interchangeable:

  inertial rate  — motion against the background stars.  This is what "the star
                   stood still" refers to for a star-like object, and what the
                   lunar-rate bound (~0.5 deg/hr) constrains.

  ground rate    — motion relative to the observer's horizon, i.e. including
                   diurnal rotation.  Standing still in this frame requires
                   cancelling Earth's ~15 deg/hr spin, which only a very close,
                   fast object could do.
"""

import numpy as np
import sob_frames as sf

LAT, LON = 31.7, 35.2
JD_MIDNIGHT_UT = 1719708.5          # 22 April 5 BCE (Julian), 00:00 UT
STOP_THRESH_INERTIAL = 0.5          # deg/hr, lunar-rate upper bound
STOP_THRESH_GROUND = 0.30           # deg/hr, apparent rest w.r.t. horizon


def sun_altitude(jd_ut1):
    """Geocentric Sun altitude at Bethlehem. Sun direction = -Earth heliocentric."""
    jd_tt = jd_ut1 + sf.DELTA_T_5BCE
    earth = sf.earth_helio_ecl_j2000(jd_tt)
    sun_ecl = -np.atleast_2d(earth.T).T if earth.ndim == 1 else -earth
    sun_equ = sf.ecl_to_equ(sun_ecl)
    P = sf.precession_matrix(float(np.mean(jd_tt)))
    sd = P @ sun_equ
    d = np.linalg.norm(sd, axis=0)
    ra = np.arctan2(sd[1], sd[0]) % (2 * np.pi)
    dec = np.arcsin(np.clip(sd[2] / d, -1, 1))
    alt, _ = sf.altaz_from_of_date(ra, dec, np.atleast_1d(jd_ut1), LAT, LON)
    return alt


def find_crossing(jd, values, target, rising):
    for k in range(len(jd) - 1):
        a, b = values[k], values[k + 1]
        if (rising and a < target <= b) or (not rising and a > target >= b):
            f = (target - a) / (b - a)
            return jd[k] + f * (jd[k + 1] - jd[k])
    return None


def local_hours(jd):
    """Hours since local midnight (local mean solar time at LON)."""
    return (jd - (JD_MIDNIGHT_UT - LON / 15.0 / 24.0)) * 24.0


def analyse(name, elements):
    print("=" * 78)
    print(f"{name}")
    print("=" * 78)
    print(f"  e={elements['e']}  q={elements['q']} AU  i={elements['i_deg']}"
          f"  Om={elements['Om_deg']}  om={elements['om_deg']}"
          f"  T_peri(TT)={elements['T_peri_tt']}")

    # ---- closest approach over +/- 60 d ----
    jd_c = JD_MIDNIGHT_UT + np.arange(-60, 60, 0.005)
    o_c = sf.observe(jd_c, elements, lat_deg=LAT, lon_deg_east=LON, topocentric=False)
    k = int(np.argmin(o_c['dist_km']))
    print(f"\n  Closest approach : {o_c['dist_km'][k]:,.0f} km "
          f"({o_c['dist_au'][k]:.6f} AU) at {jd_c[k]-JD_MIDNIGHT_UT:+.2f} d from 22 Apr")
    print(f"                     = {o_c['dist_km'][k]/384400:.2f} lunar distances")

    # ---- the night of 22 April ----
    jd_n = JD_MIDNIGHT_UT - 0.4 + np.arange(0, 1.4, 1.0 / 1440)   # 1-min steps
    sun_alt = sun_altitude(jd_n)
    sunset = find_crossing(jd_n, sun_alt, -0.833, rising=False)
    sunrise = find_crossing(jd_n, sun_alt, -0.833, rising=True)

    o = sf.observe(jd_n, elements, lat_deg=LAT, lon_deg_east=LON, topocentric=True)
    alt, az = o['alt_deg'], o['az_deg']

    inert = sf.angular_rate_deg_per_hr(o['ra_date'], o['dec_date'], jd_n)
    altr = np.radians(alt); azr = np.radians(az)
    ground = sf.angular_rate_deg_per_hr(azr, altr, jd_n)
    jd_m = 0.5 * (jd_n[:-1] + jd_n[1:])
    alt_m = 0.5 * (alt[:-1] + alt[1:])
    az_m = 0.5 * (az[:-1] + az[1:])
    sun_m = 0.5 * (sun_alt[:-1] + sun_alt[1:])

    vis = (alt_m > 0) & (sun_m < -6)      # above horizon, sky dark enough
    print(f"\n  Sunset  {local_hours(sunset):5.2f} h local   "
          f"Sunrise {local_hours(sunrise):5.2f} h local")

    if not np.any(vis):
        print("\n  Comet NOT visible at any point during darkness on this night.")
        return

    print(f"  Visible (dark & above horizon): "
          f"{local_hours(jd_m[vis][0]):.2f} h to {local_hours(jd_m[vis][-1]):.2f} h local"
          f"  ({vis.sum()/60.0:.2f} h)")

    print(f"\n  {'local h':>8} {'alt':>7} {'az':>8} {'dist km':>12} "
          f"{'inertial':>10} {'ground':>9}")
    idx = np.where(vis)[0]
    for j in idx[:: max(1, len(idx) // 12)]:
        print(f"  {local_hours(jd_m[j]):8.2f} {alt_m[j]:7.2f} {az_m[j]:8.2f} "
              f"{o['dist_km'][j]:12,.0f} {inert[j]:10.2f} {ground[j]:9.2f}")

    iv, gv = inert[vis], ground[vis]
    print(f"\n  Inertial rate while visible: min {iv.min():.3f}  max {iv.max():.2f} deg/hr")
    print(f"  Ground   rate while visible: min {gv.min():.3f}  max {gv.max():.2f} deg/hr")

    # ---- criteria ----
    print("\n  CRITERIA")
    stop_i = iv.min() < STOP_THRESH_INERTIAL
    print(f"    stopping (inertial < {STOP_THRESH_INERTIAL} deg/hr): "
          f"{'reached' if stop_i else 'NEVER'} "
          f"(min {iv.min():.3f}, i.e. {iv.min()/STOP_THRESH_INERTIAL:.1f}x threshold)")
    stop_g = gv.min() < STOP_THRESH_GROUND
    print(f"    stopping (ground   < {STOP_THRESH_GROUND} deg/hr): "
          f"{'reached' if stop_g else 'NEVER'} "
          f"(min {gv.min():.3f}, i.e. {gv.min()/STOP_THRESH_GROUND:.1f}x threshold)")

    # change of state: was it ever fast, then slow?
    if stop_i:
        first_slow = np.where(iv < STOP_THRESH_INERTIAL)[0][0]
        was_fast = np.any(iv[:first_slow] > STOP_THRESH_INERTIAL)
        print(f"    change of state (aorist): "
              f"{'YES' if was_fast else 'NO - already slow at first visibility'}")

    corridor = (az_m[vis] >= 190) & (az_m[vis] <= 210)
    print(f"    azimuth in corridor [190,210]: "
          f"{'yes, ' + format(corridor.sum()/60.0, '.2f') + ' h' if corridor.any() else 'never'}"
          f"   (range {az_m[vis].min():.1f} to {az_m[vis].max():.1f} deg)")
    print()


if __name__ == "__main__":
    sf.validate()
    print()
    analyse("ORBIT 1  (reviewer-proposed)", sf.R2_ORBIT_1)
    analyse("ORBIT 2  (reviewer-proposed, corrected elements)", sf.R2_ORBIT_2)
