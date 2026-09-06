#!/usr/bin/env python3
"""
sob_frames.py — verified coordinate-frame utilities for the Star of Bethlehem analysis.

Written after two frame bugs were identified in earlier code:

  Bug A  A low-precision Sun-longitude series (280.460 + 0.9856474 n) returns
         ecliptic longitude referred to the MEAN EQUINOX OF DATE, but its output
         was consumed by a J2000 pipeline.  At 5 BCE this misplaces Earth by the
         accumulated general precession, ~28 deg.

  Bug B  Topocentric alt/az was formed by combining J2000 RA/Dec with a
         Greenwich Mean Sidereal Time referred to the equinox of date.  The two
         are different frames; at 5 BCE the mismatch is ~23 deg in azimuth.

Both are avoided here by keeping an explicit frame label on every quantity and
converting deliberately:

    comet elements (J2000 ecliptic)  ->  J2000 ecliptic rectangular
    Earth VSOP87                     ->  J2000 ecliptic rectangular
    difference                       ->  geocentric J2000 ecliptic
    obliquity rotation               ->  geocentric J2000 equatorial
    topocentric parallax correction  ->  topocentric J2000 equatorial
    IAU 1976 precession              ->  topocentric equinox-of-date equatorial
    GMST(UT1) + longitude            ->  local hour angle  ->  alt/az

Validated against reviewer-supplied Stellarium positions; see validate() below.

Time systems: comet perihelion epochs are TT.  Sidereal time requires UT1.
TT = UT1 + DeltaT; at 5 BCE DeltaT ~ 0.12236 d (the value Stellarium uses).
"""

import numpy as np

try:
    from pymeeus.Earth import Earth as _MeeusEarth
    from pymeeus.Epoch import Epoch as _MeeusEpoch
    _HAVE_VSOP = True
except ImportError:                                     # pragma: no cover
    _HAVE_VSOP = False

D2R = np.pi / 180.0
R2D = 180.0 / np.pi
GM_AU_DAY2 = 2.959122082855911e-4      # heliocentric gravitational parameter, AU^3/day^2
EPS_J2000 = 23.4392911 * D2R           # mean obliquity at J2000
AU_KM = 1.495978707e8
R_EARTH_KM = 6378.137                  # equatorial radius, WGS84
R_EARTH_AU = R_EARTH_KM / AU_KM
F_FLAT = 1.0 / 298.257223563           # flattening

DELTA_T_5BCE = 0.12236                 # days; TT = UT1 + DeltaT (Stellarium convention)


# ────────────────────────────────────────────────────────────────────────────
#  Rotations
# ────────────────────────────────────────────────────────────────────────────

def _Rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _Ry(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def ecl_to_equ(v, eps=EPS_J2000):
    """Ecliptic -> equatorial, same equinox. v shape (3,) or (3,N)."""
    R = np.array([[1.0, 0.0, 0.0],
                  [0.0, np.cos(eps), -np.sin(eps)],
                  [0.0, np.sin(eps),  np.cos(eps)]])
    return R @ v


def precession_matrix(jd_tt):
    """
    IAU 1976 precession, mean J2000 equator/equinox -> mean equator/equinox of date.

    Composition follows Meeus, *Astronomical Algorithms* ch. 21 (rigorous method):

        RA_date = atan2(A, B) + z,   with the intermediate rotation by +zeta
                                     about z, then -theta about y, then +z about z

    i.e.  P = Rz(z) . Ry(-theta) . Rz(zeta).

    NOTE ON SIGN CONVENTION.  Many texts write this as Rz(-z) Ry(theta) Rz(-zeta);
    that form assumes rotation OF THE AXES rather than of the vector, and using it
    here produced the exact inverse rotation.  Because T < 0 for antiquity the
    three angles are themselves negative, so the error is not visible as an
    obvious sign flip -- it shows up as a ~2x precession offset (~60 deg in RA at
    5 BCE).  The form below is validated against independent Stellarium output in
    validate(); do not "simplify" it without re-running that check.

    Returns a single (3,3) matrix; jd_tt must be scalar.  Over a single night the
    change is negligible, so evaluating once at mid-window is standard practice.
    """
    T = (jd_tt - 2451545.0) / 36525.0
    zeta = (2306.2181 * T + 0.30188 * T**2 + 0.017998 * T**3) / 3600.0 * D2R
    z = (2306.2181 * T + 1.09468 * T**2 + 0.018203 * T**3) / 3600.0 * D2R
    theta = (2004.3109 * T - 0.42665 * T**2 - 0.041833 * T**3) / 3600.0 * D2R
    return _Rz(z) @ _Ry(-theta) @ _Rz(zeta)


# ────────────────────────────────────────────────────────────────────────────
#  Earth position
# ────────────────────────────────────────────────────────────────────────────

def earth_helio_ecl_j2000(jd_tt):
    """
    Earth heliocentric rectangular position, J2000 ECLIPTIC frame, in AU.

    Uses VSOP87 (via pymeeus), which is accurate over roughly +/-4000 yr from
    J2000 -- unlike the truncated series that caused Bug A, and unlike the
    bundled DE kernels, which do not cover antiquity.

    Accepts scalar or array jd_tt; returns (3,) or (3,N).
    """
    if not _HAVE_VSOP:
        raise RuntimeError(
            "pymeeus is required for VSOP87 Earth positions "
            "(pip install pymeeus). Do NOT substitute a low-order series: "
            "see the Bug A note in this module's docstring."
        )
    scalar = np.isscalar(jd_tt) or np.asarray(jd_tt).ndim == 0
    jds = np.atleast_1d(np.asarray(jd_tt, dtype=float))
    out = np.empty((3, jds.size))
    for k, jd in enumerate(jds):
        l, b, r = _MeeusEarth.geometric_heliocentric_position_j2000(_MeeusEpoch(float(jd)))
        L = float(l) * D2R
        B = float(b) * D2R
        R = float(r)
        out[:, k] = (R * np.cos(B) * np.cos(L),
                     R * np.cos(B) * np.sin(L),
                     R * np.sin(B))
    return out[:, 0] if scalar else out


# ────────────────────────────────────────────────────────────────────────────
#  Comet position
# ────────────────────────────────────────────────────────────────────────────

def _kepler_elliptic(M, e, tol=1e-13, itmax=200):
    M = np.atleast_1d(np.asarray(M, dtype=float))
    E = np.where(e < 0.8, M, np.pi * np.ones_like(M))
    for _ in range(itmax):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))
        E = E + dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def _kepler_hyperbolic(Mh, e, tol=1e-13, itmax=300):
    Mh = np.atleast_1d(np.asarray(Mh, dtype=float))
    H = np.arcsinh(Mh / max(e, 1.0 + 1e-9))
    for _ in range(itmax):
        f = e * np.sinh(H) - H - Mh
        dH = -f / (e * np.cosh(H) - 1.0)
        H = H + dH
        if np.max(np.abs(dH)) < tol:
            break
    return H


def _barker(dt, q):
    """Parabolic (Barker's equation) -> tan(nu/2)."""
    A = 1.5 * np.sqrt(GM_AU_DAY2 / (2.0 * q**3)) * dt
    B = np.cbrt(A + np.sqrt(A * A + 1.0))
    return B - 1.0 / B


def comet_helio_ecl_j2000(jd_tt, q, e, i_deg, Om_deg, om_deg, T_peri_tt):
    """
    Comet heliocentric rectangular position, J2000 ECLIPTIC frame, in AU.

    Orbital elements are assumed referred to the J2000 ecliptic and equinox,
    which is the convention Stellarium's comet .ini entries use.
    """
    jds = np.atleast_1d(np.asarray(jd_tt, dtype=float))
    dt = jds - T_peri_tt

    i_r, Om_r, om_r = i_deg * D2R, Om_deg * D2R, om_deg * D2R
    ci, si = np.cos(i_r), np.sin(i_r)
    cO, sO = np.cos(Om_r), np.sin(Om_r)
    co, so = np.cos(om_r), np.sin(om_r)

    # Perifocal basis expressed in the ecliptic frame
    P = np.array([cO * co - sO * so * ci,
                  sO * co + cO * so * ci,
                  so * si])
    Q = np.array([-cO * so - sO * co * ci,
                  -sO * so + cO * co * ci,
                  co * si])

    if abs(e - 1.0) < 1e-9:
        tanv2 = _barker(dt, q)
        nu = 2.0 * np.arctan(tanv2)
        r = q * (1.0 + tanv2**2)
    elif e < 1.0:
        a = q / (1.0 - e)
        n = np.sqrt(GM_AU_DAY2 / a**3)
        E = _kepler_elliptic(n * dt, e)
        nu = 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(E / 2.0),
                              np.sqrt(1.0 - e) * np.cos(E / 2.0))
        r = a * (1.0 - e * np.cos(E))
    else:
        a = q / (1.0 - e)                      # negative
        n = np.sqrt(GM_AU_DAY2 / (-a)**3)
        H = _kepler_hyperbolic(n * dt, e)
        nu = 2.0 * np.arctan2(np.sqrt(e + 1.0) * np.sinh(H / 2.0),
                              np.sqrt(e - 1.0) * np.cosh(H / 2.0))
        r = a * (1.0 - e * np.cosh(H))

    pos = np.outer(P, r * np.cos(nu)) + np.outer(Q, r * np.sin(nu))
    return pos[:, 0] if pos.shape[1] == 1 else pos


# ────────────────────────────────────────────────────────────────────────────
#  Sidereal time and the observer
# ────────────────────────────────────────────────────────────────────────────

def gmst_rad(jd_ut1):
    """
    Greenwich Mean Sidereal Time, radians, referred to the EQUINOX OF DATE.

    Because the result is of-date, it must only ever be combined with of-date
    right ascension.  Pairing it with J2000 RA was Bug B.
    """
    T = (jd_ut1 - 2451545.0) / 36525.0
    deg = (280.46061837
           + 360.98564736629 * (jd_ut1 - 2451545.0)
           + 0.000387933 * T**2
           - T**3 / 38710000.0) % 360.0
    return deg * D2R


def observer_geocentric_vector(jd_ut1, lat_deg, lon_deg_east, height_m=0.0):
    """
    Observer's position relative to Earth's centre, in the equatorial frame
    OF DATE, in AU.  Needed to convert geocentric to topocentric coordinates --
    a correction of up to ~1.5 deg for an object at lunar distance, and therefore
    not optional when checking close-flyby geometry.
    """
    lat = lat_deg * D2R
    # geodetic -> geocentric parameters (Meeus ch. 11)
    u_ = np.arctan((1.0 - F_FLAT) * np.tan(lat))
    h = height_m / 1000.0 / AU_KM
    rho_sin = (1.0 - F_FLAT) * np.sin(u_) + h * np.sin(lat)
    rho_cos = np.cos(u_) + h * np.cos(lat)
    lst = np.atleast_1d(gmst_rad(np.asarray(jd_ut1, dtype=float))
                        + lon_deg_east * D2R)
    return np.array([R_EARTH_AU * rho_cos * np.cos(lst),
                     R_EARTH_AU * rho_cos * np.sin(lst),
                     R_EARTH_AU * rho_sin * np.ones_like(lst)])


def altaz_from_of_date(ra_date, dec_date, jd_ut1, lat_deg, lon_deg_east):
    """
    Alt/az from OF-DATE equatorial coordinates (radians) and UT1.
    Azimuth measured from North through East.
    """
    lat = lat_deg * D2R
    lst = (gmst_rad(jd_ut1) + lon_deg_east * D2R) % (2.0 * np.pi)
    ha = lst - ra_date
    sin_alt = np.clip(np.sin(lat) * np.sin(dec_date)
                      + np.cos(lat) * np.cos(dec_date) * np.cos(ha), -1.0, 1.0)
    alt = np.arcsin(sin_alt)
    az_y = -np.cos(dec_date) * np.sin(ha)
    az_x = (np.sin(dec_date) - np.sin(lat) * sin_alt) / np.maximum(np.cos(lat), 1e-12)
    az = np.arctan2(az_y, az_x) % (2.0 * np.pi)
    return alt * R2D, az * R2D


# ────────────────────────────────────────────────────────────────────────────
#  Full pipeline
# ────────────────────────────────────────────────────────────────────────────

def observe(jd_ut1, elements, lat_deg, lon_deg_east, height_m=0.0,
            delta_t_days=DELTA_T_5BCE, topocentric=True):
    """
    Full geometry for a comet as seen by a ground observer.

    elements: dict with q, e, i_deg, Om_deg, om_deg, T_peri_tt  (J2000 ecliptic).
    jd_ut1:   scalar or array of Julian dates in UT1.

    Returns dict of arrays:
        dist_au, dist_km   topocentric (or geocentric) distance
        ra_j2000, dec_j2000    radians, J2000 equatorial
        ra_date,  dec_date     radians, equinox of date
        alt_deg,  az_deg
    """
    jd_ut1 = np.atleast_1d(np.asarray(jd_ut1, dtype=float))
    jd_tt = jd_ut1 + delta_t_days

    comet = comet_helio_ecl_j2000(jd_tt, **elements)
    earth = earth_helio_ecl_j2000(jd_tt)
    comet = np.atleast_2d(comet.T).T if comet.ndim == 1 else comet
    earth = np.atleast_2d(earth.T).T if earth.ndim == 1 else earth

    geo_ecl = comet - earth                       # geocentric, J2000 ecliptic
    geo_equ = ecl_to_equ(geo_ecl)                 # geocentric, J2000 equatorial

    P = precession_matrix(float(np.mean(jd_tt)))
    geo_date = P @ geo_equ                        # geocentric, of-date equatorial

    d_geo = np.linalg.norm(geo_date, axis=0)

    # Topocentric shift: matters at ~1 deg for objects near lunar distance, and
    # is what a ground observer (and Stellarium's alt/az readout) actually sees.
    if topocentric:
        obs = observer_geocentric_vector(jd_ut1, lat_deg, lon_deg_east, height_m)
        topo_date = geo_date - obs
    else:
        topo_date = geo_date

    d_topo = np.linalg.norm(topo_date, axis=0)
    ra_date = np.arctan2(topo_date[1], topo_date[0]) % (2.0 * np.pi)
    dec_date = np.arcsin(np.clip(topo_date[2] / np.maximum(d_topo, 1e-15), -1.0, 1.0))
    alt, az = altaz_from_of_date(ra_date, dec_date, jd_ut1, lat_deg, lon_deg_east)

    # J2000 equivalents of the geocentric direction (for comparison with
    # catalogues and with reviewer-reported J2000 positions)
    ra_j2000 = np.arctan2(geo_equ[1], geo_equ[0]) % (2.0 * np.pi)
    dec_j2000 = np.arcsin(np.clip(geo_equ[2] / np.maximum(d_geo, 1e-15), -1.0, 1.0))

    return {
        'jd_ut1': jd_ut1, 'jd_tt': jd_tt,
        'dist_au': d_geo, 'dist_km': d_geo * AU_KM,          # geocentric
        'dist_topo_au': d_topo, 'dist_topo_km': d_topo * AU_KM,
        'ra_date': ra_date, 'dec_date': dec_date,            # topocentric, of date
        'ra_j2000': ra_j2000, 'dec_j2000': dec_j2000,        # geocentric, J2000
        'alt_deg': alt, 'az_deg': az,
    }


def angular_rate_deg_per_hr(ra, dec, jd):
    """
    Great-circle rate between consecutive samples, deg/hr.
    Feed it of-date or J2000 RA/Dec to get the INERTIAL rate (motion against the
    stars); feed it alt/az to get the GROUND-FRAME rate (which includes diurnal
    rotation).  These are different quantities and the paper uses both --
    keeping them distinct matters.
    """
    c = (np.sin(dec[:-1]) * np.sin(dec[1:])
         + np.cos(dec[:-1]) * np.cos(dec[1:]) * np.cos(ra[1:] - ra[:-1]))
    sep = np.arccos(np.clip(c, -1.0, 1.0)) * R2D
    dt = np.diff(jd) * 24.0
    return sep / dt


# ────────────────────────────────────────────────────────────────────────────
#  Validation against reviewer-supplied Stellarium output
# ────────────────────────────────────────────────────────────────────────────

R2_ORBIT_1 = dict(q=1.00284101, e=0.8070114, i_deg=10.9963,
                  Om_deg=56.1727, om_deg=193.7693, T_peri_tt=1719718.38394)

R2_ORBIT_2 = dict(q=0.9008354, e=0.9668455, i_deg=10.9153,
                  Om_deg=56.1658, om_deg=220.6361, T_peri_tt=1719735.20246)

BETHLEHEM = dict(lat_deg=31.7, lon_deg_east=35.2)

# Reviewer's Stellarium check values for Orbit 2, 22 April 5 BCE
_R2_CHECK = [
    # jd_ut1,        az,            el,           dist_km, ra_date_h,   dec_date
    (1719708.56888, 211 + 43 / 60, 36 + 52 / 60, 271000, 16 + 8 / 60,  -(15 + 19 / 60)),
    (1719708.65221, 207 + 28 / 60, 22 + 51 / 60, 248000, 17 + 55 / 60, -(29 + 29 / 60)),
    (1719708.73554, 203 + 49 / 60, 16 + 35 / 60, 290000, 19 + 57 / 60, -(36 + 38 / 60)),
]


def validate(verbose=True):
    """Compare this module against the reviewer's independent Stellarium run."""
    rows = []
    for jd, az_r2, el_r2, d_r2, ra_r2_h, dec_r2 in _R2_CHECK:
        o = observe(jd, R2_ORBIT_2, **BETHLEHEM)
        rows.append(dict(
            jd=jd,
            d_ours=float(o['dist_km'][0]), d_r2=d_r2,
            az_ours=float(o['az_deg'][0]), az_r2=az_r2,
            el_ours=float(o['alt_deg'][0]), el_r2=el_r2,
            ra_ours_h=float(o['ra_date'][0]) * R2D / 15.0, ra_r2_h=ra_r2_h,
            dec_ours=float(o['dec_date'][0]) * R2D, dec_r2=dec_r2,
        ))
    if verbose:
        print("=" * 78)
        print("VALIDATION vs reviewer Stellarium output — Orbit 2, 22 Apr 5 BCE")
        print("=" * 78)
        for r in rows:
            print(f"\n  JD(UT1) {r['jd']:.5f}")
            print(f"    distance  ours {r['d_ours']:>11,.0f} km   ref {r['d_r2']:>9,} km"
                  f"   ratio {r['d_ours']/r['d_r2']:.3f}")
            print(f"    RA(date)  ours {r['ra_ours_h']:>11.3f} h    ref {r['ra_r2_h']:>9.3f} h"
                  f"   d {(r['ra_ours_h']-r['ra_r2_h'])*15:+.2f} deg")
            print(f"    Dec(date) ours {r['dec_ours']:>11.2f} d    ref {r['dec_r2']:>9.2f} d"
                  f"   d {r['dec_ours']-r['dec_r2']:+.2f} deg")
            print(f"    azimuth   ours {r['az_ours']:>11.2f} d    ref {r['az_r2']:>9.2f} d"
                  f"   d {r['az_ours']-r['az_r2']:+.2f} deg")
            print(f"    altitude  ours {r['el_ours']:>11.2f} d    ref {r['el_r2']:>9.2f} d"
                  f"   d {r['el_ours']-r['el_r2']:+.2f} deg")
    return rows


if __name__ == "__main__":
    validate()
