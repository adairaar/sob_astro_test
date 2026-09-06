#!/usr/bin/env python3
"""
sob_precession.py — shared J2000 → equinox-of-date precession.

Why this module exists
----------------------
Comet and Earth positions in this project are built in the J2000 (ICRS)
frame, but Greenwich Mean Sidereal Time is referred to the mean equinox OF
DATE.  Forming an hour angle from a J2000 right ascension and an of-date
sidereal time silently omits the general precession accumulated since the
epoch of interest — roughly 28° for 5 BCE, which propagates to about 23° of
error in azimuth and 13° in altitude.

Every script that converts a J2000 direction into alt/az must precess first.

Sign convention
---------------
Composition follows Meeus, *Astronomical Algorithms*, ch. 21 (rigorous method):

    P = Rz(z) · Ry(−θ) · Rz(ζ)

This is NOT the Rz(−z)·Ry(θ)·Rz(−ζ) form quoted in many references; that form
describes rotation of the coordinate axes rather than of the vector, and using
it here yields the exact inverse rotation.  Because ζ, z and θ are themselves
negative for epochs before J2000, the mistake does not appear as an obvious
sign flip — it appears as a doubled offset (~60° in RA at 5 BCE).

Validated against independent Stellarium output in sob_frames.validate().
"""

import numpy as np

_D2R = np.pi / 180.0


def _rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _ry(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def precession_matrix(jd_tt):
    """Mean J2000 equator/equinox -> mean equator/equinox of date. jd_tt scalar."""
    T = (jd_tt - 2451545.0) / 36525.0
    zeta = (2306.2181 * T + 0.30188 * T**2 + 0.017998 * T**3) / 3600.0 * _D2R
    z = (2306.2181 * T + 1.09468 * T**2 + 0.018203 * T**3) / 3600.0 * _D2R
    theta = (2004.3109 * T - 0.42665 * T**2 - 0.041833 * T**3) / 3600.0 * _D2R
    return _rz(z) @ _ry(-theta) @ _rz(zeta)


def precess_radec(ra_r, dec_r, jd_tt):
    """
    Precess J2000 (ra, dec) in radians to the mean equinox of date.
    Accepts scalars or arrays; returns (ra_date, dec_date) in radians.
    """
    P = precession_matrix(float(np.mean(jd_tt)))
    v = np.array([np.cos(dec_r) * np.cos(ra_r),
                  np.cos(dec_r) * np.sin(ra_r),
                  np.sin(dec_r)])
    w = P @ v
    n = np.linalg.norm(w, axis=0)
    return (np.arctan2(w[1], w[0]) % (2 * np.pi),
            np.arcsin(np.clip(w[2] / np.maximum(n, 1e-15), -1.0, 1.0)))


def precess_unit_vectors(rhat, jd_tt):
    """Precess an array of J2000 unit vectors, shape (3,) or (3,N), to of date."""
    return precession_matrix(float(np.mean(jd_tt))) @ rhat
