#!/usr/bin/env python3
"""
comet_orbit.py

Determine orbital elements of a comet from Earth-bound astrometric
observations (RA, Dec, Julian Date).

Two modes
---------
Parabolic (default, e = 1 fixed):
    Five free parameters: q, i, Ω, ω, T.

Free-eccentricity (--fit-eccentricity):
    Six free parameters: q, e, i, Ω, ω, T.
    Covers elliptic (0 ≤ e < 1), parabolic (e = 1), and hyperbolic (e > 1).
    Orbit type is labelled automatically from the fitted value.

Algorithm
---------
1. Convert each (RA, Dec) observation to a geocentric unit direction vector.
2. Retrieve Earth's heliocentric equatorial (ICRS J2000) position from the
   JPL ephemeris via astropy at each epoch.
3. Propagate the comet with the appropriate Kepler equation:
     Elliptic   : M = E − e·sin E  (Newton-Raphson)
     Parabolic  : Barker's equation (closed-form Cardano root)
     Hyperbolic : M = e·sinh F − F  (Newton-Raphson)
4. Rotate perifocal → ecliptic → equatorial, subtract Earth's position,
   compare predicted sky direction with observed (RA, Dec).
5. Minimise sum-of-squared angular residuals with global differential
   evolution + Nelder-Mead polish.
6. Accept if RMS residual < 60 arcsec; otherwise report no good fit.

Orbital elements returned (ecliptic J2000):
    e     - eccentricity  (fixed 1.0 unless --fit-eccentricity)
    q     - perihelion distance (AU)
    i     - inclination (deg)
    Omega - longitude of ascending node (deg)
    omega - argument of perihelion (deg)
    T     - time of perihelion passage (JD TDB)

Minimum observations
--------------------
    Parabolic  (e fixed)  : 3 recommended (2 accepted with warning)
    Free e                : 4 recommended (3 accepted with warning)

Requirements
------------
    pip install astropy scipy numpy
"""

import argparse
import warnings
from typing import List, Optional, Tuple

import numpy as np
from scipy.optimize import differential_evolution, minimize

try:
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric
    import astropy.units as u
except ImportError as exc:
    raise SystemExit("astropy is required: pip install astropy\n" + str(exc))

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

K_GAUSS = 0.01720209895          # Gaussian gravitational constant (AU^3/2 day^-1)
MU      = K_GAUSS ** 2           # GM_sun  [AU^3 day^-2]

EPS_J2000 = np.radians(23.439291111)   # obliquity of ecliptic at J2000 (IAU 1976)

_ce, _se = np.cos(EPS_J2000), np.sin(EPS_J2000)
R_ECL_TO_EQU = np.array([
    [1.0,  0.0,   0.0],
    [0.0,  _ce,  -_se],
    [0.0,  _se,   _ce],
])

# Eccentricity tolerance for parabolic branch
_E_PARA_TOL = 1e-4

# Fit quality threshold
RMS_THRESHOLD_ARCSEC = 60.0

# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------

def ra_dec_to_unit_vector(ra_deg: float, dec_deg: float) -> np.ndarray:
    """Geocentric ICRS unit vector for (RA, Dec) in degrees."""
    ra  = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.array([np.cos(dec) * np.cos(ra),
                     np.cos(dec) * np.sin(ra),
                     np.sin(dec)])


def get_earth_heliocentric(jd: float) -> np.ndarray:
    """Heliocentric ICRS J2000 position of Earth (AU) at JD (TDB)."""
    t     = Time(jd, format="jd", scale="tdb")
    earth = get_body_barycentric("earth", t)
    sun   = get_body_barycentric("sun",   t)
    return (earth - sun).xyz.to(u.AU).value


# ---------------------------------------------------------------------------
# Kepler equation solvers  →  true anomaly
# ---------------------------------------------------------------------------

def _solve_barker(dt: float, q: float) -> float:
    """
    Barker's equation (parabolic, e = 1).
    Returns D = tan(ν/2).  dt = t − T_peri [days].
    """
    W         = 3.0 * np.sqrt(MU / (2.0 * q ** 3)) * dt
    sqrt_disc = np.sqrt((W / 2.0) ** 2 + 1.0)
    return np.cbrt(W / 2.0 + sqrt_disc) + np.cbrt(W / 2.0 - sqrt_disc)


def _solve_kepler_elliptic(M: float, e: float,
                           tol: float = 1e-13, max_iter: int = 200) -> float:
    """
    Solve M = E − e·sin E for eccentric anomaly E (radians).
    Uses Newton-Raphson with Kepler's classic initial guess.
    """
    E = M if e < 0.8 else np.sign(np.sin(M)) * np.pi
    for _ in range(max_iter):
        f  = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return E


def _solve_kepler_hyperbolic(M_h: float, e: float,
                             tol: float = 1e-13, max_iter: int = 200) -> float:
    """
    Solve e·sinh F − F = M_h for hyperbolic anomaly F.
    Uses Newton-Raphson with Barker-inspired initial guess.
    """
    # Initial guess that works well for all M_h
    F = np.sign(M_h) * np.log(2.0 * abs(M_h) / e + 1.8)
    for _ in range(max_iter):
        f  = e * np.sinh(F) - F - M_h
        fp = e * np.cosh(F) - 1.0
        dF = -f / fp
        F += dF
        if abs(dF) < tol:
            break
    return F


def true_anomaly_from_time(dt: float, q: float, e: float) -> float:
    """
    Compute true anomaly ν (radians) from time past perihelion dt [days],
    perihelion distance q [AU], and eccentricity e.

    Dispatches to the correct Kepler solver:
        |e − 1| < 1e-4  →  Barker's equation
        e < 1           →  elliptic Kepler
        e > 1           →  hyperbolic Kepler
    """
    if abs(e - 1.0) < _E_PARA_TOL:
        # Parabolic
        D  = _solve_barker(dt, q)
        return 2.0 * np.arctan(D)

    elif e < 1.0:
        # Elliptic
        a  = q / (1.0 - e)
        n  = np.sqrt(MU / a ** 3)          # mean motion [rad/day]
        M  = (n * dt) % (2.0 * np.pi)
        E  = _solve_kepler_elliptic(M, e)
        return 2.0 * np.arctan2(
            np.sqrt(1.0 + e) * np.sin(E / 2.0),
            np.sqrt(1.0 - e) * np.cos(E / 2.0),
        )

    else:
        # Hyperbolic (e > 1)
        # a is negative; |a| = q / (e - 1)
        a_abs = q / (e - 1.0)
        n_h   = np.sqrt(MU / a_abs ** 3)   # hyperbolic mean motion [rad/day]
        M_h   = n_h * dt
        F     = _solve_kepler_hyperbolic(M_h, e)
        nu    = 2.0 * np.arctan2(
            np.sqrt(e + 1.0) * np.sinh(F / 2.0),
            np.sqrt(e - 1.0) * np.cosh(F / 2.0),
        )
        # Clamp to the physical range |nu| < arccos(-1/e)
        nu_max = np.arccos(-1.0 / e) - 1e-10
        return float(np.clip(nu, -nu_max, nu_max))


# ---------------------------------------------------------------------------
# Orbit propagation:  elements  →  heliocentric equatorial position
# ---------------------------------------------------------------------------

def comet_heliocentric_equatorial(
    t_jd: float,
    q: float, e: float,
    i: float, Omega: float, omega: float,
    T: float,
) -> np.ndarray:
    """
    Heliocentric equatorial ICRS J2000 position of the comet (AU).

    Parameters
    ----------
    t_jd  : epoch (JD TDB)
    q     : perihelion distance (AU)
    e     : eccentricity
    i     : inclination (deg, ecliptic J2000)
    Omega : longitude of ascending node (deg, ecliptic J2000)
    omega : argument of perihelion (deg, ecliptic J2000)
    T     : time of perihelion passage (JD TDB)
    """
    i_r = np.radians(i)
    O_r = np.radians(Omega)
    o_r = np.radians(omega)

    nu    = true_anomaly_from_time(t_jd - T, q, e)
    p     = q * (1.0 + e) if abs(e - 1.0) >= _E_PARA_TOL else 2.0 * q
    r_mag = p / (1.0 + e * np.cos(nu))

    # Perifocal (PQW) unit vectors in ecliptic J2000 frame
    ci, si = np.cos(i_r), np.sin(i_r)
    cO, sO = np.cos(O_r), np.sin(O_r)
    co, so = np.cos(o_r), np.sin(o_r)

    P_hat = np.array([ cO*co - sO*so*ci,
                        sO*co + cO*so*ci,
                        so*si ])
    Q_hat = np.array([-cO*so - sO*co*ci,
                      -sO*so + cO*co*ci,
                       co*si ])

    pos_ecl = r_mag * (np.cos(nu) * P_hat + np.sin(nu) * Q_hat)
    return R_ECL_TO_EQU @ pos_ecl


# ---------------------------------------------------------------------------
# Residual / cost functions
# ---------------------------------------------------------------------------

Observation = Tuple[float, float, float, np.ndarray]   # ra, dec, jd, R_earth


def angular_residuals(
    params: np.ndarray,
    observations: List[Observation],
    fit_eccentricity: bool,
) -> np.ndarray:
    """
    Angular residuals (radians) in (RA·cos δ, δ) for every observation.
    Returns array of length 2 * N_obs.

    params layout
    -------------
    fit_eccentricity=False : [q, i, Omega, omega, T]          (e = 1 fixed)
    fit_eccentricity=True  : [q, e, i, Omega, omega, T]
    """
    if fit_eccentricity:
        q, e, i, Omega, omega, T = params
    else:
        q, i, Omega, omega, T = params
        e = 1.0

    if q <= 0.0 or e < 0.0:
        return np.full(2 * len(observations), 1e6)

    residuals: List[float] = []
    for ra_deg, dec_deg, jd, R_earth in observations:
        try:
            r_comet = comet_heliocentric_equatorial(jd, q, e, i, Omega, omega, T)
        except Exception:
            residuals += [1e6, 1e6]
            continue

        rho_vec = r_comet - R_earth
        rho_mag = np.linalg.norm(rho_vec)
        if rho_mag < 1e-10:
            residuals += [1e6, 1e6]
            continue

        rho_hat  = rho_vec / rho_mag
        pred_dec = np.arcsin(np.clip(rho_hat[2], -1.0, 1.0))
        pred_ra  = np.arctan2(rho_hat[1], rho_hat[0]) % (2.0 * np.pi)

        obs_dec  = np.radians(dec_deg)
        obs_ra   = np.radians(ra_deg) % (2.0 * np.pi)

        dra  = (pred_ra - obs_ra) * np.cos(obs_dec)
        ddec =  pred_dec - obs_dec
        dra  = (dra + np.pi) % (2.0 * np.pi) - np.pi   # wrap to [−π, π]

        residuals += [dra, ddec]

    return np.array(residuals)


def cost_function(
    params: np.ndarray,
    observations: List[Observation],
    fit_eccentricity: bool,
) -> float:
    return float(np.sum(angular_residuals(params, observations, fit_eccentricity) ** 2))


# ---------------------------------------------------------------------------
# Main fitting routine
# ---------------------------------------------------------------------------

def fit_orbit(
    observations_input: List[Tuple[float, float, float]],
    fit_eccentricity: bool = False,
    verbose: bool = True,
) -> dict:
    """
    Fit a cometary orbit to astrometric observations.

    Parameters
    ----------
    observations_input : list of (ra_deg, dec_deg, jd)
        Minimum 2 observations.  RA/Dec in decimal degrees; JD in TDB.
    fit_eccentricity : bool
        False  →  assume parabolic (e = 1, 5 free parameters).
        True   →  fit eccentricity as a free parameter (6 free parameters).
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict with keys:
        success      : bool
        q            : perihelion distance (AU)
        e            : eccentricity
        i            : inclination (deg, ecliptic J2000)
        Omega        : longitude of ascending node (deg)
        omega        : argument of perihelion (deg)
        T            : perihelion time (JD TDB)
        orbit_type   : 'elliptic' | 'parabolic' | 'hyperbolic'
        rms_arcsec   : RMS angular residual (arcsec)
        message      : human-readable status string
    """
    n_obs = len(observations_input)
    if n_obs < 2:
        return {"success": False, "message": "At least 2 observations required."}

    # Underdetermined-system warnings
    n_params   = 6 if fit_eccentricity else 5
    n_equations = 2 * n_obs
    if n_equations < n_params and verbose:
        print(
            f"WARNING: {n_obs} observation(s) give {n_equations} equations "
            f"for {n_params} unknowns — system is underdetermined.  "
            f"The solution may not be unique."
        )
    elif n_equations == n_params and verbose:
        print(
            f"WARNING: {n_obs} observation(s) give exactly {n_equations} "
            f"equations for {n_params} unknowns — solution is barely "
            f"determined.  Consider adding more observations."
        )

    # Pre-fetch Earth positions (expensive; do once)
    if verbose:
        print("Fetching Earth ephemeris for each epoch…")
    observations: List[Observation] = []
    for ra, dec, jd in observations_input:
        try:
            R_earth = get_earth_heliocentric(jd)
        except Exception as exc:
            return {"success": False,
                    "message": f"Failed to retrieve Earth ephemeris: {exc}"}
        observations.append((ra, dec, jd, R_earth))

    jds    = [o[2] for o in observations_input]
    T_mid  = (min(jds) + max(jds)) / 2.0
    T_half = max(max(jds) - min(jds), 30.0)

    if fit_eccentricity:
        bounds = [
            (0.01,  10.0),                                   # q [AU]
            (0.0,   20.0),                                   # e
            (0.0,  180.0),                                   # i [deg]
            (0.0,  360.0),                                   # Omega [deg]
            (0.0,  360.0),                                   # omega [deg]
            (T_mid - 3.0*T_half, T_mid + 3.0*T_half),       # T [JD]
        ]
    else:
        bounds = [
            (0.01,  10.0),
            (0.0,  180.0),
            (0.0,  360.0),
            (0.0,  360.0),
            (T_mid - 3.0*T_half, T_mid + 3.0*T_half),
        ]

    if verbose:
        mode = "free eccentricity" if fit_eccentricity else "parabolic (e = 1)"
        print(f"Running global optimisation ({mode})…")

    args = (observations, fit_eccentricity)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        de_result = differential_evolution(
            cost_function,
            bounds,
            args=args,
            maxiter=3000,
            tol=1e-12,
            seed=42,
            workers=1,
            polish=True,
            popsize=25,
            mutation=(0.5, 1.5),
            recombination=0.9,
        )

    if verbose:
        print("Polishing with Nelder-Mead…")

    nm_result = minimize(
        cost_function,
        de_result.x,
        args=args,
        method="Nelder-Mead",
        options={"maxiter": 100_000, "xatol": 1e-12, "fatol": 1e-14},
    )

    best = nm_result.x if nm_result.fun < de_result.fun else de_result.x

    # Unpack and normalise
    if fit_eccentricity:
        q, e, i, Omega, omega, T = best
        e = float(abs(e))           # eccentricity must be non-negative
    else:
        q, i, Omega, omega, T = best
        e = 1.0

    i     = float(abs(i)) % 180.0
    Omega = float(Omega)  % 360.0
    omega = float(omega)  % 360.0
    q     = float(abs(q))
    T     = float(T)

    # RMS in arcseconds
    best_norm = np.array([q, e, i, Omega, omega, T] if fit_eccentricity
                         else [q, i, Omega, omega, T])
    res_rad    = angular_residuals(best_norm, observations, fit_eccentricity)
    rms_arcsec = float(np.degrees(np.sqrt(np.mean(res_rad ** 2))) * 3600.0)

    orbit_type = _orbit_type_label(e)

    base = dict(q=q, e=e, i=i, Omega=Omega, omega=omega, T=T,
                orbit_type=orbit_type, rms_arcsec=rms_arcsec)

    if rms_arcsec > RMS_THRESHOLD_ARCSEC:
        return {
            "success": False,
            **base,
            "message": (
                f"No good fit found.  Best RMS = {rms_arcsec:.1f}\".  "
                f"Threshold = {RMS_THRESHOLD_ARCSEC}\".  "
                "Try adding more observations or a wider arc."
            ),
        }

    return {
        "success": True,
        **base,
        "message": f"Solution converged.  RMS = {rms_arcsec:.3f}\".",
    }


def fit_parabolic_orbit(
    observations_input: List[Tuple[float, float, float]],
    verbose: bool = True,
) -> dict:
    """Backward-compatible wrapper: fit_orbit with e = 1 fixed."""
    return fit_orbit(observations_input, fit_eccentricity=False, verbose=verbose)


# ---------------------------------------------------------------------------
# Orbit-type labelling
# ---------------------------------------------------------------------------

def _orbit_type_label(e: float) -> str:
    if abs(e - 1.0) < _E_PARA_TOL:
        return "parabolic"
    elif e < 1.0:
        return "elliptic"
    else:
        return "hyperbolic"


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_result(result: dict) -> None:
    orbit_type = result.get("orbit_type", "parabolic")
    title = f"Orbit Determination  [{orbit_type}]"

    print("\n" + "=" * 54)
    print(f"  {title}")
    print("=" * 54)

    if not result["success"]:
        print(f"STATUS : FAILED")
        print(f"REASON : {result['message']}")
        if "q" in result:
            print("\nBest (non-accepted) solution:")
            _print_elements(result)
    else:
        print(f"STATUS : SUCCESS")
        print(f"        {result['message']}")
        print(f"\nOrbital Elements  (ecliptic J2000):")
        _print_elements(result)

    print("=" * 54)


def _print_elements(result: dict) -> None:
    e          = result["e"]
    orbit_type = result.get("orbit_type", _orbit_type_label(e))
    t_peri     = Time(result["T"], format="jd", scale="tdb")

    print(f"  Orbit type           = {orbit_type}")
    print(f"  Eccentricity     e   = {e:.6f}")
    print(f"  Perihelion dist  q   = {result['q']:.6f}  AU")

    if orbit_type == "elliptic":
        a = result["q"] / (1.0 - e)
        print(f"  Semi-major axis  a   = {a:.6f}  AU")

    print(f"  Inclination      i   = {result['i']:.4f}  deg")
    print(f"  Ascending node   Ω   = {result['Omega']:.4f}  deg")
    print(f"  Arg. perihelion  ω   = {result['omega']:.4f}  deg")
    print(f"  Perihelion time  T   = {result['T']:.4f}  JD")
    print(f"                       ({t_peri.iso} TDB)")
    print(f"  Fit RMS              = {result['rms_arcsec']:.3f}  arcsec")


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _parse_observations(args) -> List[Tuple[float, float, float]]:
    observations = []
    for triplet in (args.obs or []):
        if len(triplet) != 3:
            raise ValueError(
                f"Each --obs must have exactly 3 values (RA Dec JD), got: {triplet}"
            )
        observations.append((float(triplet[0]), float(triplet[1]), float(triplet[2])))
    return observations


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Fit a cometary orbit to Earth-based RA/Dec/JD observations.\n"
            "Assumes a parabolic orbit by default (e = 1); use\n"
            "--fit-eccentricity to fit e as a free parameter."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
# Parabolic fit (default):
python comet_orbit.py \\
    --obs 157.5 -15.0 2446470.5 \\
    --obs 168.0  -8.3 2446500.5 \\
    --obs 180.0   2.1 2446530.5

# Free-eccentricity fit:
python comet_orbit.py --fit-eccentricity \\
    --obs 157.5 -15.0 2446470.5 \\
    --obs 168.0  -8.3 2446500.5 \\
    --obs 180.0   2.1 2446530.5 \\
    --obs 190.0   8.5 2446560.5

# RA is in decimal degrees (not hours).  JD is Julian Date (TDB).
""",
    )
    parser.add_argument(
        "--obs",
        metavar=("RA", "DEC", "JD"),
        nargs=3,
        action="append",
        help="One observation: RA (deg), Dec (deg), JD (TDB).  Repeat for each.",
    )
    parser.add_argument(
        "--fit-eccentricity", "--fit-ecc",
        dest="fit_eccentricity",
        action="store_true",
        default=False,
        help=(
            "Fit eccentricity as a free parameter instead of fixing e = 1.  "
            "Allows elliptic, parabolic, or hyperbolic solutions.  "
            "Requires at least 4 observations for a well-determined system."
        ),
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress messages.",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help=(
            "Built-in demo.  Without --fit-eccentricity: synthetic parabolic "
            "comet (e=1, q=1.0 AU, i=45°).  With --fit-eccentricity: "
            "synthetic elliptic comet (e=0.7, q=0.6 AU, i=30°)."
        ),
    )

    args = parser.parse_args()

    if args.demo or not args.obs:
        print("=== Demo mode ===\n")

        if args.fit_eccentricity:
            # Elliptic demo comet
            q_true, e_true = 0.6, 0.7
            i_true, O_true, o_true = 30.0, 80.0, 150.0
        else:
            # Parabolic demo comet
            q_true, e_true = 1.0, 1.0
            i_true, O_true, o_true = 45.0, 120.0, 200.0

        T_true = 2460310.0
        test_jds = [T_true - 60, T_true - 30, T_true + 15, T_true + 45]

        print("True elements:")
        print(f"  e={e_true}  q={q_true} AU  i={i_true}°  "
              f"Ω={O_true}°  ω={o_true}°  T={T_true} JD")
        print("\nSynthetic observations (noise-free):")

        observations_input = []
        for jd in test_jds:
            R_earth = get_earth_heliocentric(jd)
            r_comet = comet_heliocentric_equatorial(
                jd, q_true, e_true, i_true, O_true, o_true, T_true
            )
            rho     = r_comet - R_earth
            rho_hat = rho / np.linalg.norm(rho)
            dec = float(np.degrees(np.arcsin(np.clip(rho_hat[2], -1.0, 1.0))))
            ra  = float(np.degrees(np.arctan2(rho_hat[1], rho_hat[0])) % 360.0)
            observations_input.append((ra, dec, jd))
            t = Time(jd, format="jd")
            print(f"  JD {jd:.1f}  ({t.iso})  RA={ra:.4f}°  Dec={dec:+.4f}°")

    else:
        observations_input = _parse_observations(args)

    mode_str = "free-eccentricity" if args.fit_eccentricity else "parabolic"
    print(f"\nFitting {mode_str} orbit to "
          f"{len(observations_input)} observation(s)…\n")

    result = fit_orbit(
        observations_input,
        fit_eccentricity=args.fit_eccentricity,
        verbose=not args.quiet,
    )
    print_result(result)


if __name__ == "__main__":
    main()
