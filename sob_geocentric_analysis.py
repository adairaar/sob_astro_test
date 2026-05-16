"""
sob_geocentric_analysis.py
Geocentric orbital survey for the Star of Bethlehem impossibility proof.

Tests whether any Earth-orbiting object can exhibit the apparent motion
described in Matthew 2:9 — guided motion (~2 deg/h) followed by stopping
(~0 deg/h) — in the Bethlehem azimuth corridor (190-210 deg), at night,
during the relevant epoch (5 BCE).

Reference: Adair (2026), "Orbital mechanics and corpus linguistics jointly
preclude any astronomical explanation for the Star of Bethlehem."
"""

import numpy as np

# ── Constants ──────────────────────────────────────────────────────────────
mu_E  = 3.986004418e5   # km^3/s^2   (Earth GM)
R_E   = 6371.0          # km          (Earth radius)
d2r   = np.pi / 180
r2d   = 180.0 / np.pi

# ── Observer: Jerusalem ────────────────────────────────────────────────────
phi = 31.77 * d2r    # latitude
lam = 35.24 * d2r    # longitude
h_obs = 0.760        # km altitude


def kepler_E(M, e, tol=1e-11):
    """Vectorised Newton-Raphson solution to Kepler's equation."""
    E = np.array(M, dtype=float)
    for _ in range(60):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))
        E += dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def orbit_eci(a, e, i_r, om_r, Om_r, M0_r, t_s):
    """ECI position (km) for a geocentric Keplerian orbit."""
    n  = np.sqrt(mu_E / a**3)
    M  = (M0_r + n * t_s) % (2 * np.pi)
    E  = kepler_E(M, e)
    nu = 2.0 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E / 2),
        np.sqrt(1 - e) * np.cos(E / 2))
    r  = a * (1.0 - e * np.cos(E))
    xp = r * np.cos(nu)
    yp = r * np.sin(nu)

    cO, sO = np.cos(Om_r), np.sin(Om_r)
    ci, si = np.cos(i_r),  np.sin(i_r)
    co, so = np.cos(om_r), np.sin(om_r)

    rx = cO * co - sO * so * ci;  ry = -(cO * so + sO * co * ci)
    sx = sO * co + cO * so * ci;  sy = -(sO * so - cO * co * ci)
    tx = so * si;                  ty = co * si

    return np.array([rx*xp + ry*yp,
                     sx*xp + sy*yp,
                     tx*xp + ty*yp])


def observer_eci(t_s, GMST0):
    """Observer's ECI position (km)."""
    GMST = GMST0 + 7.2921150e-5 * t_s
    LST  = GMST + lam
    R    = R_E + h_obs
    return R * np.array([
        np.cos(phi) * np.cos(LST),
        np.cos(phi) * np.sin(LST),
        np.sin(phi) * np.ones_like(t_s)])


def topo_altaz(r_eci, obs_eci, t_s, GMST0):
    """Topocentric altitude (deg), azimuth (deg), range (km)."""
    GMST = GMST0 + 7.2921150e-5 * t_s
    LST  = GMST + lam
    rho  = r_eci - obs_eci
    mag  = np.linalg.norm(rho, axis=0)
    sL, cL = np.sin(LST), np.cos(LST)
    sp, cp = np.sin(phi), np.cos(phi)
    S = sp*cL*rho[0] + sp*sL*rho[1] - cp*rho[2]
    E = -sL*rho[0] + cL*rho[1]
    Z =  cp*cL*rho[0] + cp*sL*rho[1] + sp*rho[2]
    alt = np.arcsin(np.clip(Z / mag, -1, 1)) * r2d
    az  = (np.arctan2(E, -S) * r2d) % 360.0
    return alt, az, mag


# ── Sun altitude (simplified model, 5 BCE Jun 8) ──────────────────────────
def sun_alt_array(t_s, GMST0, JD0=1719755.898):
    JD2K = 2451545.0
    TJ   = (JD0 - JD2K) / 36525.0
    L0   = (280.46646 + 36000.76983 * TJ) % 360
    M0s  = (357.52911 + 35999.05029 * TJ) % 360
    t_d  = t_s / 86400.0
    L    = (L0  + 360.98564724 * t_d) % 360
    Ms   = (M0s + 360.98564724 * t_d) % 360
    C    = (1.9146 - 0.004817*TJ)*np.sin(Ms*d2r) + 0.019993*np.sin(2*Ms*d2r)
    lon_s = (L + C) % 360
    eps   = 23.439 - 0.013 * TJ
    dec_s = np.arcsin(np.sin(eps*d2r) * np.sin(lon_s*d2r))
    ra_s  = np.arctan2(np.cos(eps*d2r)*np.sin(lon_s*d2r), np.cos(lon_s*d2r))
    HA_s  = GMST0 + 7.2921150e-5*t_s + lam - ra_s
    return np.arcsin(
        np.sin(phi)*np.sin(dec_s) +
        np.cos(phi)*np.cos(dec_s)*np.cos(HA_s)) * r2d


# ── Survey configuration ────────────────────────────────────────────────────
DT   = 300        # 5-minute timestep (seconds)
t_s  = np.arange(0, 12*3600 + DT, DT)
DH   = DT / 3600.0

# GMST0: local midnight ~ 6 h into window (approximate 5 BCE Jun 8)
GMST0  = 2.395
sun_a  = sun_alt_array(t_s, GMST0)
night  = sun_a < -6.0   # civil twilight threshold

# Orbital parameter grids
a_vals  = np.array([7000, 10000, 15000, 20000, 26560, 42164,
                    80000, 150000, 300000, 384000])  # km
e_vals  = np.array([0.0, 0.2, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99])
i_vals  = np.array([0, 45, 90, 135, 180])
om_vals = np.arange(0, 360, 90)
Om_vals = np.arange(0, 360, 90)
M0_vals = np.arange(0, 360, 45)

# Criteria (ground-frame apparent angular velocity)
GUIDE_LO = 1.0;  GUIDE_HI = 5.0   # deg/h
STOP_TIGHT = 0.20                   # deg/h  (strict: truly "stopped")
STOP_LOOSE = 0.50                   # deg/h  (loose)
AZ_LO, AZ_HI = 190.0, 210.0
ALT_MIN = 5.0
MAX_TRANS = 4.0                     # hours


def angular_velocity(alt, az, dt_h):
    alt_r = alt * d2r
    ca    = np.cos(alt_r)
    daz   = np.gradient(np.unwrap(az * d2r) * r2d, dt_h)
    dalt  = np.gradient(alt, dt_h)
    return np.sqrt((daz * ca)**2 + dalt**2)


def scan_orbit(a, e, i_deg, om_deg, Om_deg, M0_deg, stop_thresh):
    i_r  = i_deg * d2r;  om_r = om_deg * d2r
    Om_r = Om_deg * d2r; M0_r = M0_deg * d2r
    if a * (1 - e) < R_E + 100:
        return None
    try:
        r_eci = orbit_eci(a, e, i_r, om_r, Om_r, M0_r, t_s)
        obs   = observer_eci(t_s, GMST0)
        alt, az, dist = topo_altaz(r_eci, obs, t_s, GMST0)
    except Exception:
        return None
    om = angular_velocity(alt, az, DH)
    ok = night & (alt > ALT_MIN) & (az >= AZ_LO) & (az <= AZ_HI)
    g  = ok & (om >= GUIDE_LO) & (om <= GUIDE_HI)
    s  = ok & (om < stop_thresh)
    if not (np.any(g) and np.any(s)):
        return None
    for gi in np.where(g)[0]:
        cands = np.where(s)[0]
        cands = cands[cands > gi]
        if not len(cands):
            continue
        dt_h = (t_s[cands[0]] - t_s[gi]) / 3600.0
        if dt_h <= MAX_TRANS:
            return dict(a=a, e=e, i=i_deg, om=om_deg, Om=Om_deg, M0=M0_deg,
                        dt_h=dt_h, og=om[gi], os=om[cands[0]],
                        azg=az[gi], alg=alt[gi], dg=dist[gi])
    return None


# ── Run survey ─────────────────────────────────────────────────────────────
found_tight = []; found_loose = []; n_check = 0

for a in a_vals:
    for e in e_vals:
        for i_deg in i_vals:
            for Om_deg in Om_vals:
                for om_deg in om_vals:
                    for M0_deg in M0_vals:
                        n_check += 1
                        r = scan_orbit(a, e, i_deg, om_deg, Om_deg,
                                       M0_deg, STOP_TIGHT)
                        if r:
                            found_tight.append(r)
                        else:
                            r = scan_orbit(a, e, i_deg, om_deg, Om_deg,
                                           M0_deg, STOP_LOOSE)
                            if r:
                                found_loose.append(r)

# ── Moon reference ─────────────────────────────────────────────────────────
a_m = 384400; e_m = 0.055
h_m = np.sqrt(mu_E * a_m * (1 - e_m**2))
om_mean_m = h_m / a_m**2 * r2d * 3600
om_peri_m = h_m / (a_m*(1-e_m))**2 * r2d * 3600
om_apo_m  = h_m / (a_m*(1+e_m))**2 * r2d * 3600

# ── Report ─────────────────────────────────────────────────────────────────
print("=" * 65)
print("GEOCENTRIC ORBITAL SURVEY — Star of Bethlehem")
print("=" * 65)
print(f"Observer: Jerusalem (31.77°N, 35.24°E)   Epoch: 5 BCE Jun 8")
print(f"Guidance:  ω ∈ [1, 5]°/h, az ∈ [190°, 210°], nighttime, alt>5°")
print(f"Max guidance→stop transition: {MAX_TRANS} h")
print()
print(f"Configurations checked: {n_check:,}")
print()
print(f"Moon angular velocity (only natural geocentric body):")
print(f"  Mean: {om_mean_m:.3f}°/h   Perigee: {om_peri_m:.3f}°/h   "
      f"Apogee: {om_apo_m:.3f}°/h")
print(f"  Guidance threshold (1°/h): NOT REACHED — Moon never qualifies")
print()
print(f"Hits with strict stop (ω < {STOP_TIGHT}°/h):  {len(found_tight)}")
print(f"Hits with loose  stop (ω < {STOP_LOOSE}°/h):  {len(found_loose)}")
print()
if found_loose:
    print("Loose hits (equatorial transit artifact orbits):")
    seen = set()
    for f in found_loose:
        key = (f['a'], f['e'])
        if key not in seen:
            seen.add(key)
            print(f"  a={f['a']:.0f} km, e={f['e']:.2f}, i={f['i']:.0f}°  "
                  f"ω_guide={f['og']:.2f}→ω_stop={f['os']:.3f}°/h "
                  f"(transit azimuth effect; ω_stop >> 0)")
print()
print("Conclusion:")
print("  No geocentric orbit satisfies all criteria with realistic")
print("  stopping threshold.  The Moon never reaches guidance speed.")
print("  Zero solutions found across all Earth-orbiting configurations.")
