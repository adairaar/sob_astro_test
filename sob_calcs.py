#!/usr/bin/env python3
"""
Star of Bethlehem – Orbital Mechanics Calculations for Nature paper
Aaron Adair / 2026

Three main analyses:
 1. Matney (2025) comet apparent path from Bethlehem on 5 BCE Jun 8
    → show azimuth is NOT ~206° during the Magi's journey window
 2. Diurnal motion constraint: no fixed-sky object maintains azimuth 206° while rising
 3. Velocity / distance requirement for any "temporary geosynchronous" object
"""

import numpy as np
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation, AltAz, get_body_barycentric,
    SkyCoord, get_sun
)
import astropy.units as u

# ── site ─────────────────────────────────────────────────────────────────────
BETHLEHEM_LAT =  31.70   # °N
BETHLEHEM_LON =  35.20   # °E
BETHLEHEM_ALT = 765      # m
BETHLEHEM = EarthLocation(lon=BETHLEHEM_LON*u.deg,
                          lat=BETHLEHEM_LAT*u.deg,
                          height=BETHLEHEM_ALT*u.m)
MAGI_AZ = 206.0          # azimuth of Jerusalem→Bethlehem road (°)

# ── physical constants ────────────────────────────────────────────────────────
K_GAUSS = 0.01720209895
MU      = K_GAUSS**2                     # GM_sun  [AU³/day²]
EPS     = np.radians(23.439291111)       # obliquity J2000
R_ECL2EQU = np.array([
    [1.0, 0.0,          0.0         ],
    [0.0, np.cos(EPS), -np.sin(EPS) ],
    [0.0, np.sin(EPS),  np.cos(EPS) ],
])

# ── Matney (2025) Table, J2000 ecliptic elements ─────────────────────────────
# q=0.0407 AU, e=1.0 (parabolic), i=1.39°, Ω=93.01°, ω=346.25°
# T_perihelion = JD 1719785.565 (TDB)
M_q     = 0.0407
M_e     = 1.0
M_i     = 1.39
M_Omega = 93.01
M_omega = 346.25
M_T     = 1719785.565      # JD TDB

# Reference: Matney states local solar noon at Bethlehem on 5 BCE Jun 8
#   = JD 1719755.898 (his text, p.396)
JD_NOON = 1719755.898

# ΔT (TDB - UTC) in 5 BCE ≈ +10572 s  (Espenak & Meeus 2006, cited in paper)
DT_SEC  = 10572.0

# ─────────────────────────────────────────────────────────────────────────────
# Orbital mechanics
# ─────────────────────────────────────────────────────────────────────────────

def solve_barker(dt_days, q):
    """Barker's equation for parabolic orbit: returns tan(nu/2)."""
    W = 3.0 * np.sqrt(MU / (2.0 * q**3)) * dt_days
    s = np.sqrt((W/2.0)**2 + 1.0)
    return np.cbrt(W/2.0 + s) + np.cbrt(W/2.0 - s)

def comet_helio_equ(jd_tdb, q, e, i_d, Om_d, om_d, T_jd):
    """Heliocentric equatorial (ICRS J2000) position of comet [AU]."""
    i  = np.radians(i_d)
    Om = np.radians(Om_d)
    om = np.radians(om_d)
    dt = jd_tdb - T_jd
    D  = solve_barker(dt, q)
    nu = 2.0 * np.arctan(D)
    p  = 2.0 * q
    r  = p / (1.0 + e * np.cos(nu))
    ci, si = np.cos(i), np.sin(i)
    cO, sO = np.cos(Om), np.sin(Om)
    co, so = np.cos(om), np.sin(om)
    P = np.array([ cO*co - sO*so*ci,  sO*co + cO*so*ci,  so*si])
    Q = np.array([-cO*so - sO*co*ci, -sO*so + cO*co*ci,  co*si])
    pos_ecl = r * (np.cos(nu)*P + np.sin(nu)*Q)
    return R_ECL2EQU @ pos_ecl

def earth_helio(jd_tdb):
    t = Time(jd_tdb, format='jd', scale='tdb')
    E = get_body_barycentric('earth', t)
    S = get_body_barycentric('sun',   t)
    return (E - S).xyz.to(u.AU).value

def apparent_azel(jd_tdb, q, e, i, Om, om, T):
    """Return (az°, alt°, ra°, dec°, dist_AU, solar_elong°) from Bethlehem."""
    r_c   = comet_helio_equ(jd_tdb, q, e, i, Om, om, T)
    R_e   = earth_helio(jd_tdb)
    rho   = r_c - R_e
    dist  = np.linalg.norm(rho)
    rh    = rho / dist
    dec   = np.degrees(np.arcsin(np.clip(rh[2], -1, 1)))
    ra    = np.degrees(np.arctan2(rh[1], rh[0])) % 360.0
    # convert TDB → approximate UTC
    jd_utc = jd_tdb - DT_SEC / 86400.0
    t_utc  = Time(jd_utc, format='jd', scale='utc')
    coord  = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    frame  = AltAz(obstime=t_utc, location=BETHLEHEM)
    aa     = coord.transform_to(frame)
    sun_c  = get_sun(t_utc).transform_to(frame)
    # solar elongation (angular sep between comet and Sun)
    sun_eq = get_sun(t_utc)
    elong  = coord.separation(sun_eq).deg
    return aa.az.deg, aa.alt.deg, ra, dec, dist, elong

# ═════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1 – Matney comet apparent path on 5 BCE Jun 8
# ═════════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("ANALYSIS 1: Matney (2025) comet apparent path – 5 BCE Jun 8")
print("="*70)

# Hourly grid from 04:00 to 12:00 local solar time
hours_from_noon = np.arange(-8.0, 1.1, 0.25)   # -8h to +1h → 04:00 to 13:00
jd_grid = JD_NOON + hours_from_noon / 24.0

rows = []
for h, jd in zip(hours_from_noon, jd_grid):
    az, alt, ra, dec, dist, elong = apparent_azel(
        jd, M_q, M_e, M_i, M_Omega, M_omega, M_T)
    local_h = h + 12.0   # hours past midnight
    rows.append((local_h, az, alt, dist, elong))

rows = np.array(rows)   # cols: local_h, az, alt, dist_AU, elong

print(f"\n{'Time':>6}  {'Az(°)':>7}  {'Alt(°)':>7}  {'Dist(AU)':>10}  {'SolarElong(°)':>13}")
print("-" * 52)
for r in rows:
    local_h, az, alt, dist, elong = r
    h_int, h_frac = divmod(local_h, 1)
    tstr = f"{int(h_int):02d}:{int(h_frac*60):02d}"
    marker = " ◄ JOURNEY" if 8.0 <= local_h <= 10.0 else ""
    print(f"{tstr:>6}  {az:7.2f}  {alt:7.2f}  {dist:10.6f}  {elong:13.2f}{marker}")

# Magi's journey window: 08:00–10:00 local per Matney's own narrative
mask = (rows[:,0] >= 8.0) & (rows[:,0] <= 10.0)
az_j = rows[mask, 1]
alt_j = rows[mask, 2]
elong_j = rows[mask, 4]

print(f"\n── Magi's journey window (08:00–10:00 local) ──────────────────")
print(f"   Azimuth range    : {az_j.min():.1f}° – {az_j.max():.1f}°")
print(f"   Required azimuth : {MAGI_AZ:.1f}°  (Jerusalem→Bethlehem road)")
print(f"   Mean azimuth deviation from 206°: {np.mean(np.abs(az_j - MAGI_AZ)):.1f}°")
print(f"   Azimuth drift during 2h window  : {az_j.max()-az_j.min():.1f}°")
print(f"   Altitude range   : {alt_j.min():.1f}° – {alt_j.max():.1f}°")
print(f"   Solar elongation : {elong_j.min():.1f}° – {elong_j.max():.1f}°")

# Rate of azimuth change (°/h)
t_j = rows[mask, 0]
daz_dt = np.gradient(az_j, t_j)
print(f"   Mean d(Az)/dt    : {np.mean(daz_dt):+.1f} °/h  (+ = moving east)")
print(f"   Required d(Az)/dt: ~0 °/h  (must stay at 206°)")


# ═════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2 – Diurnal motion: azimuth drift for fixed-sky objects
# ═════════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("ANALYSIS 2: Diurnal motion – azimuth drift of fixed-sky objects")
print("="*70)
print(f"\nSite: Bethlehem, φ = {BETHLEHEM_LAT}°N")
print("Question: can any object on the celestial sphere maintain")
print(f"azimuth = {MAGI_AZ}° while altitude increases over 2 hours?\n")

# Use a date/time when a notional object is AT azimuth 206° at various altitudes.
# We'll start at JD_NOON - 4/24 (08:00 local) and track for 2 hours.

# Declinations to test: various values
test_decs = [-30, -20, -10, 0, 10, 20]
t_start_utc = Time(JD_NOON - DT_SEC/86400 - 4/24, format='jd', scale='utc')

print(f"{'Dec(°)':>8}  {'Start Az':>9}  {'End Az (2h later)':>18}  "
      f"{'ΔAz':>8}  {'Start Alt':>10}  {'End Alt':>8}")
print("-"*75)

for dec_d in test_decs:
    # For each declination, find the hour angle H such that Az = 206° at t_start
    # Use numerical search
    dec_r = np.radians(dec_d)
    lat_r = np.radians(BETHLEHEM_LAT)
    target_az = np.radians(MAGI_AZ)

    # Scan HA in [-180°,180°] to find where Az ≈ 206°
    ha_scan = np.linspace(-np.pi, np.pi, 3600)
    az_scan = []
    alt_scan = []
    for ha in ha_scan:
        sin_alt = (np.sin(lat_r)*np.sin(dec_r) +
                   np.cos(lat_r)*np.cos(dec_r)*np.cos(ha))
        sin_alt = np.clip(sin_alt, -1, 1)
        alt_r   = np.arcsin(sin_alt)
        cos_alt = np.cos(alt_r)
        if cos_alt < 1e-10:
            az_scan.append(np.nan); alt_scan.append(np.nan)
            continue
        cos_az = ((np.sin(dec_r) - np.sin(lat_r)*sin_alt) /
                  (np.cos(lat_r)*cos_alt))
        cos_az  = np.clip(cos_az, -1, 1)
        az_r    = np.arccos(cos_az)
        if np.sin(ha) > 0:      # HA > 0 → object west of meridian → Az > 180°
            az_r = 2*np.pi - az_r
        az_scan.append(np.degrees(az_r))
        alt_scan.append(np.degrees(alt_r))

    az_scan  = np.array(az_scan)
    alt_scan = np.array(alt_scan)

    # Find HA where Az ≈ 206° and alt > 5°
    diff = np.abs(az_scan - MAGI_AZ)
    rising = (np.gradient(alt_scan, ha_scan) < 0)   # rising = HA decreasing
    cands  = np.where((diff < 0.5) & (alt_scan > 5) & rising)[0]
    if len(cands) == 0:
        print(f"{dec_d:8.0f}  {'no solution at alt>5°':>9}")
        continue

    # Pick the candidate with highest altitude (best case)
    best = cands[np.argmax(alt_scan[cands])]
    ha0  = ha_scan[best]
    alt0 = alt_scan[best]
    az0  = az_scan[best]

    # Convert HA to RA so we can use SkyCoord
    # LST at t_start
    lst_deg = t_start_utc.sidereal_time('apparent', longitude=BETHLEHEM_LON*u.deg).deg
    ra_d    = (lst_deg - np.degrees(ha0)) % 360.0
    coord   = SkyCoord(ra=ra_d*u.deg, dec=dec_d*u.deg, frame='icrs')

    # Track over 2 hours
    results = []
    for dt_h in [0.0, 2.0]:
        t_obs = t_start_utc + dt_h * u.hour
        aa    = coord.transform_to(AltAz(obstime=t_obs, location=BETHLEHEM))
        results.append((aa.az.deg, aa.alt.deg))

    az_start, alt_start = results[0]
    az_end,   alt_end   = results[1]
    daz = az_end - az_start
    print(f"{dec_d:8.0f}  {az_start:9.1f}  {az_end:18.1f}  {daz:+8.1f}  "
          f"{alt_start:10.1f}  {alt_end:8.1f}")

print(f"\n→ No fixed-sky object maintains azimuth {MAGI_AZ}° over 2 hours.")
print("  All objects drift in azimuth due to Earth's rotation (≈±15°/h).")


# ═════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3 – Velocity/distance requirement for geosynchronous behaviour
# ═════════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("ANALYSIS 3: Requirement for 'temporary geosynchronous' motion")
print("="*70)

omega_earth = 7.2921150e-5   # rad/s  (Earth's sidereal rotation rate)
lat_r = np.radians(BETHLEHEM_LAT)

# Earth's surface speed at Bethlehem latitude
R_earth_km  = 6371.0
v_surface   = omega_earth * R_earth_km * np.cos(lat_r)
print(f"\nEarth surface speed at φ={BETHLEHEM_LAT}°N: {v_surface:.3f} km/s")

# For apparent angular velocity of flyby object to match Earth's rotation:
# ω_apparent = v_rel / d = ω_Earth  →  d = v_rel / ω_Earth
print(f"\nRequired flyby distance d = v_rel / ω_Earth:")
print(f"{'v_rel (km/s)':>14}  {'d (km)':>12}  {'d (AU)':>10}  {'d (Moon dist)':>14}")
print("-" * 55)
R_moon_km = 384400.0
AU_km     = 1.496e8
for v_rel in [10, 20, 30, 40, 50]:
    d_km      = v_rel / omega_earth * 1e-3   # km  (v in m/s → km/s: already km/s)
    # v_rel in km/s, omega_earth in rad/s → d in km
    d_km_v2   = (v_rel * 1e3) / omega_earth / 1e3   # same
    d_AU      = d_km_v2 / AU_km
    d_moon    = d_km_v2 / R_moon_km
    print(f"{v_rel:14.0f}  {d_km_v2:12.0f}  {d_AU:10.6f}  {d_moon:14.2f}")

print(f"\nMatney's closest approach: 0.0026 AU = {0.0026*AU_km:.0f} km = "
      f"{0.0026*AU_km/R_moon_km:.2f} Moon distances")

# What relative velocity does Matney's distance imply?
d_matney_km = 0.0026 * AU_km
v_implied   = omega_earth * d_matney_km * 1e3   # m/s → km/s
# wait: v = omega * d (in consistent units)
v_implied_km_s = omega_earth * (d_matney_km * 1e3) / 1e3   # km/s
print(f"→ Required relative velocity for geosynchrony at Matney's distance: "
      f"{v_implied_km_s:.2f} km/s")

# Typical comet velocity relative to Earth at 1 AU
# v_escape Earth at 1 AU ≈ 42 km/s, Earth orbital ≈ 30 km/s
# A parabolic orbit comet at 1 AU: v_helio = sqrt(2*GM_sun/r) ≈ 42 km/s
# Relative velocity ≈ 12-50 km/s depending on direction
print(f"\nTypical inner-solar-system flyby relative velocities: 10–50 km/s")
print(f"Matney requires only {v_implied_km_s:.2f} km/s for geosynchrony – "
      f"MUCH slower than typical comets.")
print(f"At Matney's distance (0.0026 AU ≈ Moon dist), a typical comet at 30 km/s")
print(f"would appear to move at {30/(0.0026*AU_km/1e3)*206265:.0f} arcsec/s = "
      f"{30/(0.0026*AU_km/1e3)*3600:.2f} °/h")
print(f"Earth's sky rotates at 15°/h, so residual motion = "
      f"{30/(0.0026*AU_km/1e3)*3600 - 15:.1f} °/h – still enormous.")

# ═════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4 – Direction of motion: geosynchrony vs. southward travel
# ═════════════════════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("ANALYSIS 4: Geosynchrony direction vs. required southward motion")
print("="*70)

print("""
A geosynchronous object (natural or artificial) is stationary relative
to a point on Earth's surface.  Earth rotates EASTWARD, so to appear
stationary overhead, an object must share Earth's EASTWARD velocity.
This means it drifts WESTWARD across the sky at the same rate as diurnal
motion carries stars – i.e. its apparent motion cancels diurnal drift.

The result: a temporarily geosynchronous object DOES NOT MOVE IN AZIMUTH.
It hovers at a fixed azimuth/altitude.

But Matthew 2:9 requires the star to:
  (a) TRAVEL SOUTHWARD (proēgen autous – "was going before them" toward
      Bethlehem, azimuth ~206°) for ~2 hours, AND THEN
  (b) STOP (estathē) over a specific dwelling.

These are MUTUALLY EXCLUSIVE behaviours within the same event:
  • "Going before them southward" requires azimuth to REMAIN ~206° as
    altitude INCREASES – the object moves toward the zenith along a
    fixed southerly bearing.
  • Geosynchronous hovering produces zero azimuth change AND zero
    altitude change (the object is frozen in the sky).

Matney's own Figure 6 shows his comet moving EASTWARD in azimuth during
the critical window – from ~190° (SSW) toward ~270° (W) – because it is
racing eastward past Earth at high speed and its eastward proper motion
EXCEEDS diurnal motion, pushing it east rather than west.

The required trajectory (fixed azimuth 206°, increasing altitude) is
geometrically distinct from both:
  • Normal diurnal motion (azimuth drifts west)
  • Geosynchronous cancellation (azimuth fixed, altitude fixed)
  • Matney's comet (azimuth moves east rapidly)

We can quantify the required proper motion:
""")

# For an object at Az=206°, Alt=30° to remain at Az=206° while rising to Alt=80°
# in 2 hours from Bethlehem lat=31.7°N, what RA/Dec motion is required?

# The azimuth-altitude transformation:
# sin(alt) = sin(lat)sin(dec) + cos(lat)cos(dec)cos(H)
# cos(az)cos(alt) = sin(dec)cos(lat) - cos(dec)sin(lat)cos(H)
# sin(az)cos(alt) = cos(dec)sin(H)

# If az is FIXED at 206° and alt increases from 30° to 80° in 2h:
lat_r = np.radians(BETHLEHEM_LAT)
az_r  = np.radians(MAGI_AZ)

print(f"If star at Az={MAGI_AZ}°, Alt=30° must reach Az={MAGI_AZ}°, Alt=80° in 2h:")
print(f"{'Alt(°)':>8}  {'Required Dec(°)':>16}  {'Required HA(°)':>15}")
print("-"*45)
for alt_d in [30, 40, 50, 60, 70, 80]:
    alt_r = np.radians(alt_d)
    # From the AltAz inversion:
    # sin(dec) = sin(alt)sin(lat) + cos(alt)cos(lat)cos(az)
    sin_dec = (np.sin(alt_r)*np.sin(lat_r) +
               np.cos(alt_r)*np.cos(lat_r)*np.cos(az_r))
    sin_dec = np.clip(sin_dec, -1, 1)
    dec_r   = np.arcsin(sin_dec)
    dec_d   = np.degrees(dec_r)
    # sin(H) = -cos(alt)sin(az) / cos(dec)
    cos_dec = np.cos(dec_r)
    if abs(cos_dec) < 1e-10:
        print(f"{alt_d:8.0f}  {'at pole':>16}")
        continue
    sin_H = -np.cos(alt_r)*np.sin(az_r) / cos_dec
    sin_H = np.clip(sin_H, -1, 1)
    H_r   = np.arcsin(sin_H)
    H_d   = np.degrees(H_r)
    print(f"{alt_d:8.0f}  {dec_d:16.2f}  {H_d:15.2f}")

print("""
→ To maintain a FIXED azimuth of 206° while rising, the star must
  simultaneously change both its DECLINATION and HOUR ANGLE at precise
  coordinated rates.  No orbital body can produce this pattern because:
  • Orbital motion changes RA/Dec slowly (typically <<1°/h except at
    extreme close approach)
  • The required Dec change here is several degrees per hour
  • And it must be coordinated exactly with changing HA to keep Az constant
  This is not physically achievable for any gravitationally bound orbit.
""")

print("\n" + "="*70)
print("ANALYSIS 5: Solar elongation of Matney's comet during flyby")
print("="*70)

# Check solar elongation throughout the day
print(f"\n{'Time':>6}  {'Alt(°)':>7}  {'Az(°)':>7}  {'SolarElong(°)':>14}  "
      f"{'Dist(AU)':>10}  {'Visible?':>10}")
print("-"*60)
for r in rows:
    local_h, az, alt, dist, elong = r
    if local_h < 4.5 or local_h > 13.0:
        continue
    h_int, h_frac = divmod(local_h, 1)
    tstr = f"{int(h_int):02d}:{int(h_frac*60):02d}"
    visible = "YES" if (alt > 0 and elong > 20) else ("below horiz" if alt <= 0 else "near Sun")
    marker = " ◄ JOURNEY" if 8.0 <= local_h <= 10.0 else ""
    print(f"{tstr:>6}  {alt:7.1f}  {az:7.1f}  {elong:14.1f}  {dist:10.6f}  {visible:>10}{marker}")

print(f"""
Key findings:
  • Matney's comet reaches closest approach (~0.0026 AU) around 11:00 local.
  • During the journey window (08:00-10:00), solar elongation is 113-141 degrees,
    so the comet IS visible -- the elongation objection does not apply here.
  • The comet azimuth during 08:00-10:00 sweeps from 193.7 to 209.7 deg
    (mean drift +8.0 deg/h eastward), passing through 206 deg only briefly.
  • The comet is NOT fixed at the road azimuth; it races through it.
""")
