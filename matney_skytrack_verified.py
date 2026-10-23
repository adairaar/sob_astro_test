#!/usr/bin/env python3
"""
matney_skytrack_verified.py — sky track of Matney's 5 BCE comet, verified frames.

Regenerates the table and figure in the appendix of matney_response.tex
(the response to Matney 2025 submitted to *The Observatory*).

Supersedes the sky-track portion of sob_calcs.py, which had two defects that
matter at this object's extreme proximity (geocentric distance ~0.0024 AU,
i.e. inside the lunar orbit):

  1. Earth's position came from astropy's bundled ephemeris, which is only
     claimed valid 1900-2100 and is off by ~43" at 5 BCE.  At 3.6e5 km that
     displaces the sightline by several degrees.

  2. The comet's direction was passed to astropy as a SkyCoord with no
     distance, so no topocentric parallax was applied.  At this range the
     observer's offset from Earth's centre moves the object by 0.1-0.43 deg.

Both are fixed here: Earth from VSOP87, and the full geocentric->topocentric
chain via sob_frames, which agrees with astropy's rigorous AltAz transform
(given the same Earth position and the distance supplied) to 0.12 arcmin.

Outputs
    matney_skytrack_table.tex     - LaTeX table body, ready to paste
    matney_fig6_clean.png/.pdf    - two-panel figure
    matney_skytrack_verified.csv  - full track
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sob_frames as sf

# ── Matney (2025) orbital elements, J2000 ecliptic ───────────────────────────
MATNEY = dict(q=0.0407, e=1.0, i_deg=1.39, Om_deg=93.01, om_deg=346.25,
              T_peri_tt=1719785.565)

JD_NOON = 1719755.898        # local solar noon at Jerusalem, 8 June 5 BCE
DT_SEC = 10572.0             # TT - UT1 at this epoch, seconds
DT_DAYS = DT_SEC / 86400.0

LAT, LON = 31.70, 35.20      # Jerusalem
ROAD_AZ = 206.0              # Jerusalem -> Bethlehem bearing
JOURNEY = (8.0, 10.0)        # local hours, per Matney's own narrative
ZENITH_CUT = 80.0            # 'near-zenith' definition used in the paper


def track(local_hours):
    """Sky track at the given local solar hours. Returns dict of arrays."""
    jd_tdb = JD_NOON + (np.asarray(local_hours, dtype=float) - 12.0) / 24.0
    jd_ut1 = jd_tdb - DT_DAYS
    o = sf.observe(jd_ut1, MATNEY, lat_deg=LAT, lon_deg_east=LON,
                   delta_t_days=DT_DAYS, topocentric=True)
    o['local_h'] = np.asarray(local_hours, dtype=float)
    return o


def main():
    # ── fine track across the day ────────────────────────────────────────────
    hours = np.arange(3.0, 13.0 + 1e-9, 1.0 / 120.0)     # 30-second steps
    t = track(hours)
    alt, az = t['alt_deg'], t['az_deg']
    dist = t['dist_au']

    k_min = int(np.argmin(dist))
    h_min, d_min = hours[k_min], dist[k_min]

    # ── table rows ───────────────────────────────────────────────────────────
    rows = []
    for h in (4.0, 8.0, 10.0, h_min, 12.0):
        r = track([h])
        rows.append((h,
                     float(r['alt_deg'][0]), float(r['az_deg'][0]),
                     float(r['dist_au'][0]),
                     float(r['ra_date'][0]) * sf.R2D / 15.0,
                     float(r['dec_date'][0]) * sf.R2D))

    def compass(a):
        pts = ['N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE',
               'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW']
        return pts[int((a % 360) / 22.5 + 0.5) % 16]

    lines = []
    for h, a, z, d, ra, dec in rows:
        label = (f"${{\\approx}}${int(h):02d}:{int(round((h % 1) * 60)):02d}"
                 if abs(h - round(h)) > 1e-6
                 else f"{int(h):02d}:00")
        zen = r"\textbf{Yes}" if a > ZENITH_CUT else "No"
        lines.append(f"{label} & {a:.1f} & {z:.0f} ({compass(z)}) & "
                     f"{d:.5f} & {ra:05.2f} & ${dec:+.1f}$ & {zen} \\\\")
    table_body = "\n".join(lines)

    with open('matney_skytrack_table.tex', 'w') as fh:
        fh.write(table_body + "\n")

    # ── journey-window azimuth behaviour ─────────────────────────────────────
    a8 = float(track([JOURNEY[0]])['az_deg'][0])
    a10 = float(track([JOURNEY[1]])['az_deg'][0])
    drift = a10 - a8
    rate = drift / (JOURNEY[1] - JOURNEY[0])

    in_win = (hours >= JOURNEY[0]) & (hours <= JOURNEY[1])
    cross = np.where(np.diff(np.sign(az[in_win] - ROAD_AZ)) != 0)[0]
    h_cross = hours[in_win][cross[0] + 1] if len(cross) else None

    hi = alt > ZENITH_CUT
    h_zen = hours[np.where(hi)[0][0]] if hi.any() else None

    print("=" * 78)
    print("MATNEY COMET SKY TRACK — verified frames (VSOP87 + precession + parallax)")
    print("=" * 78)
    print(f"  {'local':>7} {'alt':>7} {'az':>8} {'Delta(AU)':>10} {'RA(h)':>7} {'Dec':>7}")
    for h, a, z, d, ra, dec in rows:
        print(f"  {h:7.2f} {a:7.1f} {z:8.1f} {d:10.5f} {ra:7.2f} {dec:+7.1f}")
    print()
    print(f"  closest approach : {d_min:.5f} AU "
          f"({d_min * sf.AU_KM:,.0f} km) at local {h_min:.2f} h")
    print(f"  journey window   : az {a8:.1f} -> {a10:.1f} deg "
          f"(drift {drift:+.1f} deg, {rate:+.1f} deg/h)")
    print(f"  crosses {ROAD_AZ:.0f} deg  : "
          + (f"once, at local {h_cross:.2f} h" if h_cross else "not in window"))
    print(f"  exceeds {ZENITH_CUT:.0f} deg alt: "
          + (f"from local {h_zen:.2f} h (max {alt.max():.1f} deg)" if h_zen else "never"))
    print()
    print("  LaTeX table body written to matney_skytrack_table.tex")

    # ── figure ───────────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))

    ax1.plot(hours, alt, color='#1f77b4', lw=2.2, zorder=3)
    ax1.axvspan(*JOURNEY, color='#ff7f0e', alpha=0.20, zorder=1,
                label="Magi's journey (08:00–10:00)")
    ax1.axvline(h_min, color='#d62728', ls=':', lw=1.8, zorder=2,
                label=fr"closest approach $\Delta_{{\min}}={d_min:.5f}$ au")
    ax1.axhline(ZENITH_CUT, color='grey', ls='--', lw=1.0, zorder=2,
                label=f'near-zenith ({ZENITH_CUT:.0f}°)')
    ax1.set_xlabel('Local solar time (h)')
    ax1.set_ylabel('Altitude (°)')
    ax1.set_title('Altitude vs local time')
    ax1.set_xlim(hours[0], hours[-1])
    ax1.set_ylim(0, 90)
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8, loc='upper left')

    vis = alt > 0
    sc = ax2.scatter(az[vis], alt[vis], c=hours[vis], cmap='viridis', s=7, zorder=3)
    ax2.axvline(ROAD_AZ, color='green', ls='--', lw=1.6, zorder=2,
                label=f'Jerusalem→Bethlehem bearing ({ROAD_AZ:.0f}°)')
    for h in JOURNEY:
        r = track([h])
        ax2.plot(float(r['az_deg'][0]), float(r['alt_deg'][0]), 'o',
                 ms=9, mfc='none', mec='#d62728', mew=2, zorder=4)
        ax2.annotate(f'{int(h):02d}:00',
                     (float(r['az_deg'][0]), float(r['alt_deg'][0])),
                     textcoords='offset points', xytext=(8, -12),
                     fontsize=9, color='#d62728')
    plt.colorbar(sc, ax=ax2, label='Local solar time (h)')
    ax2.set_xlabel('Azimuth (°)')
    ax2.set_ylabel('Altitude (°)')
    ax2.set_title('Alt–az sky track')
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8, loc='upper left')

    fig.suptitle("Matney's 8 June comet from Jerusalem, 8 June 5 BCE "
                 "(VSOP87 ephemeris, topocentric)", fontsize=11)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(f'matney_fig6_clean.{ext}', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print("  figure written to matney_fig6_clean.png / .pdf")

    # ── csv ──────────────────────────────────────────────────────────────────
    np.savetxt('matney_skytrack_verified.csv',
               np.column_stack([hours, alt, az, dist,
                                t['ra_date'] * sf.R2D / 15.0,
                                t['dec_date'] * sf.R2D]),
               delimiter=',', header='local_h,alt_deg,az_deg,dist_au,ra_h,dec_deg',
               comments='', fmt='%.6f')
    print("  full track written to matney_skytrack_verified.csv")


if __name__ == '__main__':
    main()
