#!/usr/bin/env python3
"""
sob_flyby_figure.py  —  Visualisation of the best close-flyby orbit
Aaron Adair / 2026

Three-panel figure:
  (a) Stereographic sky chart (planetarium view, N up / E left):
      full night track of the comet, colour-coded by ground-frame apparent
      velocity ω_app, with road corridor [190°,210°] and guidance/stopping
      windows marked.
  (b) Azimuth and altitude time series, with corridor and twilight bands.
  (c) ω_app and ω_ICRS time series, with threshold lines and ω_sid.

Saves: figure_flyby_bestorbit.pdf  (and .png for quick view)
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.colorbar import ColorbarBase
from astropy.time import Time
from astropy.coordinates import get_body_barycentric, get_sun
import astropy.units as u

# ── Best orbit parameters ─────────────────────────────────────────────────────
PERI    = 0.998156
ECC     = 0.709313
I_R     = np.radians(2.172)
OM_R    = np.radians(95.196)
OM_r    = np.radians(204.638)
DT_DAYS = 13.723          # T_peri = JD_PRIMARY + DT_DAYS

# ── Constants ─────────────────────────────────────────────────────────────────
K_GAUSS = 0.01720209895
MU      = K_GAUSS**2
EPS     = np.radians(23.4392911)
AU_KM   = 1.495978707e8
R_EQ    = np.array([[1,0,0],
                    [0, np.cos(EPS), -np.sin(EPS)],
                    [0, np.sin(EPS),  np.cos(EPS)]])

LAT_DEG = 31.70;  LON_DEG = 35.20
LAT     = np.radians(LAT_DEG)
JD_PRIMARY   = 1719755.898
DT_TDB_UTC   = 10572.0

OMEGA_SID        = 15.041   # °/h
AZ_LO, AZ_HI    = 190.0, 210.0
AZ_DRIFT_THRESH  = 1.5
ALT_RATE_MIN     = 0.3
STOP_THRESH      = 0.30
STEP_MIN         = 1
WINDOW_HRS       = 14


# ── Orbital mechanics ─────────────────────────────────────────────────────────

def kepler_elliptic(Mv, ecc, tol=1e-12):
    Ev = Mv.copy()
    for _ in range(60):
        dE = (Mv - Ev + ecc*np.sin(Ev)) / (1 - ecc*np.cos(Ev))
        Ev += dE
        if np.max(np.abs(dE)) < tol:
            break
    return Ev

def comet_pos(dt_days):
    dt   = np.atleast_1d(np.asarray(dt_days, float))
    ci, si = np.cos(I_R), np.sin(I_R)
    cO, sO = np.cos(OM_R), np.sin(OM_R)
    co, so = np.cos(OM_r), np.sin(OM_r)
    Pecl = np.array([cO*co - sO*so*ci,  sO*co + cO*so*ci,  so*si])
    Qecl = np.array([-cO*so - sO*co*ci, -sO*so + cO*co*ci, co*si])
    Pvec = R_EQ @ Pecl
    Qvec = R_EQ @ Qecl
    a    = PERI / (1 - ECC)
    n    = np.sqrt(MU / a**3)
    Mv   = (n * dt) % (2*np.pi)
    Ev   = kepler_elliptic(Mv, ECC)
    nu   = 2*np.arctan2(np.sqrt(1+ECC)*np.sin(Ev/2), np.sqrt(1-ECC)*np.cos(Ev/2))
    r    = a * (1 - ECC*np.cos(Ev))
    return np.outer(Pvec, r*np.cos(nu)) + np.outer(Qvec, r*np.sin(nu))

def gmst_rad(jd_utc):
    T = (jd_utc - 2451545.0) / 36525.0
    return np.radians((280.46061837 + 360.98564736629*(jd_utc-2451545.0)
                       + 0.000387933*T**2 - T**3/38710000.0) % 360.0)

def lst_rad(jd_utc):
    return (gmst_rad(jd_utc) + np.radians(LON_DEG)) % (2*np.pi)

def altaz(ra_r, dec_r, lst_r):
    ha      = lst_r - ra_r
    sin_alt = (np.sin(LAT)*np.sin(dec_r)
               + np.cos(LAT)*np.cos(dec_r)*np.cos(ha))
    sin_alt = np.clip(sin_alt, -1, 1)
    alt     = np.arcsin(sin_alt)
    cos_alt = np.cos(alt)
    az_sin  = -np.cos(dec_r)*np.sin(ha)
    az_cos  = ((np.sin(dec_r) - np.sin(LAT)*sin_alt)
               / (np.cos(LAT)*np.maximum(cos_alt, 1e-9)))
    az      = np.arctan2(az_sin, az_cos) % (2*np.pi)
    return np.degrees(alt), np.degrees(az)


# ── Compute orbit data ────────────────────────────────────────────────────────

def compute_orbit():
    T_peri_jd = JD_PRIMARY + DT_DAYS
    jd_utc_midnight = (JD_PRIMARY - DT_TDB_UTC/86400.0) - 0.5
    dt_h   = np.arange(-WINDOW_HRS/2, WINDOW_HRS/2, STEP_MIN/60.0)
    jd_utc = jd_utc_midnight + dt_h / 24.0
    jd_tdb = jd_utc + DT_TDB_UTC / 86400.0

    times_tdb = Time(jd_tdb, format='jd', scale='tdb')
    times_utc = Time(jd_utc, format='jd', scale='utc')
    eb        = get_body_barycentric('earth', times_tdb)
    sb        = get_body_barycentric('sun',   times_tdb)
    earth_pos = (eb - sb).xyz.to(u.AU).value
    sun_icrs  = get_sun(times_utc)
    sun_ra    = sun_icrs.ra.rad
    sun_dec   = sun_icrs.dec.rad
    lst       = lst_rad(jd_utc)
    sun_alt, sun_az = altaz(sun_ra, sun_dec, lst)

    dt_from_peri = jd_tdb - T_peri_jd
    r_c  = comet_pos(dt_from_peri)
    rho  = r_c - earth_pos
    d    = np.linalg.norm(rho, axis=0)
    rhat = rho / d
    ra_r  = np.arctan2(rhat[1], rhat[0]) % (2*np.pi)
    dec_r = np.arcsin(np.clip(rhat[2], -1, 1))
    alt, az = altaz(ra_r, dec_r, lst)

    h_rad = STEP_MIN / 60.0
    dalt  = np.gradient(alt, h_rad)
    daz   = np.zeros_like(az)
    for k in range(1, len(az)-1):
        da = az[k+1] - az[k-1]
        if da > 180:  da -= 360
        if da < -180: da += 360
        daz[k] = da / (2*h_rad)
    daz[0] = daz[1];  daz[-1] = daz[-2]

    omega_app = np.sqrt(dalt**2 + (daz*np.cos(np.radians(alt)))**2)

    # ω_ICRS via finite differences of geocentric vector
    rho_km   = rho * AU_KM
    dt_sec   = STEP_MIN * 60.0
    v_rho    = np.gradient(rho_km, dt_sec, axis=1)
    v_r_rad  = np.sum(v_rho * rhat, axis=0)
    v_transv = np.sqrt(np.maximum(np.sum(v_rho**2, axis=0) - v_r_rad**2, 0))
    omega_icrs = np.degrees(v_transv / (d * AU_KM)) * 3600.0

    t_event_h = (jd_utc - (JD_PRIMARY - DT_TDB_UTC/86400.0)) * 24.0

    # Masks
    night   = sun_alt < -6.0
    visible = (alt > 10.0) & night
    in_az   = (az >= AZ_LO) & (az <= AZ_HI)
    guidance_mask = visible & in_az & (np.abs(daz) < AZ_DRIFT_THRESH) & (np.abs(dalt) > ALT_RATE_MIN)
    stopping_mask = visible & (omega_app < STOP_THRESH)

    return dict(
        t=t_event_h, d=d, alt=alt, az=az, dalt=dalt, daz=daz,
        omega_app=omega_app, omega_icrs=omega_icrs,
        sun_alt=sun_alt, night=night, visible=visible,
        guidance_mask=guidance_mask, stopping_mask=stopping_mask,
    )


# ── Stereographic sky-chart helpers ──────────────────────────────────────────

def sky_xy(alt_deg, az_deg):
    """
    Zenith-centred stereographic projection, N-up / E-left (standard star chart).
    r = 90 - alt,  azimuth → angle with North at top, clockwise = West→East map.
    """
    r  = 90.0 - np.asarray(alt_deg)
    az = np.radians(np.asarray(az_deg))
    # N at top: y = r cos(az); E at LEFT: x = -r sin(az)
    x = -r * np.sin(az)
    y =  r * np.cos(az)
    return x, y


# ── Main figure ───────────────────────────────────────────────────────────────

def make_figure(data):
    fig = plt.figure(figsize=(15, 10))
    gs  = fig.add_gridspec(2, 2,
                           left=0.06, right=0.97,
                           top=0.93,  bottom=0.08,
                           wspace=0.32, hspace=0.40,
                           width_ratios=[1, 1.35])

    ax_sky  = fig.add_subplot(gs[:, 0])    # left: sky chart (spans both rows)
    ax_az   = fig.add_subplot(gs[0, 1])    # top right: azimuth / altitude
    ax_om   = fig.add_subplot(gs[1, 1])    # bottom right: ω traces

    cmap_speed = plt.cm.plasma_r
    vmin, vmax = 0.0, 15.0

    # ─── Panel (a): Planetarium sky chart ────────────────────────────────────
    ax = ax_sky
    ax.set_aspect('equal')
    ax.set_xlim(-95, 95)
    ax.set_ylim(-95, 95)
    ax.axis('off')

    # Horizon and altitude rings
    for r_ring, label, ls in [(90, '0° (horizon)', '-'), (60, '30°', '--'), (30, '60°', ':')]:
        circle = plt.Circle((0,0), r_ring, color='#aaaaaa',
                             fill=False, lw=0.8, ls=ls, zorder=1)
        ax.add_patch(circle)
        if r_ring < 90:
            ax.text(0, r_ring+2, label, ha='center', va='bottom',
                    fontsize=7, color='#888888')

    # Cardinal direction labels
    for az_label, label_txt in [(0,'N'), (90,'E'), (180,'S'), (270,'W')]:
        xc, yc = sky_xy(0, az_label)
        ax.text(xc*1.06, yc*1.06, label_txt, ha='center', va='center',
                fontsize=9, fontweight='bold', color='#555555')

    # Road corridor wedge [190°,210°] — fill from zenith to horizon
    # Build wedge boundary at horizon and fill
    az_wedge  = np.linspace(AZ_LO, AZ_HI, 80)
    x_lo, y_lo = sky_xy(np.zeros_like(az_wedge), az_wedge)
    x_wedge = np.concatenate([[0], x_lo, [0]])
    y_wedge = np.concatenate([[0], y_lo, [0]])
    ax.fill(x_wedge, y_wedge, color='gold', alpha=0.25, zorder=2,
            label='Road corridor [190°–210°]')
    # Corridor boundary lines
    for az_edge in [AZ_LO, AZ_HI]:
        xe, ye = sky_xy([10, 0], [az_edge, az_edge])
        ax.plot(xe, ye, color='goldenrod', lw=1.2, ls='--', zorder=3)

    # Below-horizon region (hatched arc) to show Sun direction
    # Draw a hatched semicircle below the horizon (S side, since Sun is near S in daylight)
    theta_horiz = np.linspace(np.radians(90), np.radians(270), 200)  # E through S to W
    x_h = 90*np.cos(theta_horiz)
    y_h = 90*np.sin(theta_horiz)
    ax.fill(np.concatenate([[0], x_h, [0]]),
            np.concatenate([[0], y_h, [0]]),
            color='#fffbe6', alpha=0.0, zorder=0)   # just a region reference

    # Full visible track coloured by ω_app
    vis = data['visible']
    t_v   = data['t'][vis]
    alt_v = data['alt'][vis]
    az_v  = data['az'][vis]
    om_v  = data['omega_app'][vis]

    x_v, y_v = sky_xy(alt_v, az_v)

    # Draw coloured segments
    pts  = np.array([x_v, y_v]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    lc   = LineCollection(segs, cmap=cmap_speed, norm=norm,
                          linewidth=3.5, zorder=5, capstyle='round')
    lc.set_array(om_v[:-1])
    ax.add_collection(lc)

    # Comet direction arrows every ~90 min (roughly every 90 steps at 1-min resolution)
    stride = 90
    for k in range(stride//2, len(x_v)-1, stride):
        dx = x_v[k+1] - x_v[k-1]
        dy = y_v[k+1] - y_v[k-1]
        norm_d = np.hypot(dx, dy)
        if norm_d > 0:
            ax.annotate('', xy=(x_v[k]+dx/norm_d*4, y_v[k]+dy/norm_d*4),
                        xytext=(x_v[k], y_v[k]),
                        arrowprops=dict(arrowstyle='->', color='#555555',
                                        lw=1.2),
                        zorder=6)

    # Guidance segments
    guide = data['guidance_mask'][vis]
    if np.any(guide):
        ax.scatter(x_v[guide], y_v[guide], s=18, color='limegreen',
                   zorder=7, label='Guidance window', edgecolors='none')

    # Stopping point (minimum ω_app in visible window)
    stop = data['stopping_mask'][vis]
    if np.any(stop):
        idx_stop = np.argmin(om_v[stop])
        xs = x_v[stop][idx_stop]
        ys = y_v[stop][idx_stop]
        ax.scatter([xs], [ys], s=120, marker='*', color='red',
                   zorder=8, label='Stopping point (min ω_app)')
        ax.annotate('Stopping\n(ω_app = 0.21°/h)', xy=(xs, ys),
                    xytext=(xs-10, ys-12),
                    fontsize=7.5, color='red',
                    arrowprops=dict(arrowstyle='->', color='red', lw=1.0),
                    zorder=9)

    # Rise label
    ax.text(x_v[0]+3, y_v[0]-4, 'Rise\n(t=−11.5h)', fontsize=7,
            color='steelblue', ha='left', va='top', zorder=9)
    # Closest approach label
    idx_ca = np.argmin(data['d'][vis])
    ax.text(x_v[idx_ca]-3, y_v[idx_ca]+4, f"Closest\napproach\n(d={data['d'][vis][idx_ca]*AU_KM/384400:.2f}×lunar)",
            fontsize=7, color='darkred', ha='right', va='bottom', zorder=9)

    # Colourbar for ω_app
    sm = plt.cm.ScalarMappable(cmap=cmap_speed, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.03, 0.12, 0.016, 0.30])
    cb = fig.colorbar(sm, cax=cbar_ax)
    cb.set_label('ω_app (°/h, ground frame)', fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # Legend for sky chart
    from matplotlib.lines import Line2D
    legend_elements = [
        mpatches.Patch(facecolor='gold', alpha=0.5, label='Road corridor\n[190°–210°]'),
        Line2D([0], [0], color=cmap_speed(0.15), lw=3, label='Comet track (slow)'),
        Line2D([0], [0], color=cmap_speed(0.80), lw=3, label='Comet track (fast)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='limegreen',
               markersize=7, label='Guidance window'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='red',
               markersize=10, label='Stopping point'),
    ]
    ax.legend(handles=legend_elements, loc='lower left',
              fontsize=7.5, framealpha=0.9, borderpad=0.8)

    ax.set_title('(a)  Sky track: Bethlehem, 5 BCE\n'
                 'Best close-flyby orbit — N up, E left',
                 fontsize=10, pad=6)

    # Zenith dot
    ax.plot(0, 0, '+', color='#888888', ms=8, zorder=4)
    ax.text(2, 2, 'Zenith', fontsize=7, color='#888888')

    # ─── Panel (b): Azimuth and altitude vs time ──────────────────────────────
    ax = ax_az
    t  = data['t']
    vis_t = t[data['visible']]

    # Background shading: twilight / daylight
    ax.axhspan(AZ_LO, AZ_HI, color='gold', alpha=0.18, label='Road corridor')

    ax.plot(t[data['visible']], data['az'][data['visible']],
            color='steelblue', lw=1.8, label='Azimuth (°)', zorder=4)

    # Guidance and stopping regions
    for k, (mask, col, lbl) in enumerate([
            (data['guidance_mask'], 'limegreen', 'Guidance'),
            (data['stopping_mask'], 'red',       'Stopping'),
    ]):
        if np.any(mask):
            # Fill vertical spans
            in_mask = mask & data['visible']
            changes = np.diff(in_mask.astype(int))
            starts  = np.where(changes == 1)[0] + 1
            ends    = np.where(changes == -1)[0] + 1
            if in_mask[0]:  starts = np.insert(starts, 0, 0)
            if in_mask[-1]: ends   = np.append(ends, len(in_mask)-1)
            for s, e in zip(starts, ends):
                ax.axvspan(t[s], t[e], color=col, alpha=0.20 if k==0 else 0.35,
                           zorder=3, label=lbl if s == starts[0] else '')

    # Altitude on second y-axis
    ax2 = ax.twinx()
    ax2.plot(t[data['visible']], data['alt'][data['visible']],
             color='darkorange', lw=1.5, ls='--', label='Altitude (°)', zorder=4)
    ax2.set_ylabel('Altitude (°)', color='darkorange', fontsize=9)
    ax2.tick_params(axis='y', colors='darkorange', labelsize=8)
    ax2.set_ylim(0, 70)

    ax.set_xlim(t[data['visible']][0]-0.3, t[data['visible']][-1]+0.3)
    ax.set_ylim(140, 230)
    ax.axhline(AZ_LO, color='goldenrod', lw=0.8, ls=':')
    ax.axhline(AZ_HI, color='goldenrod', lw=0.8, ls=':')
    ax.set_xlabel('Time relative to event epoch (h)', fontsize=9)
    ax.set_ylabel('Azimuth (°)', fontsize=9)
    ax.tick_params(labelsize=8)

    # Combine legends
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=7.5, loc='upper left',
              framealpha=0.9, ncol=2)
    ax.set_title('(b)  Azimuth and altitude vs time', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Night window shading
    night_changes = np.diff(data['night'].astype(int))
    night_starts  = np.where(night_changes == 1)[0] + 1
    night_ends    = np.where(night_changes == -1)[0] + 1
    if data['night'][0]:  night_starts = np.insert(night_starts, 0, 0)
    if data['night'][-1]: night_ends   = np.append(night_ends, len(data['night'])-1)
    for s, e in zip(night_starts, night_ends):
        ax.axvspan(t[s], t[e], color='#e8e8f8', alpha=0.5, zorder=0)

    # ─── Panel (c): ω_app and ω_ICRS vs time ─────────────────────────────────
    ax = ax_om

    ax.fill_between(t[data['visible']],
                    data['omega_app'][data['visible']],
                    alpha=0.15, color='steelblue')
    ax.plot(t[data['visible']], data['omega_app'][data['visible']],
            color='steelblue', lw=2.0, label=r'$\omega_\mathrm{app}$ (ground frame)', zorder=4)
    ax.plot(t[data['visible']], data['omega_icrs'][data['visible']],
            color='darkorange', lw=1.8, ls='--',
            label=r'$\omega_\mathrm{ICRS}$ (stellar frame)', zorder=4)

    # Thresholds
    ax.axhline(OMEGA_SID, color='black', lw=1.0, ls=':', zorder=3,
               label=f'Sidereal rate ({OMEGA_SID:.1f}°/h)')
    ax.axhline(STOP_THRESH, color='red', lw=1.2, ls='--', zorder=3,
               label=f'Stopping threshold ({STOP_THRESH}°/h)')
    ax.axhline(AZ_DRIFT_THRESH, color='limegreen', lw=1.2, ls='--', zorder=3,
               label=f'Drift threshold ({AZ_DRIFT_THRESH}°/h)')

    # Guidance and stopping windows
    for mask, col in [(data['guidance_mask'], 'limegreen'),
                      (data['stopping_mask'], 'red')]:
        in_mask = mask & data['visible']
        changes = np.diff(in_mask.astype(int))
        starts  = np.where(changes == 1)[0] + 1
        ends    = np.where(changes == -1)[0] + 1
        if in_mask[0]:  starts = np.insert(starts, 0, 0)
        if in_mask[-1]: ends   = np.append(ends, len(in_mask)-1)
        for s, e in zip(starts, ends):
            ax.axvspan(t[s], t[e], color=col, alpha=0.20, zorder=2)

    # Annotate the minimum ω_app
    vis = data['visible']
    idx_min = np.argmin(data['omega_app'][vis])
    t_min   = t[vis][idx_min]
    om_min  = data['omega_app'][vis][idx_min]
    ax.annotate(f'min ω_app = {om_min:.2f}°/h\n(≈ stopping point)',
                xy=(t_min, om_min),
                xytext=(t_min-1.5, om_min+1.5),
                fontsize=7.5, color='red',
                arrowprops=dict(arrowstyle='->', color='red', lw=1.0))

    ax.set_xlim(t[data['visible']][0]-0.3, t[data['visible']][-1]+0.3)
    ax.set_ylim(0, 18)
    ax.set_xlabel('Time relative to event epoch (h)', fontsize=9)
    ax.set_ylabel('Angular velocity (°/h)', fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7.5, loc='upper left', framealpha=0.9, ncol=2)
    ax.set_title('(c)  Apparent angular velocities vs time', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Night shading on panel c too
    for s, e in zip(night_starts, night_ends):
        ax.axvspan(t[s], t[e], color='#e8e8f8', alpha=0.5, zorder=0)

    # ── Overall title ─────────────────────────────────────────────────────────
    fig.suptitle(
        'Best close-flyby comet satisfying the literal Matt 2:9 criteria\n'
        r'$q=0.998$ AU, $e=0.709$, $i=2.2°$ — guidance: 0.58 h (58% of 1 h minimum)',
        fontsize=11, y=0.98)

    return fig


if __name__ == '__main__':
    import os
    print("Computing orbit data...")
    data = compute_orbit()
    print(f"  Visible: {data['visible'].sum()} min  |  "
          f"Guidance: {data['guidance_mask'].sum()} min  |  "
          f"Stopping: {data['stopping_mask'].sum()} min")

    print("Building figure...")
    fig = make_figure(data)

    out_dir = os.path.dirname(os.path.abspath(__file__))
    pdf_path = os.path.join(out_dir, 'figure_flyby_bestorbit.pdf')
    png_path = os.path.join(out_dir, 'figure_flyby_bestorbit.png')

    fig.savefig(pdf_path, dpi=200, bbox_inches='tight')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {os.path.basename(pdf_path)}")
    print(f"  Saved: {os.path.basename(png_path)}")
    plt.close(fig)
    print("Done.")
