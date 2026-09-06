#!/usr/bin/env python3
"""
sob_falsepositive_figure.py — anatomy of the one orbit that "passed".

Of 6,116 close-flyby orbits drawn from 10^9 samples, exactly one satisfied a
displacement stopping criterion of 5 deg per hour, and only under a corridor
widened to 40-60 deg, an azimuth-drift tolerance raised tenfold, a guidance
leg cut to fifteen minutes, and a stopping window of eight hours.

This figure shows that it is a false positive, on three counts, all verified
against the track rather than asserted:

  (i)  The object never stops. Its total apparent motion reaches a minimum of
       3.85 deg/h and never falls below it at any point in the night.
  (ii) What vanishes is only the altitude component. At the best-displacement
       epoch 98% of the residual motion is azimuthal (|dalt/dt| = 0.80,
       |daz/dt| cos alt = 3.92 deg/h); the altitude turning point proper lies
       half an hour later still, and even there the total rate is 3.85 deg/h.
  (iii) It happens in the wrong place at the wrong time: azimuth 233 deg,
       24 deg west of the guidance corridor, 3.8 h after guidance ended.

A criterion that scores this track as a "stop" is measuring the moment a
descending object levels out, not a body coming to rest over a house.
"""

import os
import pickle
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sob_close_flyby_mc as M                                 # noqa: E402
import sob_criteria_sweep as S                                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
IDX, AZ_LO, K = 1178, 170.0, 12


def main():
    sky = M.precompute_sky(M.JD_PRIMARY)
    jd = sky[0]
    dt = np.gradient(jd) * 24.0
    hr = (jd - jd[0]) * 24.0
    c = pickle.load(open(os.path.join(
        HERE, "_big_stratified_1000000000_2026.pkl"), "rb"))[IDX]
    q, e, i_r, Om, om_, T = c
    alt, az, night, dm = S.track(*c, sky)
    dalt = np.gradient(alt) / dt
    azu = np.degrees(np.unwrap(np.radians(az)))
    daz = np.gradient(azu) / dt
    omega = np.sqrt(dalt**2 + (daz * np.cos(np.radians(alt)))**2)
    shift = np.full_like(omega, np.inf)
    shift[:-K] = S.great_circle(alt[:-K], azu[:-K], alt[K:], azu[K:])

    g = (night & (alt > M.ALT_MIN) & (az >= AZ_LO) & (az <= M.AZ_HI)
         & (np.abs(dalt) > M.ALT_RATE_MIN) & (np.abs(daz) < 15.0))
    gi = int(np.where(g)[0][-1])
    j = int(np.nanargmin(np.where(np.isfinite(shift), shift, np.inf)))

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.4, 5.6))

    # (a) sky track
    ax.axvspan(AZ_LO, M.AZ_HI, color="0.88", zorder=0)
    ax.text((AZ_LO + M.AZ_HI) / 2, 21, "guidance corridor", ha="center",
            fontsize=9, color="0.35")
    m = night & (alt > 0)
    ax.plot(az[m], alt[m], lw=1.2, color="0.55")
    ax.plot(az[g], alt[g], lw=3.0, color="steelblue", label="guidance window")
    ax.scatter([az[gi]], [alt[gi]], s=70, color="steelblue", zorder=5,
               label=f"guidance ends (az {az[gi]:.0f}°)")
    ax.scatter([az[j]], [alt[j]], s=140, marker="*", color="crimson", zorder=6,
               label=f'best 1 h displacement (az {az[j]:.0f}°, {shift[j]:.2f}°)')
    ax.annotate("", xy=(az[j], alt[j] + 3), xytext=(az[gi], alt[gi] + 3),
                arrowprops=dict(arrowstyle="->", color="crimson", lw=1.4))
    ax.text((az[gi] + az[j]) / 2, (alt[gi] + alt[j]) / 2 - 4.5,
            f"{az[j]-az[gi]:.0f}° further west,\n"
            f"{(hr[j]-hr[gi]):.1f} h later", ha="center", fontsize=8.5,
            color="crimson")
    ax.axhline(M.ALT_MIN, color="0.6", ls=":", lw=1)
    ax.set_xlabel("azimuth [deg]")
    ax.set_ylabel("altitude [deg]")
    ax.set_title("(a) the apparent path over Bethlehem", fontsize=10.5)
    ax.legend(fontsize=8, loc="lower left")
    ax.set_ylim(0, None)

    # (b) the two components of the motion
    w = slice(max(gi - 6, 0), min(j + 20, len(hr)))
    bx.plot(hr[w], np.abs(dalt[w]), lw=2, color="seagreen",
            label="|d(alt)/dt|  — vanishes at the turning point")
    bx.plot(hr[w], np.abs(daz[w]) * np.cos(np.radians(alt[w])), lw=2,
            color="darkorange", label="|d(az)/dt| cos(alt)  — never stops")
    bx.plot(hr[w], omega[w], lw=2.4, color="crimson", label="total apparent motion")
    om_min = float(np.nanmin(omega[night]))
    bx.axvline(hr[j], color="crimson", ls="--", lw=1.2)
    bx.text(hr[j] + 0.06, 7.6, 'best 1 h\ndisplacement', color="crimson",
            fontsize=8.2, va="top")
    bx.axvline(hr[gi], color="steelblue", ls="--", lw=1.2)
    bx.text(hr[gi] + 0.06, 9.3, "guidance ends", color="steelblue",
            fontsize=8.2, va="top")
    bx.axhspan(0, 5, color="0.9", zorder=0)
    bx.text(hr[w][0] + 0.1, 4.6, "criterion band ($\\leq$ 5)", fontsize=8.5,
            color="0.35")
    bx.axhline(om_min, color="crimson", lw=1.3, ls="-.", alpha=0.85)
    bx.text(hr[w][0] + 0.15, om_min - 0.55,
            f"total motion never falls below {om_min:.2f}°/h",
            ha="left", fontsize=8.8, color="crimson", fontweight="bold")
    bx.set_xlabel("hours from start of window")
    bx.set_ylabel("apparent angular rate [deg/h]")
    bx.set_title("(b) the altitude component vanishes; the object never stops",
                 fontsize=10.5)
    bx.legend(fontsize=8.2, loc="upper right")
    bx.grid(alpha=0.3)

    a_hyp = q / (1 - e)
    fig.text(0.5, 0.005,
             f"q = {q:.3f} au   e = {e:.2f} (hyperbolic, "
             f"v$_\\infty$ = {29.785/np.sqrt(abs(a_hyp)):.0f} km/s)   "
             f"i = {np.degrees(i_r):.1f}°   "
             f"closest approach {dm*1.496e8:,.0f} km "
             f"({dm*1.496e8/384400:.1f} lunar distances)",
             ha="center", fontsize=9, color="0.25")

    fig.tight_layout(rect=(0, 0.035, 1, 1))
    for ext in ("png", "pdf"):
        p = os.path.join(HERE, f"fig_false_positive.{ext}")
        fig.savefig(p, dpi=190 if ext == "png" else None, bbox_inches="tight")
        print("wrote", p)


if __name__ == "__main__":
    main()
