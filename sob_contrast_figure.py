#!/usr/bin/env python3
"""
WITHDRAWN 2026-09-05 — DO NOT USE. Kept only as a record of the defect.

This figure was built for the round-4 response and pulled before submission.
Its headline claim, "no orbit slows down", is false, and the accompanying
claim that no displacement threshold separates the two phases is also false.

The defect. The comparison conditions the guidance phase on low azimuth drift
and then contrasts it with an unconditioned later window, so the later window
is systematically the faster one by construction. Tested without any guidance
definition -- simply asking whether an orbit's hourly displacement ever falls
sharply anywhere in the night -- the picture reverses:

    orbits slowing by >=2x in some hour                 521 of 3,142
    orbits slowing by >=5x                               74
    largest slowdown found                              9.70x
    slowing >=2x to <=5 deg/h INSIDE the corridor       226
    tightest such displacement                          1.95 deg/h

So slow epochs are common, and 226 orbits reach the referee's threshold in the
corridor. What does not occur -- and what the sweep in sob_bigsweep.py actually
tests -- is a slow epoch standing in the required relation to a PRECEDING
guidance leg. The published result (0 of 525 parameter combinations) is
unaffected. Only this figure's framing of it was wrong.

The lesson worth keeping: a contrast between a selected phase and an
unselected one is not evidence about the object's behaviour.

---- original docstring below ----

sob_contrast_figure.py — the guidance/stopping contrast, without a threshold.

Reviewer 2 observes that instantaneous motion perception requires rates far
above any at issue here (~0.03 deg/s with a fixed reference), and concludes
that a finite angular displacement over a finite time is the better metric.
That is right, and it has a consequence he does not draw.

Whatever displacement threshold D is adopted, it must apply to BOTH phases of
Matthew's narrative. The Star is described as going before the Magi
(proegen, imperfect: continuous, ongoing motion) and then as standing still
(estathe, aorist: a punctual change of state). For that contrast to be
perceptible at all, the object's displacement per unit time must fall from
above D during guidance to below D afterwards. A threshold D separating the
two phases must therefore exist.

This figure tests whether any orbit produces such a separation, at any value
of D. For each close-flyby candidate admitting a guidance window under the
loosest settings in the sweep, it plots

    x = median 1 h displacement DURING the guidance leg
    y = minimum 1 h displacement WHILE STOPPING, within the corridor

The comparison is deliberately generous to the hypothesis: a typical value
during guidance against the best value during stopping. Any orbit exhibiting
a perceptible halt must fall well BELOW the diagonal. None does. The median
orbit falls above it -- these objects move faster after guidance than during
it -- and the largest slowdown achieved anywhere is a factor of 1.63.

The result is threshold-free: it cannot be defeated by loosening D, because
no value of D separates the two phases for any orbit in the sample.
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
AZ_LO, DRIFT, GAP = 150.0, 30.0, 8.0        # loosest settings in the sweep


def collect(sky):
    jd = sky[0]
    step = (jd[1] - jd[0]) * 24.0
    dt = np.gradient(jd) * 24.0
    K = max(1, int(round(1.0 / step)))
    cands = pickle.load(open(os.path.join(
        HERE, "_big_stratified_1000000000_2026.pkl"), "rb"))
    out = []
    for c in cands:
        alt, az, night, _ = S.track(*c, sky)
        if alt is None:
            continue
        dalt = np.gradient(alt) / dt
        azu = np.degrees(np.unwrap(np.radians(az)))
        daz = np.gradient(azu) / dt
        sh = np.full_like(dalt, np.inf)
        if len(alt) > K:
            sh[:-K] = S.great_circle(alt[:-K], azu[:-K], alt[K:], azu[K:])
        g = (night & (alt > M.ALT_MIN) & (az >= AZ_LO) & (az <= M.AZ_HI)
             & (np.abs(dalt) > M.ALT_RATE_MIN) & (np.abs(daz) < DRIFT))
        if not g.any():
            continue
        gi = int(np.where(g)[0][-1])
        dg = sh[g]
        dg = dg[np.isfinite(dg)]
        if not dg.size:
            continue
        sl = slice(gi, min(gi + int(np.ceil(GAP / step)) + 1, len(sh)))
        cm = night[sl] & (az[sl] >= AZ_LO) & (az[sl] <= M.AZ_HI)
        ds = sh[sl][cm]
        ds = ds[np.isfinite(ds)]
        if not ds.size:
            continue
        out.append((float(np.median(dg)), float(ds.min())))
    return np.array(out)


def main():
    sky = M.precompute_sky(M.JD_PRIMARY)
    a = collect(sky)
    x, y = a[:, 0], a[:, 1]
    ratio = x / np.maximum(y, 1e-9)

    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    lim = [0, max(x.max(), y.max()) * 1.06]

    # region in which a perceptible halt would lie
    ax.fill_between(lim, 0, [l * 0.5 for l in lim], color="#DFF0DA", zorder=0)
    ax.text(lim[1] * 0.62, lim[1] * 0.135,
            "a perceptible halt\nwould lie in here\n(stopping motion less\n"
            "than half the guidance motion)",
            fontsize=9, color="#2E6B34", ha="center", va="center")

    ax.plot(lim, lim, color="0.35", lw=1.4, ls="--", zorder=2,
            label="no change in speed")
    ax.scatter(x, y, s=13, color="steelblue", alpha=0.45, edgecolors="none",
               zorder=3, label=f"close-flyby orbits (n = {len(a):,})")

    b = int(np.argmax(x - y))
    ax.scatter([x[b]], [y[b]], s=110, marker="*", color="crimson", zorder=5,
               label=f"largest slowdown found ({ratio[b]:.2f}$\\times$)")

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("median 1 h displacement during guidance  [deg]", fontsize=10)
    ax.set_ylabel("minimum 1 h displacement while stopping  [deg]", fontsize=10)
    ax.set_title("No orbit slows down\n"
                 "Every candidate lies on or above the diagonal: the object "
                 "moves\nno more slowly after guidance than during it",
                 fontsize=11.5)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)

    fig.text(0.5, 0.005,
             f"median contrast ratio {np.median(ratio):.2f}$\\times$   |   "
             f"orbits achieving even a 2$\\times$ slowdown: "
             f"{int((ratio > 2).sum())} of {len(a):,}   |   "
             f"comparison is generous: typical guidance motion against "
             f"best-case stopping motion",
             ha="center", fontsize=8.6, color="0.25")

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "pdf"):
        p = os.path.join(HERE, f"fig_contrast.{ext}")
        fig.savefig(p, dpi=190 if ext == "png" else None, bbox_inches="tight")
        print("wrote", p)

    print(f"\n  n = {len(a)}")
    print(f"  guidance median {np.median(x):.2f} deg, stopping median "
          f"{np.median(y):.2f} deg")
    print(f"  contrast ratio: median {np.median(ratio):.2f}x, "
          f"best {ratio.max():.2f}x")
    print(f"  orbits with >=2x slowdown: {int((ratio > 2).sum())}")


if __name__ == "__main__":
    main()
