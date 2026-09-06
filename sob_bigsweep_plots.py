#!/usr/bin/env python3
"""
sob_bigsweep_plots.py — figures for the criteria sweep.

Four panels, addressing the two objections a referee can raise against a
negative result of this kind.

  (a) Where guidance is possible at all, as a function of the corridor and the
      azimuth-drift tolerance.  Shows that the drift tolerance, not the
      stopping criterion, is what binds.
  (b) The tightest stopping metric reached, as a function of drift tolerance
      and the guidance duration demanded.  Shows how far the result has to be
      pushed before anything approaches the referee's range.
  (c) Cumulative distributions of both stopping metrics over all close-flyby
      candidates.  Shows the floor is a property of the population, not of a
      handful of orbits.
  (d) Published priors against flat priors.  Shows the conclusion is not an
      artefact of the stratified sampling.
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
REF_AZ_LO, REF_DRIFT, REF_GUIDE, REF_GAP = 150.0, 30.0, 0.25, 8.0


def candidate_metrics(prior, sky):
    """Per-candidate stopping metrics at the loosest configuration."""
    pk = [f for f in os.listdir(HERE)
          if f.startswith(f"_big_{prior}_") and f.endswith(".pkl")]
    if not pk:
        return np.array([]), np.array([])
    pk.sort(key=lambda f: os.path.getsize(os.path.join(HERE, f)))
    cands = pickle.load(open(os.path.join(HERE, pk[-1]), "rb"))
    jd = sky[0]
    step_h = (jd[1] - jd[0]) * 24.0
    dt_h = np.gradient(jd) * 24.0
    k = max(1, int(round(1.0 / step_h)))
    rate, shift = [], []
    for c in cands:
        alt, az, night, _ = S.track(*c, sky)
        if alt is None:
            continue
        dalt = np.gradient(alt) / dt_h
        az_u = np.degrees(np.unwrap(np.radians(az)))
        daz = np.gradient(az_u) / dt_h
        om = np.sqrt(dalt**2 + (daz * np.cos(np.radians(alt)))**2)
        sh = np.full_like(om, np.inf)
        if len(alt) > k:
            sh[:-k] = S.great_circle(alt[:-k], az_u[:-k], alt[k:], az_u[k:])
        g = (night & (alt > M.ALT_MIN) & (az >= REF_AZ_LO) & (az <= M.AZ_HI)
             & (np.abs(dalt) > M.ALT_RATE_MIN) & (np.abs(daz) < REF_DRIFT))
        if not g.any() or g.sum() * step_h < REF_GUIDE:
            continue
        gi = int(np.where(g)[0][-1])
        sl = slice(gi, min(gi + int(np.ceil(REF_GAP / step_h)) + 1, len(om)))
        nm = night[sl]
        if not nm.any():
            continue
        s = sh[sl][nm]
        s = s[np.isfinite(s)]
        rate.append(float(np.min(om[sl][nm])))
        shift.append(float(np.min(s)) if s.size else np.inf)
    return np.array(rate), np.array(shift)


def main():
    f_s = os.path.join(HERE, "bigsweep_stratified.npz")
    if not os.path.exists(f_s):
        raise SystemExit("run sob_bigsweep.py first")
    d = np.load(f_s, allow_pickle=True)
    az_lo, drift = d["az_lo"], d["drift"]
    guide, stop = d["guide"], d["stop"]
    n_guide, min_rate = d["n_guide"], d["min_rate"]
    min_shift_corr = d["min_shift_corr"]
    n_cand, n_mc = int(d["n_cand"]), int(d["n_mc"])

    sky = M.precompute_sky(M.JD_PRIMARY)
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 10.4))

    # (a) guidance counts over corridor x drift, summed over guide/gap
    ax = axes[0, 0]
    g2 = n_guide.sum(axis=(2, 3))
    im = ax.imshow(g2, cmap="viridis", aspect="auto", origin="lower")
    ax.set_xticks(range(len(drift)))
    ax.set_xticklabels([f"{x:g}" for x in drift])
    ax.set_yticks(range(len(az_lo)))
    ax.set_yticklabels([f"{a:.0f}-{M.AZ_HI:.0f}" for a in az_lo])
    ax.set_xlabel("azimuth-drift tolerance [deg/h]")
    ax.set_ylabel("guidance corridor [deg]")
    ax.set_title("(a) orbit-configurations admitting guidance", fontsize=10.5)
    for i in range(g2.shape[0]):
        for j in range(g2.shape[1]):
            ax.text(j, i, f"{g2[i, j]:,}", ha="center", va="center", fontsize=7.5,
                    color="w" if g2[i, j] < g2.max() * 0.6 else "k")
    fig.colorbar(im, ax=ax, pad=0.02)
    ax.axvline(np.argmin(np.abs(drift - 1.5)) + 0.5, color="crimson", lw=1.6)
    ax.text(np.argmin(np.abs(drift - 1.5)) - 0.35, len(az_lo) - 0.4,
            "published", color="crimson", fontsize=8, rotation=90, va="top")

    # (b) tightest rate over drift x guidance duration, at loosest corridor
    ax = axes[0, 1]
    mr = min_shift_corr[-1].min(axis=2)    # loosest corridor, min over gap
    mr = np.where(np.isfinite(mr), mr, np.nan)
    im = ax.imshow(mr.T, cmap="magma_r", aspect="auto", origin="lower",
                   vmin=0, vmax=max(12, np.nanmax(mr)))
    ax.set_xticks(range(len(drift)))
    ax.set_xticklabels([f"{x:g}" for x in drift])
    ax.set_yticks(range(len(guide)))
    ax.set_yticklabels([f"{g:g}" for g in guide])
    ax.set_xlabel("azimuth-drift tolerance [deg/h]")
    ax.set_ylabel("guidance duration demanded [h]")
    ax.set_title(f"(b) tightest 1 h displacement with the stop required\n"
                 f"in the corridor ({az_lo[-1]:.0f}-{M.AZ_HI:.0f} deg); "
                 f"criterion is 5 deg", fontsize=10.5)
    for i in range(mr.shape[0]):
        for j in range(mr.shape[1]):
            v = mr[i, j]
            if np.isfinite(v):
                ax.text(i, j, f"{v:.1f}", ha="center", va="center",
                        fontsize=7.5, color="w" if v > 9 else "k")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("deg per hour of arc")
    ax.axhspan(-0.5, len(guide) - 0.5, xmin=0, xmax=0, color="none")

    # (c) cumulative distributions of both metrics
    ax = axes[1, 0]
    rate, shift = candidate_metrics("stratified", sky)
    for v, lab, col in ((rate, "rate [deg/h]", "steelblue"),
                        (shift, "shift over 1 h [deg]", "darkred")):
        v = v[np.isfinite(v)]
        if v.size:
            xs = np.sort(v)
            ax.plot(xs, np.arange(1, xs.size + 1) / xs.size, lw=2,
                    color=col, label=f"{lab}  (n={xs.size:,})")
    ax.axvspan(1, 5, color="0.85", zorder=0)
    ax.text(3, 0.04, "referee's range", ha="center", fontsize=8.5, color="0.3")
    ax.set_xlim(0, 30)
    ax.set_xlabel("tightest value achieved by an orbit")
    ax.set_ylabel("cumulative fraction of candidates")
    ax.set_title("(c) distribution over all close-flyby candidates",
                 fontsize=10.5)
    ax.legend(loc="lower right", fontsize=8.5)
    ax.grid(alpha=0.3)

    # (d) prior comparison
    ax = axes[1, 1]
    f_f = os.path.join(HERE, "bigsweep_flat.npz")
    txt = [f"published priors: {n_mc:,} orbits, {n_cand:,} close flybys"]
    fin = min_shift_corr[np.isfinite(min_shift_corr)]
    txt.append(f"   tightest in-corridor displacement  {fin.min():.2f} deg"
               if fin.size else "")
    if os.path.exists(f_f):
        e = np.load(f_f, allow_pickle=True)
        r2, s2 = candidate_metrics("flat", sky)
        for v, lab, col, ls in ((shift, "published prior", "steelblue", "-"),
                                (s2, "flat prior", "seagreen", "--")):
            v = v[np.isfinite(v)]
            if v.size:
                xs = np.sort(v)
                ax.plot(xs, np.arange(1, xs.size + 1) / xs.size, lw=2, ls=ls,
                        color=col, label=f"{lab}  (n={xs.size:,})")
        fin2 = e["min_shift_corr"][np.isfinite(e["min_shift_corr"])]
        txt.append(f"flat priors: {int(e['n_mc']):,} orbits, "
                   f"{int(e['n_cand']):,} close flybys")
        if fin2.size:
            txt.append(f"   tightest in-corridor displacement  "
                       f"{fin2.min():.2f} deg")
        ax.axvspan(1, 5, color="0.85", zorder=0)
        ax.set_xlim(0, 30)
        ax.set_xlabel("tightest 1 h displacement [deg]")
        ax.set_ylabel("cumulative fraction")
        ax.legend(loc="lower right", fontsize=8.5)
        ax.grid(alpha=0.3)
    else:
        ax.axis("off")
        txt.append("\nflat-prior run not present")
    ax.set_title("(d) sensitivity to the sampling prior\n(displacement metric, as in panel b)", fontsize=10.5)
    fig.text(0.52, 0.015, "   |   ".join(t for t in txt if t), fontsize=8.5,
             ha="center", color="0.25")

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    for ext in ("png", "pdf"):
        p = os.path.join(HERE, f"fig_criteria_sweep.{ext}")
        fig.savefig(p, dpi=190 if ext == "png" else None, bbox_inches="tight")
        print("wrote", p)


if __name__ == "__main__":
    main()
