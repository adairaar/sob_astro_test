#!/usr/bin/env python3
"""
sob_bigsweep.py — the criteria sweep at survey scale, with a prior control.

The exploratory sweeps established the shape of the answer on 24 million
orbits and 148 close flybys.  That is too thin to put to a referee who is
being asked to accept a negative result, and it is open to two objections:

    (i)  small-sample.  148 candidates cannot characterise the tail of the
         stopping-metric distribution, and the qualification found in the
         itinerary sweep (4 orbits reaching 3.7-5.0 deg/h) rests on four.

    (ii) biased sampling.  The published priors are stratified: 60 % of q in
         [0.80, 1.20] au and 60 % of inclinations below 20 deg, on the
         argument that close flybys able to cancel diurnal drift must lie
         near the ecliptic and near 1 au.  The stratification concentrates
         effort; it does not exclude anything, since 15 % of draws are
         isotropic.  But a referee is entitled to see that demonstrated
         rather than asserted.

This script answers both.  It runs to ~1e9 orbits under the published priors,
and repeats a smaller run under FLAT priors -- isotropic inclination, q
log-uniform over the whole range -- so the two candidate populations can be
compared directly.  All derived quantities are precomputed once per track, so
the grid evaluation costs seconds rather than tens of minutes.

Outputs `bigsweep_<prior>.npz` for the plotting script.
"""

import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sob_close_flyby_mc as M                                 # noqa: E402
import sob_criteria_sweep as S                                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
N_MC = int(os.environ.get("SOB_N_MC", 1_000_000_000))
# The host has ~3 GB of RAM. A batch of n draws holds six element arrays plus
# a (3, n) position array plus Barker intermediates, so 4e6 peaks near 2 GB and
# is killed. 1e6 is safe with room to spare.
BATCH = int(os.environ.get("SOB_BATCH", 1_000_000))
SEED = int(os.environ.get("SOB_SEED", 2026))
PRIOR = os.environ.get("SOB_PRIOR", "stratified")     # or "flat"

AZ_LO_GRID = np.array([190.0, 180.0, 170.0, 160.0, 150.0])
DRIFT_GRID = np.array([1.5, 3.0, 5.0, 10.0, 15.0, 20.0, 30.0])
GUIDE_GRID = np.array([0.25, 0.5, 1.0, 1.5, 2.0])
GAP_GRID = np.array([2.0, 4.0, 8.0])
STOP_GRID = np.array([0.3, 1.0, 2.0, 3.0, 4.0, 5.0])
SHIFT_WIN = 1.0


def draw(n, rng, prior):
    """Orbital elements under the published or a flat prior."""
    if prior == "flat":
        q = np.exp(rng.uniform(np.log(0.30), np.log(1.40), n))
        e = rng.uniform(0.70, 5.00, n)
        i_r = np.arccos(rng.uniform(-1.0, 1.0, n))            # isotropic
    else:
        n0, n1 = int(0.60 * n), int(0.25 * n)
        n2 = n - n0 - n1
        q = np.concatenate([
            rng.uniform(0.80, 1.20, n0),
            np.where(rng.integers(0, 2, n1).astype(bool),
                     rng.uniform(0.50, 0.80, n1), rng.uniform(1.20, 1.40, n1)),
            rng.uniform(0.30, 0.50, n2)])
        e = np.concatenate([
            rng.uniform(0.70, 1.00, n // 3), rng.uniform(1.00, 1.01, n // 6),
            rng.uniform(1.01, 5.00, n - n // 3 - n // 6)])
        m0, m1 = int(0.60 * n), int(0.25 * n)
        m2 = n - m0 - m1
        i_r = np.concatenate([
            np.arccos(rng.uniform(np.cos(np.radians(20)), 1.0, m0)),
            np.arccos(rng.uniform(np.cos(np.radians(45)),
                                  np.cos(np.radians(20)), m1)),
            np.arccos(rng.uniform(-1.0, np.cos(np.radians(45)), m2))])
    return (q, e, i_r, rng.uniform(0, 2 * np.pi, n), rng.uniform(0, 2 * np.pi, n),
            M.JD_PRIMARY + rng.uniform(-90.0, 90.0, n))


def harvest(n_mc, sky, seed, prior):
    """Accumulate close flybys, batching to bound memory."""
    cache = os.path.join(HERE, f"_big_{prior}_{n_mc}_{seed}.pkl")
    if os.path.exists(cache):
        c = pickle.load(open(cache, "rb"))
        print(f"  loaded {len(c):,} cached flybys ({prior} prior)", flush=True)
        return c
    # Partial checkpoint, so an interrupted run resumes rather than restarts.
    part = cache + ".part"
    jd_utc, earth_pos = sky[0], sky[1]
    emid = earth_pos[:, len(jd_utc) // 2]
    out, k = [], 0
    if os.path.exists(part):
        k, out = pickle.load(open(part, "rb"))
        print(f"  resuming at batch {k} with {len(out):,} flybys", flush=True)
    done, t0 = k * BATCH, time.time()
    while done < n_mc:
        n = min(BATCH, n_mc - done)
        rng = np.random.default_rng(seed + 7919 * k)
        q, e, i_r, Om, om, T = draw(n, rng, prior)
        pos = M.comet_pos_parabolic_batch(M.JD_PRIMARY - T, q, i_r, Om, om)
        keep = np.linalg.norm(pos - emid[:, None], axis=0) < M.D_REJECT
        del pos
        for j in np.where(keep)[0]:
            _, _, _, dm = S.track(q[j], e[j], i_r[j], Om[j], om[j], T[j], sky)
            if dm <= M.D_FLYBY:
                out.append((q[j], e[j], i_r[j], Om[j], om[j], T[j]))
        del q, e, i_r, Om, om, T, keep
        done += n
        k += 1
        if k % 100 == 0 or done >= n_mc:
            el = max(time.time() - t0, 1e-9)
            rate = (done - (k - 1) * 0) / el
            pickle.dump((k, out), open(part, "wb"))
            print(f"    {done:>13,}/{n_mc:,}  flybys {len(out):>6,}  "
                  f"{done/el/1e6:.2f} M/s  eta "
                  f"{max(n_mc-done,0)/max(done/el,1)/60:.1f} min", flush=True)
    pickle.dump(out, open(cache, "wb"))
    try:                       # the workspace mount may forbid unlink
        os.remove(part)
    except OSError:
        pass
    return out


def derive(cands, sky):
    """Per-track quantities, computed once."""
    jd_utc = sky[0]
    step_h = (jd_utc[1] - jd_utc[0]) * 24.0
    dt_h = np.gradient(jd_utc) * 24.0
    k = max(1, int(round(SHIFT_WIN / step_h)))
    D = []
    for c in cands:
        alt, az, night, _ = S.track(*c, sky)
        if alt is None:
            continue
        dalt = np.gradient(alt) / dt_h
        az_u = np.degrees(np.unwrap(np.radians(az)))
        daz = np.gradient(az_u) / dt_h
        omega = np.sqrt(dalt**2 + (daz * np.cos(np.radians(alt)))**2)
        shift = np.full_like(omega, np.inf)
        if len(alt) > k:
            shift[:-k] = S.great_circle(alt[:-k], az_u[:-k], alt[k:], az_u[k:])
        D.append((alt, az, night, dalt, daz, omega, shift))
    return D, step_h


def main():
    print("=" * 78)
    print(f"BIG SWEEP — {PRIOR} prior, N = {N_MC:,}")
    print("=" * 78)
    sky = M.precompute_sky(M.JD_PRIMARY)
    cands = harvest(N_MC, sky, SEED, PRIOR)
    print(f"\n  close flybys: {len(cands):,}")
    if not cands:
        return
    D, step_h = derive(cands, sky)
    print(f"  usable tracks: {len(D):,}\n")

    shape = (len(AZ_LO_GRID), len(DRIFT_GRID), len(GUIDE_GRID), len(GAP_GRID))
    n_guide = np.zeros(shape, dtype=np.int64)
    min_rate = np.full(shape, np.inf)
    min_shift = np.full(shape, np.inf)
    # Corridor-constrained displacement: the stop must occur at a bearing
    # inside the guidance corridor, because Matthew has the star stand over
    # the house it has led them to. This is the quantity the claim rests on.
    min_shift_corr = np.full(shape, np.inf)
    n_pass = np.zeros(shape + (len(STOP_GRID), 2), dtype=np.int64)
    n_pass_corr = np.zeros(shape + (len(STOP_GRID),), dtype=np.int64)

    t0 = time.time()
    for alt, az, night, dalt, daz, omega, shift in D:
        high = night & (alt > M.ALT_MIN)
        moving = np.abs(dalt) > M.ALT_RATE_MIN
        adz = np.abs(daz)
        for ia, az_lo in enumerate(AZ_LO_GRID):
            corr = high & (az >= az_lo) & (az <= M.AZ_HI) & moving
            if not corr.any():
                continue
            for idr, drift in enumerate(DRIFT_GRID):
                g = corr & (adz < drift)
                ng = int(g.sum())
                if not ng:
                    continue
                gi = int(np.where(g)[0][-1])
                dur = ng * step_h
                for ig, gd in enumerate(GUIDE_GRID):
                    if dur < gd:
                        continue
                    for ip, gap in enumerate(GAP_GRID):
                        sl = slice(gi, min(gi + int(np.ceil(gap / step_h)) + 1,
                                           len(omega)))
                        nm = night[sl]
                        if not nm.any():
                            continue
                        om = float(np.min(omega[sl][nm]))
                        sh = shift[sl][nm]
                        sh = sh[np.isfinite(sh)]
                        sh = float(np.min(sh)) if sh.size else np.inf
                        # stop restricted to the guidance corridor
                        cm = nm & (az[sl] >= az_lo) & (az[sl] <= M.AZ_HI)
                        shc = shift[sl][cm] if cm.any() else np.array([])
                        shc = shc[np.isfinite(shc)]
                        shc = float(np.min(shc)) if shc.size else np.inf
                        idx = (ia, idr, ig, ip)
                        n_guide[idx] += 1
                        min_rate[idx] = min(min_rate[idx], om)
                        min_shift[idx] = min(min_shift[idx], sh)
                        min_shift_corr[idx] = min(min_shift_corr[idx], shc)
                        n_pass[idx + (slice(None), 0)] += (om <= STOP_GRID)
                        n_pass[idx + (slice(None), 1)] += (sh <= STOP_GRID)
                        n_pass_corr[idx + (slice(None),)] += (shc <= STOP_GRID)
    print(f"  grid evaluated in {time.time()-t0:.1f}s")

    out = os.path.join(HERE, f"bigsweep_{PRIOR}.npz")
    np.savez_compressed(
        out, n_guide=n_guide, min_rate=min_rate, min_shift=min_shift,
        min_shift_corr=min_shift_corr, n_pass_corr=n_pass_corr,
        n_pass=n_pass, az_lo=AZ_LO_GRID, drift=DRIFT_GRID, guide=GUIDE_GRID,
        gap=GAP_GRID, stop=STOP_GRID, shift_win=SHIFT_WIN,
        n_mc=N_MC, n_cand=len(D), prior=PRIOR)

    finite = min_rate[np.isfinite(min_rate)]
    fs = min_shift[np.isfinite(min_shift)]
    print("\n" + "=" * 78)
    print("HEADLINE")
    print("=" * 78)
    print(f"  orbits sampled            : {N_MC:,}")
    print(f"  close flybys              : {len(D):,}")
    print(f"  parameter sets with guidance : "
          f"{int((n_guide > 0).sum())} of {n_guide.size}")
    if finite.size:
        print(f"  tightest stopping rate    : {finite.min():.3f} deg/h")
        print(f"  tightest angular shift    : {fs.min():.3f} deg "
              f"(over {SHIFT_WIN:g} h)")
    tot = n_pass[..., 0].sum(axis=(0, 1, 2, 3))
    tos = n_pass[..., 1].sum(axis=(0, 1, 2, 3))
    toc = n_pass_corr.sum(axis=(0, 1, 2, 3))
    fc = min_shift_corr[np.isfinite(min_shift_corr)]
    if fc.size:
        print(f"  tightest shift IN CORRIDOR   : {fc.min():.3f} deg")
    print(f"  grid cells with any in-corridor pass at <=5 deg: "
          f"{int((n_pass_corr[..., -1] > 0).sum())} of {n_pass_corr[..., -1].size}")
    print("\n  orbit-configurations passing, summed over the grid:")
    print(f"    {'threshold':>16}" + "".join(f"{s:>8g}" for s in STOP_GRID))
    print(f"    {'rate':>16}" + "".join(f"{v:>8d}" for v in tot))
    print(f"    {'shift (anywhere)':>16}" + "".join(f"{v:>8d}" for v in tos))
    print(f"    {'shift (corridor)':>16}" + "".join(f"{v:>8d}" for v in toc))
    print(f"\n  saved -> {os.path.basename(out)}")


if __name__ == "__main__":
    main()
