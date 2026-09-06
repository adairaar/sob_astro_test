#!/usr/bin/env python3
"""
sob_itinerary_sweep.py — loosening the itinerary, as Reviewer 2 asks.

R2, round 4:

    "the 'itinerary' of the hypothetical Magi in the author's construction is
     rather rigid – with the requirement that the object remain fixed in the
     sky for an extended period (a ~2 hour trip with limited azimuth range
     plus another ~2 hours of 'stopping'). However, if the Magi completed
     making the ~2 hour trip from Jerusalem to Bethlehem, how much longer
     would the Star have needed to pause over the house? It does not seem
     likely that the Magi would have stopped a few meters before the entrance
     of the house for an hour or two just to make sure the Star had come to a
     halt. Here again, I believe some flexibility is warranted."

⚠ ONE POINT OF FACT TO PUT TO THE REFEREE FIRST.

The implementation does not require the object to remain at rest for any
period at all.  The stopping test takes the MINIMUM apparent motion over the
window following guidance and asks whether that single instant falls below
threshold; no duration is demanded.  The "~2 hours of stopping" R2 describes
is a feature of the paper's prose, not of the test.  On the duration of
stopping the model is therefore already at its loosest possible setting, and
the paper's wording should be corrected to say so.

What can meaningfully be loosened is:

    GUIDE_DUR   how long guidance must be sustained before the stop
                (published 1.0 h; a shorter journey demands less)
    STOP_GAP    how long after guidance the stop may occur
                (published 4.0 h; a longer window gives more opportunity)
    ALT_RATE_MIN whether the object must be seen to move at all during
                guidance (published 0.3 deg/h; dropping it is a loosening)
    drift       the azimuth-drift tolerance, which the companion sweep
                identified as the binding constraint

Every one of these is swept in the direction that makes a passing orbit
EASIER to find.  If none is found, the negative result is robust to the
whole itinerary question.
"""

import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sob_close_flyby_mc as M                                 # noqa: E402
import sob_criteria_sweep as S                                 # noqa: E402

N_MC = int(os.environ.get("SOB_N_MC", 24_000_000))
SEED = int(os.environ.get("SOB_SEED", 2026))
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     f"_flyby_cache_{N_MC}_{SEED}.pkl")

AZ_LO_GRID = [190.0, 150.0]                 # published, and loosest
DRIFT_GRID = [1.5, 3.0, 5.0, 10.0, 15.0]    # deg/h
GUIDE_GRID = [0.25, 0.5, 1.0, 2.0]          # h, published 1.0
GAP_GRID = [2.0, 4.0, 6.0, 8.0]             # h, published 4.0
STOP_GRID = [0.3, 1.0, 2.0, 3.0, 4.0, 5.0]  # deg/h and deg
SHIFT_WIN = 1.0                             # h
REQUIRE_MOVING = [True, False]              # ALT_RATE_MIN on / off


def load_candidates(sky):
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            c = pickle.load(f)
        print(f"  loaded {len(c)} cached close flybys from {os.path.basename(CACHE)}")
        return c
    c = S.sample_close_flybys(N_MC, sky, seed=SEED)
    with open(CACHE, "wb") as f:
        pickle.dump(c, f)
    print(f"  cached {len(c)} close flybys")
    return c


def best_metrics(alt, az, night, jd_utc, az_lo, drift, guide_dur, gap, moving_req):
    """Smallest stopping metrics reachable under one parameter set."""
    step_h = (jd_utc[1] - jd_utc[0]) * 24.0
    dt_h = np.gradient(jd_utc) * 24.0
    dalt = np.gradient(alt) / dt_h
    az_u = np.degrees(np.unwrap(np.radians(az)))
    daz = np.gradient(az_u) / dt_h
    omega = np.sqrt(dalt**2 + (daz * np.cos(np.radians(alt)))**2)

    k = max(1, int(round(SHIFT_WIN / step_h)))
    shift = np.full_like(omega, np.inf)
    if len(alt) > k:
        shift[:-k] = S.great_circle(alt[:-k], az_u[:-k], alt[k:], az_u[k:])

    guide = (night & (alt > M.ALT_MIN)
             & (az >= az_lo) & (az <= M.AZ_HI)
             & (np.abs(daz) < drift))
    if moving_req:
        guide &= np.abs(dalt) > M.ALT_RATE_MIN
    if (not np.any(guide)) or float(np.sum(guide)) * step_h < guide_dur:
        return None

    gi = int(np.where(guide)[0][-1])
    sl = slice(gi, min(gi + int(np.ceil(gap / step_h)) + 1, len(omega)))
    nm = night[sl]
    if not np.any(nm):
        return None
    om = omega[sl][nm]
    sh = shift[sl][nm]
    sh = sh[np.isfinite(sh)]
    return float(np.min(om)), (float(np.min(sh)) if sh.size else np.inf)


def main():
    print("=" * 78)
    print("ITINERARY SENSITIVITY SWEEP — Reviewer 2, round 4")
    print("=" * 78)
    print(f"  N_MC = {N_MC:,}")
    print(f"  corridors {AZ_LO_GRID} .. {M.AZ_HI:.0f} deg")
    print(f"  drift {DRIFT_GRID} deg/h")
    print(f"  guidance duration {GUIDE_GRID} h   (published {M.GUIDE_DUR})")
    print(f"  stop window {GAP_GRID} h           (published {M.STOP_GAP})")
    print(f"  'must be moving' during guidance: {REQUIRE_MOVING}")
    print()
    print("  NOTE: no stopping DURATION is required by the test, published or")
    print("        here; the stop is evaluated at its single best instant.")
    print()

    sky = M.precompute_sky(M.JD_PRIMARY)
    jd_utc = sky[0]
    cands = load_candidates(sky)
    if not cands:
        print("  no candidates")
        return

    tracks = []
    for c in cands:
        alt, az, night, _ = S.track(*c, sky)
        if alt is not None:
            tracks.append((alt, az, night))
    print(f"  usable tracks: {len(tracks)}\n")

    best_om, best_sh = np.inf, np.inf
    best_cfg = None
    n_guide_total = 0
    passes = 0
    rows = []

    for mv in REQUIRE_MOVING:
        for az_lo in AZ_LO_GRID:
            for drift in DRIFT_GRID:
                for gd in GUIDE_GRID:
                    for gap in GAP_GRID:
                        ng = 0
                        om_min, sh_min = np.inf, np.inf
                        for alt, az, night in tracks:
                            r = best_metrics(alt, az, night, jd_utc,
                                             az_lo, drift, gd, gap, mv)
                            if r is None:
                                continue
                            ng += 1
                            om_min = min(om_min, r[0])
                            sh_min = min(sh_min, r[1])
                            if r[0] <= max(STOP_GRID) or r[1] <= max(STOP_GRID):
                                passes += 1
                        n_guide_total += ng
                        if ng:
                            rows.append((mv, az_lo, drift, gd, gap, ng,
                                         om_min, sh_min))
                            if om_min < best_om:
                                best_om = om_min
                                best_cfg = (mv, az_lo, drift, gd, gap)
                            best_sh = min(best_sh, sh_min)

    print("=" * 78)
    print("PARAMETER SETS ADMITTING A GUIDANCE WINDOW")
    print("=" * 78)
    if not rows:
        print("  none, under any combination tested")
    else:
        print(f"  {'moving':>7} {'corr':>6} {'drift':>6} {'guide':>6} {'gap':>5} "
              f"{'n':>5} {'min rate':>9} {'min shift':>10}")
        print("  " + "-" * 68)
        for mv, a, d, gd, gap, ng, om, sh in sorted(rows, key=lambda r: r[6]):
            print(f"  {str(mv):>7} {a:>6.0f} {d:>6g} {gd:>6g} {gap:>5g} "
                  f"{ng:>5d} {om:>9.3f} {sh:>10.3f}")

    print("\n" + "=" * 78)
    print("RESULT")
    print("=" * 78)
    print(f"  parameter sets giving guidance : {len(rows)} of "
          f"{len(REQUIRE_MOVING)*len(AZ_LO_GRID)*len(DRIFT_GRID)*len(GUIDE_GRID)*len(GAP_GRID)}")
    print(f"  orbits passing any stopping test <= {max(STOP_GRID):g}: {passes}")
    print(f"\n  tightest stopping rate  reached anywhere : {best_om:.3f} deg/h")
    print(f"  tightest angular shift  reached anywhere : {best_sh:.3f} deg"
          f"  (over {SHIFT_WIN:g} h)")
    if best_cfg:
        mv, a, d, gd, gap = best_cfg
        print(f"    achieved with moving={mv}, corridor {a:.0f}-{M.AZ_HI:.0f}, "
              f"drift {d:g} deg/h, guidance {gd:g} h, window {gap:g} h")
    print(f"""
==============================================================================
READING THIS
==============================================================================
  Every parameter here has been moved in the direction that makes a passing
  orbit easier to obtain: shorter journeys, longer windows in which to stop,
  wider corridors, larger drift tolerances, and no requirement that the object
  be seen to move at all.

  If the tightest stopping value reached anywhere still exceeds the range the
  referee proposes (1-5 deg/h), then no amount of flexibility about the Magi's
  itinerary rescues the hypothesis, and that is the honest form of the reply.""")


if __name__ == "__main__":
    main()
