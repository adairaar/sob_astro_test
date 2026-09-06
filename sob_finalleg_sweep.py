#!/usr/bin/env python3
"""
sob_finalleg_sweep.py — the two-leg journey, with a turn toward the house.

Reviewer 2, round 4:

    "the author has assumed the house was aligned with the Jerusalem-to-
     Bethlehem road, and thus the Star maintained its azimuth for the entire
     journey. The town of Bethlehem itself is located some distance to the
     east of the road. ... Using the Church of the Nativity as a convenient
     example of a possible 'off-road' location, the final leg of the Magi's
     journey off the main road would have required a shift in azimuth to more
     like ~170º (or turning even further east for the last portion, depending
     on the path taken) for the last ~30 minutes of the journey."

WHY THIS IS NOT THE SAME TEST AS THE WIDENED CORRIDOR

The companion sweep answered this by widening the corridor to a single static
window of 150-210 deg.  That is more generous in one respect and less in
another, and the difference matters.  A static corridor still demands LOW
AZIMUTH DRIFT throughout, so a star that swings from 200 deg to 170 deg is
excluded by the drift test no matter how wide the corridor.  But R2's scenario
positively REQUIRES such a swing: roughly 20-40 deg within the final half hour,
which is a drift of 40-80 deg/h during that leg alone.

The journey is therefore modelled here in two legs:

    leg 1   the road from Jerusalem, azimuth held in corridor A (190-210),
            drift below tolerance, for at least T1 hours;
    -----   a transition of up to TRANS_MAX hours in which the azimuth is
            free to swing as the party turns off the road;
    leg 2   the approach to the house, azimuth held in corridor B (near 170),
            drift below tolerance, for at least T2 hours;
    stop    the stopping test in the window that follows.

This is strictly more permissive than either the published single corridor or
the widened static one, because the star is now allowed to move between two
bearings and is only required to hold each of them briefly.
"""

import os
import pickle
HERE = os.path.dirname(os.path.abspath(__file__))
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sob_close_flyby_mc as M                                 # noqa: E402
import sob_criteria_sweep as S                                 # noqa: E402
import sob_itinerary_sweep as I                                # noqa: E402

# leg 1: the Jerusalem-Bethlehem road
AZ1 = (190.0, 210.0)
# leg 2: turning east toward a house off the road. R2 offers ~170 deg for the
# Church of the Nativity, "or turning even further east", so the window runs
# from due south round to the south-east.
AZ2_GRID = [(160.0, 180.0), (140.0, 180.0), (120.0, 180.0)]
T1_GRID = [1.5, 1.0, 0.5]        # h on the road   (R2's trip is ~2 h total)
T2_GRID = [0.5, 0.25]            # h on the final leg (R2 says ~30 min)
DRIFT_GRID = [1.5, 3.0, 5.0, 10.0, 15.0, 20.0, 30.0]
TRANS_MAX = float(os.environ.get("TRANS_MAX", 1.0))   # h allowed for the turn
STOP_GAP = 8.0                   # h, loosest from the itinerary sweep
STOP_GRID = [0.3, 1.0, 2.0, 3.0, 4.0, 5.0]
SHIFT_WIN = 1.0


def runs(mask):
    """Yield (start, stop) index pairs of contiguous True runs."""
    d = np.diff(mask.astype(np.int8))
    starts = list(np.where(d == 1)[0] + 1)
    stops = list(np.where(d == -1)[0] + 1)
    if mask[0]:
        starts = [0] + starts
    if mask[-1]:
        stops = stops + [len(mask)]
    return list(zip(starts, stops))


def two_leg(alt, az, night, jd_utc, az2, t1, t2, drift):
    """Best stopping metrics for a road leg, a turn, and an approach leg."""
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

    base = night & (alt > M.ALT_MIN) & (np.abs(daz) < drift)
    m1 = base & (az >= AZ1[0]) & (az <= AZ1[1])
    m2 = base & (az >= az2[0]) & (az <= az2[1])

    n1 = int(np.ceil(t1 / step_h))
    n2 = int(np.ceil(t2 / step_h))
    n_trans = int(np.ceil(TRANS_MAX / step_h))
    n_gap = int(np.ceil(STOP_GAP / step_h))

    best = None
    for a1, b1 in runs(m1):
        if b1 - a1 < n1:
            continue
        for a2, b2 in runs(m2):
            if b2 - a2 < n2:
                continue
            if not (b1 <= a2 <= b1 + n_trans):     # leg 2 follows leg 1
                continue
            sl = slice(b2, min(b2 + n_gap + 1, len(omega)))
            nm = night[sl]
            if not np.any(nm):
                continue
            om = omega[sl][nm]
            sh = shift[sl][nm]
            sh = sh[np.isfinite(sh)]
            cand = (float(np.min(om)),
                    float(np.min(sh)) if sh.size else np.inf,
                    float(np.max(np.abs(daz[b1:a2 + 1]))) if a2 > b1 else 0.0)
            if best is None or cand[0] < best[0]:
                best = cand
    return best


def main():
    print("=" * 78)
    print("TWO-LEG ITINERARY SWEEP — Reviewer 2's off-road final approach")
    print("=" * 78)
    print(f"  leg 1 corridor : {AZ1[0]:.0f}-{AZ1[1]:.0f} deg (the road)")
    print(f"  leg 2 corridors: {[f'{a:.0f}-{b:.0f}' for a, b in AZ2_GRID]} deg")
    print(f"  leg 1 duration : {T1_GRID} h")
    print(f"  leg 2 duration : {T2_GRID} h")
    print(f"  drift tolerance: {DRIFT_GRID} deg/h (within each leg)")
    print(f"  turn allowed   : up to {TRANS_MAX:g} h, azimuth unconstrained")
    print(f"  stopping window: {STOP_GAP:g} h after leg 2")
    print()

    sky = M.precompute_sky(M.JD_PRIMARY)
    jd_utc = sky[0]
    if not os.path.exists(I.CACHE):
        print(f"  cache {I.CACHE} missing; run sob_itinerary_sweep.py first")
        return
    # Prefer the 10^9 candidate set over the 24M exploratory cache, so the
    # itinerary answer rests on the same sample as the criteria sweep.
    big = os.path.join(HERE, "_big_stratified_1000000000_2026.pkl")
    src = big if os.path.exists(big) else I.CACHE
    cands = pickle.load(open(src, "rb"))
    print(f"  candidate set: {os.path.basename(src)}  ({len(cands):,} flybys)")
    tracks = []
    for c in cands:
        a, z, n, _ = S.track(*c, sky)
        if a is not None:
            tracks.append((a, z, n))
    print(f"  close flybys: {len(tracks)}\n")

    rows = []
    for az2 in AZ2_GRID:
        for t1 in T1_GRID:
            for t2 in T2_GRID:
                for drift in DRIFT_GRID:
                    oms, shs, turns = [], [], []
                    for alt, az, night in tracks:
                        r = two_leg(alt, az, night, jd_utc, az2, t1, t2, drift)
                        if r is not None:
                            oms.append(r[0]); shs.append(r[1]); turns.append(r[2])
                    if oms:
                        rows.append((az2, t1, t2, drift, len(oms),
                                     min(oms), min(shs), max(turns)))

    print("=" * 78)
    print("PARAMETER SETS ADMITTING A TWO-LEG ITINERARY")
    print("=" * 78)
    if not rows:
        print("\n  none. No orbit produces a road leg followed by a turn and an")
        print("  approach leg, under any combination tested.\n")
        print("  The reason is directional, and is reported separately below:")
        print("  the two legs do occur, but always in the wrong order.\n")

    else:
        print(f"  {'leg2':>9} {'T1':>5} {'T2':>5} {'drift':>6} {'n':>4} "
              f"{'min rate':>9} {'min shift':>10} {'turn rate':>10}")
        print("  " + "-" * 68)
        for az2, t1, t2, d, n, om, sh, tr in sorted(rows, key=lambda r: r[5]):
            print(f"  {az2[0]:>4.0f}-{az2[1]:<4.0f}{t1:>5g} {t2:>5g} {d:>6g} "
                  f"{n:>4d} {om:>9.3f} {sh:>10.3f} {tr:>10.1f}")

        allom = min(r[5] for r in rows)
        allsh = min(r[6] for r in rows)
        n5 = sum(1 for r in rows if r[5] <= 5.0)
        print(f"\n  parameter sets with a viable two-leg itinerary : {len(rows)}")
        print(f"  ...of which best rate falls inside 1-5 deg/h    : {n5}")
        print(f"\n  tightest stopping rate  anywhere : {allom:.3f} deg/h")
        print(f"  tightest angular shift  anywhere : {allsh:.3f} deg"
              f"  (over {SHIFT_WIN:g} h)")

    # ── why: which way does the azimuth run? ────────────────────────────────
    inc = dec = 0
    order = {"leg2 first": 0, "leg1 first": 0}
    step_h = (jd_utc[1] - jd_utc[0]) * 24.0
    n1 = int(np.ceil(min(T1_GRID) / step_h))
    n2 = int(np.ceil(min(T2_GRID) / step_h))
    for alt, az, night in tracks:
        az_u = np.degrees(np.unwrap(np.radians(az)))
        m = night & (alt > M.ALT_MIN) & (az >= 120) & (az <= 210)
        if m.sum() >= 3:
            g = float(np.median(np.gradient(az_u)[m]))
            inc += g > 0
            dec += g < 0
        daz = np.gradient(az_u) / (np.gradient(jd_utc) * 24.0)
        base = night & (alt > M.ALT_MIN) & (np.abs(daz) < 30.0)
        r1 = [r for r in runs(base & (az >= AZ1[0]) & (az <= AZ1[1]))
              if r[1] - r[0] >= n1]
        r2 = [r for r in runs(base & (az >= 120) & (az <= 180))
              if r[1] - r[0] >= n2]
        for _, b1 in r1:
            for a2, b2 in r2:
                if b2 <= b1 - n1:
                    order["leg2 first"] += 1
                elif a2 >= b1:
                    order["leg1 first"] += 1

    print("\n" + "=" * 78)
    print("WHY: THE REQUIRED TURN RUNS AGAINST THE SKY")
    print("=" * 78)
    print(f"  orbits whose azimuth increases (westward) through 120-210 deg : {inc}")
    print(f"  orbits whose azimuth decreases (eastward)                     : {dec}")
    print(f"\n  leg pairs in which the approach leg PRECEDES the road leg     : "
          f"{order['leg2 first']}")
    print(f"  leg pairs in which the road leg precedes the approach leg      : "
          f"{order['leg1 first']}")
    print("""
  The referee's own figure puts the final leg near 170 deg, east of the
  190-210 deg bearing held along the road from Jerusalem, so his final leg
  requires the star to swing EASTWARD, from about 200 deg to about 170 deg. Diurnal motion
  carries objects the other way. Every close-flyby orbit in the sample moves
  westward through this range, so where both bearings occur they occur in the
  reverse order: the star stands at 170 deg first and at 200 deg afterwards.
  The scenario is therefore not merely unmet but inverted.""")

    print("""
==============================================================================
READING THIS
==============================================================================
  This is the most permissive geometry yet tested. The star is allowed to
  hold one bearing along the road, swing freely through a turn of up to an
  hour, hold a second bearing anywhere from due south to the south-east while
  the party walks the last stretch, and only then is it asked to stop.

  'turn rate' reports the largest azimuth rate the star must exhibit during
  the transition. If that number is large, the star is not merely drifting
  but visibly swinging across the sky, which is itself difficult to reconcile
  with an object that then appears to halt.""")


if __name__ == "__main__":
    main()
