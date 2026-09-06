> **Record updated 5 September 2026.** Three things in this file were true when
> written and are now out of date in their description, though not in their
> numbers:
>
> 1. The corridor-constrained displacement is no longer a side calculation. It
>    is computed inside `sob_bigsweep.py` (`min_shift_corr`, `n_pass_corr`) and
>    regenerates with the sweep. The headline — tightest 6.938°, 0 of 525 cells
>    — comes straight out of that script.
> 2. The final-leg (two-leg itinerary) sweep now runs on the 10⁹ candidate set,
>    6,116 close flybys, not the 148 of the exploratory run. Result unchanged:
>    no parameter set admits the journey, and 2,024 orbits drift westward
>    through the relevant azimuths against 0 eastward.
> 3. `fig_criteria_sweep` panel (b) now plots the corridor-constrained
>    displacement rather than the stopping rate, and panel (d) compares priors
>    on the same metric. The rate remains a diagnostic only.
>
> The open question at the foot of this file about modelling the final leg as
> an explicit azimuth change has been closed: it was modelled explicitly, and
> that is where the directional inversion was found.

# ★★★ THE CLAIM TO MAKE

**Across 10⁹ sampled orbits and 6,116 close flybys, no orbit satisfies the displacement stopping criterion at any of the 525 parameter combinations tested — 0 of 525 cells.** The tightest displacement achieved anywhere is **6.938°** per hour, against the 5° at the loosest end of the referee's proposed range: a margin of **1.39×**.

The claim rests on exactly two commitments, and everything else can be loosened to the extreme:

1. **Stopping is measured as a finite angular displacement, not an instantaneous rate.** This is the referee's own proposal, on his own perceptual argument.
2. **The stop must occur at the bearing along which the star guided.** Matthew has the star stand over the house it has led them to; guidance and stopping were never linked in azimuth in the original test, and that was an oversight.

With those two in place, the corridor may be widened from 20° to 60°, the azimuth-drift tolerance raised from 1.5 to 30 °/h, the guidance leg cut from an hour to fifteen minutes, and the stopping window stretched to eight hours — and nothing appears, under either stratified or flat sampling priors.

**No altitude cap is required.** An earlier draft of this argument used one; it is unnecessary and should be dropped, since the displacement criterion is a necessary condition and its failure alone forecloses any conjunction containing it.

## Two things to state plainly rather than bury

- **The corridor link is a tightening.** Without it, one orbit in 6,116 reaches 3.88°. Give the before-and-after numbers and justify the change from the text, not from the outcome.
- **The margin is 1.39×, not an order of magnitude.** A referee who argued for a 7° threshold would overturn it. The answer is not to hide the margin but to observe that 7° per hour is some fourteen lunar diameters, which no observer would describe as having stopped.

## The rate metric is a diagnostic, not part of the claim

Reported separately, it shows *why* a rate criterion misleads. It admits 172 orbits, and inspection shows two artefact classes: altitude turning points, where the vertical motion vanishes while the object slides on (the single orbit examined below stops 24° off-bearing, 3.8 h late, at *v*∞ = 52 km s⁻¹), and near-zenith geometry at a median altitude of 79.5°, where cos(alt) suppresses a 14.5 °/h azimuth sweep into a small angular speed while the object indicates no direction at all.

---

# ⚠⚠ SUPERSEDED IN PART — read this first

Everything below the first horizontal rule was computed on **24 million orbits and 148 close flybys**. That sample was too small, and two of its headline numbers were artefacts of it. The campaign at **10⁹ orbits and 6,116 close flybys** (`sob_bigsweep.py`, figure `fig_criteria_sweep.pdf`) supersedes them.

| quantity | at 148 candidates | at 6,116 candidates |
|---|---|---|
| tightest stopping **rate** anywhere | 3.661 °/h | **0.680 °/h** |
| tightest angular **shift** anywhere | 9.422° | **3.880°** |
| orbits passing displacement ≤ 5° | **0, claimed "nowhere"** | **1** |

**The claim that no orbit satisfies the displacement metric anywhere was a small-sample result and is withdrawn.** At 10⁹ one orbit reaches 3.88°, inside the referee's 1–5° range.

## ★★ The single passing orbit is a false positive

Candidate **#1178** of 6,116 — the only orbit satisfying a displacement criterion of 5° per hour. Examined directly (`sob_falsepositive_figure.py`, figure `fig_false_positive.pdf`), it fails on three independent grounds.

**1. It is not a comet.** *q* = 1.018 au, **e = 4.058**, *i* = 17.3°, perihelion 2.84 days before the event, closest approach 620,499 km (1.6 lunar distances). A hyperbolic excess velocity of **52 km s⁻¹** — against 26 km s⁻¹ for 1I/ʻOumuamua and 32 for 2I/Borisov. This is not a Solar System body but an interstellar interloper faster than either object ever observed.

**2. It is not stopping.** The "stop" is an **altitude turning point**. Around the minimum the altitude rate passes through zero (−0.80 → −0.00 → +0.23 ° h⁻¹) while the **azimuth rate never falls below ~3.97 ° h⁻¹** and shows no inflexion at all. The total apparent motion bottoms out at 3.85 ° h⁻¹, and that residue is *entirely* the azimuth term: 3.97 × cos 14.3° = 3.85. The 3.88° "displacement" is simply four degrees per hour of uninterrupted westward sliding, which over the two hours the Magi supposedly watched would carry the object some eight degrees across the sky.

**3. It is not over the house.** Guidance ends at azimuth 209.5°. The "stop" occurs at **azimuth 233.5°** and altitude 14.3° — **24° west of the corridor edge, and 3.8 hours later.** The star has crossed most of the south-western sky before its altitude happens to bottom out.

### The methodological correction this exposes

The stopping test as written permits the stop to occur **anywhere**, at any bearing. But Matthew has the star stand over the house it has led them to, so the stop must occur at the bearing it guided along. Adding that requirement:

| configuration | stop anywhere | stop **in corridor** |
|---|---|---|
| corr 170–210, drift 15, guide 0.25 h | shift ≤5: **1**, min 3.88° | shift ≤5: **0**, min **7.81°** |
| corr 150–210, drift 30, guide 0.25 h | shift ≤5: **1**, min 3.88° | shift ≤5: **0**, min **8.08°** |

**With the stop required to occur in the corridor, no orbit in 10⁹ satisfies the displacement criterion at any setting tested.**

> ⚠ This is a *tightening*, and must be declared as such rather than slipped in. It is a correction of an oversight in the original test — the guidance and stopping conditions were never linked in azimuth — and it is justified by the text rather than by the result it produces. State it plainly, give the before-and-after numbers above, and let the referee judge. Concealing it would be far worse than the qualification it removes.

## What does survive, and survives robustly

1. **The published configuration yields nothing.** At corridor 190–210°, drift 1.5 °/h, guidance 1.0 h, window 4.0 h: **n_guide = 0** at 10⁹ orbits. No orbit obtains a guidance window at all, so no orbit reaches the stopping test.
2. **Zero guidance windows anywhere at drift ≤ 5 °/h**, in any of the five corridors down to 150–210°. The entire left third of panel (a) is zero.
3. **Therefore the sweep the referee actually asked for — loosening the stopping criterion over 1–5 °/h — changes nothing**, because the constraint that binds is the guidance drift tolerance, not the stopping threshold.
4. **The result is not an artefact of the stratified priors.** A control run of 10⁹ orbits under **flat priors** — isotropic inclination, *q* log-uniform over 0.30–1.40 au — gives 450 close flybys, **n_guide = 0 in the published cell**, and **zero guidance windows at drift ≤ 5 °/h**. Its tightest rate is 2.55 °/h against 0.68 °/h for the published priors, i.e. slightly *worse* for the hypothesis. The stratification concentrates effort on the productive region; it does not manufacture the conclusion.

## Where the passes actually live

Passes appear only once the azimuth-drift tolerance is raised to **≥ 10 °/h** — a 6.7-fold loosening of a criterion the referee did *not* ask to loosen. Even there the two metrics behave very differently. In the single most permissive cell (corridor 170–210°, drift 15 °/h, guidance 0.25 h, window ≥ 4 h, n = 1,270):

| metric | orbits passing ≤ 5 |
|---|---|
| rate | **173** |
| displacement over 1 h | **1** |

Across the whole grid the displacement passes are essentially one orbit, at 3.88°, and it requires a stopping window of at least 4 hours; at a 2-hour window the same orbit gives 5.81° and fails.

> **Defensible statement.** No orbit in 10⁹ obtains a guidance window under the published criteria, or under any corridor, at an azimuth-drift tolerance of 5 °/h or below; this holds under both stratified and flat sampling. Raising the drift tolerance to 10–30 °/h — permitting 20–60° of azimuth wander across a two-hour journey — admits orbits, of which up to 173 satisfy a rate criterion of 5 °/h but exactly one satisfies a displacement criterion of 5° per hour.

The displacement metric remains far more discriminating, by a factor of ~170 in the same cell, and the reason is the one already given: an apparent path can be instantaneously stationary at a turning point while still carrying the object degrees across the sky within the hour.

---

# Criteria sensitivity sweep — result for Reviewer 2, round 4

Run 4 September 2026 with `sob_criteria_sweep.py`. **24,000,000 orbits sampled**, on the same stratified priors as the published survey; **148 close flybys** found at *d*_min < 0.020 au, consistent with the 66 reported from 10⁷.

---

## What the referee asked

> "Rather than choosing a rigid screening criterion and seeing if orbits match, find out what are the tightest criteria that still allow such an orbit to exist or not… Maybe try stopping criteria of 1º/hr, 2º/hr, 3º/hr, 4º/hr, and 5º/hr… That way, it will become clear how dependent the orbit matching is on the detailed criterion chosen."

Two of his other points are folded into the same run: the stopping test is now expressed as a **finite angular shift** rather than a rate, following his argument that the eye's motion-detection thresholds lie far above any rate at issue; and the azimuth corridor is **swept downward to 150°** to accommodate his observation that a house such as the Church of the Nativity lies off the road, on the eastern side of the town, and would require ~170° or further east on the final leg.

---

## ⚠ The first finding: the stopping criterion is not the binding constraint

Sweeping the stopping threshold alone is uninformative, because **no orbit reaches the stopping test at all**. The constraint that binds is the *guidance* azimuth-drift tolerance.

**Orbits admitting a guidance window** (night, altitude > 10°, in corridor, moving in altitude, azimuth drift below tolerance, sustained ≥ 1 h):

| corridor | 1.5°/h | 3°/h | 5°/h | 10°/h | 15°/h |
|---|---|---|---|---|---|
| 190–210 | 0 | 0 | 0 | 1 | 9 |
| 180–210 | 0 | 0 | 0 | 1 | 10 |
| 170–210 | 0 | 0 | 0 | 1 | 13 |
| 160–210 | 0 | 0 | 0 | 1 | 16 |
| 150–210 | 0 | 0 | 0 | 1 | 16 |

At the published tolerance of 1.5°/h — and at 3 and 5°/h — **no orbit in 24 million produces a guidance window in any corridor.** Widening the corridor by forty degrees does not change this.

Guidance only becomes possible at a drift tolerance of 10–15°/h. Over a two-hour journey that permits 20–30° of azimuth wander, at which point the object can hardly be said to have "gone before" anyone; but the sweep concedes it anyway.

## The second finding: even then, nothing stops

Among the orbits that *do* obtain a guidance window at 10–15°/h, in every corridor from 190–210 down to 150–210:

| stopping metric | tightest value any orbit achieves | R2's proposed range |
|---|---|---|
| angular rate | **8.317 °/h** | 1–5 °/h |
| angular shift over 1 h | **9.576 °** | 1–5 ° |

**Every cell of the pass table is zero**, for all six stopping thresholds (0.3, 1, 2, 3, 4, 5) crossed with all five corridors and both drift tolerances that admit guidance.

---

## What this licenses the paper to say

> The conclusion is insensitive to the stopping criterion across the entire range the referee proposes, and beyond it. No orbit satisfies the stopping test at any threshold below **8.3 °/h** on a rate metric or **9.6°** per hour on a displacement metric — roughly twice the loosest value suggested. The result is likewise insensitive to the azimuth corridor: relaxing the lower edge from 190° to 150°, which more than covers a final leg toward a house east of the road's bearing, changes no cell of the table.

This is a stronger statement than the original paper made, and it is the form the referee asked for: the reader can see the dependence rather than being asked to accept a threshold.

## Points conceded explicitly

1. **The rate-based stopping criterion is replaced by a finite angular shift**, per the referee's own reasoning about detection thresholds. Both are reported; the conclusion does not depend on which is used.
2. **The corridor assumption is relaxed**, per his point about the destination lying off the road, east of its bearing.
3. **The binding constraint is identified honestly** as the guidance drift tolerance rather than the stopping criterion, and that is swept too.

---

# Itinerary sweep — and a genuine qualification

Run 4 September 2026 with `sob_itinerary_sweep.py` on the same 148 close flybys. Every parameter moved in the direction that makes a passing orbit *easier*: shorter journeys, longer windows in which to stop, wider corridors, larger drift tolerances, and the requirement that the object be seen to move at all switched off. 320 parameter combinations; 112 admit a guidance window.

## ⚠ First, a point of fact to put to the referee

R2 writes that the construction requires "another ~2 hours of 'stopping'". **The test requires no stopping duration at all.** It takes the *minimum* apparent motion over the window following guidance and asks whether that single instant falls below threshold. On the duration of stopping the model is already at its loosest possible setting, and the paper's prose — not the test — is what gave the contrary impression. That wording should be corrected.

## The qualification

| configuration | orbits with guidance | rate ≤ 5°/h | rate ≤ 3°/h | shift ≤ 5° | min rate | min shift |
|---|---|---|---|---|---|---|
| **published**: 190–210, drift 1.5, guide 1.0 h | **0** | 0 | 0 | 0 | — | — |
| 190–210, drift 15, guide 0.25 h | 10 | 0 | 0 | 0 | 9.374 | 9.576 |
| 150–210, drift 10, guide 0.25 h | 3 | 0 | 0 | 0 | 9.374 | 9.576 |
| 150–210, drift 15, guide 1.0 h | 16 | 0 | 0 | 0 | 8.317 | 9.576 |
| 150–210, drift 15, guide 0.5 h | 25 | 1 | 0 | 0 | 4.964 | 9.528 |
| **150–210, drift 15, guide 0.25 h** | 32 | **4** | 0 | 0 | **3.661** | 9.422 |

**Under a rate criterion the hypothesis is not excluded at the loosest end of the referee's range.** With the corridor widened from 20° to 60°, the drift tolerance raised tenfold, and the guidance requirement cut to fifteen minutes, four orbits of 148 reach apparent rates between 3.7 and 5.0 °/h, which a 4 or 5 °/h threshold would admit. This must be stated. The paper can no longer claim that no orbit passes under any criterion.

## But the referee's own metric is not satisfied anywhere

Under the finite angular displacement he advocates, **no orbit in any of the 112 viable configurations achieves better than 9.42° over an hour**, against his proposed range of 1–5°. Zero of 320 parameter sets.

> **The discrepancy between the two metrics is itself the answer.** An object may be instantaneously slow — at a turning point in its apparent path, where the altitude and azimuth rates momentarily cancel — while still traversing nine degrees of sky in the following hour. A rate criterion evaluated at its best instant cannot distinguish that from stopping; a displacement criterion can, which is precisely why the referee proposed one. Adopting his metric resolves the ambiguity, and it resolves it against the hypothesis.

This is a better outcome than a flat negative would have been. The paper now reports where the result is fragile, why, and which criterion settles it — and the criterion that settles it is the referee's own.

## Suggested form of words

> Under a rate criterion the result is not robust at the extreme of the plausible range: with the guidance corridor widened to 60°, the drift tolerance raised to 15° h⁻¹ and the guidance requirement reduced to fifteen minutes, four of 148 close-flyby orbits attain apparent rates of 3.7–5.0° h⁻¹. Under the displacement criterion, however, none attains better than 9.4° per hour against a threshold of 1–5°, and this holds across all 112 parameter combinations admitting guidance. The difference is diagnostic rather than incidental: an apparent path may be instantaneously stationary at a turning point while still carrying the object nine degrees across the sky within the hour, and it is this that a displacement metric excludes and a rate metric cannot.

---

# Final-leg sweep — the turn runs against the sky

Run 4 September 2026 with `sob_finalleg_sweep.py`, on the same 148 close flybys.

R2's scenario is modelled explicitly, rather than approximated by a widened corridor: a **road leg** holding 190–210°, a **transition of up to an hour** in which the azimuth is entirely unconstrained, an **approach leg** holding a second bearing, and only then the stopping test. Grids: leg-2 corridor 160–180, 140–180 and 120–180; leg 1 of 1.5, 1.0 or 0.5 h; leg 2 of 0.5 or 0.25 h; drift tolerance 1.5 to **30 °/h**. 168 combinations.

**This matters because the widened static corridor did not actually test his scenario.** A static corridor still demands low azimuth drift throughout, so a star swinging from 200° to 170° is excluded however wide the window. R2's final leg positively *requires* that swing — some 20–40° within half an hour, a drift of 40–80 °/h during the turn.

## Result: no orbit, under any of the 168 combinations

| | count |
|---|---|
| parameter sets admitting a two-leg itinerary | **0 of 168** |

## The reason is directional, and it is worth putting to him

| | |
|---|---|
| orbits whose azimuth **increases** (westward) through 120–210° | **48** |
| orbits whose azimuth **decreases** (eastward) | **0** |
| leg pairs where the approach leg **precedes** the road leg | **19** |
| leg pairs where the road leg precedes the approach leg | **0** |

At a drift tolerance of 30 °/h, 23 orbits produce a valid road leg and 36 a valid approach leg, and 18 produce both. In **every** such case the approach bearing occurs *first* and the road bearing afterwards.

> The referee's own figure puts the final leg near 170°, east of the 190–210° bearing held along the road from Jerusalem, so the final leg requires the star to swing **eastward**, from roughly 200° to roughly 170°. Diurnal motion carries objects the other way, and every close-flyby orbit in the sample moves westward through this range. Where both bearings occur they therefore occur in the reverse order: the star stands at 170° first and at 200° afterwards. The referee's scenario is not merely unmet but inverted.

This is a physical argument rather than a tuned threshold, which is the kind of answer the exchange has been converging on. It also disposes of his "many more possible permutations" concern: the permutations exist, but they run backwards.

**Caveat to state honestly.** A close flyby can in principle produce parallactic motion large enough to reverse the apparent azimuth drift. None of the 148 does so within this azimuth range, but the claim is empirical over the sample, not a theorem, and should be worded that way.

## Still to do before resubmission

- [ ] Re-run at higher *N* to firm the tail. 148 flybys is enough to establish that the minimum sits near 8–10°/h, but a larger sample would tighten it.
- [ ] Address the referee's point that the Magi need not have waited one to two hours at the door: sweep `GUIDE_DUR` and `STOP_GAP` as well. This is the one loosening not yet tested.
- [ ] Decide whether to model the final leg as an explicit azimuth *change* (corridor for the first ~90 min, then a shift to ~170°) rather than a widened static corridor. The present treatment is more generous than his scenario, so it is safe, but the explicit version answers him more directly.
- [ ] ⚠ Tell the editor that the Matney criticisms R2 suggests expanding are already a standalone paper under review at *The Observatory*, so they cannot be folded into this one.
