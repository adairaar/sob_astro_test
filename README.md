# Star of Bethlehem: Orbital Mechanics and Corpus Linguistics

Code repository for:

**"The Star of Bethlehem belongs to the guiding-star legendary narrative tradition: evidence from orbital mechanics and Greek corpus linguistics"**
Aaron Adair (2026), submitted to *Nature*

The paper demonstrates that no Keplerian orbit of any type can simultaneously satisfy Matthew 2:9's guidance and stopping requirements, and that Matthew's motion vocabulary belongs to the legendary narrative register rather than the astronomical one.

---

## Requirements

Python 3.9 or later. Install dependencies with:

```bash
pip install numpy astropy scipy matplotlib
```

All scripts use only these four third-party libraries. The JPL planetary ephemeris (DE432) is fetched automatically by `astropy` on first use and cached locally; an internet connection is required for the initial download only.

---

## Repository contents

### Submission documents

| File | Description |
|------|-------------|
| `sob_paper.tex` | Main article (Nature format) |
| `sob_si.tex` | Supplementary Information |
| `sob_refs.bib` | BibTeX reference database |
| `cover_letter_draft.docx` | Cover letter |
| `sn-jnl.cls`, `sn-nature.bst`, `breakurl.sty` | Nature LaTeX class and style files |

### Figures (pre-generated)

| File | Used in |
|------|---------|
| `figure_paper1_sky_geometry.png` | Paper Fig. 1 — observational constraints from Matt 2:9 |
| `figure_paper2_az_motion.png` | Paper Fig. 2 — direction + guidance motion (two-criterion) |
| `figure_paper2_impossibility_scatter.png` | Paper Fig. 3 — guidance vs. stopping motion scatter |
| `figure6_impossibility_proof.png` | Extended Data Fig. 1 — orbit inversion proof |
| `figure_flyby_bestorbit.png` / `.pdf` | SI — best close-flyby orbit detail |
| `figure9_corpus_linguistics.png` | SI — corpus linguistics feature heatmap |
| `figure_mc_diagnostic.png` | SI Fig. S1 — Monte Carlo sweep diagnostic |

---

## Scripts

### Core orbital survey pipeline

These scripts establish the central impossibility result and should be run in the order shown if you want to reproduce all numerical results from scratch.

---

#### 1. `sob_proof_fine.py` — Exhaustive heliocentric orbital survey

The primary heliocentric survey. Tests 529,079,040 orbital configurations spanning elliptic, parabolic, and hyperbolic orbits at 5° resolution in Ω and ω.

**What it does:** For each configuration, evaluates whether the object simultaneously (a) lies within the road corridor [190°,210°] azimuth throughout the 2-hour guidance window, (b) moves faster than 2°/h during guidance, and (c) moves slower than 2°/h during stopping.

**Outputs:** `sob_proof_fine_progress.log`, `sob_proof_fine_results.json`

**Runtime:** Several hours on a modern multicore machine (uses `multiprocessing`). Progress is logged to `sob_proof_fine_progress.log`.

```bash
python3 sob_proof_fine.py
```

---

#### 2. `sob_mc_refined.py` — Monte Carlo refinement

Performs two dense Monte Carlo sweeps around the most promising grid orbits identified by `sob_proof_fine.py`, confirming the null result is not a grid-spacing artefact.

- **Sweep 1:** 1.5 million draws near the best-azimuth grid orbit (e=1.1, q=0.03 AU, i=10°)
- **Sweep 2:** 750,000 draws near Matney's (2025) proposed 5 BCE comet

**Outputs:** `sob_mc_refined.log`, `sob_mc_results.json`

**Runtime:** 30–90 minutes depending on hardware.

```bash
python3 sob_mc_refined.py
```

Optional flags:
```bash
python3 sob_mc_refined.py --sweep1-only   # run sweep 1 only
python3 sob_mc_refined.py --sweep2-only   # run sweep 2 only
python3 sob_mc_refined.py --n1 2000000    # change sweep 1 sample count
```

---

#### 3. `sob_fig_sky.py` — Paper Figure 1 (sky geometry)

Generates the two-panel sky-geometry figure showing the observational constraints from Matthew 2:9. Independent of all other figure scripts; can be run at any time.

**Panel (a):** Full polar sky chart (altazimuth projection, North up, azimuth clockwise) showing the azimuth corridor [190°, 210°], the guidance-phase track (arrow), and the stopping-point star symbol.

**Panel (b):** Azimuth corridor diagram from Bethlehem — road bearing 203°, corridor sweep 190°–210°, with Jerusalem marked to the NNE.

**Output:** `figure_paper1_sky_geometry.png` — Paper Fig. 1 (600 DPI)

**Runtime:** < 5 seconds.

```bash
python3 sob_fig_sky.py
```

---

#### 4. `sob_fig_update.py` — Regenerate main paper figures

Produces the main paper figures and caches the full scatter dataset. **Must be run before `sob_fig_proof.py` and `sob_fig_scatter_v2.py`.**

**Outputs:**
- `figure_paper2_az_motion.png` — Paper Fig. 2
- `figure_paper2_impossibility_scatter.png` — Paper Fig. 3 (original; retained as backup)
- `figure9_corpus_linguistics.png` — corpus linguistics heatmap
- `figure_mc_diagnostic.png` — MC diagnostic (SI Fig. S1)
- `scatter_cache.npz` — cached scatter data for `sob_fig_proof.py` and `sob_fig_scatter_v2.py`

**Runtime:** 20–60 minutes (re-runs orbit computations for the scatter plots).

```bash
python3 sob_fig_update.py
```

---

#### 4b. `sob_fig_scatter_v2.py` — Paper Figure 3 (revised version)

Loads `scatter_cache.npz` and produces the revised impossibility scatter plot used in the submitted manuscript. Key differences from the original: non-corridor points are coloured by eccentricity (plasma colourmap), and the purple Keplerian constraint band (empirical 1st–99th percentile of stopping/guidance speed ratio: ×0.995–×1.033) makes the structural nature of the exclusion visually explicit.

**Requires:** `scatter_cache.npz` must exist (produced by `sob_fig_update.py`).

**Output:** `figure_paper2_impossibility_scatter_v2.png` — Paper Fig. 3 (600 DPI)

**Runtime:** < 30 seconds.

```bash
python3 sob_fig_scatter_v2.py
```

---

#### 4. `sob_fig_proof.py` — Extended Data Figure 1

Regenerates the four-panel orbit-inversion proof figure using the scatter cache from `sob_fig_update.py`.

**Requires:** `scatter_cache.npz` must exist (produced by `sob_fig_update.py`).

**Output:** `figure6_impossibility_proof.png`

```bash
python3 sob_fig_proof.py
```

---

### Supplementary analyses

These scripts are independent of each other and of the core pipeline. They can be run in any order.

---

#### `sob_proof.py` — Kinematic analysis and preliminary survey

An earlier, self-contained script that establishes the kinematic argument (Proof A: required ICRS velocities are unachievable for any heliocentric body) and runs a coarser preliminary orbital survey. Prints full results to stdout. This script produced the foundational kinematic argument in §1.3 of the paper; `sob_proof_fine.py` supersedes its survey results with a 16× denser grid.

```bash
python3 sob_proof.py
```

---

#### `sob_epoch_robustness.py` — Multi-epoch robustness check

Tests the null result across all major proposed Star of Bethlehem dates: Halley's Comet (12 BCE), the Jupiter–Saturn triple conjunction (7 BCE), the reference epoch (5 BCE), and lunar-eclipse proposals (1 BCE). Confirms zero configurations satisfy all criteria at any proposed date within ±40 days of each epoch.

```bash
python3 sob_epoch_robustness.py
```

---

#### `sob_geocentric_survey.py` — Geocentric orbit survey

Surveys 51,200 hypothetical Earth-orbiting configurations (semi-major axes from 7,000 km to 500,000 km, all eccentricities and orientations). Identifies the 28,800 physically admissible subset (perigee > 100 km) and confirms zero satisfy all criteria.

```bash
python3 sob_geocentric_survey.py
```

---

#### `sob_geocentric_analysis.py` — Detailed geocentric analysis

Companion to `sob_geocentric_survey.py`. Provides a more detailed breakdown of the geocentric results, including the Moon's peak apparent motion at perigee and the meridian-transit azimuth artefact in geosynchronous orbits.

```bash
python3 sob_geocentric_analysis.py
```

---

#### `sob_close_flyby_mc.py` — Close-flyby Monte Carlo survey

Surveys 10 million randomly sampled comet orbits for the close-flyby scenario (Earth passage within 0.01 AU). Tests whether any close-flyby orbit can be visible at night from Bethlehem and simultaneously satisfy the guidance and stopping criteria. Finds zero of 66 genuine close-flyby orbits satisfy all criteria.

**Runtime:** 1–2 hours.

```bash
python3 sob_close_flyby_mc.py
```

---

#### `sob_close_flyby_opt.py` — Close-flyby numerical optimizer

Complements the Monte Carlo survey by finding the single best close-flyby orbit using differential evolution (global) followed by BFGS (local polish). Demonstrates that even the optimal orbit achieves only 0.58 h of the required 1.0 h guidance duration.

```bash
python3 sob_close_flyby_opt.py
```

---

#### `sob_flyby_detail.py` — Best close-flyby orbit detail

Prints a time-series table (JD, geocentric distance, ICRS angular velocity, ground-frame angular velocity, altitude, azimuth, sun altitude, azimuth and altitude rates) for the best orbit found by `sob_close_flyby_opt.py`. The orbital parameters are hardcoded from that optimizer run.

```bash
python3 sob_flyby_detail.py
```

---

#### `sob_flyby_figure.py` — Best close-flyby orbit visualisation

Generates the three-panel figure of the best close-flyby orbit: a stereographic sky chart, azimuth/altitude time series, and angular velocity time series.

**Output:** `figure_flyby_bestorbit.pdf`, `figure_flyby_bestorbit.png`

```bash
python3 sob_flyby_figure.py
```

---

#### `sob_calcs.py` — Supporting orbital mechanics calculations

Three auxiliary calculations: (1) the Matney (2025) comet's actual azimuth track from Bethlehem on 5 BCE Jun 8, showing it does not lie in the road corridor; (2) the diurnal motion constraint showing no fixed sky object maintains azimuth 206° while rising; (3) the velocity and distance requirements for any "temporary geosynchronous" interpretation. Prints results to stdout.

```bash
python3 sob_calcs.py
```

---

#### `sob_figures.py` — Early figure generation

An earlier version of the figure-generation pipeline, retained for reproducibility of the preliminary figures (figures 1–5 and 7–8). The current paper figures are produced by `sob_fig_update.py` and `sob_fig_proof.py`.

```bash
python3 sob_figures.py
```

---

### Statistical scripts

---

#### `guiding_star_bayes.py` — Corpus and within-Matthew Bayes factor computation

Computes two Bayes factors: BF_corpus (comparing Matthew 2:9's four-feature TLG vocabulary profile under the guiding-star narrative hypothesis vs. the astronomical hypothesis, using Laplace smoothing over N=4 guiding-star traditions with sensitivity analyses over the smoothing parameter k and over extended corpus definitions) and BF_Matthew (within-Matthew register analysis).

**Key output:** BF_corpus ≈ 6,012 (baseline, k=1); range 4,546–6,818 for k ∈ [0.5, 3]; conservative lower bound 1,040 when all eight omen texts are pooled.

```bash
python3 guiding_star_bayes.py
```

---

#### `corpus_pvalue.py` — Fisher combined p-value and Gaussian copula independence robustness test

Computes the Fisher combined p-value for the four TLG corpus features using TLG-verified hit counts. Implements a one-sided Fisher exact test (p = 1/(N+1)) for each feature and combines them using the Fisher chi-squared method. Also runs a Gaussian copula independence-assumption robustness analysis (Brown's method with variance-corrected chi-squared combination) over a range of equicorrelation coefficients ρ.

**Key output:** Combined p ≈ 0.006.

```bash
python3 corpus_pvalue.py
```

---

### Utility

---

#### `comet_orbit.py` — Orbit determination from observations

Fits Keplerian orbital elements to a set of astrometric observations (RA, Dec, Julian Date). Supports parabolic fits (e=1 fixed) and free-eccentricity fits (elliptic, parabolic, or hyperbolic). Uses differential evolution for global minimisation followed by Nelder-Mead polishing. Useful for independently verifying proposed comet orbital elements against historical observation records.

```bash
# Parabolic fit (default)
python3 comet_orbit.py observations.txt

# Free-eccentricity fit
python3 comet_orbit.py observations.txt --fit-eccentricity
```

The input file should be whitespace-separated columns: `JD  RA_deg  Dec_deg`.

---

## Full reproduction pipeline

To reproduce all numerical results reported in the paper from scratch:

```bash
# Step 1 — Main heliocentric survey (hours; uses all CPU cores)
python3 sob_proof_fine.py

# Step 2 — Monte Carlo refinement (30–90 min)
python3 sob_mc_refined.py

# Step 3 — Supplementary surveys (run in any order; each takes minutes to ~2 hours)
python3 sob_proof.py
python3 sob_epoch_robustness.py
python3 sob_geocentric_survey.py
python3 sob_geocentric_analysis.py
python3 sob_close_flyby_mc.py        # ~1–2 hours
python3 sob_close_flyby_opt.py
python3 sob_flyby_detail.py
python3 sob_calcs.py

# Step 4 — Statistical analyses (seconds each)
python3 guiding_star_bayes.py
python3 corpus_pvalue.py

# Step 5 — Figures (run in order; sob_fig_update.py must precede sob_fig_proof.py)
python3 sob_fig_sky.py               # Paper Fig. 1 — independent, < 5 s
python3 sob_fig_update.py            # ~20–60 min; writes scatter_cache.npz
python3 sob_fig_proof.py             # requires scatter_cache.npz
python3 sob_flyby_figure.py
python3 sob_figures.py               # early figures only
```

---

## Data files

| File | Description |
|------|-------------|
| `scatter_cache.npz` | Cached scatter data (~18 MB) from `sob_fig_update.py`; loaded by `sob_fig_proof.py` |
| `sob_mc_results.json` | Summary statistics from Monte Carlo refinement sweeps |
| `sob_mc_refined.log` | Full log of Monte Carlo refinement run |
| `sob_proof_fine_results.json` | Summary results from the fine-grid heliocentric survey |
| `sob_proof_fine_progress.log` | Progress log from fine-grid survey run |
| `tlg_corpus_results.md` | TLG corpus search results (manually verified, 2026-05-10) |
| `TLG_search_protocol.md` | Step-by-step protocol for reproducing TLG corpus searches |

---

## Building the paper

Both files use **XeLaTeX** (not pdflatex). `sob_paper.tex` uses the `xr` package to cross-reference labels defined in `sob_si.tex`, so **`sob_si.tex` must be compiled first** to produce `sob_si.aux`.

```bash
# 1. Compile the SI first (generates sob_si.aux, needed by sob_paper.tex)
xelatex sob_si.tex
bibtex sob_si
xelatex sob_si.tex
xelatex sob_si.tex

# 2. Compile the main paper (reads sob_si.aux for \ref{subsec:mc_refined})
xelatex sob_paper.tex
bibtex sob_paper
xelatex sob_paper.tex
xelatex sob_paper.tex
```

Or with `latexmk` (configuration in `.latexmkrc`):

```bash
latexmk -pdf -xelatex sob_si.tex
latexmk -pdf -xelatex sob_paper.tex
```

On **Overleaf**, set the compiler to XeLaTeX (Menu → Compiler). Overleaf compiles all project files together so the ordering is handled automatically.

> **Note on cross-references:** Section references to the SI (e.g., `§S1.7`) will appear as `??` in `sob_paper.pdf` if `sob_si.aux` does not exist or is stale. Re-run the SI compilation step to fix this.

---

## Contact

Aaron Adair — adairaar@gmail.com
