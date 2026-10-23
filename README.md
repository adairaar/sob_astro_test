# Star of Bethlehem: orbital mechanics and corpus linguistics

Code, data and manuscript files for:

**"Matthew's Star of Bethlehem is inconsistent with natural objects: evidence from orbital mechanics and corpus linguistics"**
Aaron M. Adair (Department of Physics, Massachusetts Institute of Technology), under review at *npj Heritage Science* (2026).
Preprint (earlier version): [doi:10.21203/rs.3.rs-10944256/v1](https://doi.org/10.21203/rs.3.rs-10944256/v1)

The paper tests whether the motion described in Matthew 2:9, a star that goes before travelers and then stands over a place, can be produced by any natural object, by surveying orbital parameter space against explicit criteria drawn from the text. It also tests whether Matthew's four motion terms belong to the Greek astronomical register. Within the readings and parameter space tested, no orbit satisfies the description, and the vocabulary belongs to guiding-star legendary narrative.

---

## Requirements

Python 3.10 or later.

```bash
pip install numpy scipy matplotlib astropy pymeeus
```

**Ephemeris.** The AU-scale surveys take Earth's position from the ephemeris built into `astropy` (ERFA `epv00`). At 5 BCE its error displaces sightlines by ≤0.024° beyond 0.5 AU and ≤0.12° at 0.1 AU, negligible against the 20° road corridor. Every analysis of objects within 0.05 AU (the close-flyby campaigns, and Matney's proposed comet) takes Earth's position from VSOP87 via `pymeeus`, with precession to the equinox of date and topocentric parallax, through `sob_frames.py`. No modern JPL ephemeris extends to 5 BCE.

---

## Repository contents

### Manuscript and Supplementary Information (revised version, October 2026)

| File | Description |
|------|-------------|
| `sob_paper_v4.tex`, `sob_paper_v4.pdf` | Revised manuscript |
| `sob_si_v4.tex`, `sob_si_v4.pdf` | Revised Supplementary Information (Supplementary Notes 1–30, Tables 1–20, Figures 1–5) |
| `sob_paper.tex`, `sob_si.tex` (and `.pdf`) | Earlier version, as first reviewed |
| `sob_refs.bib` | Bibliography |
| `sn-jnl.cls`, `sn-nature.bst`, `breakurl.sty` | Springer Nature LaTeX class and style files |

Build with XeLaTeX (or pdfLaTeX with `texlive-lang-greek`). The revised manuscript cites the SI by number and needs no cross-document references:

```bash
xelatex sob_si_v4.tex && bibtex sob_si_v4 && xelatex sob_si_v4.tex && xelatex sob_si_v4.tex
xelatex sob_paper_v4.tex && bibtex sob_paper_v4 && xelatex sob_paper_v4.tex && xelatex sob_paper_v4.tex
```

### Figures

Final figure files for the revised manuscript are in `npj_figures_2026-10/` (PNG, TIFF and EPS, 600 dpi, without embedded titles).

| Figure | Final file | Source image | Script |
|--------|------------|--------------|--------|
| Fig. 1 | `Figure1.*` | `figure_paper1_sky_geometry.png` | `sob_fig_sky.py` |
| Fig. 2 | `Figure2.*` | `figure_paper2_az_motion.png` | `sob_fig_update.py` |
| Fig. 3 | `Figure3.*` | `figure_paper2_impossibility_scatter_v2.png` | `sob_fig_scatter_v2.py` (reads `scatter_cache.npz`) |
| Fig. 4 | `Figure4.*` | `figure9_corpus_linguistics.png` | `sob_fig_update.py` |
| Suppl. Fig. 1 | — | `figure_mc_diagnostic.png` | `sob_fig_update.py` |
| Suppl. Fig. 2 | — | `figure_flyby_bestorbit.png` / `.pdf` | `sob_flyby_figure.py` |
| Suppl. Fig. 3 | — | `fig_criteria_sweep.png` | `sob_bigsweep_plots.py` |
| Suppl. Fig. 4 | — | `fig_false_positive.png` | `sob_falsepositive_figure.py` |
| Suppl. Fig. 5 | — | `figure6_impossibility_proof.png` | `sob_proof.py` |

Figures 2 and 3 in `npj_figures_2026-10/` are the source images with their title bands removed; Figure 4 was regenerated without its title.

---

## Scripts, by analysis

Table 1 of the manuscript lists the frame, windows and criteria of each analysis.

| Analysis (manuscript section) | Script(s) | Output |
|---|---|---|
| Validated frame pipeline (VSOP87, precession, parallax) | `sob_frames.py`, `sob_precession.py` | used by the close-flyby scripts |
| Kinematic argument and demonstration grid ("Kinematic demonstration of incompatibility") | `sob_proof.py` | Table 2; Suppl. Fig. 5 |
| Heliocentric grid, 529,079,040 configurations ("Exhaustive survey of heliocentric orbits") | `sob_proof_fine.py` | `sob_proof_fine_results.json` |
| Monte Carlo refinement, 2.25×10⁶ | `sob_mc_refined.py` | `sob_mc_results.json`, `sob_mc_refined.log` |
| Seven epochs, 12–1 BCE | `sob_epoch_robustness.py` | log output |
| Geocentric survey, 28,800 admissible ("Survey of Earth-orbiting objects") | `sob_geocentric_survey.py`, `sob_geocentric_analysis.py` | log output |
| Close flybys, 10⁷ orbits, and optimizer ("Survey of close flybys") | `sob_close_flyby_mc.py`, `sob_close_flyby_opt.py`, `sob_flyby_detail.py` | Suppl. Note 27 |
| Close flybys, 10⁹ orbits, 525 parameter combinations | `sob_bigsweep.py`, `sob_criteria_sweep.py`, `sob_itinerary_sweep.py` | `bigsweep_*.npz`; Table 3; Suppl. Note 29 |
| Two-leg itineraries (east) | `sob_finalleg_sweep.py` | Suppl. Note 29 |
| Two-leg itineraries (west and far east; revision) | `sob_finalleg_sweep_allaz.py` | `finalleg_allaz_results_2026-10-15.txt` |
| Reviewer-proposed orbits | `r2_orbits_verified.py` | Suppl. Note 28 |
| Matney's comet, corrected topocentric sky track | `matney_skytrack_verified.py` | `matney_skytrack_verified.csv`; Suppl. Note 9 |
| Corpus Fisher test and copula robustness | `corpus_pvalue.py` | Table 4; Suppl. Table 11 |
| Bayes factors and sensitivity | `guiding_star_bayes.py`, `sob_bf_sensitivity.py` | Suppl. Tables 12 and 14 |

The 10⁹-orbit candidate sets are archived as `_big_stratified_1000000000_2026.pkl` and `_big_flat_1000000000_2026.pkl` (6,116 and 450 close flybys), so the sweeps can be rerun without repeating the sampling.

### Corpus records

| File | Description |
|------|-------------|
| `TLG_search_protocol.md` | Search procedure (exact wordform, diacritics-sensitive, full paradigm) |
| `TLG_HISTEMI_RESULTS_2026-09-04.md`, `HISTEMI_CLASSIFICATION_2026-09-04.md` | ἵστημι hit record and classification |
| `TLG_AUDIT_2026-09-04.md` | Audit of the corpus counts |
| `tlg_corpus_results.md` | Earlier (May 2026) counts, superseded by the 2026-09-04 records |

### Superseded or withdrawn scripts (kept as a record)

`sob_calcs.py` (sky track superseded by `matney_skytrack_verified.py`), `sob_figures.py` (early figures), `sob_contrast_figure.py` (withdrawn before submission; see its header), `comet_orbit.py` (utility).

---

## Citation and license

See `CITATION.cff`. Code is released under the MIT License (`LICENSE`).

## Contact

Aaron M. Adair — adairaar@gmail.com — ORCID [0000-0003-3806-1047](https://orcid.org/0000-0003-3806-1047)
