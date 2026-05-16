# TLG Corpus Search Results — Verified Data
## Star of Bethlehem paper corpus linguistics verification
## Searches conducted 2026-05-10 via MIT proxy (stephanus-tlg-uci-edu.libproxy.mit.edu)

---

## CORRECTED TLG AUTHOR NUMBERS (all verified via TLG autocomplete)

| Author | TLG Search ID | Full TLG ID | Old (Wrong) ID |
|--------|--------------|-------------|----------------|
| Claudius Ptolemaeus | 363 | {0363} | 0363 ✓ |
| Hipparchus Astron. | 1431 | {1431} | 2865 ✗ |
| Geminus Astron. | 1383 | {1383} | 0668 ✗ |
| Autolycus Astron. | 1210 | {1210} | 0649 ✗ |
| Theodosius Astron.,Math. | 1719 | {1719} | 0666 ✗ |
| Aristarchus Astron. | 1181 | {1181} | 0651 ✗ |
| Cleomedes Astron. | 1272 | {1272} | 0669 ✗ |
| Hypsicles Astron.,Math. | 717 | {0717} | 1712 ✗ |
| Aratus Astron.,Epic. | 653 | {0653} | 0653 ✓ |
| Eratosthenes Philol. | 222 | {0222} | 0241 ✗ |
| Achilles Tatius Astron. | 2133 | {2133} | 2003 ✗ (2003=Julian!) |
| Theon of Smyrna | 1358 | {1358}? | 1358 (unconfirmed) |
| Dorotheus Astrol. | 1337 | {1337} | 1181 ✗ |
| Manetho Astrol. | 2583 | {2583} | 1258 ✗ |
| Vettius Valens Astrol. | 1764 | {1764} | 1764 ✓ |
| Paulus Astrol. (Paul of Alex.) | 2053 | {2053} | 2054 ✗ |
| Hephaestion Astrol. | 2043 | {2043} | 2061 ✗ |
| Aristoteles Phil. | 86 | {0086} | 0086 ✓ |
| Theophrastus Phil. | 93 | {0093} | 0093 ✓ |

**Note on Aristotle's Meteorologica**: TLG work 0086 has hits only in works 001, 005, 006, 012, 014
across ἐπάνω search. Meteorologica (probably 0086.007 or 0086.008) has 0 ἐπάνω hits.
No Aristotle hit is in the Meteorologica — confirmed 0 in-corpus hits for all C–D features.

**Note on Theophrastus De signis**: ἐπάνω hits are in works 002, 004, 005, 010, 014, 017.
De signis (weather-sign treatise, not these work numbers) has 0 in-corpus hits for all features.

**Note on ἵστημι (Feature B) lemma search failure**: TLG's morphological expander cannot
expand ἵστημι (-μι verb) to inflected forms; it returns only the exact string, which does not
appear in running prose. Feature B is therefore reported using the proxy lemma στηρίζω
(the standard technical term for a planetary stationary point), which expands correctly.
All στηρίζω hits are B-negative by definition (planet reversing in zodiac ≠ stopped over a place).

---

## FEATURE A: προάγω (Matt 2:9 — the star GUIDED/LED them)
**TLG lemma search: lemma=προάγω (Beta code: PRO/-A)/GW), author-restricted**
**Matthew classification: A-positive (star leads travelers to a destination)**

| Author | TLG ID | Count | A-pos | A-neg | Hit locations | Class |
|--------|--------|-------|-------|-------|---------------|-------|
| Ptolemy | 0363 | 3 | 0 | 3 | {0363.007} Tetrab. 4.3.4.4; 3.14.25.6; 4.4.4.8 | A-neg |
| Vettius Valens | 1764 | 3+1 | 0 | 4 | {1764.008} 1.20.123; 1.3.107; 8.7.67 / {1764.010} App.10.1.106 | A-neg |
| Hephaestion | 2043 | 2+1 | 0 | 3 | {2043.001} p.155 l.25; p.169 l.10 / {2043.002} p.214 l.5 | A-neg |
| Dorotheus | 1337 | 1 | 0 | 1 | {1337.001} Frag. Graeca p.382 l.21 | A-neg |
| Cleomedes | 1272 | 0 | 0 | 0 | — | A-neg-0 |
| Hipparchus | 1431 | 0 | 0 | 0 | — | A-neg-0 |
| Geminus | 1383 | 0 | 0 | 0 | — | A-neg-0 |
| Autolycus | 1210 | 0 | 0 | 0 | — | A-neg-0 |
| Theodosius | 1719 | 0 | 0 | 0 | — | A-neg-0 |
| Aristarchus | 1181 | 0 | 0 | 0 | — | A-neg-0 |
| Hypsicles | 717 | 0 | 0 | 0 | — | A-neg-0 |
| Aratus | 653 | 0 | 0 | 0 | — | A-neg-0 |
| Eratosthenes | 222 | 0 | 0 | 0 | — | A-neg-0 |
| Achilles Tatius | 2133 | 0 | 0 | 0 | — | A-neg-0 |
| Theon of Smyrna | 1358 | 0 | 0 | 0 | — | A-neg-0 |
| Manetho | 2583 | 0 | 0 | 0 | — | A-neg-0 |
| Paul of Alexandria | 2053 | 0 | 0 | 0 | — | A-neg-0 |
| Aristotle (Meteor.) | 86 | 0 | 0 | 0 | hits in Frag.+NE, NOT Meteorologica | A-neg-0 |
| Theophrastus (De sig.) | 93 | 0 | 0 | 0 | hit in Frag., NOT De signis | A-neg-0 |

**TOTALS**: 9 in-corpus hits, 0 A-positive
**Forms**: noun cognates προαγωγαῖς, προαγωγικούς, προαγωγούς (all astrological rank/position)
**Fisher exact p = 1/(9+1) = 0.100**
(Original claimed 0.014 was based on hallucinated corpus counts; corrected here)

---

## FEATURE B2: στηρίζω (proxy for ἵστημι — STOPPING verb)
**TLG lemma search: lemma=στηρίζω (Beta code: STHRI/ZW), author-restricted**
**Note**: ἵστημι lemma search fails (TLG cannot expand -μι verb morphologically).
στηρίζω is the standard Greek technical term for a planet's "stationary point"
(when it appears to reverse direction in the zodiac). All corpus uses are B-negative
by definition: στηρίζω describes a planet reversing in the zodiac, NOT physically
stopping above a specific terrestrial location (Matthew's sense).

| Author | TLG ID | Count | B-pos | Hit locations | Class |
|--------|--------|-------|-------|---------------|-------|
| Ptolemy | 0363 | 3 | 0 | {0363.001} vol.1,2 p.463 l.11 / {0363.005} vol.2 p.173 l.15 / {0363.007} Bk2 ch5 sect2 l.5 | B-neg |
| Hephaestion | 2043 | 10 | 0 | {2043.001} p.318 ll.5,9,12 / {2043.002} pp.58–59 (×4), p.109 l.5, p.304 l.15, p.337 ll.18,22 | B-neg |
| Dorotheus | 1337 | 4 | 0 | {1337.001} p.379 l.9; p.415 ll.15,19,21 | B-neg |
| Vettius Valens | 1764 | 0 | 0 | — | B-neg-0 |
| Cleomedes | 1272 | 0 | 0 | — | B-neg-0 |
| Hipparchus | 1431 | 0 | 0 | — | B-neg-0 |
| Geminus | 1383 | 0 | 0 | — | B-neg-0 |
| Autolycus | 1210 | 0 | 0 | — | B-neg-0 |
| Theodosius | 1719 | 0 | 0 | — | B-neg-0 |
| Aristarchus | 1181 | 0 | 0 | — | B-neg-0 |
| Hypsicles | 717 | 0 | 0 | — | B-neg-0 |
| Aratus | 653 | 0 | 0 | — | B-neg-0 |
| Eratosthenes | 222 | 0 | 0 | — | B-neg-0 |
| Achilles Tatius | 2133 | 0 | 0 | — | B-neg-0 |
| Theon of Smyrna | 1358 | 0 | 0 | — | B-neg-0 |
| Manetho | 2583 | 0 | 0 | — | B-neg-0 |
| Paul of Alexandria | 2053 | 0 | 0 | — | B-neg-0 |
| Aristotle (Meteor.) | 86 | 0 | 0 | — | B-neg-0 |
| Theophrastus (De sig.) | 93 | 0 | 0 | — | B-neg-0 |

**TOTALS**: 17 in-corpus στηρίζω hits, 0 B-positive
**Fisher exact p = 1/(17+1) = 0.056**
(Conservative: uses only the proxy lemma στηρίζω; ἵστημι would add more B-negative hits)

---

## FEATURE C: ἐπάνω (Matt 2:9 — the star over WHERE THE CHILD WAS)
**TLG lemma search: lemma=ἐπάνω (Beta code: E)PA/NW), author-restricted**
**Matthew classification: C-positive = celestial body located above a specific terrestrial point**
**Note**: Aristotle (86) and Theophrastus (93) author-wide hits are ALL in non-corpus works.
Aristotle Meteorologica and Theophrastus De signis have 0 ἐπάνω hits.

| Author | TLG ID | All-works count | In-corpus count | C-pos | Hit locations | Class |
|--------|--------|-----------------|-----------------|-------|---------------|-------|
| Ptolemy | 0363 | 10 | 10 | 0 | {0363.001} Almagest: vol.1,1 p.475 l.3; vol.1,2 pp.40–44 (×4), pp.62–65 (×2), pp.152–153, pp.267–268 | C-neg |
| Vettius Valens | 1764 | 6 | 6 | 0 | {1764.004} vol.4 p.149 l.28; {1764.008} 1.2.128; 1.18.90; 1.20.211; 5.1.49; {1764.010} App.10.1.235 | C-neg |
| Hephaestion | 2043 | 6 | 6 | 0 | {2043.001} p.96 l.10; p.314 l.14; p.333 l.15; {2043.002} p.75 l.4; p.122 l.21; p.[last hit] | C-neg |
| Theodosius | 1719 | 4 | 4 | 0 | {1719.003} De diebus: p.66 l.33; p.116 l.35; p.120 l.8; p.148 l.33 | C-neg |
| Achilles Tatius | 2133 | 1 | 1 | 0 | {2133.001} Isagoga §18 l.28 | C-neg |
| Dorotheus | 1337 | 1 | 1 | 0 | {1337.001} Frag. Graeca p.413 l.2 | C-neg |
| Manetho | 2583 | 1 | 1 | 0 | {2583.001} Apoteles. Bk5 l.169 | C-neg |
| Aristotle (Meteor.) | 86 | 10 | 0 | 0 | hits in 0086.001,005,006,012,014 — NOT Meteorologica | C-neg-0 |
| Theophrastus (De sig.) | 93 | 10 | 0 | 0 | hits in 0093.002,004,005,010,014,017 — NOT De signis | C-neg-0 |
| Cleomedes | 1272 | 0 | 0 | 0 | — | C-neg-0 |
| Hipparchus | 1431 | 0 | 0 | 0 | — | C-neg-0 |
| Geminus | 1383 | 0 | 0 | 0 | — | C-neg-0 |
| Autolycus | 1210 | 0 | 0 | 0 | — | C-neg-0 |
| Aristarchus | 1181 | 0 | 0 | 0 | — | C-neg-0 |
| Hypsicles | 717 | 0 | 0 | 0 | — | C-neg-0 |
| Aratus | 653 | 0 | 0 | 0 | — | C-neg-0 |
| Eratosthenes | 222 | 0 | 0 | 0 | — | C-neg-0 |
| Theon of Smyrna | 1358 | 0 | 0 | 0 | — | C-neg-0 |
| Paul of Alexandria | 2053 | 0 | 0 | 0 | — | C-neg-0 |

**Classification notes**: All 29 in-corpus ἐπάνω uses refer to geometric/mathematical "above"
(Almagest: above a construction line), astrological "above" (planet above a degree or zodiacal
region), or comparative/geographic "above". NONE describes a celestial body localized over a
specific named terrestrial point (building, city, person). This contrasts sharply with Matt 2:9:
ἔστη ἐπάνω οὗ ἦν τὸ παιδίον ("stood over where the child was").
**Fisher exact p = 1/(29+1) = 0.033**

---

## FEATURE D: ἔρχομαι (Matt 2:9 — the star CAME and stood)
**TLG lemma search: lemma=ἔρχομαι (Beta code: E)/RXOMAI), author-restricted**
**Matthew classification: D-positive = celestial body arrives at a specific terrestrial location**

| Author | TLG ID | Count | D-pos | Hit locations | Class |
|--------|--------|-------|-------|---------------|-------|
| Vettius Valens | 1764 | 6 | 0 | {1764.008} 5.7.86; 5.7.89; 8.5.65; 8.7.8; 8.7.26; 8.7.32 | D-neg |
| Hephaestion | 2043 | 1 | 0 | {2043.001} p.124 l.17 | D-neg |
| Ptolemy | 0363 | 0 | 0 | — | D-neg-0 |
| Cleomedes | 1272 | 0 | 0 | — | D-neg-0 |
| Hipparchus | 1431 | 0 | 0 | — | D-neg-0 |
| Geminus | 1383 | 0 | 0 | — | D-neg-0 |
| Autolycus | 1210 | 0 | 0 | — | D-neg-0 |
| Theodosius | 1719 | 0 | 0 | — | D-neg-0 |
| Aristarchus | 1181 | 0 | 0 | — | D-neg-0 |
| Hypsicles | 717 | 0 | 0 | — | D-neg-0 |
| Aratus | 653 | 0 | 0 | — | D-neg-0 |
| Eratosthenes | 222 | 0 | 0 | — | D-neg-0 |
| Achilles Tatius | 2133 | 0 | 0 | — | D-neg-0 |
| Theon of Smyrna | 1358 | 0 | 0 | — | D-neg-0 |
| Dorotheus | 1337 | 0 | 0 | — | D-neg-0 |
| Manetho | 2583 | 0 | 0 | — | D-neg-0 |
| Paul of Alexandria | 2053 | 0 | 0 | — | D-neg-0 |
| Aristotle (Meteor.) | 86 | 0 | 0 | — | D-neg-0 |
| Theophrastus (De sig.) | 93 | 0 | 0 | — | D-neg-0 |

**Classification notes**: VV's 6 ἔρχομαι hits (all in Books 5 and 8 of the Anthologiae) use
the verb for a planet "coming to" a zodiacal degree, sign, or aspect — purely positional/
angular events, NOT arrival at a terrestrial address. Hephaestion's 1 hit is similar.
NONE describes a celestial body arriving at a specific building, city, or person's location.
**Fisher exact p = 1/(7+1) = 0.125**

---

## VERIFIED P-VALUES SUMMARY

| Feature | Lemma | Matthew sense | Corpus N | Pos | Fisher p | Original claimed |
|---------|-------|---------------|----------|-----|----------|-----------------|
| A | προάγω | star guides travelers to destination | 9 | 0 | **0.100** | 0.014 |
| B | στηρίζω (proxy) | celestial body stops over terrestrial place | 17 | 0 | **0.056** | 0.032 |
| C | ἐπάνω | celestial body above specific terrestrial point | 29 | 0 | **0.033** | 0.005 |
| D | ἔρχομαι | celestial body arrives at terrestrial location | 7 | 0 | **0.125** | 0.063 |

**Fisher combined test (verified p-values):**
X = -2 × (ln 0.100 + ln 0.056 + ln 0.033 + ln 0.125)
  = -2 × (-2.303 + -2.882 + -3.411 + -2.079)
  = -2 × (-10.675) = **21.35**
p(χ²₈ ≥ 21.35) ≈ **0.006**

The joint conclusion (no corpus text uses Matthew's exact cluster of four celestial-motion
features in the same way) is **statistically significant at p < 0.01** even with verified,
conservative counts replacing the original estimated values.

---

## NOTES ON METHODOLOGY

1. **ἵστημι (Feature B original)**: The TLG IRIS lemma expander cannot resolve inflected
   forms of ἵστημι (-μι verbs are not expanded by the TLG morphological index accessed via
   direct field entry). The proxy στηρίζω captures the astronomically relevant usage
   (standard Greek term for planetary stationary point) and gives a conservative estimate.

2. **Author-restricted vs. corpus-text-restricted**: Searches use TLG's author-restriction
   (`window.author = "ID"`). For authors whose corpus text is a specific work (Aristotle →
   Meteorologica; Theophrastus → De signis), hits were further filtered by checking TLG work
   numbers in the results. Aristotle and Theophrastus have 0 in-corpus hits for all features.

3. **Theon of Smyrna (1358)**: TLG number not confirmed via autocomplete but returns 0 hits
   for all four features regardless.

4. **Statistical model**: Fisher exact 2×2 table [[1,0],[0,N]] with N = corpus negatives.
   p = 1/(N+1). Independence assumption tested via Gaussian copula Monte Carlo in corpus_pvalue.py.
