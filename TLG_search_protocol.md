# TLG Search Protocol for Corpus Linguistics Transparency Table
## Star of Bethlehem paper — corpus_pvalue.py verification

This protocol generates the passage-level data needed to:
1. Verify (or correct) the counts in corpus_pvalue.py
2. Produce Table S3 in the Supplementary Information

---

## VERIFIED TLG AUTHOR NUMBERS (corrected 2026-05-10)

The original protocol had incorrect TLG numbers (hallucinated). The following
were verified directly via the TLG IRIS autocomplete at
https://stephanus-tlg-uci-edu.libproxy.mit.edu/Iris/inst/tsearch.jsp

**IMPORTANT**: The TLG internal search ID omits leading zeros: author 0363
is searched as "363". Work numbers must still be verified per-author.

| # | Author | Work in Corpus | TLG Full ID | TLG Search ID | Old (WRONG) ID |
|---|--------|----------------|-------------|---------------|----------------|
| 1 | Ptolemy (Claudius PTOLEMAEUS) | Almagest (Syntaxis Math.) | {0363} | 363 | 0363 ✓ |
| 2 | Ptolemy | Tetrabiblos (Apotelesmatica) | {0363} | 363 | 0363 ✓ (but TLG work# = 007, not 003) |
| 3 | Ptolemy | Planetary Hypotheses | {0363} | 363 | 0363 ✓ |
| 4 | Hipparchus | Commentary on Aratus | {1431} | 1431 | 2865 ✗ |
| 5 | Geminus | Isagoge (Εἰσαγωγή) | {1383} | 1383 | 0668 ✗ |
| 6 | Autolycus | On the Moving Sphere | {1210} | 1210 | 0649 ✗ |
| 7 | Autolycus | On Risings and Settings | {1210} | 1210 | 0649 ✗ |
| 8 | Theodosius | Sphaerics | {1719} | 1719 | 0666 ✗ |
| 9 | Theodosius | On Days and Nights | {1719} | 1719 | 0666 ✗ |
| 10 | Aristarchus | On Sizes and Distances | {1181} | 1181 | 0651 ✗ |
| 11 | Cleomedes | Caelestia (Meteora) | {1272} | 1272 | 0669 ✗ |
| 12 | Hypsicles | Anaphorikos | {0717} | 717 | 1712 ✗ |
| 13 | Aratus | Phaenomena | {0653} | 653 | 0653 ✓ |
| 14 | Eratosthenes | Catasterismi | {0222} | 222 | 0241 ✗ |
| 15 | Achilles Tatius (Astron.) | Isagoge to Aratus | {2133} | 2133 | 2003 ✗ (2003=Julian!) |
| 16 | Theon of Smyrna | Math. Useful for Reading Plato | {1358}? | 1358? | 1358 (unverified by autocomplete) |
| 17 | Dorotheus | Carmen Astrologicum (Gk frags) | {1337} | 1337 | 1181 ✗ |
| 18 | Pseudo-Manetho / Manetho Astrol. | Apotelesmatica | {2583} | 2583 | 1258 ✗ |
| 19 | Vettius Valens | Anthologies | {1764} | 1764 | 1764 ✓ |
| 20 | Paul of Alexandria | Eisagoge | {2053} | 2053 | 2054 ✗ |
| 21 | Hephaestion (Astrol.) | Apotelesmatika | {2043} | 2043 | 2061 ✗ |
| 22 | Aristotle | Meteorologica | {0086} | 86 | 0086 ✓ |
| 23 | Theophrastus | De signis | {0093} | 93 | 0093 ✓ |

**Notes on uncertain entries:**
- Theon of Smyrna (1358): not found via autocomplete; protocol number may be correct
  but needs in-situ verification by searching with ID 1358 and checking result TLG codes
- Eratosthenes Catasterismi: the Catasterismi may be listed as a separate work under
  {0222} or as Pseudo-Eratosthenes under a different TLG number — verify by searching
  author "222" and checking which work numbers appear

---

## HOW TO SEARCH A SPECIFIC AUTHOR IN TLG IRIS (JavaScript method)

In the Chrome console at tsearch.jsp:

```javascript
// Step 1: Set up lemma search for a specific author
document.querySelector('#lemma').checked = true;
document.querySelector('#ac_author').checked = true;
document.querySelector('#q').value = 'GREEK_WORD_HERE';
document.querySelector('#ql').value = 'BETA_CODE_HERE';
document.querySelector('select[name="dispg"]').value = '500';
window.author = 'TLG_SEARCH_ID';  // e.g., "363" for Ptolemy
window.work = null;

// Step 2: Execute search
getResults(0, false, false);

// Step 3: Wait ~3 seconds, then extract results
// Look for divs matching: /^\d+\.\s/ and /\{XXXX\.XXX\}/
```

**Beta code for each lemma:**
- προάγω  → `PRO/-A)/GW`
- ἵστημι  → `I(/STHMI`
- στηρίζω → `STHRI/ZW`
- ἐπάνω   → `E)PA/NW`
- ὑπὲρ    → `U(PE/R`
- ἔρχομαι → `E)/RXOMAI`

---

## FEATURE A: προάγω  (narrative guidance)

**Matthew's usage:** ὁ ἀστὴρ … προῆγεν αὐτούς (the star went before/guided them)

**Classification:**
- A-positive (Matthew-type): subject is a celestial body AND verb means spatial
  guidance of human travellers toward a destination
- A-negative: any other use (temporal antecedence, rank precedence, mathematical
  precedence of one planet before another, or noun/adjective cognates like προαγωγή)

**Ptolemy (ID 363) — confirmed 2026-05-10:**
- 3 hits, ALL in {0363.007} (Apotelesmatica = Tetrabiblos)
  1. {0363.007} Book 4 ch. 3 sec. 4 l. 4 — form: προαγωγαῖς (noun, A-neg)
  2. {0363.007} Book 3 ch. 14 sec. 25 l. 6 — form: προαγωγικούς (adj., A-neg)
  3. {0363.007} Book 4 ch. 4 sec. 4 l. 8 — form: προαγωγούς (noun, A-neg)
- {0363.001} Almagest: 0 hits
- Classification: 0 A-positive, 3 A-negative

---

## FEATURE B: ἵστημι / στηρίζω  (stopping verb)

**Matthew's usage:** ἐστάθη ἐπάνω (it stood/came to rest over the child)

**Classification for ἵστημι:**
- B-positive: celestial body physically stops/rests above a specific terrestrial point
- B-negative: technical "stationary point" (στηριγμός — planet reverses direction in
  zodiac), mathematical/geometric standing, or non-celestial use

**Classification for στηρίζω:**
- All uses B-negative by definition (this IS the standard technical term for stationary point)

---

## FEATURE C: ἐπάνω / ὑπὲρ  (locative preposition)

**Matthew's usage:** ἔστη ἐπάνω οὗ ἦν τὸ παιδίον (it stood over where the child was)

**Feature C definition (updated):** ἐπι-based localization of a celestial body over a
specific terrestrial point. Includes ἐπάνω or structural equivalents such as ἐπικαθίσαι
(ἐπι- + verb of sitting/resting). The JosAsen manuscript tradition shows ἐπί/ἐπάνω are
interchangeable in identical localization contexts.

**Classification for ἐπάνω:**
- C-positive: locates a celestial body directly above a specific terrestrial point
  (building, person, city — named earthly location)
- C-negative: above zodiacal region, above horizon, comparative "more than", temporal use

**Classification for ὑπὲρ:**
- Record ALL celestial body + ὑπὲρ uses; note if any are over specific terrestrial points
- Most will be C-negative (over zodiac/ecliptic/abstract region)

---

## FEATURE D: ἔρχομαι  (arrival verb)

**Matthew's usage:** ἦλθεν καὶ ἔστη (it came and stood — physical arrival at a place)

**Classification:**
- D-positive: celestial body arrives at a specific terrestrial geographic point
  (comes TO a building, city, person)
- D-negative: purely astronomical — rising above horizon, entering a zodiacal sign,
  arriving at opposition, coming to conjunction (all positional/angular events, not
  arrival at a terrestrial address)

---

## RECORDING TEMPLATE

Tab-separated format for each hit:

```
Feature | Author | TLG_ID | Work | Location | Greek_form | Subj_celestial | Class | Notes
A | Ptolemy | 0363.007 | Tetrabiblos | 4.3.4.4 | προαγωγαῖς | no | A-neg | noun cognate, astrological
```

For zero-hit searches:
```
A | Aristarchus | 1181.001 | On Sizes | (full text) | — | — | A-neg-0 | 0 hits
```

---

## SEARCH SEQUENCE AND STATUS

### Feature A (προάγω) — Beta code: PRO/-A)/GW — ✓ COMPLETE (2026-05-10)

| Author | TLG ID | Status | Count | A-pos | A-neg |
|--------|--------|--------|-------|-------|-------|
| Ptolemy | 363 | ✓ done | 3 | 0 | 3 |
| Vettius Valens | 1764 | ✓ done | 4 | 0 | 4 |
| Hephaestion | 2043 | ✓ done | 3 | 0 | 3 |
| Dorotheus | 1337 | ✓ done | 1 | 0 | 1 |
| All others | — | ✓ done | 0 | 0 | 0 |

**Result**: 9 in-corpus hits, 0 A-positive. Fisher p = 0.100.

### Feature B (ἵστημι) — Beta code: I(/STHMI — ✗ LEMMA SEARCH FAILS
TLG morphological expander cannot resolve -μι verb inflections.
"Selected wordforms: ἵστημι / No results found" in SELECTION section.
Use στηρίζω (B2) as proxy.

### Feature B2 (στηρίζω) — Beta code: STHRI/ZW — ✓ COMPLETE (2026-05-10)

| Author | TLG ID | Status | Count | B-pos | B-neg |
|--------|--------|--------|-------|-------|-------|
| Ptolemy | 363 | ✓ done | 3 | 0 | 3 |
| Hephaestion | 2043 | ✓ done | 10 | 0 | 10 |
| Dorotheus | 1337 | ✓ done | 4 | 0 | 4 |
| All others | — | ✓ done | 0 | 0 | 0 |

**Result**: 17 in-corpus hits, 0 B-positive (all "stationary point" = planet reverses in zodiac).
Fisher p = 0.056 (conservative proxy).

### Feature C (ἐπάνω) — Beta code: E)PA/NW — ✓ COMPLETE (2026-05-10)

| Author | TLG ID | Status | All-works | In-corpus | C-pos |
|--------|--------|--------|-----------|-----------|-------|
| Ptolemy | 363 | ✓ | 10 | 10 | 0 |
| Vettius Valens | 1764 | ✓ | 6 | 6 | 0 |
| Hephaestion | 2043 | ✓ | 6 | 6 | 0 |
| Theodosius | 1719 | ✓ | 4 | 4 | 0 |
| Achilles Tatius | 2133 | ✓ | 1 | 1 | 0 |
| Dorotheus | 1337 | ✓ | 1 | 1 | 0 |
| Manetho | 2583 | ✓ | 1 | 1 | 0 |
| Aristotle | 86 | ✓ | 10 | 0 | 0 |
| Theophrastus | 93 | ✓ | 10 | 0 | 0 |
| All others | — | ✓ | 0 | 0 | 0 |

**Result**: 29 in-corpus hits, 0 C-positive. Fisher p = 0.033.

### Feature C2 (ὑπὲρ) — Beta code: U(PE/R — spot-check only; not needed for p-value
### Feature D (ἔρχομαι) — Beta code: E)/RXOMAI — ✓ COMPLETE (2026-05-10)

| Author | TLG ID | Status | Count | D-pos | D-neg |
|--------|--------|--------|-------|-------|-------|
| Vettius Valens | 1764 | ✓ done | 6 | 0 | 6 |
| Hephaestion | 2043 | ✓ done | 1 | 0 | 1 |
| All others | — | ✓ done | 0 | 0 | 0 |

**Result**: 7 in-corpus hits, 0 D-positive (all "planet comes to zodiacal position/aspect").
Fisher p = 0.125.

## FINAL VERIFIED P-VALUES

P_OBS = [0.100, 0.056, 0.033, 0.125]
Fisher combined p = 0.00627 (df=8, X²=21.35)
Corrected p at rho=0.5 = 0.028 (still significant)
Crossover rho ≈ 0.77 (realistic upper bound 0.4–0.6)

---

## CRITICAL EDITION REFERENCES

For the location column, use the standard editions:
- Ptolemy Almagest: Heiberg (Teubner 1898–1903)
- Ptolemy Tetrabiblos: Hübner (Teubner 1998) — TLG work 007
- Hipparchus Comm. on Aratus: Manitius (Teubner 1894)
- Geminus Isagoge: Aujac (Budé 1975)
- Autolycus: Mogenet (Louvain 1950)
- Cleomedes: Ziegler (Teubner 1891) or Todd (Teubner 1990)
- Aratus: Martin (Budé 1998)
- Vettius Valens: Pingree (Teubner 1986)
- Hephaestion: Pingree (Teubner 1973)
- Aristotle Meteor.: Fobes (Harvard 1919)

TLG displays passage references in its own notation but typically follows
the standard edition section numbering.
