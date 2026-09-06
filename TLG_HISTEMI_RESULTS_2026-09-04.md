# ἵστημι in the 19-author corpus — exact wordform search, 4 September 2026

Method as Aaron proposed: paradigm from Wiktionary, then one exact wordform
search per form. 164 forms tested, all 19 corpus authors loaded as a single
search selection (125 works), 500 results per page.

**Settings that matter: `Exact Match` ON and `Diacritics sensitive` ON.**
Both are OFF by default and both were off in the original work. With them off,
στῇ returns 589 instead of 18, ἔστη returns 81 instead of 5, ἱστάς returns 38
instead of 0. Any count taken with the defaults is inflated by accent-blind
matching.

---

## Result

| | |
|---|---|
| forms tested | 164 |
| forms attested | **69** |
| raw occurrences (all works of the 19 authors) | **395** |
| **in-corpus occurrences** (Aristotle → Meteorologica only, Theophrastus → De signis only) | **121** |

The SI states that ἵστημι appears with a celestial subject four times and could
not be lemma-searched. It appears **121 times in-corpus**, in 16 of the 19
authors, and it searches without difficulty one form at a time.

### By author (in-corpus)

| author | hits |
|---|---|
| Ptolemy | 30 |
| Vettius Valens | 13 |
| Dorotheus | 11 |
| Hephaestion | 11 |
| Achilles Tatius | 10 |
| Cleomedes | 8 |
| Theon of Smyrna {1724} | 6 |
| Manetho | 6 |
| Eratosthenes | 6 |
| Theodosius | 4 |
| Paul of Alexandria | 4 |
| Aratus | 4 |
| Hipparchus | 3 |
| Geminus | 2 |
| Aristotle (Meteorologica) | 2 |
| Theophrastus (De signis) | 1 |
| Aristarchus, Autolycus, Hypsicles | 0 |

Theophrastus *De signis* is TLG work **{0093.008}** and does contain one hit —
worth noting, since the earlier notes recorded De signis as having no hits on
any feature.

### Attested forms and raw counts

ἵσταται 54, ἵστασθαι 33, ἑστάναι 32, στῆναι 21, ἕστηκεν 20, στήσεται 19,
στῇ 18, ἑστῶτες 18, ἵστανται 12, ἑστώτων 11, ἱστάμενον 10, στῆσαι 8,
ἕστηκε 8, ἑστῶτος 7, ἵστησι 6, ἱσταμένου 6, ἑστώς 6, ἑστῶτα 6, ἑστηκώς 6,
ἵστησιν 5, ἔστη 5, ἱστάναι 4, ἔστησαν 4, στάντος 4, ἱστᾶσι 3, ἱστάντες 3,
σταθέντες 3, ἑστηκός 3, ἵσταμεν 2, ἱστᾶσιν 2, ἱστάν 2, ἱστῆται 2,
ἱσταμένη 2, ἱσταμένων 2, ἱσταμένοις 2, ἵστατο 2, ἵσταντο 2, ἔστησεν 2,
στήσας 2, στῶσι 2, σταθῆναι 2, σταθείς 2, ἑστᾶσιν 2, ἑστηκότι 2,
ἑστηκότα 2, and 24 forms at 1 each.

**Matthew's own form ἐστάθη: 0 occurrences in the corpus.**
That is a real and reportable result — but it is a fact about one wordform, not
about the verb, and it cannot carry Feature B on its own.

---

## What this changes

1. The stated reason for using στηρίζω as a proxy is void. ἵστημι is fully
   searchable; it just cannot be reached by typing the headword into a field
   that does prefix matching.
2. Feature B's denominator is not 17. On the verb ἵστημι alone it is 121; on
   στηρίζω it is far more than the 17 recorded; on both together, more again.
3. 121 passages now need reading and classifying. Until that is done there is
   no Feature B p-value. If all 121 are negative, p = 1/122 = 0.0082, against
   0.056 recorded.
4. Features A, C and D need the same treatment. Their lemma searches had the
   same defect **and** ran with accent-blind matching, so their counts are
   wrong in both directions at once — too narrow on inflection, too broad on
   diacritics.

## Reproducibility note for the SI

The procedure that works, and that should be written into the methods:

- load all corpus authors into My Search Selection (author autocomplete),
- set results per page to 500,
- tick **Exact Match** and **Diacritics sensitive**,
- use **Word Index** mode, one inflected form per search,
- enumerate forms from a printed paradigm, not from the search box,
- filter hits by TLG work number for authors represented in the corpus by a
  single work (Aristotle {0086.026}, Theophrastus {0093.008}).

Do not use the Lemma radio. In this deployment it prefix-matches the headword:
λέγω returns λέγω, λέγωμεν, λέγων, λέγωσι and not λέγει, ἔλεγεν or εἶπεν.
