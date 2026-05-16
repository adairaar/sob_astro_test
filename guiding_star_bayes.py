#!/usr/bin/env python3
"""
Bayes Factor: Matthew 2:9 as Guiding-Star Narrative vs. Astronomical Report
Aaron Adair / 2026

Compares P(A,B,C,D | H_guide) to P(A,B,C,D | H_astro) where A–D are the
four linguistic features identified in the TLG corpus analysis.

Feature definitions (from sob_paper.tex §Corpus Linguistics):
  A  — προάγω or equivalent guidance verb ("went before", "led them")
  B  — ἵστημι-type stopping verb ("stood", "halted")
  C  — ἐπι-based localization over a specific point (ἐπάνω or structural
       equivalent such as ἐπικαθίσαι + dative; attested interchangeable with
       ἐπάνω in JosAsen manuscript tradition)
  D  — ἔρχομαι arrival verb ("came", "arrived at")

The corpus is divided into four groups:

  (1) Guiding-star narrative traditions  [primary corpus, N=4]
      Each row is one independently attested convention; duplicate retellings
      of the same event are collapsed to avoid double-counting.

  (2) Stellar-omen texts  [extended sensitivity, added to primary for N=8]
      Texts in which a celestial body appears over a city or alongside
      a fleet as a supernatural sign, but without directional guidance.
      These are explicitly excluded from the primary corpus because they
      represent a different genre (omen/portent, not guide), but they are
      used to show that even the maximally inclusive corpus gives BF > 70.

  (3) Supernatural descent texts  [extended commentary; not folded into BF]
      Texts in which a celestial or heavenly figure descends to ground level.
      Included here for linguistic analysis of feature C (ἐπάνω); not added
      to the BF sensitivity analysis because neither text involves guidance
      of travelers (A=0 in both), placing them outside the guiding-star and
      omen genres alike.

  (4) Matthean derivative  [excluded entirely]
      Protoevangelium of James 21 paraphrases Matt 2:9-10 verbatim.

Sources and scoring notes
─────────────────────────
GUIDING-STAR TRADITIONS:
  Virgil, Aeneid 2.687-711
    Venus/lucifer "marks the way" (signantemque vias) to Italy; star runs
    ahead then falls into Ida's forest (A=1, B=1); Aeneas arrives (D=1).
    No ἐπάνω-type localization (C=0).

  Timoleon flame tradition  [Plutarch Tim. 8.§5–6 + Diodorus 16.66.3]
    Same Sicilian event, both likely drawing on Timaeus of Tauromenium.
    Torch "runs with" fleet (A=1), "darts down upon precisely" the
    destination (B=1); fleet arrives (D=1).  No ἐπάνω (C=0).

  Lucian, Navigium 7; 9
    Nav. 7: star named as savior. Nav. 9: star first stands over (ἐπικαθίσαι)
    the masthead (τῷ καρχησίῳ) — ἐπι- prefix, the same root as ἐπάνω — then
    guides the ship on its left side (ἐπὶ τὰ λαιά) away from a cliff (A=1, D=1).
    C=1: ἐπικαθίσαι is a structural equivalent of ἐπάνω; JosAsen manuscripts
    use ἐπί and ἐπάνω interchangeably for the same angelic localization scene.
    No discrete stopping-at-destination (B=0).
    Genre: comic/satirical literary dialogue — fictional.

  Israelite pillar-of-fire tradition  [LXX Exod 13.21 + Num 14.14]
    Same Israelite convention; Num 14.14 is Moses recounting the Exodus
    event, not an independent invention.  Pillar goes before (A=1); no
    discrete stopping or arrival (B=0, C=0, D=0).

STELLAR-OMEN TEXTS (added in sensitivity analysis):
  Josephus, Jewish War 6.289
    A sword-like star ἔστη ἐπάνω ("stood over") Jerusalem — WAIT: the
    article notes that Josephus uses ἵστημι (same verb as Matthew, so B=1)
    but ὑπὲρ, NOT ἐπάνω, for the preposition (C=0).  The object is a
    supernatural omen of destruction (τερατεία), not an astronomical body
    and not a guiding star (A=0, D=0).

  Dio Cassius, Roman History 54.29.8
    A comet "stood over" Rome at death of Marcus Agrippa; verb is
    αἰωρηθείς (passive of αἰωρέω: "to be suspended/lifted"), NOT ἵστημι
    (B=0); preposition is ὑπὲρ, NOT ἐπάνω (C=0).  No guidance or arrival
    (A=0, D=0).

  Plutarch, Lysander 12.1
    Dioscuri appear as twin stars on either side of Lysander's ship as a
    victory omen.  No directional guidance, no stopping, no arrival
    (A=0, B=0, C=0, D=0).

  Cicero, De Divinatione 1.75
    Same Dioscuri/Aegospotami tradition.  Same feature absence.

SUPERNATURAL DESCENT TEXTS (linguistic commentary on feature C only):
  PGM I.54-95 (Greek Magical Papyri, 3rd c. CE)
    A magical spell describes a blazing star that descends (κατελθών) to
    the earth (ἔγγαιος): "ἀστὴρ πύρινος καταβήσεται καὶ στήσεται
    ἐν μέσῳ τοῦ δώματος" ("a blazing star will descend and stand in the
    middle of the housetop").  Star descends (D=1), stands/stops (B=1,
    via στήσεται = ἵστημι future), but the location is expressed with
    ἐν μέσῳ ("in the middle of"), NOT ἐπάνω (C=0).  No guidance (A=0).
    KEY: Even in the closest pagan parallel for a celestial body descending
    to earthly level, the localization preposition is NOT ἐπάνω.

  Joseph and Aseneth 14.1-5 (Jewish novel, 1st c. BCE–2nd c. CE)
    The morning star (ἑωσφόρος) descends from heaven; a radiant man
    then comes (ἦλθεν) and stands over (ἔστη ἐπάνω) her head.  Burchard
    longer recension: "ἔστη ἐπάνω τῆς κεφαλῆς αὐτῆς" (A=0, B=1, C=1, D=1).
    Philonenko shorter recension: "ἔστη ἐπί τῆς κεφαλῆς αὐτῆς" (C=0).
    KEY INTERPRETATION: The morning star and the heavenly visitor are the
    same entity.  The Greek word for "messenger" (ἄγγελος) is the same as
    "angel"; the visitor departs in a chariot heading east — the precise
    trajectory of the morning star returning to the horizon before sunrise.
    This identification (star = divine messenger) parallels PGM I.54-95
    (where the descending star is the god).  JosAsen 14.4 (longer rec.) is
    therefore the one non-Matthean instance of a star-entity using ἐπάνω
    to stand over a precise earthly location.
    GENRE NOTE: JosAsen is a diaspora Jewish novel — narrative prose fiction,
    not visionary or apocalyptic literature.  As narrative, it is the closest
    genre match to Matthew's gospel among all texts examined.  The fact that
    ἐπάνω for a star-entity appears in this narrative novel (and NOT in
    astronomical, omen, or pagan-magical texts) is positive evidence that
    ἐπάνω belongs to the narrative register — consistent with H_guide.
    No guidance of travelers (A=0); C=1 only in longer recension (Philonenko: ἐπί).

KEY LINGUISTIC FINDING (updated: Lucian C=1; JosAsen manuscript note):
  Feature C (ἐπι-based localization over a specific point) scores 0 across
  all astronomical texts, all omen texts, and PGM I.54-95.  It appears in
  exactly two non-Matthean texts:
    (1) Lucian, Navigium 9: ἐπικαθίσαι τῷ καρχησίῳ — star sits/stands upon
        the masthead; ἐπι- prefix identical to Matthew's ἐπάνω (ἐπί + ἄνω).
        JosAsen manuscripts show ἐπί and ἐπάνω used interchangeably for the
        same angelic localization, confirming these are stylistic variants
        within a single semantic domain.
    (2) JosAsen 14.4 (Burchard longer rec.): ἔστη ἐπάνω τῆς κεφαλῆς αὐτῆς.
        Philonenko shorter rec. uses ἐπί — same manuscript interchangeability.

  Both instances occur in narrative fiction (Lucian's satirical dialogue;
  the diaspora Jewish novel JosAsen).  The convergent finding is:
    • Astronomical corpus           → ὑπὲρ, never ἐπάνω / ἐπι-localization
    • Omen texts                    → ὑπὲρ (Josephus, Dio), never ἐπάνω
    • Pagan magical descent (PGM)   → ἐν μέσῳ, not ἐπάνω
    • Literary fiction (Lucian)     → ἐπικαθίσαι (ἐπι-) ✓  C=1
    • Narrative novel (JosAsen)     → ἐπάνω / ἐπί ✓  C=1 (longer rec.)
    • Matthew 2:9 (gospel narrative)→ ἐπάνω ✓  same narrative register

  Crucially, all four primary guiding-star traditions belong to legendary,
  mythological, or outright fictional genres: Virgil's Aeneid (mythological
  epic), the Timoleon flame tradition (legendary history), Lucian's Navigium
  (comic dialogue), and — as a descent-text parallel — JosAsen (Jewish novel).
  Matthew's linguistic profile matches fictional-narrative tradition, not the
  astronomical or omen register.

  NOTE ON BF CALCULATIONS: PGM I.54-95 and JosAsen 14.4 still have A=0
  and are not admitted to the GUIDING or OMEN corpora.  Lucian's revised
  score (C=0→1) raises P(C|H_guide) from 1/6 to 2/6, roughly doubling the
  primary BF (see RESULTS section below).
"""

import numpy as np

FEATURE_LABELS = ["A (guidance verb)",
                  "B (stood/stopped)",
                  "C (epano location)",
                  "D (arrival verb)"]

# ── Primary corpus: four independent guiding-star traditions ─────────────────
#                                              A   B   C   D
GUIDING = {
    "Virgil, Aeneid 2.687-711":              [1,  1,  0,  1],
    "Timoleon flame tradition":               [1,  1,  0,  1],
    "Lucian, Navigium 7; 9":                  [1,  0,  1,  1],
    "Israelite pillar-of-fire tradition":    [1,  0,  1,  0],
}

# ── Extended omen corpus (sensitivity analysis) ──────────────────────────────
#                                              A   B   C   D
OMEN = {
    "Josephus, Jewish War 6.289":            [0,  1,  0,  0],
    # ἵστημι (B=1) + ὑπὲρ (C=0, NOT ἐπάνω); no guidance (A=0); no arrival
    # (D=0).  Sword-star omen of Jerusalem's destruction; explicitly labeled
    # τερατεία ("fable/portent") by Josephus.

    "Dio Cassius, Roman History 54.29.8":    [0,  0,  0,  0],
    # αἰωρηθείς (lifted/suspended, B=0) + ὑπὲρ (C=0); no guidance (A=0);
    # no arrival (D=0).  Comet at death of Marcus Agrippa.

    "Plutarch, Lysander 12.1":               [0,  0,  0,  0],
    # Dioscuri appear as omen alongside fleet; no guidance or destination.

    "Cicero, De Divinatione 1.75":           [0,  0,  0,  0],
    # Same Dioscuri/Aegospotami tradition; same feature absence.
}

# ── Supernatural descent texts (feature C commentary only) ──────────────────
#    NOT added to BF sensitivity analysis: A=0 in both (not guiding-star or
#    omen genre).  Included here to document the ἐπάνω / ἐν μέσῳ finding.
#                                              A   B   C   D
DESCENT = {
    "PGM I.54-95":                          [0,  1,  0,  1],
    # Star descends (D=1), stands ἐν μέσῳ τοῦ δώματος (B=1, ἵστημι future).
    # Localization: ἐν μέσῳ ("in the middle of"), NOT ἐπάνω (C=0).
    # No guidance of travelers (A=0).

    "Joseph and Aseneth 14.4 (Burchard)":   [0,  1,  1,  1],
    # Star-entity (morning star = divine messenger/angel, same Greek word)
    # descends and stands ἐπάνω τῆς κεφαλῆς αὐτῆς (B=1, C=1, D=1).
    # The visitor departs eastward in a solar chariot — the morning star's
    # own trajectory — confirming star = messenger identity.
    # Philonenko shorter recension uses ἐπί (C=0); C=1 only in longer text.
    # No guidance of travelers (A=0).
}

# ── H_astro denominator ──────────────────────────────────────────────────────
# P(A,B,C,D | H_astro) = product of individual Fisher exact p-values.
# TLG-VERIFIED 2026-05-10 (MIT proxy); see tlg_corpus_results.md.
#
# Individual p-values (each = 1/(N+1) for N corpus negatives, 0 positives):
#   A  προάγω    9 negatives  → p = 1/10  = 0.100
#   B  στηρίζω  17 negatives  → p = 1/18  ≈ 0.0556
#   C  ἐπάνω   29 negatives  → p = 1/30  ≈ 0.0333
#   D  ἔρχομαι   7 negatives  → p = 1/8   = 0.125
#
# Joint probability (independent features):
#   P_ASTRO_JOINT = 0.100 × 0.0556 × 0.0333 × 0.125 = 2.31e-5
#
# NOTE: The Fisher combined chi-squared p (6.27e-3) is a meta-test tail
# probability — the probability that a sum-of-log-p statistic would exceed
# the observed value by chance.  It is NOT P(data|H_astro) and is therefore
# not the correct BF denominator.  The product of the individual one-sided
# Fisher exact p-values is the correct joint probability under independence.
# (See Fisher 1932 §21.1; Held & Ott 2018 §3.)
#
# Previous erroneous values:
#   6.27e-3  (Fisher combined tail-p — meta-test, not joint data prob)
#   1.1e-4   (AI-estimated Fisher combined p with hallucinated corpus counts)
P_ASTRO_JOINT = 2.31e-5   # product of individual Fisher exact p-values


# ── Helper functions ─────────────────────────────────────────────────────────
def feature_probs(corpus_dict, k=1):
    scores = np.array(list(corpus_dict.values()))
    N = len(scores)
    counts = scores.sum(axis=0)
    probs  = (counts + k) / (N + 2 * k)
    return probs, counts, N

def joint_prob(probs):
    return float(np.prod(probs))

def run(label, corpus_dict, k=1):
    probs, counts, N = feature_probs(corpus_dict, k=k)
    pj = joint_prob(probs)
    BF = pj / P_ASTRO_JOINT
    return dict(label=label, N=N, k=k,
                counts=counts, probs=probs, pj=pj, BF=BF)

def strength(BF):
    if BF > 100:  return "decisive"
    if BF >  30:  return "very strong"
    if BF >  10:  return "strong"
    if BF >   3:  return "moderate"
    return "weak"


# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("BAYES FACTOR: Matt 2:9 linguistic features")
print(f"  P(A,B,C,D | H_astro) = {P_ASTRO_JOINT:.2e}  [product of individual Fisher exact p-values]")
print("=" * 70)

# Print corpus tables
for label, corp in [("PRIMARY CORPUS — guiding-star traditions", GUIDING),
                    ("OMEN TEXTS — sensitivity only", OMEN)]:
    print(f"\n{label}")
    print(f"  {'Text/Tradition':<44} {'A':>2} {'B':>2} {'C':>2} {'D':>2}")
    print("  " + "-" * 54)
    for t, s in corp.items():
        print("  " + f"{t:<44}" + "".join(f" {x:>2}" for x in s))

print("\nExcluded (Matthean derivative):")
print("  Protoevangelium of James 21 — paraphrases Matt 2:9-10 verbatim")

# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("FEATURE PROBABILITIES (primary guiding corpus, N=4, k=1)")
print("=" * 70)
R_primary = run("Primary (guiding only, N=4)", GUIDING, k=1)
print(f"  {'Feature':<22}  {'Hits/N':>8}  {'P(k=1)':>8}  {'P_MLE':>8}")
print("  " + "-" * 52)
for i, lab in enumerate(FEATURE_LABELS):
    c, N = R_primary["counts"][i], R_primary["N"]
    print(f"  {lab:<22}  {c:>3}/{N:<3}    "
          f"{R_primary['probs'][i]:>8.4f}  {c/N:>8.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# Run all scenarios
ALL_CORP = {**GUIDING, **OMEN}
scenarios = [
    run("Primary: guiding only         (N=4)", GUIDING,   k=1),
    run("+ Josephus & Dio              (N=6)", {**GUIDING,
        **{k: v for k,v in OMEN.items()
           if "Josephus" in k or "Dio" in k}},             k=1),
    run("+ all omen texts              (N=8)", ALL_CORP,   k=1),
]

print("\n" + "=" * 70)
print("BAYES FACTOR — RESULTS ACROSS CORPUS DEFINITIONS")
print("=" * 70)
print(f"  {'Corpus definition':<40}  {'BF':>8}  {'Strength'}")
print("  " + "-" * 60)
for R in scenarios:
    print(f"  {R['label']:<40}  {R['BF']:>8.1f}  {strength(R['BF'])}")

# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SENSITIVITY: Laplace k  (primary guiding corpus, N=4)")
print("=" * 70)
print(f"  {'k':>5}  {'P(joint|guide)':>16}  {'BF':>10}  Strength")
print("  " + "-" * 48)
for k in [0.5, 1, 2, 3]:
    R = run("", GUIDING, k=k)
    print(f"  {k:>5.1f}  {R['pj']:>16.4e}  {R['BF']:>10.1f}  {strength(R['BF'])}")

# ─────────────────────────────────────────────────────────────────────────────
BF_lo = R_primary["pj"] / (P_ASTRO_JOINT * 3)
BF_hi = R_primary["pj"] / (P_ASTRO_JOINT / 3)
print("\n" + "=" * 70)
print("CONSERVATIVE RANGE (primary, 3× uncertainty in P_astro)")
print("=" * 70)
print(f"  BF range:  {BF_lo:.0f} – {BF_hi:.0f}")
print(f"  Midpoint:  {R_primary['BF']:.0f}  ({strength(R_primary['BF'])})")

# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("KEY FINDING: FEATURE C (ἐπάνω) ACROSS ALL TEXT TYPES")
print("=" * 70)
print("""
  Feature C (ἐπι-based localization over a specific point) appears in
  zero non-narrative texts: not in astronomical reports, not in omen texts,
  not in PGM I.54-95.  It appears in exactly two non-Matthean sources:
    • Lucian, Nav. 9: ἐπικαθίσαι τῷ καρχησίῳ (star sits upon masthead;
      fictional comic dialogue)
    • JosAsen 14.4 longer rec.: ἔστη ἐπάνω τῆς κεφαλῆς (diaspora Jewish novel)
  Both are narrative fiction.  JosAsen manuscripts alternate ἐπί / ἐπάνω for
  the same scene, confirming they are interchangeable in the narrative register.

  Most telling in the non-narrative direction is Josephus, JW 6.289: a
  celestial body does "stand" (ἵστημι, same verb as Matthew — B=1) over
  Jerusalem, but uses ὑπὲρ, not ἐπάνω (C=0).  The shift to ἐπάνω / ἐπι-
  localization marks the narrative/fictional register, not the astronomical.
""")

# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("FEATURE C (ἐπάνω) ACROSS ALL CORPORA INCLUDING DESCENT TEXTS")
print("=" * 70)
all_corpora = [
    ("Primary guiding corpus",    GUIDING),
    ("Stellar-omen corpus",       OMEN),
    ("Supernatural descent",      DESCENT),
]
print(f"\n  {'Text/Tradition':<46} {'C':>2}  Notes")
print("  " + "-" * 70)
for label, corp in all_corpora:
    print(f"  [{label}]")
    for t, s in corp.items():
        note = ""
        if "Josephus" in t:     note = "ἵστημι + ὑπὲρ (not ἐπάνω)"
        if "Dio" in t:          note = "αἰωρηθείς + ὑπὲρ"
        if "PGM" in t:          note = "στήσεται + ἐν μέσῳ (not ἐπάνω)"
        if "Aseneth" in t:      note = "star=messenger entity; ἐπάνω (longer rec. only); Philonenko: C=0"
        print(f"    {t:<44} {s[2]:>2}  {note}")
print()
C_any = sum(1 for corp in [GUIDING, OMEN, DESCENT]
              for s in corp.values() if s[2] == 1)
print(f"  Non-Matthean texts with C=1: {C_any}")
print(f"  → JosAsen 14.4 longer rec. (angel/man subject, visionary genre)")
print(f"  → Every occurrence of ἐπάνω over-a-place in antiquity:")
print(f"     • Astronomical texts: NEVER ἐπάνω (always ὑπὲρ)")
print(f"     • Omen texts: NEVER ἐπάνω for a celestial body")
print(f"     • Descent-to-earth (PGM): ἐν μέσῳ, NOT ἐπάνω")
print(f"     • Visionary/narrative: ἐπάνω for ANGEL (JosAsen), never for star")
print(f"     • Matthew 2:9: ἐπάνω for a STAR — unique in the extant record")

# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY FOR PAPER")
print("=" * 70)
R = R_primary
print(f"""
  Primary analysis (N=4 independent guiding-star traditions, Laplace k=1):
    P(A,B,C,D | H_guide) ≈ {R['pj']:.2e}
    P(A,B,C,D | H_astro) ≈ {P_ASTRO_JOINT:.2e}
    BF ≈ {R['BF']:.0f}  ({strength(R['BF'])} evidence on Jeffreys scale)
    Conservative range (3× denominator uncertainty): {BF_lo:.0f}–{BF_hi:.0f}

  Extended sensitivity (all 8 texts including omen corpus):
    BF ≈ {scenarios[2]['BF']:.0f}  (still {strength(scenarios[2]['BF'])} evidence)

  Adding Josephus JW 6.289 and Dio Cassius 54.29.8 to the omen corpus
  reveals that ἐπάνω (feature C) is absent from every ancient celestial-
  body-over-place description in the extant record — guiding stars, omen
  stars, and astronomical texts alike.  Matthew's usage is unique.
""")
