#!/usr/bin/env python3
"""
sob_bf_sensitivity.py — robustness of the corpus Bayes factor and the Fisher
combination to the two dependence assumptions embedded in them.

The reported figures rest on a product over four features. Two independence
assumptions sit inside that product, and neither is exactly true.

  rho_w  WITHIN a feature. p = 1/(N+1) is a rule of succession over N trials
         with zero successes, which is exact only if the N occurrences are
         independent opportunities. They are clustered by author: 98 of the
         356 occurrences of erchomai are Vettius Valens. Effective N is
         computed with the Kish design effect for unequal clusters,
         N_eff = N / (1 + (m_bar - 1) rho_w), where m_bar is the
         size-weighted mean cluster size.

         Against this: the absence is replicated across sixteen authors,
         several centuries and two sub-genres. Were it an idiolect one would
         expect between-author heterogeneity, and there is none. That bounds
         rho_w from above and is the reason the tabulated range stops at 0.2.

  rho_b  BETWEEN features. The four features are four properties of one clause
         of one sentence. The manuscript already concedes this and applies a
         Gaussian copula to the Fisher combination; the same concession is
         applied here to the Bayes factor, which is the correction the
         published figure omits. The effective feature count is taken from
         Brown's variance-corrected degrees of freedom, k_eff = df/2, and the
         independence-product BF is downweighted geometrically to
         BF^(k_eff/4). This is a heuristic, and is labelled as one.

The point of the table is not the central value but its range: across every
combination tested the Bayes factor remains decisive on the Jeffreys scale and
the Fisher combination remains significant at 0.05.
"""

import math
from scipy.stats import chi2

# In-corpus occurrences by author, from the re-verified hit tables (2026-09-04)
A = {'Ptolemy': 3, 'VettiusValens': 4, 'Hephaestion': 3, 'Dorotheus': 1}
B = {'Ptolemy': 30, 'VettiusValens': 13, 'Dorotheus': 11, 'Hephaestion': 11,
     'AchillesTatius': 10, 'Cleomedes': 8, 'TheonSmyrna': 6, 'Manetho': 6,
     'Eratosthenes': 6, 'Theodosius': 4, 'PaulAlex': 4, 'Aratus': 4,
     'Hipparchus': 3, 'Geminus': 2, 'Aristotle': 2, 'Theophrastus': 1}
C = {'Ptolemy': 14, 'TheonSmyrna': 4, 'VettiusValens': 6, 'Hephaestion': 6,
     'Theodosius': 4, 'Aristotle': 3, 'Dorotheus': 1, 'AchillesTatius': 1,
     'Manetho': 1}
D = {'VettiusValens': 98, 'Hephaestion': 65, 'Dorotheus': 42, 'Ptolemy': 28,
     'Theodosius': 24, 'Manetho': 19, 'Aratus': 18, 'Aristotle': 15,
     'TheonSmyrna': 14, 'Eratosthenes': 14, 'Autolycus': 6, 'Hipparchus': 5,
     'Cleomedes': 3, 'Theophrastus': 3, 'Geminus': 1, 'AchillesTatius': 1}
FEATURES = [('A', A), ('B', B), ('C', C), ('D', D)]

# Guiding-star side: Laplace-smoothed over the four attested traditions
P_GUIDE = dict(A=5/6, B=3/6, C=3/6, D=4/6)
BF_MATTHEW = 6 * 21 * 8 * 87            # within-Matthew Bayes factor, 87,696

# Brown effective degrees of freedom at each rho_b, from corpus_pvalue.py
BROWN_DF = {0.0: 8.00, 0.2: 7.17, 0.4: 5.45, 0.5: 4.62, 0.7: 3.26, 0.95: 2.16}


def n_effective(counts, rho_w):
    """Kish effective sample size for unequal cluster sizes."""
    n = sum(counts.values())
    m_bar = sum(c * c for c in counts.values()) / n
    return n / (1.0 + (m_bar - 1.0) * rho_w)


def cell(rho_w, rho_b):
    ps = {k: 1.0 / (n_effective(cd, rho_w) + 1.0) for k, cd in FEATURES}
    bf_independent = math.prod(P_GUIDE.values()) / math.prod(ps.values())
    x = -2.0 * sum(math.log(p) for p in ps.values())
    df = BROWN_DF[rho_b]
    k_eff = df / 2.0
    bf_corpus = bf_independent ** (k_eff / 4.0)
    bf_joint = bf_corpus * BF_MATTHEW ** (k_eff / 4.0)
    p_fisher = chi2.sf(x * df / 8.0, df)
    return ps, bf_corpus, bf_joint, p_fisher


def main():
    print("=" * 72)
    print("SENSITIVITY OF THE CORPUS BAYES FACTOR AND FISHER COMBINATION")
    print("=" * 72)
    print("\nEffective sample size per feature:\n")
    print(f"  {'rho_w':>6}" + "".join(f"{k:>12}" for k, _ in FEATURES))
    for rw in (0.0, 0.1, 0.2):
        print(f"  {rw:6.2f}" + "".join(
            f"{n_effective(cd, rw):12.1f}" for _, cd in FEATURES))

    print("\n\nBayes factor and Fisher p across both dependence parameters:\n")
    print(f"  {'rho_w':>6} {'rho_b':>6} | {'BF_corpus':>11} {'BF_joint':>11}"
          f" | {'Fisher p':>10}")
    print("  " + "-" * 54)
    rows = []
    for rw in (0.0, 0.1, 0.2):
        for rb in (0.0, 0.4, 0.5, 0.7):
            _, bc, bj, pf = cell(rw, rb)
            rows.append((rw, rb, bc, bj, pf))
            print(f"  {rw:6.2f} {rb:6.2f} | {bc:11.2e} {bj:11.2e} | {pf:10.2e}")

    lo = min(r[3] for r in rows)
    hi = max(r[3] for r in rows)
    worst_p = max(r[4] for r in rows)
    print("  " + "-" * 54)
    print(f"\n  BF_joint spans {lo:.1e} to {hi:.1e} across the grid")
    print(f"  worst-case Fisher p: {worst_p:.2e}")
    print(f"  every cell decisive on the Jeffreys scale (BF > 100): "
          f"{all(r[3] > 100 for r in rows)}")
    print(f"  every cell significant at alpha = 0.05: "
          f"{all(r[4] < 0.05 for r in rows)}")
    print("\n  Central estimate for reporting: rho_w in [0, 0.2], "
          "rho_b in [0.4, 0.5]")
    central = [r[3] for r in rows if r[0] <= 0.2 and r[1] in (0.4, 0.5)]
    print(f"  -> BF_joint between {min(central):.1e} and {max(central):.1e}")


if __name__ == "__main__":
    main()
