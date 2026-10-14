"""Descriptive Route B direction views; frozen inference and cluster inputs are untouched."""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.pathway_centric import prepare_proteins
from src.pathway_comparison import pathway_key

FIELDS = ['observed_direction', 'null_tail', 'empirical_lower_tail_probability',
          'empirical_upper_tail_probability', 'n_positive_members', 'n_negative_members',
          'n_zero_members', 'positive_fraction', 'negative_fraction', 'dominant_direction_fraction']
KEYS = ['pathway_library', 'pathway']


def describe_pathways(results, proteins, definitions, *, permutations=100000, seed=42):
    """Replay the recorded null only to recover missing tails; never replace p/q or select candidates.

    All coherence counts use unique eligible gene means, not original protein rows.
    Balanced tails have an explicit tie label. Tail choice does not imply significance.
    """
    stats, audit = prepare_proteins(proteins)
    values = stats.to_numpy(dtype=float)
    nulls = {}
    rows = []
    for row in results.itertuples():
        genes = sorted({str(g).strip().upper() for g in definitions[row.pathway_library][row.pathway]} & set(stats.index))
        member_stats = stats.loc[genes]
        n = len(genes)
        if n != int(row.pathway_size) or n == 0:
            raise ValueError('Saved measured pathway size disagrees with the frozen input')
        mean = float(member_stats.mean())
        if not np.isclose(mean, row.mean_t, rtol=1e-12, atol=1e-14):
            raise ValueError('Saved pathway mean disagrees with frozen statistics')
        if n not in nulls:
            rng = np.random.default_rng(np.random.SeedSequence([seed, n]))
            nulls[n] = np.array([values[rng.choice(len(values), n, replace=False)].mean()
                                for _ in range(permutations)])
        null = nulls[n]
        lower = (1 + np.count_nonzero(null <= mean)) / (permutations + 1)
        upper = (1 + np.count_nonzero(null >= mean)) / (permutations + 1)
        reproduced_p = min(1., 2 * min(lower, upper))
        if not np.isclose(reproduced_p, row.p_value, rtol=1e-12, atol=1e-14):
            raise ValueError('Replayed empirical tails do not reproduce the frozen p-value')
        positive, negative, zero = (int((member_stats > 0).sum()), int((member_stats < 0).sum()), int((member_stats == 0).sum()))
        rows.append(dict(pathway_library=row.pathway_library, pathway=row.pathway,
            observed_direction='positive' if mean > 0 else 'negative' if mean < 0 else 'zero',
            null_tail='lower_than_random' if lower < upper else 'higher_than_random' if upper < lower else 'balanced_or_tied',
            empirical_lower_tail_probability=float(lower), empirical_upper_tail_probability=float(upper),
            n_positive_members=positive, n_negative_members=negative, n_zero_members=zero,
            positive_fraction=positive/n, negative_fraction=negative/n,
            dominant_direction_fraction=max(positive,negative)/n))
    annotations = pd.DataFrame(rows, columns=KEYS+FIELDS)
    return attach_fields(results, annotations), annotations, stats, audit


def attach_fields(table, annotations):
    """Validated one-to-one diagnostic join preserving saved row order and all core columns."""
    return table.drop(columns=[c for c in FIELDS if c in table],errors='ignore').merge(annotations, on=KEYS, how='left', validate='one_to_one', sort=False)


def enrich_comparison(table, annotations):
    """Attach B diagnostics by the comparison's existing exact library/ID matching rule."""
    lookup = {}
    for row in annotations.to_dict('records'):
        key = pathway_key(row['pathway_library'], row['pathway'])
        if key in lookup:
            raise ValueError('Ambiguous Route B identity while joining comparison diagnostics')
        lookup[key] = row
    result = table.copy()
    for field in FIELDS:
        result['route_b_'+field] = [lookup.get((r.pathway_library,r.pathway_identifier),{}).get(field,np.nan)
                                   for r in table.itertuples()]
    return result
