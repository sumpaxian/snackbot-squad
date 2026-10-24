"""Gruppenvergleiche für Populationsanteile.

Nichtparametrisch, wie in der bisherigen Auswertung: Mann-Whitney-U (exakt bei kleinen
Gruppen), Kruskal-Wallis für mehr als zwei Gruppen, paarweise gegen eine Referenzgruppe
mit wählbarer Korrektur. Dazu die Zähl-Unsicherheit seltener Populationen (Poisson,
CV = 100/√n; Cossarizza et al., Eur J Immunol 2021).
"""
from __future__ import annotations

import numpy as np
from scipy import stats as st


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * p[i])
        adj[i] = min(run, 1.0)
    return adj


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)[::-1]
    adj = np.empty(m)
    run = 1.0
    for k, i in enumerate(order):
        rank = m - k
        run = min(run, p[i] * m / rank)
        adj[i] = run
    return adj


def correct(p, method):
    if method in (None, "", "none"):
        return np.asarray(p, float)
    return holm(p) if method == "holm" else bh(p)


def hodges_lehmann(a, b):
    """Verschiebungsschätzer (Median aller paarweisen Differenzen b − a)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.median(b[:, None] - a[None, :])) if a.size and b.size else float("nan")


def count_cv(n):
    """Zähl-Variationskoeffizient in % (Poisson)."""
    return float("inf") if n <= 0 else 100.0 / np.sqrt(n)


def compare_groups(values: dict, reference=None, correction="holm"):
    """values: {Gruppe: [Werte je Probe]}. Gibt Kennzahlen je Gruppe und Tests zurück."""
    groups = [g for g in values if len(values[g])]
    out = {"groups": [], "tests": [], "global": None}
    for g in groups:
        v = np.asarray(values[g], float)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        out["groups"].append({"group": g, "n": int(v.size), "median": float(med), "q1": float(q1), "q3": float(q3),
                              "mean": float(v.mean()), "values": v.tolist()})
    if len(groups) < 2:
        return out
    if len(groups) > 2:
        h, p = st.kruskal(*[values[g] for g in groups])
        out["global"] = {"test": "Kruskal-Wallis", "statistic": float(h), "p": float(p)}
    ref = reference if reference in groups else groups[0]
    pairs = [(ref, g) for g in groups if g != ref]
    raw = []
    for a, b in pairs:
        va, vb = np.asarray(values[a], float), np.asarray(values[b], float)
        method = "exact" if va.size + vb.size <= 20 else "asymptotic"
        u, p = st.mannwhitneyu(va, vb, alternative="two-sided", method=method)
        raw.append(p)
        out["tests"].append({"a": a, "b": b, "test": f"Mann-Whitney-U ({'exakt' if method == 'exact' else 'asymptotisch'})",
                             "U": float(u), "p": float(p), "diff_median": float(np.median(vb) - np.median(va)),
                             "hodges_lehmann": hodges_lehmann(va, vb), "n_a": int(va.size), "n_b": int(vb.size),
                             "min_p_possible": _min_exact_p(va.size, vb.size)})
    adj = correct(raw, correction)
    for t, pa in zip(out["tests"], adj):
        t["p_adj"] = float(pa)
        t["correction"] = correction or "keine"
    return out


def _min_exact_p(n1, n2):
    from math import comb
    return 2.0 / comb(n1 + n2, n1) if n1 and n2 else float("nan")
