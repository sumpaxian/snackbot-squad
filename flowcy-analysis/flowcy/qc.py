"""Qualitätskontrolle je Probe: Flussrate über die Zeit, Randereignisse, negatives Streulicht.

Die Fluss-QC folgt dem Gedanken von flowAI/flowClean in einfacher, nachvollziehbarer Form:
Ereignisse werden in Zeitfenster geteilt, die Rate je Fenster mit Median und MAD robust
standardisiert, auffällige Fenster markiert und zu Abschnitten zusammengefasst.
"""
from __future__ import annotations

import numpy as np

from .transforms import is_scatter

MAD_K = 1.4826


def flow_rate(sample, bins: int = 100, z_limit: float = 3.5, rel_dev: float = 0.25):
    tch = sample.time_channel
    if tch is None:
        return {"available": False, "reason": "kein Zeitkanal"}
    t = sample.raw[:, sample.index(tch)]
    if not np.all(np.isfinite(t)) or np.ptp(t) <= 0:
        return {"available": False, "reason": "Zeitkanal ohne Spannweite"}
    edges = np.linspace(t.min(), t.max(), bins + 1)
    counts, _ = np.histogram(t, edges)
    width = np.diff(edges)
    rate = counts / width
    med = np.median(rate)
    mad = np.median(np.abs(rate - med)) * MAD_K
    z = (rate - med) / mad if mad > 0 else np.zeros_like(rate)
    # Einbrüche (Verstopfung, Luft) und Spitzen; dazu eine Mindestabweichung vom Median,
    # damit sehr gleichmäßige Läufe nicht wegen winziger MAD auffallen
    bad = (np.abs(z) > z_limit) & (np.abs(rate - med) > rel_dev * med)
    segments = []
    i = 0
    while i < bins:
        if bad[i]:
            j = i
            while j + 1 < bins and bad[j + 1]:
                j += 1
            segments.append({"start": float(edges[i]), "end": float(edges[j + 1]),
                             "events": int(counts[i:j + 1].sum())})
            i = j + 1
        else:
            i += 1
    affected = sum(s["events"] for s in segments)
    span = float(t.max() - t.min())
    time_frac = sum(s["end"] - s["start"] for s in segments) / span if span else 0.0
    excluded = np.zeros(sample.n_events, bool)
    for s in segments:
        excluded |= (t >= s["start"]) & (t <= s["end"])
    return {
        "available": True, "channel": tch,
        "time": ((edges[:-1] + edges[1:]) / 2).tolist(), "rate": rate.tolist(), "z": z.tolist(),
        "median_rate": float(med), "segments": segments,
        "affected_fraction": affected / max(sample.n_events, 1),
        "affected_time_fraction": time_frac,
        "excluded_mask": excluded,
    }


def margin_events(sample):
    """Anteil der Ereignisse am oberen Messbereichsende je Streulichtkanal (Sättigung)."""
    out = {}
    for c in sample.channels:
        if not is_scatter(c):
            continue
        v = sample.raw[:, sample.index(c)]
        top = sample.max_range(c)
        out[c] = float(np.mean(v >= top * 0.999)) if top else 0.0
    return out


def negative_scatter(sample):
    out = {}
    for c in sample.channels:
        if is_scatter(c):
            out[c] = float(np.mean(sample.raw[:, sample.index(c)] <= 0))
    return out


def run_qc(sample, margin_limit=0.02, flow_limit=0.05):
    fr = flow_rate(sample)
    mg = margin_events(sample)
    ng = negative_scatter(sample)
    issues = []
    if fr.get("available") and max(fr["affected_fraction"], fr["affected_time_fraction"]) > flow_limit:
        segs = ", ".join(f"{s['start']:.4g}–{s['end']:.4g}" for s in fr["segments"][:3])
        issues.append(f"Flussrate instabil ({fr['affected_time_fraction']:.0%} der Messzeit, "
                      f"{fr['affected_fraction']:.1%} der Ereignisse; Zeit {segs})")
    for c, f in mg.items():
        if f > margin_limit:
            issues.append(f"Sättigung in {c}: {f:.1%} am Bereichsende")
    if not fr.get("available"):
        issues.append(f"Fluss-QC nicht möglich: {fr.get('reason')}")
    real = [i for i in issues if "nicht möglich" not in i]
    status = "warn" if real else ("info" if issues else "ok")
    fr_public = {k: v for k, v in fr.items() if k != "excluded_mask"}
    return {"sample": sample.name, "events": sample.n_events, "status": status, "issues": issues,
            "flow": fr_public, "margin": mg, "negative_scatter": ng,
            "_excluded_mask": fr.get("excluded_mask")}
