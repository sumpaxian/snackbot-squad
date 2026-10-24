"""Dichte in der Anzeigeskala.

Gezählt wird in Klassen der transformierten Achse. Bei logarithmischer oder biexponentieller
Achse sind die Klassen damit automatisch logarithmisch – lineare Klassen auf einer Log-Achse
erzeugen grobe Blöcke im unteren Bereich.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

LO, HI = -0.02, 1.02


def hist2d(dx, dy, bins=256, smooth=1.0, lo=LO, hi=HI):
    ok = np.isfinite(dx) & np.isfinite(dy)
    h, xe, ye = np.histogram2d(dx[ok], dy[ok], bins=bins, range=[[lo, hi], [lo, hi]])
    if smooth:
        h = ndi.gaussian_filter(h, smooth)
    return h, xe, ye


def hist1d(dx, bins=256, smooth=1.0, lo=LO, hi=HI):
    ok = np.isfinite(dx)
    h, e = np.histogram(dx[ok], bins=bins, range=(lo, hi))
    h = h.astype(float)
    if smooth:
        h = ndi.gaussian_filter1d(h, smooth)
    return h, e


def auto_range(d, lo_q=0.1, hi_q=99.9, pad=0.04):
    """Datengetriebene Achsengrenzen in der Anzeigeskala (die Hauptpopulation soll groß erscheinen)."""
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return LO, HI
    a, b = np.percentile(d, [lo_q, hi_q])
    if b <= a:
        a, b = a - 0.05, b + 0.05
    span = b - a
    return float(a - pad * span), float(b + pad * span)


def encode_grid(h, mode="log"):
    """Dichtegitter als 0–255 für die Oberfläche (Zeilen = y von unten nach oben)."""
    g = np.asarray(h, float)
    if mode == "log":
        g = np.log1p(g)
    top = g.max()
    if top <= 0:
        return np.zeros(g.shape, np.uint8)
    q = np.clip(g / top * 255.0, 0, 255)
    q[(h > 0) & (q < 1)] = 1           # schwach besetzte Klassen bleiben sichtbar
    return q.astype(np.uint8)
