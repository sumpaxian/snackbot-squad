"""Gate-Vorschläge. Sie sind Vorschläge zum Annehmen oder Überzeichnen, keine gesetzten Gates.

- Schwelle aus dem Dichtetal zwischen zwei Gipfeln (z. B. Lebend/Tot)
- Einzelzell-Band um die FSC-H/FSC-A-Diagonale (Median ± k robuste SD)
- Hauptpopulation als Dichtekontur (z. B. Zellen gegen Debris)
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.signal import find_peaks

from .density import auto_range, hist1d, hist2d
from .gates import PolygonGate, RectangleGate

MAD_K = 1.4826


def _rows(sample, rows):
    return np.arange(sample.n_events) if rows is None else np.asarray(rows)


def valley_threshold(sample, channel, transform, rows=None, bins=256, smooth=3.0, keep="below"):
    """Tiefster Punkt zwischen den zwei stärksten Gipfeln. keep='below' behält die Ereignisse darunter."""
    r = _rows(sample, rows)
    x = sample.events()[r, sample.index(channel)]
    d = transform.forward(x)
    lo, hi = auto_range(d)
    h, e = hist1d(d, bins=bins, smooth=smooth, lo=lo, hi=hi)
    centers = (e[:-1] + e[1:]) / 2
    peaks, props = find_peaks(h, prominence=h.max() * 0.02)
    info = {"modes": int(len(peaks))}
    if len(peaks) < 2:
        cut_d = float(np.percentile(d[np.isfinite(d)], 95))
        info["fallback"] = "nur ein Gipfel – 95. Perzentil"
    else:
        top2 = sorted(peaks[np.argsort(props["prominences"])[-2:]])
        seg = h[top2[0]:top2[1] + 1]
        cut_d = float(centers[top2[0] + int(np.argmin(seg))])
    cut = float(transform.inverse(cut_d))
    gate = RectangleGate(channel, None, xmax=cut, xform=transform) if keep == "below" \
        else RectangleGate(channel, None, xmin=cut, xform=transform)
    frac = float(np.mean(gate.contains(d, None)))
    info.update(cut=cut, cut_display=cut_d, fraction=frac, hist=h.tolist(), edges=e.tolist())
    return gate, info


def singlet_band(sample, area="FSC-A", height="FSC-H", rows=None, k=3.0, xform=None, yform=None):
    """Band um die Diagonale: Verhältnis Höhe/Fläche im Bereich Median ± k·robuste SD."""
    from .transforms import default_transform
    r = _rows(sample, rows)
    ev = sample.events()
    a = ev[r, sample.index(area)]
    hgt = ev[r, sample.index(height)]
    ok = (a > 0) & (hgt > 0)
    ratio = hgt[ok] / a[ok]
    med = float(np.median(ratio))
    rsd = float(np.median(np.abs(ratio - med)) * MAD_K)
    lo, hi = med - k * rsd, med + k * rsd
    x0, x1 = np.percentile(a[ok], [0.2, 99.8])
    x1 = max(x1, sample.max_range(area))
    verts = [(x0, lo * x0), (x1, lo * x1), (x1, hi * x1), (x0, hi * x0)]
    xform = xform or default_transform(area, sample.max_range(area))
    yform = yform or default_transform(height, sample.max_range(height))
    gate = PolygonGate(area, height, verts, xform=xform, yform=yform)
    frac = float(np.mean((ratio >= lo) & (ratio <= hi)))
    return gate, {"median_ratio": med, "robust_sd": rsd, "lower": lo, "upper": hi, "fraction": frac}


def simplify_ring(pts, eps):
    """Douglas-Peucker für einen geschlossenen Umriss.

    Bei identischem Anfangs- und Endpunkt ist die Bezugsgerade entartet; der Ring wird
    deshalb an den zwei am weitesten auseinanderliegenden Punkten in zwei Bögen geteilt.
    """
    P = np.asarray(pts, float)
    if np.allclose(P[0], P[-1]):
        P = P[:-1]
    if len(P) <= 4:
        return P
    i0 = int(np.argmax(((P - P.mean(axis=0)) ** 2).sum(axis=1)))
    P = np.roll(P, -i0, axis=0)
    i1 = int(np.argmax(((P - P[0]) ** 2).sum(axis=1)))

    def dp(seg):
        if len(seg) < 3:
            return seg
        a, b = seg[0], seg[-1]
        ab = b - a
        n = np.hypot(*ab)
        if n == 0:
            d = np.hypot(*(seg - a).T)
        else:
            d = np.abs(ab[0] * (seg[:, 1] - a[1]) - ab[1] * (seg[:, 0] - a[0])) / n
        i = int(np.argmax(d))
        if d[i] > eps:
            return np.vstack([dp(seg[:i + 1])[:-1], dp(seg[i:])])
        return np.vstack([a, b])

    arc1 = dp(np.vstack([P[:i1 + 1]]))
    arc2 = dp(np.vstack([P[i1:], P[:1]]))
    return np.vstack([arc1[:-1], arc2[:-1]])


def main_population(sample, x, y, xform, yform, rows=None, coverage=0.9, bins=200, smooth=2.0,
                    eps=0.008, seed=None, min_x=0.1):
    """Kontur einer Dichtepopulation, die `coverage` der Ereignisse ihres Einzugsgebiets umfasst.

    Die Dichte wird per Wasserscheide in Einzugsgebiete der Gipfel geteilt, damit eine
    Population nicht in benachbarte (z. B. Debris) ausläuft. Gewählt wird das Gebiet um
    `seed` (Anzeigekoordinaten) oder sonst das massereichste Gebiet, dessen Gipfel nicht im
    untersten Bereich der x-Achse liegt (dort sitzt in der Regel Debris).
    """
    import contourpy
    from .gates import points_in_polygon
    r = _rows(sample, rows)
    ev = sample.events()
    dx = xform.forward(ev[r, sample.index(x)])
    dy = yform.forward(ev[r, sample.index(y)])
    h, xe, ye = hist2d(dx, dy, bins=bins, smooth=smooth)
    xc = (xe[:-1] + xe[1:]) / 2
    yc = (ye[:-1] + ye[1:]) / 2
    # Gipfel auf stärker geglätteter Dichte, damit ein verrauschtes Plateau ein Gipfel bleibt
    hs = ndi.gaussian_filter(h, 3)
    pk_mask = (hs == ndi.maximum_filter(hs, size=15)) & (hs > 0.03 * hs.max())
    peaks = [tuple(p) for p in np.argwhere(pk_mask)]
    levels = np.unique(np.quantile(h[h > 0], np.linspace(1.0, 0.0, 240)))[::-1]

    def grow(pk):
        """Region um den Gipfel bei sinkender Schwelle, bis sie einen anderen Gipfel erreicht."""
        hist = []
        for lev in levels:
            if h[pk] < lev:
                continue
            lab, _ = ndi.label(h >= lev)
            comp = lab == lab[pk]
            if any(comp[q] for q in peaks if q != pk):
                break
            hist.append((lev, float(h[comp].sum()), comp))
        return hist

    if seed is not None:
        i = int(np.clip(np.searchsorted(xc, seed[0]), 0, bins - 1))
        j = int(np.clip(np.searchsorted(yc, seed[1]), 0, bins - 1))
        pk = min(peaks, key=lambda q: (q[0] - i) ** 2 + (q[1] - j) ** 2)
        hist = grow(pk)
    else:
        grown = [(q, grow(q)) for q in peaks]
        grown = [(q, g) for q, g in grown if g]
        cand = [(q, g) for q, g in grown if xc[q[0]] >= min_x] or grown
        pk, hist = max(cand, key=lambda t: t[1][-1][1])
    final_mass = hist[-1][1]
    level, _, comp = next(t for t in hist if t[1] >= coverage * final_mass)
    region = ndi.binary_fill_holes(comp).astype(float)
    n = len(peaks)
    region = np.pad(region, 1)
    xcp = np.concatenate([[xc[0] - (xc[1] - xc[0])], xc, [xc[-1] + (xc[1] - xc[0])]])
    ycp = np.concatenate([[yc[0] - (yc[1] - yc[0])], yc, [yc[-1] + (yc[1] - yc[0])]])
    lines = contourpy.contour_generator(xcp, ycp, region.T).lines(0.5)
    px, py = xc[pk[0]], yc[pk[1]]
    best = None
    for ln in lines:
        if len(ln) >= 4 and points_in_polygon(np.array([px]), np.array([py]), ln[:, 0], ln[:, 1])[0]:
            area = 0.5 * abs(np.dot(ln[:, 0], np.roll(ln[:, 1], 1)) - np.dot(ln[:, 1], np.roll(ln[:, 0], 1)))
            if best is None or area > best[0]:
                best = (area, ln)
    if best is None:
        raise ValueError("Keine geschlossene Kontur um den Dichtegipfel gefunden")
    ring = simplify_ring(best[1], eps)
    verts = np.column_stack([xform.inverse(ring[:, 0]), yform.inverse(ring[:, 1])])
    gate = PolygonGate(x, y, verts.tolist(), xform=xform, yform=yform)
    frac = float(np.mean(gate.contains(dx, dy)))
    return gate, {"level": level, "fraction": frac, "vertices": len(verts), "populations": n}
