"""Achsentransformationen.

Jede Transformation bildet Rohwerte auf eine Anzeigeskala von ungefähr 0 bis 1 ab.
Gates werden in dieser Anzeigeskala ausgewertet: eine gerade Polygonkante auf dem
Bildschirm ist auch eine gerade Kante bei der Zugehörigkeitsprüfung. So verhalten sich
FlowJo und Gating-ML, wenn ein Gate auf transformierten Achsen gezeichnet wurde.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from flowutils import transforms as _fu

KINDS = ("linear", "log", "logicle", "asinh", "arcsinh", "flowjo_biex", "external")


def _col(x):
    return np.ascontiguousarray(np.asarray(x, dtype=np.float64).reshape(-1, 1))


@dataclass
class Transform:
    kind: str = "linear"
    params: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError(f"Unbekannte Transformation: {self.kind}")
        p = dict(self.params)
        if self.kind == "linear":
            p.setdefault("min", 0.0)
            p.setdefault("max", 262144.0)
        elif self.kind == "log":
            p.setdefault("T", 262144.0)
            p.setdefault("M", 4.5)
        elif self.kind == "logicle":
            p.setdefault("T", 262144.0)
            p.setdefault("W", 0.5)
            p.setdefault("M", 4.5)
            p.setdefault("A", 0.0)
        elif self.kind == "asinh":
            p.setdefault("T", 262144.0)
            p.setdefault("M", 4.5)
            p.setdefault("A", 0.0)
        elif self.kind == "arcsinh":
            # asinh(x / Kofaktor), auf asinh(T / Kofaktor) normiert
            p.setdefault("cofactor", 5.0)
            p.setdefault("T", 262144.0)
        elif self.kind == "flowjo_biex":
            p.setdefault("negative", 0.0)
            p.setdefault("width", -10.0)
            p.setdefault("positive", 4.41854)
            p.setdefault("max_value", 262144.000029)
        self.params = p
        self._fk = None

    # FlowJos Biex-Variante weicht numerisch von Logicle ab; FlowKit bildet sie nach.
    @classmethod
    def external(cls, obj, name=""):
        """Beliebiges Objekt mit apply()/inverse() auf Spaltenarrays, z. B. eine FlowKit-Transformation."""
        t = cls("external", {"name": name or type(obj).__name__})
        t._fk = obj
        return t

    def _flowkit(self):
        if self.kind == "external":
            if self._fk is None:
                raise ValueError("Externe Transformation ohne Objekt (nicht speicherbar)")
            return self._fk
        if self._fk is None:
            import flowkit as fk
            p = self.params
            self._fk = fk.transforms.WSPBiexTransform(negative=p["negative"], width=p["width"],
                                        positive=p["positive"], max_value=p["max_value"])
        return self._fk

    def forward(self, x):
        x = np.asarray(x, dtype=np.float64)
        shape = x.shape
        p = self.params
        if self.kind == "linear":
            y = (x - p["min"]) / (p["max"] - p["min"])
        elif self.kind == "log":
            with np.errstate(divide="ignore", invalid="ignore"):
                y = (np.log10(np.clip(x, 1e-300, None)) - np.log10(p["T"])) / p["M"] + 1.0
            y = np.where(x > 0, y, -np.inf)
        elif self.kind == "logicle":
            y = _fu.logicle(_col(x), [0], t=p["T"], m=p["M"], w=p["W"], a=p["A"])[:, 0]
        elif self.kind == "asinh":
            y = _fu.asinh(_col(x), [0], t=p["T"], m=p["M"], a=p["A"])[:, 0]
        elif self.kind == "arcsinh":
            c = p["cofactor"]
            y = np.arcsinh(x / c) / np.arcsinh(p["T"] / c)
        else:
            y = self._flowkit().apply(_col(x))[:, 0]
        return np.asarray(y, dtype=np.float64).reshape(shape)

    def inverse(self, y):
        y = np.asarray(y, dtype=np.float64)
        shape = y.shape
        p = self.params
        if self.kind == "linear":
            x = y * (p["max"] - p["min"]) + p["min"]
        elif self.kind == "log":
            x = 10 ** ((y - 1.0) * p["M"] + np.log10(p["T"]))
        elif self.kind == "logicle":
            x = _fu.logicle_inverse(_col(y), [0], t=p["T"], m=p["M"], w=p["W"], a=p["A"])[:, 0]
        elif self.kind == "asinh":
            x = _fu.asinh_inverse(_col(y), [0], t=p["T"], m=p["M"], a=p["A"])[:, 0]
        elif self.kind == "arcsinh":
            c = p["cofactor"]
            x = np.sinh(y * np.arcsinh(p["T"] / c)) * c
        else:
            x = self._flowkit().inverse(_col(y))[:, 0]
        return np.asarray(x, dtype=np.float64).reshape(shape)

    def to_dict(self):
        return {"kind": self.kind, "params": dict(self.params)}

    @classmethod
    def from_dict(cls, d):
        if d is None:
            return cls()
        return cls(d["kind"], dict(d.get("params", {})))

    def ticks(self, lo=None, hi=None):
        """Beschriftete Achsenstriche im sichtbaren Bereich als (Rohwert, Anzeigewert, Text)."""
        lo = -0.02 if lo is None else lo
        hi = 1.02 if hi is None else hi
        r0, r1 = sorted(float(v) for v in self.inverse(np.array([lo, hi])))
        if self.kind == "linear" or not np.isfinite(r0 * r1):
            span = r1 - r0
            if span <= 0:
                return []
            step = 10 ** np.floor(np.log10(span / 5))
            for m in (1, 2, 5, 10):
                if span / (step * m) <= 6:
                    step *= m
                    break
            raw = np.arange(np.ceil(r0 / step) * step, r1 + step * 1e-9, step) + 0.0
            raw[np.abs(raw) < step * 1e-9] = 0.0
            return [(float(r), float(self.forward(r)), _short(r, step)) for r in raw]
        top = max(abs(r0), abs(r1), 10.0)
        decades = range(0, int(np.floor(np.log10(top))) + 1)
        if self.kind == "log":
            cand = [m * 10.0 ** e for e in range(-3, decades.stop) for m in (1, 2, 5)]
            prio = sorted(cand, key=lambda v: (0 if np.isclose(np.log10(v) % 1, 0) else 1, -v))
        else:
            # Vorrang: 0, dann große positive Dekaden, dann negative – dichte Striche fallen weg
            prio = [0.0] + [10.0 ** e for e in reversed(decades)] + [-(10.0 ** e) for e in reversed(decades)]
        min_gap = 0.07 * (hi - lo)
        chosen = []
        for v in prio:
            if not (r0 <= v <= r1) or (self.kind == "log" and v <= 0):
                continue
            pos = float(self.forward(v))
            if np.isfinite(pos) and all(abs(pos - c[1]) >= min_gap for c in chosen):
                chosen.append((v, pos, _pow10(v)))
        if self.kind == "log":
            pure = [c for c in chosen if np.isclose(np.log10(c[0]) % 1, 0) or np.isclose(np.log10(c[0]) % 1, 1)]
            if len(pure) >= 3:
                chosen = pure
        return sorted(chosen, key=lambda t: t[1])


def _short(v, step=1.0):
    v = float(v)
    if abs(v) >= 1000 and step >= 1000:
        return f"{v / 1000:.0f}K"
    if step < 1:
        return f"{v:.{int(np.ceil(-np.log10(step)))}f}".replace(".", ",")
    return f"{v:.0f}"


def _pow10(v):
    if v == 0:
        return "0"
    e = np.log10(abs(v))
    if abs(e - round(e)) > 1e-9:          # 2·10³, 5·10³ auf Log-Achsen
        return f"{v:g}"
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    e = int(round(np.log10(abs(v))))
    s = "−" if v < 0 else ""
    return f"{s}10{str(e).translate(sup)}"


SCATTER_PREFIXES = ("FSC", "SSC", "FS ", "SS ", "FS-", "SS-")
TIME_NAMES = ("time", "hdr-t")


def is_scatter(name: str) -> bool:
    n = name.upper()
    return n.startswith(SCATTER_PREFIXES) or n in ("FS", "SS")


def is_time(name: str) -> bool:
    return name.strip().lower() in TIME_NAMES


def estimate_logicle_w(values, T=262144.0, M=4.5):
    """Breite des linearen Bereichs nach Parks et al. 2006 aus dem 5. Perzentil der negativen Werte."""
    v = np.asarray(values, dtype=np.float64)
    neg = v[v < 0]
    if neg.size < 50:
        return 0.5
    r = abs(np.percentile(neg, 5))
    if r <= 0:
        return 0.5
    w = (M - np.log10(T / r)) / 2.0
    return float(np.clip(w, 0.0, 2.0))


def estimate_cofactor(values):
    """Kofaktor für asinh(x/c): robuste SD der Rauschverteilung um null.

    Die untere Hälfte wird am Median gespiegelt, davon MAD·1,4826. Folgt der Idee aus der
    Triple-C-Auswertung (01_kofaktoren.py); deren genaue Formel ist hier nicht überprüft.
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size < 50:
        return 5.0
    med = np.median(v)
    c = 1.4826 * np.median(med - v[v <= med])
    return float(c) if c > 0 else 5.0


def default_transform(name: str, max_range: float = 262144.0, values=None) -> Transform:
    """Streulicht und Zeit linear, Fluoreszenz Logicle mit geschätzter Breite."""
    if is_time(name):
        lo = 0.0
        hi = float(np.nanmax(values)) if values is not None and len(values) else max_range
        return Transform("linear", {"min": lo, "max": hi or 1.0})
    if is_scatter(name):
        return Transform("linear", {"min": 0.0, "max": float(max_range)})
    w = estimate_logicle_w(values, max_range) if values is not None else 0.5
    return Transform("logicle", {"T": float(max_range), "W": w, "M": 4.5, "A": 0.0})
