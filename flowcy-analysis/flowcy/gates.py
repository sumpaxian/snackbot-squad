"""Gate-Typen.

Koordinaten werden in Rohwerten gespeichert (lesbar, austauschbar), die Zugehörigkeit
wird in der Anzeigeskala der jeweiligen Achsentransformation geprüft.
"""
from __future__ import annotations

import uuid
import numpy as np

from .transforms import Transform


def _new_id():
    return uuid.uuid4().hex[:10]


class Gate:
    kind = "gate"

    def __init__(self, x, y=None, xform=None, yform=None, id=None):
        self.x = x
        self.y = y
        self.xform = xform if isinstance(xform, Transform) else Transform.from_dict(xform)
        self.yform = (yform if isinstance(yform, Transform) else Transform.from_dict(yform)) if y else None
        self.id = id or _new_id()

    @property
    def dims(self):
        return [self.x] if self.y is None else [self.x, self.y]

    def display_coords(self, sample, compensated=True, rows=None):
        ev = sample.events(compensated)
        xs = ev[:, sample.index(self.x)]
        if rows is not None:
            xs = xs[rows]
        dx = self.xform.forward(xs)
        dy = None
        if self.y is not None:
            ys = ev[:, sample.index(self.y)]
            if rows is not None:
                ys = ys[rows]
            dy = self.yform.forward(ys)
        return dx, dy

    def evaluate(self, sample, compensated=True, rows=None) -> np.ndarray:
        dx, dy = self.display_coords(sample, compensated, rows)
        return self.contains(dx, dy)

    def contains(self, dx, dy):  # pragma: no cover - abstrakt
        raise NotImplementedError

    def base_dict(self):
        return {"kind": self.kind, "id": self.id, "x": self.x, "y": self.y,
                "xform": self.xform.to_dict(), "yform": self.yform.to_dict() if self.yform else None}


def points_in_polygon(px, py, vx, vy):
    """Gerade-Ungerade-Regel, vektorisiert über alle Punkte."""
    px = np.asarray(px, float)
    py = np.asarray(py, float)
    inside = np.zeros(px.shape, dtype=bool)
    n = len(vx)
    j = n - 1
    for i in range(n):
        xi, yi, xj, yj = vx[i], vy[i], vx[j], vy[j]
        cond = (yi > py) != (yj > py)
        with np.errstate(divide="ignore", invalid="ignore"):
            xcross = (xj - xi) * (py - yi) / (yj - yi) + xi
        inside ^= cond & (px < xcross)
        j = i
    return inside


class PolygonGate(Gate):
    kind = "polygon"

    def __init__(self, x, y, vertices, **kw):
        super().__init__(x, y, **kw)
        self.vertices = [(float(a), float(b)) for a, b in vertices]
        if len(self.vertices) < 3:
            raise ValueError("Ein Polygon braucht mindestens drei Eckpunkte")

    def display_vertices(self):
        v = np.array(self.vertices)
        return self.xform.forward(v[:, 0]), self.yform.forward(v[:, 1])

    def contains(self, dx, dy):
        vx, vy = self.display_vertices()
        ok = np.isfinite(dx) & np.isfinite(dy)
        out = np.zeros(dx.shape, bool)
        out[ok] = points_in_polygon(dx[ok], dy[ok], vx, vy)
        return out

    def to_dict(self):
        d = self.base_dict()
        d["vertices"] = [list(v) for v in self.vertices]
        return d


class RectangleGate(Gate):
    """Rechteck, Bereich oder Schwelle. Fehlende Grenzen (None) sind offen."""
    kind = "rectangle"

    def __init__(self, x, y=None, xmin=None, xmax=None, ymin=None, ymax=None, **kw):
        super().__init__(x, y, **kw)
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax

    def contains(self, dx, dy):
        out = np.isfinite(dx)
        if self.xmin is not None:
            out &= dx >= self.xform.forward(self.xmin)
        if self.xmax is not None:
            out &= dx < self.xform.forward(self.xmax)
        if self.y is not None:
            out &= np.isfinite(dy)
            if self.ymin is not None:
                out &= dy >= self.yform.forward(self.ymin)
            if self.ymax is not None:
                out &= dy < self.yform.forward(self.ymax)
        return out

    def to_dict(self):
        d = self.base_dict()
        d.update(xmin=self.xmin, xmax=self.xmax, ymin=self.ymin, ymax=self.ymax)
        return d


class EllipseGate(Gate):
    """Ellipse in der Anzeigeskala: Mittelpunkt in Rohwerten, Halbachsen in Anzeigeeinheiten."""
    kind = "ellipse"

    def __init__(self, x, y, center, rx, ry, angle=0.0, **kw):
        super().__init__(x, y, **kw)
        self.center = (float(center[0]), float(center[1]))
        self.rx, self.ry, self.angle = float(rx), float(ry), float(angle)

    def contains(self, dx, dy):
        cx = self.xform.forward(self.center[0])
        cy = self.yform.forward(self.center[1])
        a = np.deg2rad(self.angle)
        u = (dx - cx) * np.cos(a) + (dy - cy) * np.sin(a)
        v = -(dx - cx) * np.sin(a) + (dy - cy) * np.cos(a)
        with np.errstate(invalid="ignore"):
            return np.nan_to_num((u / self.rx) ** 2 + (v / self.ry) ** 2, nan=np.inf) <= 1.0

    def to_dict(self):
        d = self.base_dict()
        d.update(center=list(self.center), rx=self.rx, ry=self.ry, angle=self.angle)
        return d


QUADRANTS = {"+-": "x⁺ y⁻", "++": "x⁺ y⁺", "-+": "x⁻ y⁺", "--": "x⁻ y⁻"}


class QuadrantGate(Gate):
    """Ein Quadrant. Vier Knoten mit gleicher Gruppe teilen sich die Teilungslinien."""
    kind = "quadrant"

    def __init__(self, x, y, xdiv, ydiv, quadrant="++", group=None, **kw):
        super().__init__(x, y, **kw)
        if quadrant not in QUADRANTS:
            raise ValueError(f"Quadrant muss einer von {list(QUADRANTS)} sein")
        self.xdiv, self.ydiv, self.quadrant = float(xdiv), float(ydiv), quadrant
        self.group = group or _new_id()

    def contains(self, dx, dy):
        xd = self.xform.forward(self.xdiv)
        yd = self.yform.forward(self.ydiv)
        ok = np.isfinite(dx) & np.isfinite(dy)
        xs = dx >= xd if self.quadrant[0] == "+" else dx < xd
        ys = dy >= yd if self.quadrant[1] == "+" else dy < yd
        return ok & xs & ys

    def to_dict(self):
        d = self.base_dict()
        d.update(xdiv=self.xdiv, ydiv=self.ydiv, quadrant=self.quadrant, group=self.group)
        return d


def quadrant_set(x, y, xdiv, ydiv, xform=None, yform=None, names=None):
    """Die vier Quadranten einer Teilung, gemeinsame Gruppe."""
    group = _new_id()
    order = ["-+", "++", "--", "+-"]
    gates = [QuadrantGate(x, y, xdiv, ydiv, q, group=group, xform=xform, yform=yform) for q in order]
    names = names or [f"{x}{q[0]} {y}{q[1]}" for q in order]
    return list(zip(names, gates))


class BooleanGate(Gate):
    """Verknüpfung anderer Knoten: and / or / not. Ausgewertet in der Hierarchie (tree.masks)."""
    kind = "boolean"

    def __init__(self, op, refs, complements=None, id=None, **_):
        if op not in ("and", "or", "not"):
            raise ValueError("op muss and, or oder not sein")
        self.op = op
        self.refs = list(refs)
        self.complements = list(complements or [False] * len(self.refs))
        self.x = self.y = None
        self.xform = self.yform = None
        self.id = id or _new_id()

    @property
    def dims(self):
        return []

    def combine(self, masks):
        ms = [~masks[r] if c else masks[r] for r, c in zip(self.refs, self.complements)]
        if self.op == "and":
            return np.logical_and.reduce(ms)
        if self.op == "or":
            return np.logical_or.reduce(ms)
        return ~ms[0]

    def evaluate(self, sample, compensated=True, rows=None):
        raise RuntimeError("Boolesche Gates werden in der Hierarchie ausgewertet")

    def to_dict(self):
        return {"kind": "boolean", "id": self.id, "op": self.op, "refs": self.refs,
                "complements": self.complements}


_KINDS = {"polygon": PolygonGate, "boolean": BooleanGate, "rectangle": RectangleGate, "ellipse": EllipseGate, "quadrant": QuadrantGate}


def gate_from_dict(d: dict) -> Gate:
    d = dict(d)
    kind = d.pop("kind")
    cls = _KINDS[kind]
    if kind == "boolean":
        return cls(d["op"], d["refs"], d.get("complements"), id=d.get("id"))
    common = dict(x=d.pop("x"), y=d.pop("y", None), xform=d.pop("xform", None),
                  yform=d.pop("yform", None), id=d.pop("id", None))
    if kind == "polygon":
        return cls(common.pop("x"), common.pop("y"), d["vertices"], **common)
    if kind == "rectangle":
        return cls(common.pop("x"), common.pop("y"), d.get("xmin"), d.get("xmax"),
                   d.get("ymin"), d.get("ymax"), **common)
    if kind == "ellipse":
        return cls(common.pop("x"), common.pop("y"), d["center"], d["rx"], d["ry"], d.get("angle", 0.0), **common)
    return cls(common.pop("x"), common.pop("y"), d["xdiv"], d["ydiv"], d.get("quadrant", "++"),
               group=d.get("group"), **common)
