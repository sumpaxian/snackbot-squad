"""Gates aus einfachen Beschreibungen übernehmen, z. B. aus Bildern ausgelesene Polygone.

Format (Liste, Reihenfolge = Hierarchie von oben nach unten):

    [
      {"name": "Zellen", "parent": null, "x": "FSC-A", "y": "SSC-A",
       "scale": "log10", "vertices": [[143.8, 14.1], [125.0, 16.2], ...]},
      {"name": "Lebend", "parent": "Einzelzellen", "x": "FL5-A",
       "scale": "arcsinh", "cofactor": 8.45, "max": 184.8}
    ]

scale: linear | log10 | logicle | arcsinh. Ein Eintrag mit vertices wird ein Polygon, mit
min/max (und optional ymin/ymax) ein Rechteck bzw. eine Schwelle.
"""
from __future__ import annotations

from .gates import PolygonGate, RectangleGate
from .transforms import Transform
from .tree import ROOT, GatingTree


def scale_transform(scale: str, spec: dict, axis: str = "x") -> Transform:
    scale = (scale or "linear").lower()
    T = float(spec.get("T", 262144.0))
    if scale in ("log", "log10"):
        # affine Abbildung von log10 – die Zugehörigkeit ist identisch mit der Prüfung in log10
        return Transform("log", {"T": T, "M": float(spec.get("decades", 6))})
    if scale == "logicle":
        return Transform("logicle", {"T": T, "W": float(spec.get("W", 0.5)), "M": 4.5, "A": 0.0})
    if scale in ("arcsinh", "asinh"):
        key = "cofactor" if axis == "x" else "ycofactor"
        return Transform("arcsinh", {"cofactor": float(spec.get(key, spec.get("cofactor", 5.0))), "T": T})
    return Transform("linear", {"min": 0.0, "max": T})


def tree_from_spec(spec: list, name: str = "Importierte Gates") -> GatingTree:
    tree = GatingTree(name)
    ids = {None: ROOT}
    for e in spec:
        parent = ids.get(e.get("parent"))
        if parent is None:
            raise KeyError(f"Elterngate „{e.get('parent')}“ für „{e['name']}“ fehlt (Reihenfolge prüfen)")
        xs = scale_transform(e.get("scale", e.get("xscale")), e, "x")
        ys = scale_transform(e.get("yscale", e.get("scale")), e, "y") if e.get("y") else None
        if "vertices" in e:
            g = PolygonGate(e["x"], e["y"], e["vertices"], xform=xs, yform=ys)
        else:
            g = RectangleGate(e["x"], e.get("y"), e.get("min"), e.get("max"), e.get("ymin"), e.get("ymax"),
                              xform=xs, yform=ys)
        ids[e["name"]] = tree.add(e["name"], g, parent, note=e.get("note", ""))
    return tree
