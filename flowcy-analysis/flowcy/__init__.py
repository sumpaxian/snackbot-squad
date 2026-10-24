"""FlowCy Analysis – Auswertung durchflusszytometrischer Daten.

Die Engine ist ohne Oberfläche nutzbar:

    from flowcy import load_sample, GatingTree
    s = load_sample("probe.fcs")
    tree = GatingTree.load("strategie.flowcy.json")
    stats = tree.statistics(s)
"""
from .fcs import Sample, load_sample, load_folder
from .transforms import Transform, default_transform
from .gates import PolygonGate, RectangleGate, EllipseGate, QuadrantGate, gate_from_dict
from .tree import GatingTree

__version__ = "0.1.0"

__all__ = [
    "Sample", "load_sample", "load_folder",
    "Transform", "default_transform",
    "PolygonGate", "RectangleGate", "EllipseGate", "QuadrantGate", "gate_from_dict",
    "GatingTree",
]
