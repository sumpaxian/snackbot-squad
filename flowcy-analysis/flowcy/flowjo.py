"""FlowJo-Workspaces (.wsp) lesen.

Zwei Dinge werden gebraucht:
1. die Zählwerte, die FlowJo selbst in den Workspace geschrieben hat – die Referenz für
   den Konkordanztest;
2. die Gates, Transformationen und Kompensation, umgewandelt in FlowCy-Objekte. Das Parsen
   übernimmt FlowKit; Koordinaten liegen dort in der transformierten Skala und werden für
   FlowCy in Rohwerte zurückgerechnet.
"""
from __future__ import annotations

import pathlib
import warnings
import xml.etree.ElementTree as ET

import numpy as np

from .fcs import load_sample
from .gates import BooleanGate, EllipseGate, PolygonGate, RectangleGate
from .transforms import Transform
from .tree import ROOT, GatingTree

POP_TAGS = {"Population", "AndNode", "OrNode", "NotNode"}


def _local(tag):
    return tag.rsplit("}", 1)[-1]


def reference_counts(wsp_path) -> dict:
    """{Probenname: {Gate-Pfad (Tupel ohne 'root'): Ereigniszahl}} wie von FlowJo gespeichert."""
    root = ET.parse(str(wsp_path)).getroot()
    out = {}
    for sample in root.iter():
        if _local(sample.tag) != "Sample":
            continue
        node = next((c for c in sample if _local(c.tag) == "SampleNode"), None)
        if node is None:
            continue
        name = node.get("name")
        counts = {(): int(float(node.get("count", "0") or 0))}

        def walk(el, path):
            for sub in el:
                if _local(sub.tag) != "Subpopulations":
                    continue
                for pop in sub:
                    if _local(pop.tag) in POP_TAGS:
                        p = path + (pop.get("name"),)
                        counts[p] = int(float(pop.get("count", "0") or 0))
                        walk(pop, p)
        walk(node, ())
        out[name] = counts
    return out


def convert_transform(fkt, name="") -> Transform:
    """FlowKit-Transformation → FlowCy. Bekannte Typen speicherbar, sonst extern eingebunden."""
    if fkt is None:
        return Transform("linear", {"min": 0.0, "max": 1.0})
    cls = type(fkt).__name__
    v = vars(fkt)
    if cls == "LinearTransform":
        return Transform("linear", {"min": -float(v["param_a"]), "max": float(v["param_t"])})
    if cls == "LogicleTransform":
        return Transform("logicle", {"T": float(v["param_t"]), "W": float(v["param_w"]),
                                     "M": float(v["param_m"]), "A": float(v["param_a"])})
    if cls == "AsinhTransform":
        return Transform("asinh", {"T": float(v["param_t"]), "M": float(v["param_m"]), "A": float(v["param_a"])})
    return Transform.external(fkt, name or cls)


def _raw(t: Transform, val):
    return None if val is None else float(t.inverse(np.array([val], float))[0])


def convert_gate(fkg, transforms: dict, node_ids: dict):
    """FlowKit-Gate → FlowCy-Gate. Gibt (gate, None) oder (None, Grund) zurück."""
    kind = type(fkg).__name__
    if kind == "BooleanGate":
        refs, comps = [], []
        for r in fkg.gate_refs:
            key = (r["ref"], tuple(r["path"]))
            if key not in node_ids:
                return None, f"Verweis nicht gefunden: {r['ref']}"
            refs.append(node_ids[key])
            comps.append(bool(r.get("complement", False)))
        return BooleanGate(fkg.type, refs, comps), None
    dims = fkg.dimensions
    if len(dims) > 2:
        return None, f"{kind} mit {len(dims)} Dimensionen"

    def tf(d):
        return convert_transform(transforms.get(d.transformation_ref) if d.transformation_ref else None, d.id) \
            if d.transformation_ref else Transform("linear", {"min": 0.0, "max": 1.0})

    if kind == "RectangleGate":
        tx = tf(dims[0])
        ty = tf(dims[1]) if len(dims) > 1 else None
        g = RectangleGate(dims[0].id, dims[1].id if len(dims) > 1 else None,
                          _raw(tx, dims[0].min), _raw(tx, dims[0].max),
                          _raw(ty, dims[1].min) if ty else None, _raw(ty, dims[1].max) if ty else None,
                          xform=tx, yform=ty)
    elif kind == "PolygonGate":
        tx, ty = tf(dims[0]), tf(dims[1])
        v = np.asarray(fkg.vertices, float)
        g = PolygonGate(dims[0].id, dims[1].id, np.column_stack([tx.inverse(v[:, 0]), ty.inverse(v[:, 1])]),
                        xform=tx, yform=ty)
    elif kind == "EllipsoidGate":
        tx, ty = tf(dims[0]), tf(dims[1])
        cov = np.asarray(fkg.covariance_matrix, float)
        d2 = float(fkg.distance_square)
        w, vec = np.linalg.eigh(cov)
        rx, ry = np.sqrt(w * d2)
        angle = float(np.degrees(np.arctan2(vec[1, 0], vec[0, 0])))
        cx, cy = fkg.coordinates
        g = EllipseGate(dims[0].id, dims[1].id, (_raw(tx, cx), _raw(ty, cy)), rx, ry, angle, xform=tx, yform=ty)
    else:
        return None, f"Gate-Typ {kind} nicht unterstützt"
    if getattr(fkg, "use_complement", False):
        return None, "Komplement-Gate nicht unterstützt"
    return g, None


def import_workspace(wsp_path, fcs_dir=None, samples=None):
    """Alle Proben eines Workspaces: [(Sample, GatingTree, Liste nicht übernommener Gates)]."""
    import flowkit as fk
    wsp_path = pathlib.Path(wsp_path)
    fcs_dir = pathlib.Path(fcs_dir) if fcs_dir else wsp_path.parent
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ws = fk.Workspace(str(wsp_path), fcs_samples=str(fcs_dir))
    out = []
    files = {p.name: p for p in fcs_dir.rglob("*") if p.suffix.lower() == ".fcs"}
    for sid in ws.get_sample_ids():
        if samples and sid not in samples:
            continue
        if sid not in files:
            out.append((None, None, [f"FCS-Datei fehlt: {sid}"]))
            continue
        s = load_sample(files[sid])
        comp = ws.get_comp_matrix(sid)
        if comp is not None:
            s.set_compensation(list(comp.detectors), np.asarray(comp.matrix, float))
        transforms = ws.get_transforms(sid)
        tree = GatingTree(name=f"FlowJo: {wsp_path.name}")
        tree.log.clear()
        node_ids = {}
        skipped = []
        for name, path in ws.get_gate_ids(sid):
            parent_key = (path[-1], tuple(path[:-1])) if len(path) > 1 else None
            parent = node_ids.get(parent_key, ROOT) if parent_key else ROOT
            if parent_key and parent_key not in node_ids:
                skipped.append(f"{'/'.join(path[1:] + (name,))}: Elterngate fehlt")
                continue
            fkg = ws.get_gate(sid, name, path)
            g, why = convert_gate(fkg, transforms, node_ids)
            if g is None:
                skipped.append(f"{'/'.join(path[1:] + (name,))}: {why}")
                continue
            nid = tree.add(name, g, parent)
            node_ids[(name, tuple(path))] = nid
        tree.log.clear()
        tree.log.append({"time": "", "user": "", "action": f"Aus FlowJo-Workspace {wsp_path.name} übernommen"})
        out.append((s, tree, skipped))
    return out


def flowkit_counts(wsp_path, fcs_dir=None) -> dict:
    """Zweite, unabhängige Referenz: FlowKit wertet denselben Workspace selbst aus."""
    import flowkit as fk
    wsp_path = pathlib.Path(wsp_path)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ws = fk.Workspace(str(wsp_path), fcs_samples=str(fcs_dir or wsp_path.parent))
        for sid in ws.get_sample_ids():
            try:
                ws.analyze_samples(sample_id=sid, use_mp=False)
                rep = ws.get_gating_results(sid).report
            except Exception as e:  # noqa: BLE001 – Referenz ist optional
                out[sid] = {"_error": str(e)}
                continue
            out[sid] = {tuple(gp[1:]) + (gn,): int(c) for gp, gn, c in
                        zip(rep["gate_path"], rep["gate_name"], rep["count"])}
    return out


def classify(ours, ref, rel_tol=0.005, abs_tol=5):
    if ref is None or ref < 0:
        return "ohne Referenz", None
    diff = ours - ref
    rel = abs(diff) / max(ref, 1)
    if diff == 0:
        return "gleich", rel
    if rel <= rel_tol or abs(diff) <= abs_tol:
        return "≈", rel
    return "abweichend", rel


def concordance(wsp_path, fcs_dir=None, rel_tol=0.005, abs_tol=5, with_flowkit=True):
    """FlowCy-Zählung gegen FlowJo (im Workspace gespeichert) und FlowKit, je Probe und Gate.

    "gleich" = identisch, "≈" = innerhalb rel_tol oder abs_tol Ereignissen.
    """
    ref = reference_counts(wsp_path)
    fkc = flowkit_counts(wsp_path, fcs_dir) if with_flowkit else {}
    rows = []
    for s, tree, skipped in import_workspace(wsp_path, fcs_dir):
        if s is None:
            rows.append({"sample": None, "gate": None, "status": "fehlt", "note": "; ".join(skipped)})
            continue
        r = ref.get(s.name, {})
        k = fkc.get(s.name, {})
        masks = tree.masks(s)
        for nid, _ in tree.order()[1:]:
            path = tuple(tree.path(nid)[1:])
            ours = int(masks[nid].sum())
            theirs = r.get(path)
            status, rel = classify(ours, theirs, rel_tol, abs_tol)
            fk_n = k.get(path)
            rows.append({"sample": s.name, "gate": "/".join(path), "flowjo": theirs, "flowcy": ours,
                         "flowkit": fk_n, "diff": None if rel is None else ours - theirs, "rel_diff": rel,
                         "status": status,
                         "flowcy_eq_flowkit": None if fk_n is None else ours == fk_n})
        for sk in skipped:
            rows.append({"sample": s.name, "gate": sk.split(":")[0], "flowjo": None, "flowcy": None,
                         "flowkit": None, "diff": None, "rel_diff": None, "status": "nicht übernommen",
                         "note": sk})
    return rows
