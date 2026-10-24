"""Lokaler Server der Desktop-App. Läuft nur auf 127.0.0.1; die Oberfläche liegt in web/."""
from __future__ import annotations

import base64
import csv
import io
import json
import pathlib
import threading

import numpy as np
from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles

from . import __version__
from .density import encode_grid, hist1d
from .gates import QuadrantGate, quadrant_set
from .project import Project
from .proposals import main_population, singlet_band, valley_threshold
from .stats import compare_groups, count_cv
from .tree import ROOT

WEB = pathlib.Path(__file__).parent / "web"


class State:
    project = Project()
    lock = threading.RLock()
    results_cache = {}


S = State()
app = FastAPI(title="FlowCy Analysis", version=__version__)


def P() -> Project:
    return S.project


def _err(e, code=400):
    raise HTTPException(code, detail=str(e))


# --- Projekt & Proben -------------------------------------------------------------

@app.get("/api/info")
def info():
    p = P()
    return {"version": __version__, "project": p.name, "folder": p.folder, "file": p.file,
            "samples": len(p.samples), "channels": p.channels() if p.samples else []}


@app.post("/api/open")
def open_folder(body: dict = Body(...)):
    with S.lock:
        try:
            if body.get("new", True):
                S.project = Project()
            P().open_folder(body["folder"])
        except Exception as e:  # noqa: BLE001
            _err(e)
        S.results_cache.clear()
        return info()


@app.post("/api/demo")
def demo(body: dict = Body(default={})):
    from .synth import make_study
    folder = pathlib.Path(body.get("folder") or (pathlib.Path.home() / "FlowCy" / "Beispielstudie"))
    if not folder.exists() or not any(folder.glob("*.fcs")):
        make_study(folder)
    return open_folder({"folder": str(folder)})


@app.get("/api/samples")
def samples():
    return [{k: e.get(k) for k in ("name", "events", "group", "qc", "error", "cytometer", "date", "vol_nl")}
            for e in P().samples.values()]


@app.get("/api/samples/{name}")
def sample_detail(name: str):
    p = P()
    if name not in p.samples:
        _err("Probe unbekannt", 404)
    s = p.sample(name)
    d = s.summary()
    d["excluded_time"] = p.excluded_time.get(name, [])
    d["metadata"] = {k: str(v)[:200] for k, v in s.metadata.items() if not k.startswith(("p", "spill"))}
    if s.compensation:
        chans, m = s.compensation
        d["spill"] = {"channels": chans, "matrix": np.round(m, 6).tolist(),
                      "labels": [s.labels[s.index(c)] for c in chans]}
    return d


@app.put("/api/samples/{name}/group")
def set_group(name: str, body: dict = Body(...)):
    P().set_group(name, body.get("group", ""))
    S.results_cache.clear()
    return {"ok": True}


@app.post("/api/groups/regex")
def groups_regex(body: dict = Body(...)):
    try:
        n = P().groups_from_regex(body["regex"], body.get("mapping"))
    except Exception as e:  # noqa: BLE001
        _err(e)
    S.results_cache.clear()
    return {"assigned": n}


@app.post("/api/groups/table")
def groups_table(body: dict = Body(...)):
    text = body["csv"]
    dialect = csv.Sniffer().sniff(text[:2000], delimiters=",;\t")
    rows = list(csv.DictReader(io.StringIO(text), dialect=dialect))
    n = P().groups_from_table(rows, body["key_col"], body["group_col"], body.get("regex") or None)
    S.results_cache.clear()
    return {"assigned": n, "columns": list(rows[0].keys()) if rows else []}


# --- QC ---------------------------------------------------------------------------

@app.post("/api/qc/{name}")
def qc(name: str):
    with S.lock:
        return P().run_qc(name)


@app.post("/api/exclude")
def exclude(body: dict = Body(...)):
    P().exclude_segments(body["sample"], body.get("segments", []), P().user)
    S.results_cache.clear()
    return {"ok": True}


# --- Achsen -----------------------------------------------------------------------

def _axis_payload(ch, s=None):
    t, lo, hi = P().axis(ch, s)
    ticks = [{"pos": (pos - lo) / (hi - lo), "label": lab} for _, pos, lab in t.ticks(lo, hi) if lo <= pos <= hi]
    labels = {c["name"]: c["display"] for c in P().channels()}
    return {"channel": ch, "label": labels.get(ch, ch), "lo": lo, "hi": hi, "ticks": ticks,
            "transform": t.to_dict()}


@app.get("/api/axes/{channel}")
def get_axis(channel: str):
    return _axis_payload(channel)


@app.put("/api/axes/{channel}")
def put_axis(channel: str, body: dict = Body(...)):
    with S.lock:
        try:
            P().set_axis(channel, body.get("kind"), body.get("params"), body.get("lo"), body.get("hi"),
                         body.get("raw_lo"), body.get("raw_hi"))
        except Exception as e:  # noqa: BLE001
            _err(e)
        return _axis_payload(channel)


# --- Plot ---------------------------------------------------------------------------

def _shape(g, tx, ty, lo_x, hi_x, lo_y, hi_y):
    """Gate in normierten Plotkoordinaten (0–1) der aktuellen Achsen."""
    nx = lambda v: (v - lo_x) / (hi_x - lo_x)  # noqa: E731
    ny = lambda v: (v - lo_y) / (hi_y - lo_y) if ty is not None else None  # noqa: E731
    same = g.xform.to_dict() == tx.to_dict() and (ty is None or g.yform.to_dict() == ty.to_dict())

    def to_axis_x(raw):
        return float(tx.forward(np.array([raw]))[0])

    def to_axis_y(raw):
        return float(ty.forward(np.array([raw]))[0])

    if g.kind == "polygon":
        v = np.asarray(g.vertices, float)
        if not same:   # Kanten in Gate-Raum verdichten, damit sie auf geänderter Achse richtig liegen
            dv = np.column_stack([g.xform.forward(v[:, 0]), g.yform.forward(v[:, 1])])
            pts = []
            for i in range(len(dv)):
                a, b = dv[i], dv[(i + 1) % len(dv)]
                for t in np.linspace(0, 1, 16, endpoint=False):
                    pts.append(a + (b - a) * t)
            dv = np.asarray(pts)
            v = np.column_stack([g.xform.inverse(dv[:, 0]), g.yform.inverse(dv[:, 1])])
        pts = np.column_stack([nx(tx.forward(v[:, 0])), ny(ty.forward(v[:, 1]))])
        return {"type": "polygon", "points": pts.tolist(), "editable": same}
    if g.kind == "rectangle":
        def f(raw, conv, n, default):
            return default if raw is None else n(conv(raw))
        return {"type": "rect" if g.y else "range",
                "x0": f(g.xmin, to_axis_x, nx, -0.05), "x1": f(g.xmax, to_axis_x, nx, 1.05),
                "y0": f(g.ymin, to_axis_y, ny, -0.05) if g.y else None,
                "y1": f(g.ymax, to_axis_y, ny, 1.05) if g.y else None, "editable": same}
    if g.kind == "ellipse":
        cx, cy = g.xform.forward(np.array([g.center[0]]))[0], g.yform.forward(np.array([g.center[1]]))[0]
        a = np.deg2rad(g.angle)
        t = np.linspace(0, 2 * np.pi, 72, endpoint=False)
        ex = cx + g.rx * np.cos(t) * np.cos(a) - g.ry * np.sin(t) * np.sin(a)
        ey = cy + g.rx * np.cos(t) * np.sin(a) + g.ry * np.sin(t) * np.cos(a)
        rx, ry = g.xform.inverse(ex), g.yform.inverse(ey)
        pts = np.column_stack([nx(tx.forward(rx)), ny(ty.forward(ry))])
        return {"type": "polygon", "points": pts.tolist(), "editable": False, "ellipse": True}
    if g.kind == "quadrant":
        return {"type": "quad", "x": nx(to_axis_x(g.xdiv)), "y": ny(to_axis_y(g.ydiv)),
                "quadrant": g.quadrant, "group": g.group, "editable": same}
    return None


@app.post("/api/plot")
def plot(body: dict = Body(...)):
    p = P()
    with S.lock:
        try:
            s = p.sample(body["sample"])
        except Exception as e:  # noqa: BLE001
            _err(e)
        tree = p.tree
        node = body.get("node", ROOT)
        if node not in tree.nodes:
            node = ROOT
        masks = tree.masks(s)
        rows = masks[node]
        n_parent = int(rows.sum())
        x, y = body["x"], body.get("y") or None
        bins = int(body.get("bins", 256))
        tx, xlo, xhi = p.axis(x, s)
        dx = tx.forward(s.column(x)[rows])
        out = {"sample": s.name, "node": node, "n": n_parent, "x": _axis_payload(x, s), "bins": bins}
        if y:
            ty, ylo, yhi = p.axis(y, s)
            dy = ty.forward(s.column(y)[rows])
            ok = np.isfinite(dx) & np.isfinite(dy)
            h, _, _ = np.histogram2d(dx[ok], dy[ok], bins=bins, range=[[xlo, xhi], [ylo, yhi]])
            if body.get("smooth"):
                from scipy import ndimage as ndi
                h = ndi.gaussian_filter(h, float(body["smooth"]))
            grid = encode_grid(h, "log").T[::-1]          # Zeilen von oben nach unten
            out["grid"] = base64.b64encode(np.ascontiguousarray(grid).tobytes()).decode()
            out["y"] = _axis_payload(y, s)
            out["offscale"] = int(n_parent - h.sum())
        else:
            ty, ylo, yhi = None, 0, 1
            h, _ = hist1d(dx, bins=bins, smooth=float(body.get("smooth", 1.0)), lo=xlo, hi=xhi)
            out["hist"] = (h / h.max()).tolist() if h.max() > 0 else h.tolist()
        gates = []
        seen_quads = set()
        for cid in tree.children(node):
            g = tree.gate_for(cid, s.name)
            if g is None or g.kind == "boolean":
                continue
            if g.x != x or (g.y or None) != y:
                continue
            shape = _shape(g, tx, ty, xlo, xhi, ylo, yhi)
            n = int(masks[cid].sum())
            entry = {"id": cid, "name": tree.nodes[cid].name, "kind": g.kind, "shape": shape, "count": n,
                     "pct": 100.0 * n / n_parent if n_parent else 0.0,
                     "override": cid in tree.overrides.get(s.name, {})}
            if g.kind == "quadrant":
                entry["first_of_group"] = g.group not in seen_quads
                seen_quads.add(g.group)
            gates.append(entry)
        out["gates"] = gates
        return JSONResponse(out)


# --- Hierarchie -----------------------------------------------------------------------

def _tree_payload(sample=None):
    p = P()
    tree = p.tree
    stats = {}
    if sample and sample in p.samples:
        s = p.sample(sample)
        stats = {r["id"]: r for r in tree.statistics(s)}
        vol = p.samples[sample].get("vol_nl")
    else:
        vol = None
    nodes = []
    for nid, depth in tree.order():
        n = tree.nodes[nid]
        g = n.gate
        st = stats.get(nid, {})
        cnt = st.get("count")
        nodes.append({
            "id": nid, "name": n.name, "parent": n.parent, "depth": depth,
            "kind": g.kind if g else None, "x": getattr(g, "x", None), "y": getattr(g, "y", None),
            "count": cnt, "pct_parent": st.get("pct_parent"), "pct_total": st.get("pct_total"),
            "override": st.get("override", False),
            "count_cv": count_cv(cnt) if cnt is not None else None,
            "per_ul": (cnt / (vol / 1000.0)) if (cnt is not None and vol) else None,
            "path": "/".join(tree.path(nid)[1:]),
            "refs": getattr(g, "refs", None), "op": getattr(g, "op", None),
        })
    return {"name": tree.name, "nodes": nodes, "log": tree.log[-50:]}


@app.get("/api/tree")
def get_tree(sample: str | None = None):
    with S.lock:
        return _tree_payload(sample)


@app.post("/api/gates")
def add_gate(body: dict = Body(...)):
    p = P()
    with S.lock:
        try:
            s = p.sample(body["sample"]) if body.get("sample") else None
            spec = body["spec"]
            parent = body.get("parent", ROOT)
            if spec["kind"] == "quadrant":
                g0 = p.gate_from_display({**spec, "quadrant": "++"}, s)
                names = body.get("names") or None
                ids = []
                xl = _short(p, g0.x)
                yl = _short(p, g0.y)
                default = [f"{xl}⁻ {yl}⁺", f"{xl}⁺ {yl}⁺", f"{xl}⁻ {yl}⁻", f"{xl}⁺ {yl}⁻"]
                for name, g in quadrant_set(g0.x, g0.y, g0.xdiv, g0.ydiv, g0.xform, g0.yform, names or default):
                    ids.append(p.tree.add(name, g, parent, user=p.user))
                S.results_cache.clear()
                return {"ids": ids}
            if spec["kind"] == "boolean":
                from .gates import BooleanGate
                g = BooleanGate(spec["op"], spec["refs"], spec.get("complements"))
            else:
                g = p.gate_from_display(spec, s)
            nid = p.tree.add(body.get("name") or "Gate", g, parent, user=p.user)
        except Exception as e:  # noqa: BLE001
            _err(e)
        S.results_cache.clear()
        return {"id": nid}


def _short(p, ch):
    lab = next((c["label"] for c in p.channels() if c["name"] == ch), "")
    return (lab.split()[0] if lab else ch)


@app.put("/api/gates/{nid}")
def update_gate(nid: str, body: dict = Body(...)):
    p = P()
    with S.lock:
        if nid not in p.tree.nodes:
            _err("Gate unbekannt", 404)
        try:
            s = p.sample(body["sample"]) if body.get("sample") else None
            spec = body["spec"]
            only = body.get("scope") == "sample" and s is not None
            old = p.tree.gate_for(nid, s.name if s else None)
            if spec["kind"] == "quadrant" and isinstance(old, QuadrantGate):
                # alle vier Quadranten der Gruppe gemeinsam verschieben
                for qid, n in p.tree.nodes.items():
                    q = p.tree.gate_for(qid, s.name if (only and s) else None)
                    if isinstance(q, QuadrantGate) and q.group == old.group:
                        g = p.gate_from_display({**spec, "quadrant": q.quadrant, "group": q.group}, s)
                        g.id = q.id
                        p.tree.set_gate(qid, g, sample=s.name if only else None, user=p.user)
            else:
                g = p.gate_from_display(spec, s)
                g.id = old.id if old is not None else g.id
                p.tree.set_gate(nid, g, sample=s.name if only else None, user=p.user)
        except Exception as e:  # noqa: BLE001
            _err(e)
        S.results_cache.clear()
        return {"ok": True}


@app.delete("/api/gates/{nid}")
def delete_gate(nid: str):
    p = P()
    with S.lock:
        g = p.tree.nodes.get(nid).gate if nid in p.tree.nodes else None
        if isinstance(g, QuadrantGate):
            for qid in [q for q, n in p.tree.nodes.items() if isinstance(n.gate, QuadrantGate) and n.gate.group == g.group]:
                if qid in p.tree.nodes:
                    p.tree.remove(qid, user=p.user)
        else:
            p.tree.remove(nid, user=p.user)
        S.results_cache.clear()
        return {"ok": True}


@app.post("/api/gates/{nid}/rename")
def rename_gate(nid: str, body: dict = Body(...)):
    P().tree.rename(nid, body["name"], user=P().user)
    S.results_cache.clear()
    return {"ok": True}


@app.post("/api/gates/{nid}/reset")
def reset_gate(nid: str, body: dict = Body(...)):
    P().tree.reset_override(nid, body["sample"], user=P().user)
    S.results_cache.clear()
    return {"ok": True}


@app.post("/api/propose")
def propose(body: dict = Body(...)):
    """Gate-Vorschlag in Anzeigekoordinaten – die Oberfläche zeigt ihn zum Annehmen an."""
    p = P()
    with S.lock:
        s = p.sample(body["sample"])
        rows = np.flatnonzero(p.tree.masks(s)[body.get("node", ROOT)])
        x, y = body["x"], body.get("y")
        tx, xlo, xhi = p.axis(x, s)
        try:
            if body["kind"] == "valley":
                g, info_ = valley_threshold(s, x, tx, rows=rows, keep=body.get("keep", "below"))
            elif body["kind"] == "singlet":
                ty = p.axis(y, s)[0]
                g, info_ = singlet_band(s, x, y, rows=rows, k=float(body.get("k", 3.0)), xform=tx, yform=ty)
            elif body["kind"] == "population":
                ty, ylo, yhi = p.axis(y, s)
                seed = body.get("seed")
                if seed is not None:   # normierte Plotkoordinaten → Anzeigeskala
                    seed = (xlo + seed[0] * (xhi - xlo), ylo + seed[1] * (yhi - ylo))
                g, info_ = main_population(s, x, y, tx, ty, rows=rows, coverage=float(body.get("coverage", 0.9)),
                                           seed=seed)
            else:
                _err("Vorschlagsart unbekannt")
        except Exception as e:  # noqa: BLE001
            _err(e)
        ty = p.axis(y, s)[0] if y else None
        ylo, yhi = (p.axis(y, s)[1:] if y else (0, 1))
        shape = _shape(g, tx, ty, xlo, xhi, ylo, yhi)
        disp = Project.gate_to_display(g)
        info_.pop("hist", None)
        info_.pop("edges", None)
        return {"shape": shape, "spec": {"kind": g.kind, "x": g.x, "y": g.y, **disp.get("display", {})},
                "info": info_}


# --- Ergebnisse & Vergleich -------------------------------------------------------------

def _results(channels=()):
    key = (len(P().tree.log), tuple(channels), tuple((n, e["group"]) for n, e in P().samples.items()))
    if key in S.results_cache:
        return S.results_cache[key]
    rows = []
    for name, e in P().samples.items():
        if e.get("error"):
            continue
        s = P().sample(name)
        s.group = e.get("group", "")
        for r in P().tree.statistics(s, channels=[c for c in channels if c in s.channels]):
            vol = e.get("vol_nl")
            r["per_ul"] = r["count"] / (vol / 1000.0) if vol else None
            r["count_cv"] = count_cv(r["count"])
            rows.append(r)
    S.results_cache = {key: rows}
    return rows


@app.get("/api/results")
def results(channels: str = ""):
    with S.lock:
        ch = [c for c in channels.split(",") if c]
        return _results(ch)


@app.post("/api/compare")
def compare(body: dict = Body(...)):
    with S.lock:
        metric = body.get("metric", "pct_parent")
        ch = [metric.split(":", 1)[1]] if ":" in metric else []
        rows = [r for r in _results(ch) if r["id"] == body["node"]]
        vals, per_sample = {}, []
        for r in rows:
            v = r.get(metric)
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            g = r["group"] or "ohne Gruppe"
            vals.setdefault(g, []).append(float(v))
            per_sample.append({"sample": r["sample"], "group": g, "value": float(v), "count": r["count"],
                               "count_cv": r["count_cv"],
                               "qc": (P().samples[r["sample"]].get("qc") or {}).get("status")})
        res = compare_groups(vals, body.get("reference"), body.get("correction", "holm"))
        res["samples"] = per_sample
        res["metric"] = metric
        return res


@app.get("/api/export/stats.csv")
def export_stats(channels: str = ""):
    import pandas as pd
    ch = [c for c in channels.split(",") if c]
    df = pd.DataFrame(_results(ch))
    buf = io.StringIO()
    df.to_csv(buf, index=False, sep=";", decimal=",")
    return PlainTextResponse(buf.getvalue(), media_type="text/csv",
                             headers={"Content-Disposition": "attachment; filename=flowcy_statistik.csv"})


@app.get("/api/export/strategy.json")
def export_strategy():
    return Response(json.dumps(P().tree.to_dict(), ensure_ascii=False, indent=1), media_type="application/json",
                    headers={"Content-Disposition": "attachment; filename=gating.flowcy.json"})


@app.post("/api/import/strategy")
def import_strategy(body: dict = Body(...)):
    from .tree import GatingTree
    with S.lock:
        try:
            P().tree = GatingTree.from_dict(body["strategy"])
        except Exception as e:  # noqa: BLE001
            _err(e)
        S.results_cache.clear()
        return _tree_payload()


@app.post("/api/import/wsp")
def import_wsp(body: dict = Body(...)):
    from .flowjo import import_workspace
    with S.lock:
        try:
            res = import_workspace(body["path"], body.get("fcs_dir") or P().folder)
        except Exception as e:  # noqa: BLE001
            _err(e)
        first = next((r for r in res if r[0] is not None), None)
        if first is None:
            _err("Keine Probe des Workspaces gefunden")
        P().tree = first[1]
        S.results_cache.clear()
        return {"skipped": first[2], "tree": _tree_payload()}


# --- Speichern -------------------------------------------------------------------------

@app.post("/api/project/save")
def save_project(body: dict = Body(default={})):
    with S.lock:
        return {"file": P().save(body.get("path"))}


@app.post("/api/project/load")
def load_project(body: dict = Body(...)):
    with S.lock:
        try:
            S.project = Project.load(body["path"])
        except Exception as e:  # noqa: BLE001
            _err(e)
        S.results_cache.clear()
        return info()


@app.post("/api/project/user")
def set_user(body: dict = Body(...)):
    P().user = body.get("user", "")
    return {"ok": True}


# --- Oberfläche --------------------------------------------------------------------------

@app.get("/")
def index():
    return FileResponse(WEB / "index.html")


app.mount("/static", StaticFiles(directory=str(WEB)), name="static")
