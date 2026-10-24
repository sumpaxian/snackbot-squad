"""Projekt: Proben, Gruppen, Achseneinstellungen, Gating-Strategie und Ergebnisse.

Der Datenordner wird nie beschrieben. Das Projekt wird als eine JSON-Datei gespeichert,
standardmäßig unter ~/FlowCy/Projekte/.
"""
from __future__ import annotations

import collections
import datetime as dt
import json
import pathlib
import re
import threading

import numpy as np

from .density import auto_range
from .fcs import load_sample
from .gates import gate_from_dict
from .qc import run_qc
from .transforms import Transform, default_transform, estimate_cofactor, is_scatter, is_time
from .tree import ROOT, GatingTree

PROJECT_DIR = pathlib.Path.home() / "FlowCy" / "Projekte"
CACHE_DIR = pathlib.Path.home() / "FlowCy" / "Cache"


def header(path):
    """Nur den TEXT-Teil lesen – schnell, auch bei großen Dateien."""
    import flowio
    f = flowio.FlowData(str(path), only_text=True, ignore_offset_error=True, ignore_offset_discrepancy=True)
    ch = [f.channels[i] for i in sorted(f.channels)]
    return {
        "events": int(f.event_count),
        "channels": [c.get("pnn", "") for c in ch],
        "labels": [c.get("pns", "") or "" for c in ch],
        "ranges": [float(c.get("pnr", 0) or 0) for c in ch],
        "text": {k: v for k, v in f.text.items() if not k.startswith(("p", "spill", "comp")) or k in ("par",)},
    }


class Project:
    def __init__(self, name="Neues Projekt", user=""):
        self.name = name
        self.user = user
        self.folder = None
        self.samples = collections.OrderedDict()   # Name → Eintrag (Pfad, Kopf, Gruppe, QC)
        self.trees = {"Standard": GatingTree("Standard")}   # Panel → Gating-Strategie
        self.panel = "Standard"                               # aktives Panel der Oberfläche
        self.axes = {}                              # Kanal → {"transform": dict, "lo": float, "hi": float}
        self.excluded_time = {}                     # Probe → [[start, ende], …] aus der Fluss-QC
        self.file = None
        self._cache = collections.OrderedDict()
        self._lock = threading.RLock()
        self.cache_size = 3
        self.cache_dir = CACHE_DIR
        self.apply_gain = True

    # --- Panels ------------------------------------------------------------------
    @property
    def tree(self):
        return self.trees.setdefault(self.panel, GatingTree(self.panel))

    @tree.setter
    def tree(self, t):
        self.trees[self.panel] = t

    def panel_of(self, name):
        return (self.samples.get(name) or {}).get("panel") or "Standard"

    def tree_of(self, name, activate=True):
        panel = self.panel_of(name)
        if activate:
            self.panel = panel
        return self.trees.setdefault(panel, GatingTree(panel))

    def tree_with(self, nid):
        """Strategie, die den Knoten enthält – bevorzugt die aktive."""
        if nid in self.tree.nodes:
            return self.tree
        return next((t for t in self.trees.values() if nid in t.nodes), self.tree)

    def set_panel(self, name, panel):
        self.samples[name]["panel"] = panel or "Standard"
        self.trees.setdefault(panel or "Standard", GatingTree(panel or "Standard"))

    def panels_from_regex(self, regex, mapping=None):
        rx = re.compile(regex)
        n = 0
        for name in self.samples:
            m = rx.search(name)
            if m:
                key = m.group(1) if m.groups() else m.group(0)
                self.set_panel(name, (mapping or {}).get(key, key))
                n += 1
        return n

    def copy_strategy(self, src, dst, user=""):
        """Strategie eines Panels in ein anderes übernehmen (z. B. gemeinsame Streulicht-Gates)."""
        t = GatingTree.from_dict(self.trees[src].to_dict())
        t.name = dst
        t.overrides.clear()
        t._touch(f"Strategie aus Panel „{src}“ übernommen", user)
        self.trees[dst] = t
        return t

    # --- Proben -----------------------------------------------------------------
    def open_folder(self, folder, pattern=r"\.fcs$"):
        folder = pathlib.Path(folder).expanduser()
        if not folder.is_dir():
            raise FileNotFoundError(f"Ordner nicht gefunden: {folder}")
        self.folder = str(folder)
        files = sorted(p for p in folder.rglob("*") if p.is_file() and re.search(pattern, p.name, re.I))
        for p in files:
            name = p.name
            if name in self.samples:
                continue
            try:
                h = header(p)
                self.samples[name] = {"name": name, "path": str(p), "events": h["events"],
                                      "channels": h["channels"], "labels": h["labels"], "ranges": h["ranges"],
                                      "cytometer": h["text"].get("cyt", ""), "date": h["text"].get("date", ""),
                                      "vol_nl": _float(h["text"].get("vol")), "group": "", "panel": "Standard", "qc": None,
                                      "error": None}
            except Exception as e:  # noqa: BLE001 – Datei bleibt sichtbar, mit Fehler
                self.samples[name] = {"name": name, "path": str(p), "events": 0, "channels": [], "labels": [],
                                      "ranges": [], "group": "", "panel": "Standard", "qc": None,
                                      "error": str(e)}
        if not self.name or self.name == "Neues Projekt":
            self.name = folder.name
        return list(self.samples.values())

    def sample(self, name):
        with self._lock:
            if name in self._cache:
                self._cache.move_to_end(name)
                return self._cache[name]
            entry = self.samples[name]
            s = load_sample(entry["path"], checksum=False, cache_dir=self.cache_dir, apply_gain=self.apply_gain)
            s.group = entry.get("group", "")
            ex = self.excluded_time.get(name)
            if ex and s.time_channel:
                t = s.column(s.time_channel, compensated=False)
                m = np.zeros(s.n_events, bool)
                for a, b in ex:
                    m |= (t >= a) & (t <= b)
                s.excluded = m
            self._cache[name] = s
            while len(self._cache) > self.cache_size:
                self._cache.popitem(last=False)
            return s

    def forget(self, name):
        with self._lock:
            self._cache.pop(name, None)

    def channels(self):
        """Kanäle der Mehrheit der Proben, mit häufigster Beschriftung."""
        count = collections.Counter()
        labels = collections.defaultdict(collections.Counter)
        for e in self.samples.values():
            for c, l in zip(e["channels"], e["labels"]):
                count[c] += 1
                labels[c][l] += 1
        n = max(len(self.samples), 1)
        out = []
        first = next(iter(self.samples.values()), {"channels": []})["channels"]
        order = {c: i for i, c in enumerate(first)}
        for c in sorted(count, key=lambda c: order.get(c, 999)):
            if count[c] >= n / 2:
                lab = labels[c].most_common(1)[0][0]
                kind = "time" if is_time(c) or lab.upper() == "HDR-T" else ("scatter" if is_scatter(c) else "fluor")
                out.append({"name": c, "label": lab, "kind": kind,
                            "display": f"{lab} ({c})" if lab and lab != c else c})
        return out

    def run_qc(self, name):
        s = self.sample(name)
        q = run_qc(s)
        q.pop("_excluded_mask", None)
        self.samples[name]["qc"] = {k: q[k] for k in ("status", "issues", "flow", "margin")}
        return self.samples[name]["qc"]

    def exclude_segments(self, name, segments, user=""):
        self.excluded_time[name] = [[float(a), float(b)] for a, b in segments]
        self.forget(name)
        self.tree_of(name, activate=False)._touch(f"Zeitbereiche in {name} ausgeschlossen: " +
                         ", ".join(f"{a:.4g}–{b:.4g}" for a, b in segments) if segments else
                         f"Zeitausschluss in {name} aufgehoben", user)

    # --- Gruppen ----------------------------------------------------------------
    def set_group(self, name, group):
        self.samples[name]["group"] = group
        if name in self._cache:
            self._cache[name].group = group

    def groups_from_regex(self, regex, mapping=None):
        """Gruppe aus dem Dateinamen: erste Klammergruppe des regulären Ausdrucks, optional übersetzt."""
        rx = re.compile(regex)
        n = 0
        for name in self.samples:
            m = rx.search(name)
            if m:
                key = m.group(1) if m.groups() else m.group(0)
                self.set_group(name, (mapping or {}).get(key, key))
                n += 1
        return n

    def groups_from_table(self, rows, key_col, group_col, regex=None):
        """Gruppen aus einer Tabelle (z. B. Tier-ID → Arm). regex zieht die ID aus dem Dateinamen."""
        table = {str(r[key_col]).strip(): str(r[group_col]).strip() for r in rows if r.get(key_col) is not None}
        rx = re.compile(regex) if regex else None
        n = 0
        for name in self.samples:
            key = None
            if rx:
                m = rx.search(name)
                key = (m.group(1) if m and m.groups() else (m.group(0) if m else None))
            else:
                key = next((k for k in table if k and k in name), None)
            if key is not None and key in table:
                self.set_group(name, table[key])
                n += 1
        return n

    # --- Achsen -----------------------------------------------------------------
    def axis(self, channel, sample=None):
        """Transformation und sichtbarer Bereich eines Kanals – projektweit gleich für alle Proben."""
        if channel not in self.axes:
            s = sample or self.sample(next(iter(self.samples)))
            vals = s.column(channel)
            if is_scatter(channel):
                pos = vals[vals > 0]
                if pos.size > 100 and np.percentile(pos, 99) / max(np.percentile(pos, 1), 1e-9) > 100:
                    top = float(10 ** np.ceil(np.log10(np.percentile(pos, 99.99))))
                    t = Transform("log", {"T": top, "M": float(np.log10(top))})
                else:
                    t = default_transform(channel, s.max_range(channel), vals)
            elif is_time(channel):
                t = default_transform(channel, s.max_range(channel), vals)
            else:
                top = s.max_range(channel)
                fin = vals[np.isfinite(vals)]
                q = float(np.percentile(fin, 99.99)) if fin.size else top
                if 0 < q < top / 100:      # Daten nutzen nur einen kleinen Teil des Messbereichs
                    top = float(10 ** np.ceil(np.log10(q * 1.5)))
                t = Transform("logicle", default_transform(channel, top, vals).params)
            lo, hi = auto_range(t.forward(vals))
            self.axes[channel] = {"transform": t.to_dict(), "lo": lo, "hi": hi}
        a = self.axes[channel]
        return Transform.from_dict(a["transform"]), a["lo"], a["hi"]

    def set_axis(self, channel, kind=None, params=None, lo=None, hi=None, raw_lo=None, raw_hi=None, sample=None):
        t, old_lo, old_hi = self.axis(channel, sample)
        if kind:
            p = dict(params or {})
            if kind == "arcsinh" and "cofactor" not in p:
                s = sample or self.sample(next(iter(self.samples)))
                p["cofactor"] = estimate_cofactor(s.column(channel))
                p.setdefault("T", s.max_range(channel))
            t = Transform(kind, p)
            s = sample or self.sample(next(iter(self.samples)))
            old_lo, old_hi = auto_range(t.forward(s.column(channel)))
        if raw_lo is not None:
            lo = float(t.forward(raw_lo))
        if raw_hi is not None:
            hi = float(t.forward(raw_hi))
        self.axes[channel] = {"transform": t.to_dict(), "lo": old_lo if lo is None else float(lo),
                              "hi": old_hi if hi is None else float(hi)}
        return self.axes[channel]

    # --- Gates aus der Oberfläche (Anzeigekoordinaten) --------------------------
    def gate_from_display(self, spec, sample=None):
        """Gate aus Anzeigekoordinaten der Oberfläche in Rohwerte umrechnen."""
        kind = spec["kind"]
        x, y = spec["x"], spec.get("y")
        tx = self.axis(x, sample)[0]
        ty = self.axis(y, sample)[0] if y else None
        d = {"kind": kind, "x": x, "y": y, "xform": tx.to_dict(), "yform": ty.to_dict() if ty else None}
        inv_x = lambda v: None if v is None else float(tx.inverse(np.array([v]))[0])  # noqa: E731
        inv_y = lambda v: None if v is None else float(ty.inverse(np.array([v]))[0])  # noqa: E731
        if kind == "polygon":
            v = np.asarray(spec["vertices"], float)
            d["vertices"] = np.column_stack([tx.inverse(v[:, 0]), ty.inverse(v[:, 1])]).tolist()
        elif kind == "rectangle":
            d.update(xmin=inv_x(spec.get("xmin")), xmax=inv_x(spec.get("xmax")),
                     ymin=inv_y(spec.get("ymin")) if y else None, ymax=inv_y(spec.get("ymax")) if y else None)
        elif kind == "ellipse":
            d.update(center=[inv_x(spec["cx"]), inv_y(spec["cy"])], rx=spec["rx"], ry=spec["ry"],
                     angle=spec.get("angle", 0.0))
        elif kind == "quadrant":
            d.update(xdiv=inv_x(spec["xdiv"]), ydiv=inv_y(spec["ydiv"]), quadrant=spec.get("quadrant", "++"),
                     group=spec.get("group"))
        else:
            raise ValueError(f"Gate-Typ {kind}")
        return gate_from_dict(d)

    @staticmethod
    def gate_to_display(g):
        """Gate in Anzeigekoordinaten seiner eigenen Transformationen (für die Oberfläche)."""
        d = g.to_dict()
        if g.kind == "boolean":
            return d
        fx = lambda v: None if v is None else float(g.xform.forward(np.array([v]))[0])  # noqa: E731
        fy = lambda v: None if v is None else float(g.yform.forward(np.array([v]))[0])  # noqa: E731
        if g.kind == "polygon":
            v = np.asarray(g.vertices)
            d["display"] = {"vertices": np.column_stack([g.xform.forward(v[:, 0]), g.yform.forward(v[:, 1])]).tolist()}
        elif g.kind == "rectangle":
            d["display"] = {"xmin": fx(g.xmin), "xmax": fx(g.xmax),
                            "ymin": fy(g.ymin) if g.y else None, "ymax": fy(g.ymax) if g.y else None}
        elif g.kind == "ellipse":
            d["display"] = {"cx": fx(g.center[0]), "cy": fy(g.center[1]), "rx": g.rx, "ry": g.ry, "angle": g.angle}
        elif g.kind == "quadrant":
            d["display"] = {"xdiv": fx(g.xdiv), "ydiv": fy(g.ydiv)}
        return d

    # --- Speichern --------------------------------------------------------------
    def to_dict(self):
        return {"format": "flowcy-project", "version": 1, "name": self.name, "folder": self.folder,
                "samples": [{k: e.get(k) for k in ("name", "path", "group", "panel")} for e in self.samples.values()],
                "axes": self.axes, "excluded_time": self.excluded_time, "panel": self.panel,
                "apply_gain": self.apply_gain,
                "trees": {k: t.to_dict() for k, t in self.trees.items()},
                "saved": dt.datetime.now().isoformat(timespec="seconds")}

    def save(self, path=None):
        path = pathlib.Path(path or self.file or PROJECT_DIR / f"{_slug(self.name)}.flowcy-projekt.json")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=1), encoding="utf-8")
        self.file = str(path)
        return self.file

    @classmethod
    def load(cls, path):
        d = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
        p = cls(d.get("name", "Projekt"))
        if d.get("folder") and pathlib.Path(d["folder"]).is_dir():
            p.open_folder(d["folder"])
        for s in d.get("samples", []):
            if s["name"] in p.samples:
                p.samples[s["name"]]["group"] = s.get("group", "")
                p.samples[s["name"]]["panel"] = s.get("panel") or "Standard"
        p.axes = d.get("axes", {})
        p.apply_gain = d.get("apply_gain", True)
        p.excluded_time = d.get("excluded_time", {})
        if "trees" in d:
            p.trees = {k: GatingTree.from_dict(t) for k, t in d["trees"].items()}
        else:
            p.trees = {"Standard": GatingTree.from_dict(d["tree"])}
        p.panel = d.get("panel", next(iter(p.trees)))
        p.file = str(path)
        return p


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _slug(s):
    return re.sub(r"[^A-Za-z0-9ÄÖÜäöüß_-]+", "_", s).strip("_") or "projekt"
