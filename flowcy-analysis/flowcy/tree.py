"""Gating-Hierarchie: Vorlage für alle Proben, mit Anpassungen je Probe und Änderungsprotokoll."""
from __future__ import annotations

import datetime as _dt
import json
import pathlib
import uuid
from dataclasses import dataclass, field

import numpy as np

from .gates import Gate, gate_from_dict

ROOT = "root"
FORMAT_VERSION = 1


@dataclass
class Node:
    id: str
    name: str
    parent: str | None
    gate: Gate | None = None
    note: str = ""

    def to_dict(self):
        return {"id": self.id, "name": self.name, "parent": self.parent,
                "gate": self.gate.to_dict() if self.gate else None, "note": self.note}


def robust_cv(values):
    """Robuster CV wie in FlowJo: ½·(P84,13 − P15,87) / Median · 100."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size < 3:
        return float("nan")
    med = np.median(v)
    if med == 0:
        return float("nan")
    p84, p16 = np.percentile(v, [84.13, 15.87])
    return float(100.0 * 0.5 * (p84 - p16) / abs(med))


def channel_stats(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"median": float("nan"), "mean": float("nan"), "geomean": float("nan"), "rcv": float("nan")}
    pos = v[v > 0]
    return {
        "median": float(np.median(v)),
        "mean": float(np.mean(v)),
        "geomean": float(np.exp(np.mean(np.log(pos)))) if pos.size else float("nan"),
        "rcv": robust_cv(v),
    }


@dataclass
class GatingTree:
    name: str = "Gating-Strategie"
    nodes: dict = field(default_factory=dict)
    overrides: dict = field(default_factory=dict)   # Probe → {Knoten-ID: Gate}
    log: list = field(default_factory=list)
    compensated: bool = True

    def __post_init__(self):
        if ROOT not in self.nodes:
            self.nodes[ROOT] = Node(ROOT, "Alle Ereignisse", None)
        self._cache = {}

    # --- Struktur -------------------------------------------------------------
    def add(self, name, gate: Gate, parent: str = ROOT, note: str = "", user: str = "") -> str:
        if parent not in self.nodes:
            raise KeyError(f"Elternknoten fehlt: {parent}")
        nid = uuid.uuid4().hex[:10]
        self.nodes[nid] = Node(nid, name, parent, gate, note)
        self._touch(f"Gate „{name}“ angelegt unter „{self.nodes[parent].name}“", user)
        return nid

    def children(self, nid: str) -> list:
        return [n.id for n in self.nodes.values() if n.parent == nid]

    def descendants(self, nid: str) -> list:
        out, stack = [], [nid]
        while stack:
            cur = stack.pop()
            for c in self.children(cur):
                out.append(c)
                stack.append(c)
        return out

    def remove(self, nid: str, user: str = ""):
        if nid == ROOT:
            raise ValueError("Die Wurzel kann nicht gelöscht werden")
        name = self.nodes[nid].name
        for d in [nid] + self.descendants(nid):
            self.nodes.pop(d, None)
            for ov in self.overrides.values():
                ov.pop(d, None)
        self._touch(f"Gate „{name}“ gelöscht", user)

    def rename(self, nid: str, name: str, user: str = ""):
        old = self.nodes[nid].name
        self.nodes[nid].name = name
        self._touch(f"Gate „{old}“ umbenannt in „{name}“", user)

    def set_gate(self, nid: str, gate: Gate, sample: str | None = None, user: str = ""):
        """Gate der Vorlage ersetzen oder – mit sample – nur für diese Probe."""
        if sample:
            self.overrides.setdefault(sample, {})[nid] = gate
            self._touch(f"Gate „{self.nodes[nid].name}“ für {sample} angepasst", user)
        else:
            self.nodes[nid].gate = gate
            self._touch(f"Gate „{self.nodes[nid].name}“ angepasst", user)

    def reset_override(self, nid: str, sample: str, user: str = ""):
        if self.overrides.get(sample, {}).pop(nid, None) is not None:
            self._touch(f"Gate „{self.nodes[nid].name}“ für {sample} auf Vorlage zurückgesetzt", user)

    def gate_for(self, nid: str, sample: str | None):
        if sample and nid in self.overrides.get(sample, {}):
            return self.overrides[sample][nid]
        return self.nodes[nid].gate

    def path(self, nid: str) -> list:
        out = []
        while nid is not None:
            out.append(self.nodes[nid].name)
            nid = self.nodes[nid].parent
        return list(reversed(out))

    def order(self) -> list:
        """Knoten in Baumreihenfolge (Tiefensuche), mit Tiefe."""
        out = []

        def walk(nid, depth):
            out.append((nid, depth))
            for c in self.children(nid):
                walk(c, depth + 1)
        walk(ROOT, 0)
        return out

    def _eval_order(self):
        """Baumreihenfolge, boolesche Knoten erst nach den Knoten, auf die sie verweisen."""
        pending = [nid for nid, _ in self.order()[1:]]
        done, out = {ROOT}, []
        while pending:
            progressed = False
            for nid in list(pending):
                g = self.nodes[nid].gate
                need = {self.nodes[nid].parent} | (set(g.refs) if g is not None and g.kind == "boolean" else set())
                if need <= done:
                    out.append(nid)
                    done.add(nid)
                    pending.remove(nid)
                    progressed = True
            if not progressed:
                raise ValueError("Zirkuläre Verweise in booleschen Gates: " + ", ".join(pending))
        return out

    def find(self, name_or_path: str) -> str:
        for nid, n in self.nodes.items():
            if n.name == name_or_path or "/".join(self.path(nid)[1:]) == name_or_path:
                return nid
        raise KeyError(name_or_path)

    def _touch(self, text, user=""):
        self._cache.clear()
        self.log.append({"time": _dt.datetime.now().isoformat(timespec="seconds"), "user": user, "action": text})

    # --- Auswertung -----------------------------------------------------------
    def masks(self, sample) -> dict:
        key = (id(sample), sample.name, len(self.log))
        if key in self._cache:
            return self._cache[key]
        root = np.ones(sample.n_events, dtype=bool)
        if sample.excluded is not None:
            root &= ~sample.excluded
        res = {ROOT: root}
        for nid in self._eval_order():
            parent_mask = res[self.nodes[nid].parent]
            rows = np.flatnonzero(parent_mask)
            m = np.zeros(sample.n_events, dtype=bool)
            gate = self.gate_for(nid, sample.name)
            if gate is not None and gate.kind == "boolean":
                m = gate.combine(res) & parent_mask
            elif gate is not None and rows.size:
                m[rows] = gate.evaluate(sample, self.compensated, rows=rows)
            res[nid] = m
        self._cache = {key: res}
        return res

    def statistics(self, sample, channels=None) -> list:
        masks = self.masks(sample)
        total = int(masks[ROOT].sum())
        ev = sample.events(self.compensated)
        rows = []
        for nid, depth in self.order():
            m = masks[nid]
            n = int(m.sum())
            parent = self.nodes[nid].parent
            np_ = int(masks[parent].sum()) if parent else total
            row = {
                "sample": sample.name, "group": sample.group, "id": nid, "name": self.nodes[nid].name,
                "path": "/".join(self.path(nid)[1:]) or self.nodes[nid].name, "depth": depth,
                "count": n,
                "pct_parent": 100.0 * n / np_ if np_ else float("nan"),
                "pct_total": 100.0 * n / total if total else float("nan"),
                "override": nid in self.overrides.get(sample.name, {}),
            }
            for ch in channels or []:
                st = channel_stats(ev[m, sample.index(ch)])
                for k, v in st.items():
                    row[f"{k}:{ch}"] = v
            rows.append(row)
        return rows

    def batch_statistics(self, samples, channels=None):
        import pandas as pd
        rows = []
        for s in samples:
            rows.extend(self.statistics(s, channels))
        return pd.DataFrame(rows)

    # --- Speichern ------------------------------------------------------------
    def to_dict(self):
        return {
            "format": "flowcy-gating", "version": FORMAT_VERSION, "name": self.name,
            "compensated": self.compensated,
            "nodes": [n.to_dict() for n in self.nodes.values()],
            "overrides": {s: {nid: g.to_dict() for nid, g in ov.items()} for s, ov in self.overrides.items()},
            "log": self.log,
        }

    @classmethod
    def from_dict(cls, d):
        nodes = {}
        for nd in d["nodes"]:
            nodes[nd["id"]] = Node(nd["id"], nd["name"], nd["parent"],
                                   gate_from_dict(nd["gate"]) if nd.get("gate") else None, nd.get("note", ""))
        t = cls(name=d.get("name", "Gating-Strategie"), nodes=nodes,
                overrides={s: {nid: gate_from_dict(g) for nid, g in ov.items()} for s, ov in d.get("overrides", {}).items()},
                log=list(d.get("log", [])), compensated=d.get("compensated", True))
        return t

    def save(self, path):
        pathlib.Path(path).write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=1), encoding="utf-8")

    @classmethod
    def load(cls, path):
        return cls.from_dict(json.loads(pathlib.Path(path).read_text(encoding="utf-8")))
