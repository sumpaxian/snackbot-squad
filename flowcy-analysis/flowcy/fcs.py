"""FCS-Dateien lesen, Metadaten und Kompensation.

Das Einlesen übernimmt FlowKit (FCS 2.0/3.0/3.1, $PnE-Dekodierung, $PnG-Verstärkung).
FlowCy hält die Rohdaten und berechnet Kompensation und Transformationen selbst, damit
jede Rechenstufe nachvollziehbar und austauschbar bleibt.
"""
from __future__ import annotations

import hashlib
import pathlib
from dataclasses import dataclass, field

import numpy as np

from .transforms import is_scatter, is_time


def parse_spill(text: str):
    """$SPILL / $SPILLOVER: 'n,Kanal1,…,Kanaln,w11,w12,…' → (Kanäle, n×n-Matrix)."""
    parts = [p.strip() for p in str(text).split(",")]
    n = int(float(parts[0]))
    chans = parts[1:1 + n]
    vals = np.array([float(v) for v in parts[1 + n:1 + n + n * n]], dtype=np.float64)
    return chans, vals.reshape(n, n)


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class Sample:
    name: str
    path: str | None
    channels: list            # $PnN
    labels: list              # $PnS, leer wenn nicht gesetzt
    raw: np.ndarray           # Ereignisse × Kanäle, unkompensiert
    metadata: dict = field(default_factory=dict)
    ranges: list = field(default_factory=list)
    spill: tuple | None = None     # (Kanäle, Matrix) aus der Datei
    checksum: str | None = None
    group: str = ""
    excluded: np.ndarray | None = None   # bool je Ereignis, z. B. aus der Fluss-QC

    def __post_init__(self):
        self._comp_cache = None
        self._comp_override = None

    # --- Kanäle -------------------------------------------------------------
    @property
    def n_events(self) -> int:
        return int(self.raw.shape[0])

    def index(self, channel: str) -> int:
        if channel in self.channels:
            return self.channels.index(channel)
        if channel in self.labels and channel:
            return self.labels.index(channel)
        low = channel.lower()
        for i, (c, l) in enumerate(zip(self.channels, self.labels)):
            if c.lower() == low or (l and l.lower() == low):
                return i
        raise KeyError(f"Kanal nicht gefunden: {channel} in {self.name}")

    def display(self, channel: str) -> str:
        i = self.index(channel)
        lab, ch = self.labels[i], self.channels[i]
        return f"{lab} ({ch})" if lab and lab != ch else ch

    @property
    def fluorescence(self) -> list:
        return [c for c in self.channels if not is_scatter(c) and not is_time(c)]

    @property
    def time_channel(self) -> str | None:
        for c in self.channels:
            if is_time(c):
                return c
        return None

    def max_range(self, channel: str) -> float:
        i = self.index(channel)
        r = self.ranges[i] if i < len(self.ranges) else None
        return float(r) if r else float(np.nanmax(self.raw[:, i]) or 1.0)

    # --- Kompensation -------------------------------------------------------
    def set_compensation(self, channels, matrix):
        """Eigene Spillover-Matrix (Zeile = Fluorochrom, Spalte = Detektor). None = die der Datei."""
        self._comp_override = None if channels is None else (list(channels), np.asarray(matrix, float))
        self._comp_cache = None

    @property
    def compensation(self):
        return self._comp_override or self.spill

    def events(self, compensated: bool = True) -> np.ndarray:
        if not compensated or self.compensation is None:
            return self.raw
        if self._comp_cache is None:
            chans, m = self.compensation
            idx = [self.index(c) for c in chans]
            out = self.raw.copy()
            out[:, idx] = self.raw[:, idx] @ np.linalg.inv(m)
            self._comp_cache = out
        return self._comp_cache

    def column(self, channel: str, compensated: bool = True) -> np.ndarray:
        return self.events(compensated)[:, self.index(channel)]

    def summary(self) -> dict:
        md = self.metadata
        return {
            "name": self.name,
            "path": self.path,
            "events": self.n_events,
            "channels": [{"name": c, "label": l, "range": r} for c, l, r in zip(self.channels, self.labels, self.ranges)],
            "cytometer": md.get("cyt", ""),
            "date": md.get("date", ""),
            "start": md.get("btim", ""),
            "end": md.get("etim", ""),
            "fcs_version": md.get("fcsversion", md.get("version", "")),
            "has_spill": self.spill is not None,
            "sha256": self.checksum,
            "group": self.group,
        }


def load_sample(path, checksum: bool = True, name: str | None = None) -> Sample:
    import flowkit as fk

    path = pathlib.Path(path)
    s = fk.Sample(str(path), ignore_offset_error=True, ignore_offset_discrepancy=True,
                  subsample=0)
    raw = np.asarray(s.get_events(source="raw"), dtype=np.float64)
    md = {str(k).lower().lstrip("$"): v for k, v in s.get_metadata().items()}
    n = len(s.pnn_labels)
    ranges = []
    for i in range(1, n + 1):
        try:
            ranges.append(float(md.get(f"p{i}r", "nan")))
        except ValueError:
            ranges.append(float("nan"))
    spill = None
    for key in ("spill", "spillover", "comp"):
        if key in md and str(md[key]).strip():
            try:
                chans, m = parse_spill(md[key])
                # manche Geräte schreiben Kanalnummern statt Namen
                chans = [s.pnn_labels[int(c) - 1] if c.isdigit() else c for c in chans]
                spill = (chans, m)
                break
            except (ValueError, IndexError):
                continue
    labels = [l or "" for l in s.pns_labels]
    return Sample(
        name=name or path.name,
        path=str(path),
        channels=list(s.pnn_labels),
        labels=labels,
        raw=raw,
        metadata=md,
        ranges=ranges,
        spill=spill,
        checksum=sha256(path) if checksum else None,
    )


def load_folder(folder, pattern: str = "*.fcs", recursive: bool = True, **kw) -> list:
    folder = pathlib.Path(folder)
    files = sorted(folder.rglob(pattern) if recursive else folder.glob(pattern))
    files += sorted(p for p in (folder.rglob(pattern.upper()) if recursive else folder.glob(pattern.upper()))
                    if p not in files)
    return [load_sample(f, **kw) for f in files]


def from_array(name, channels, data, labels=None, ranges=None, spill=None, metadata=None) -> Sample:
    """Probe aus einem Array – für Tests und für Daten aus anderen Quellen."""
    data = np.asarray(data, dtype=np.float64)
    return Sample(name=name, path=None, channels=list(channels),
                  labels=list(labels) if labels else [""] * len(channels), raw=data,
                  metadata=dict(metadata or {}),
                  ranges=list(ranges) if ranges else [262144.0] * len(channels), spill=spill)
