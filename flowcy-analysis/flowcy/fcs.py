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
            if self.raw.dtype != np.float64:      # Memory-Map aus dem Cache: einmal in float64
                self.raw = np.array(self.raw, dtype=np.float64)
            return self.raw
        if self._comp_cache is None:
            chans, m = self.compensation
            idx = [self.index(c) for c in chans]
            out = np.array(self.raw, dtype=np.float64)
            out[:, idx] = out[:, idx] @ np.linalg.inv(m)
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


CACHE_VERSION = 1


def _cache_key(path: pathlib.Path, apply_gain=True) -> str:
    st = path.stat()
    raw = f"{path.resolve()}|{st.st_size}|{st.st_mtime_ns}|{CACHE_VERSION}|gain={apply_gain}"
    return hashlib.sha1(raw.encode()).hexdigest()[:20]


def gains(metadata, n):
    """Kanäle mit Verstärkung $PnG ≠ 1 (nur lineare Kanäle, $PnE = 0,0)."""
    out = {}
    for i in range(1, n + 1):
        try:
            g = float(metadata.get(f"p{i}g", 1) or 1)
        except ValueError:
            continue
        e = str(metadata.get(f"p{i}e", "0,0")).replace(" ", "")
        if g != 1 and e in ("0,0", "0.0,0.0", "0,0.0", "0.0,0"):
            out[i - 1] = g
    return out


def prune_cache(cache_dir, max_bytes=8e9):
    """Älteste Cache-Einträge löschen, bis der Ordner unter max_bytes liegt (Standard 8 GB)."""
    cache_dir = pathlib.Path(cache_dir)
    if not cache_dir.is_dir():
        return 0
    files = sorted(cache_dir.glob("*.npy"), key=lambda f: f.stat().st_atime)
    total = sum(f.stat().st_size for f in files)
    removed = 0
    for f in files:
        if total <= max_bytes:
            break
        total -= f.stat().st_size
        f.unlink(missing_ok=True)
        f.with_suffix(".json").unlink(missing_ok=True)
        removed += 1
    return removed


def load_sample(path, checksum: bool = True, name: str | None = None, cache_dir=None,
                apply_gain: bool = True) -> Sample:
    """FCS-Datei lesen. Mit cache_dir wird ein Abbild (.npy + .json) angelegt und beim nächsten
    Mal per Memory-Map geladen – um ein Vielfaches schneller als erneutes Parsen.

    apply_gain: Werte linearer Kanäle durch $PnG teilen (FCS-Standard, FlowKit). False liefert
    die gespeicherten Werte wie flowio ohne Vorverarbeitung – nötig, wenn Gates in diesen
    Einheiten gezeichnet wurden."""
    path = pathlib.Path(path)
    if cache_dir:
        import json
        cache_dir = pathlib.Path(cache_dir)
        key = _cache_key(path, apply_gain)
        npy, meta = cache_dir / f"{key}.npy", cache_dir / f"{key}.json"
        if npy.exists() and meta.exists():
            try:
                m = json.loads(meta.read_text(encoding="utf-8"))
                raw = np.load(npy, mmap_mode="r")
                spill = (m["spill"][0], np.asarray(m["spill"][1], float)) if m.get("spill") else None
                return Sample(name=name or path.name, path=str(path), channels=m["channels"], labels=m["labels"],
                              raw=raw, metadata=m["metadata"], ranges=m["ranges"], spill=spill,
                              checksum=m.get("sha256") if checksum else None)
            except Exception:  # noqa: BLE001 – defekter Cache: neu einlesen
                pass
        s = load_sample(path, checksum=checksum, name=name, apply_gain=apply_gain)
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            r32 = s.raw.astype(np.float32)
            store = r32 if np.array_equal(r32, s.raw) else s.raw   # verlustfrei, sonst float64
            np.save(npy, store)
            prune_cache(cache_dir)
            meta.write_text(json.dumps({
                "channels": s.channels, "labels": s.labels, "ranges": s.ranges,
                "metadata": {k: str(v) for k, v in s.metadata.items()},
                "spill": [s.spill[0], s.spill[1].tolist()] if s.spill else None, "sha256": s.checksum,
                "source": str(path)}, ensure_ascii=False), encoding="utf-8")
        except OSError:
            pass
        return s
    import flowkit as fk

    s = fk.Sample(str(path), ignore_offset_error=True, ignore_offset_discrepancy=True,
                  subsample=0)
    raw = np.asarray(s.get_events(source="raw"), dtype=np.float64)
    md = {str(k).lower().lstrip("$"): v for k, v in s.get_metadata().items()}
    n = len(s.pnn_labels)
    if not apply_gain:
        for i, g in gains(md, n).items():
            raw[:, i] *= g
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
