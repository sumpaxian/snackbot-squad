"""Test-Analyzer: prüft FlowCy an vorhandenen Daten, ohne sie zu verändern.

Läuft Datei für Datei (wenig Speicher auch bei vielen großen Dateien) und schreibt alles in
einen eigenen Ausgabeordner:

  inventar.tsv           Datei, Gerät, Ereignisse, Kanäle, Kompensation, Prüfsumme, Ladezeit
  qc.tsv                 Flussrate, Sättigung, negatives Streulicht, Status
  kanalbeschriftung.tsv  Kanäle, deren $PnS-Beschriftung von der Mehrheit abweicht
  auto_gating.tsv        automatische Plausibilitäts-Strategie (Zellen, Einzelzellen, lebend)
  statistik.tsv          Ergebnisse einer übergebenen Gating-Strategie (--gates)
  konkordanz_flowjo.tsv  FlowCy gegen FlowJo und FlowKit für jeden gefundenen .wsp
  bericht.html           Zusammenfassung zum Ansehen
  protokoll.txt          Ablauf, Fehler, Laufzeiten
"""
from __future__ import annotations

import collections
import datetime as dt
import html
import json
import pathlib
import platform
import re
import time
import traceback

import numpy as np
import pandas as pd

from . import __version__
from .fcs import gains, load_sample
from .proposals import singlet_band, valley_threshold
from .qc import run_qc
from .transforms import Transform, estimate_cofactor, is_scatter, is_time
from .tree import GatingTree

LIVE_MARKERS = ("PI", "7-AAD", "7AAD", "DAPI", "ZOMBIE", "LIVE", "VIABILITY", "FVD", "LD", "AQUA", "SYTOX")


def find_live_channel(sample):
    for c, l in zip(sample.channels, sample.labels):
        text = (l or "").upper()
        token = re.split(r"[\s\-_/]", text)[0] if text else ""
        if token in LIVE_MARKERS or any(text.startswith(m + " ") or text == m for m in LIVE_MARKERS):
            return c
    return None


def scatter_transform(values, max_range):
    v = np.asarray(values, float)
    pos = v[v > 0]
    if pos.size > 100 and np.percentile(pos, 99) / max(np.percentile(pos, 1), 1e-9) > 100:
        top = 10 ** np.ceil(np.log10(np.percentile(pos, 99.99)))
        return Transform("log", {"T": float(top), "M": float(np.log10(top))})
    return Transform("linear", {"min": 0.0, "max": float(max_range)})


def auto_gating(sample):
    """Plausibilitätsprüfung ohne Vorwissen: FSC-Tal → Einzelzell-Band → Lebend-Tal."""
    out = {"sample": sample.name}
    fsc_a = next((c for c in sample.channels if c.upper().startswith("FSC") and c.upper().endswith("A")), None)
    fsc_h = next((c for c in sample.channels if c.upper().startswith("FSC") and c.upper().endswith("H")), None)
    if fsc_a is None:
        out["note"] = "kein FSC-A"
        return out
    t = scatter_transform(sample.column(fsc_a), sample.max_range(fsc_a))
    g, info = valley_threshold(sample, fsc_a, t, keep="above")
    cells = g.evaluate(sample)
    out.update(fsc_scale=t.kind, fsc_cut=round(info["cut"], 2), fsc_modes=info["modes"],
               cells=int(cells.sum()), cells_pct=100 * cells.mean())
    rows = np.flatnonzero(cells)
    if fsc_h and rows.size > 100:
        g2, i2 = singlet_band(sample, fsc_a, fsc_h, rows=rows)
        sing = np.zeros(sample.n_events, bool)
        sing[rows] = g2.evaluate(sample, rows=rows)
        out.update(singlets=int(sing.sum()), singlets_pct_of_cells=100 * sing.sum() / max(cells.sum(), 1),
                   ratio_median=round(i2["median_ratio"], 4))
        rows = np.flatnonzero(sing)
    live_ch = find_live_channel(sample)
    if live_ch and rows.size > 100:
        vals = sample.column(live_ch)[rows]
        tl = Transform("arcsinh", {"cofactor": estimate_cofactor(vals), "T": sample.max_range(live_ch)})
        g3, i3 = valley_threshold(sample, live_ch, tl, rows=rows, keep="below")
        out.update(live_channel=live_ch, live_label=sample.labels[sample.index(live_ch)],
                   live_cut=round(i3["cut"], 2), live_modes=i3["modes"], live_pct=100 * i3["fraction"],
                   live_cofactor=round(tl.params["cofactor"], 3))
    return out


def label_consistency(inventory_channels):
    """Je $PnN die häufigste $PnS-Beschriftung; Abweichungen auflisten."""
    by_ch = collections.defaultdict(collections.Counter)
    for name, chans in inventory_channels.items():
        for c, l in chans:
            by_ch[c][l] += 1
    rows = []
    for name, chans in inventory_channels.items():
        for c, l in chans:
            modal, n = by_ch[c].most_common(1)[0]
            if l != modal:
                rows.append({"file": name, "channel": c, "label": l, "majority_label": modal, "majority_n": n})
    return rows


def _files(folders, pattern):
    seen, out = set(), []
    for f in folders:
        f = pathlib.Path(f)
        cands = [f] if f.is_file() else f.rglob("*")
        for p in cands:
            if p.is_file() and p.suffix.lower() == ".fcs" and re.search(pattern, p.name) and p not in seen:
                seen.add(p)
                out.append(p)
    return sorted(out)


def run(folders, out, gates=None, wsp=(), pattern=".", max_files=None, checksum=True,
        stats_channels=None, log=print, apply_gain=True):
    folders = [pathlib.Path(f).expanduser().resolve() for f in folders]
    out = pathlib.Path(out).expanduser().resolve()
    for f in folders:
        if f.is_dir() and (out == f or f in out.parents):
            raise SystemExit(f"Ausgabeordner {out} liegt im Datenordner {f} – bitte einen eigenen Ordner wählen")
    out.mkdir(parents=True, exist_ok=True)
    logf = open(out / "protokoll.txt", "a", encoding="utf-8")

    def say(msg):
        line = f"{dt.datetime.now():%H:%M:%S}  {msg}"
        log(line)
        logf.write(line + "\n")
        logf.flush()

    say(f"FlowCy {__version__} · Python {platform.python_version()} · {platform.platform()} · "
        f"$PnG {'angewandt' if apply_gain else 'nicht angewandt'}")
    files = _files(folders, pattern)
    if max_files:
        files = files[:max_files]
    say(f"{len(files)} FCS-Dateien gefunden")
    tree = GatingTree.load(gates) if gates else None
    inv, qcs, autos, stats, errors, chans = [], [], [], [], [], {}
    for i, p in enumerate(files, 1):
        try:
            t0 = time.perf_counter()
            s = load_sample(p, checksum=checksum, apply_gain=apply_gain)
            t_load = time.perf_counter() - t0
            sm = s.summary()
            inv.append({"file": p.name, "path": str(p), "events": s.n_events, "channels": len(s.channels),
                        "cytometer": sm["cytometer"], "date": sm["date"], "fcs_version": sm["fcs_version"],
                        "has_spill": sm["has_spill"], "sha256": sm["sha256"], "size_mb": p.stat().st_size / 1e6,
                        "gain_channels": ", ".join(f"{s.channels[i]}×{g:g}" for i, g in gains(s.metadata, len(s.channels)).items()),
                        "load_s": round(t_load, 3)})
            chans[p.name] = list(zip(s.channels, s.labels))
            q = run_qc(s)
            fl = q["flow"]
            qcs.append({"file": p.name, "status": q["status"], "issues": " | ".join(q["issues"]),
                        "flow_available": fl.get("available"), "flow_segments": len(fl.get("segments", [])),
                        "flow_time_pct": 100 * fl.get("affected_time_fraction", 0) if fl.get("available") else None,
                        "flow_event_pct": 100 * fl.get("affected_fraction", 0) if fl.get("available") else None,
                        **{f"margin_{k}": round(100 * v, 3) for k, v in q["margin"].items()},
                        "_rate": fl.get("rate")})
            t1 = time.perf_counter()
            try:
                a = auto_gating(s)
            except Exception as e:  # noqa: BLE001
                a = {"sample": p.name, "note": f"Auto-Gating fehlgeschlagen: {e}"}
            a["auto_s"] = round(time.perf_counter() - t1, 3)
            autos.append(a)
            if tree is not None:
                t2 = time.perf_counter()
                rows = tree.statistics(s, channels=[c for c in (stats_channels or []) if c in s.channels])
                for r in rows:
                    r["gating_s"] = round(time.perf_counter() - t2, 3)
                stats.extend(rows)
            say(f"[{i}/{len(files)}] {p.name}: {s.n_events:,} Ereignisse, QC {q['status']}, "
                f"geladen in {t_load:.1f}s")
            del s
        except Exception as e:  # noqa: BLE001
            errors.append({"file": p.name, "error": f"{type(e).__name__}: {e}"})
            say(f"[{i}/{len(files)}] FEHLER {p.name}: {e}")
            logf.write(traceback.format_exc() + "\n")

    conc = []
    wsp_files = [pathlib.Path(w) for w in wsp] or [w for f in folders if f.is_dir() for w in f.rglob("*.wsp")]
    for w in wsp_files:
        say(f"FlowJo-Workspace: {w}")
        try:
            from .flowjo import concordance
            for r in concordance(w, _fcs_dir_for(w, files)):
                r["workspace"] = str(w)
                conc.append(r)
        except Exception as e:  # noqa: BLE001
            errors.append({"file": str(w), "error": f"{type(e).__name__}: {e}"})
            say(f"FEHLER Workspace {w}: {e}")

    def save(rows, name):
        if rows:
            pd.DataFrame(rows).drop(columns=[c for c in ("_rate",) if c in pd.DataFrame(rows).columns]) \
                .to_csv(out / name, sep="\t", index=False)
    save(inv, "inventar.tsv")
    save(qcs, "qc.tsv")
    labels = label_consistency(chans)
    save(labels, "kanalbeschriftung.tsv")
    save(autos, "auto_gating.tsv")
    save(stats, "statistik.tsv")
    save(conc, "konkordanz_flowjo.tsv")
    save(errors, "fehler.tsv")
    summary = {"version": __version__, "files": len(files), "loaded": len(inv), "errors": len(errors),
               "qc_warn": sum(q["status"] == "warn" for q in qcs), "label_issues": len(labels),
               "workspaces": [str(w) for w in wsp_files], "finished": dt.datetime.now().isoformat(timespec="seconds")}
    if conc:
        c = pd.DataFrame(conc)
        summary["concordance"] = c["status"].value_counts().to_dict()
        if "flowcy_eq_flowkit" in c:
            summary["flowcy_eq_flowkit"] = c["flowcy_eq_flowkit"].value_counts().to_dict()
    (out / "zusammenfassung.json").write_text(json.dumps(summary, ensure_ascii=False, indent=1, default=str),
                                             encoding="utf-8")
    write_html(out / "bericht.html", summary, inv, qcs, autos, labels, conc, errors)
    say(f"Fertig: {out}")
    logf.close()
    return summary


def _fcs_dir_for(wsp_path, files):
    """FCS-Ordner zum Workspace: der Ordner, der die meisten der gelisteten Dateien enthält."""
    from .flowjo import reference_counts
    names = set(reference_counts(wsp_path))
    hits = collections.Counter(p.parent for p in files if p.name in names)
    return hits.most_common(1)[0][0] if hits else wsp_path.parent


# --- HTML-Bericht ----------------------------------------------------------------

def _spark(rate, w=160, h=28):
    if not rate:
        return ""
    r = np.asarray(rate, float)
    top = r.max() or 1
    pts = " ".join(f"{i * w / (len(r) - 1):.1f},{h - 2 - (v / top) * (h - 4):.1f}" for i, v in enumerate(r))
    return (f'<svg width="{w}" height="{h}" role="img" aria-label="Flussrate über die Zeit">'
            f'<polyline points="{pts}" fill="none" stroke="currentColor" stroke-width="1.2"/></svg>')


def _table(rows, cols, fmt=None):
    fmt = fmt or {}
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    body = []
    for r in rows:
        tds = []
        for c in cols:
            v = r.get(c)
            if c in fmt:
                v = fmt[c](v)
            elif isinstance(v, float):
                v = "" if np.isnan(v) else f"{v:,.2f}".replace(",", " ")
            elif isinstance(v, int):
                v = f"{v:,}".replace(",", " ")
            tds.append(f"<td>{v if isinstance(v, str) and v.startswith('<svg') else html.escape(str(v if v is not None else ''))}</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def write_html(path, summary, inv, qcs, autos, labels, conc, errors):
    css = """
    :root{--bg:#F6F7F9;--fg:#16181D;--muted:#5B6270;--line:#D9DCE1;--card:#fff;--warn:#92400E;--warnbg:#FEF3C7;--ok:#0F766E}
    @media (prefers-color-scheme: dark){:root{--bg:#14161A;--fg:#E8EAED;--muted:#A3A9B4;--line:#30343B;--card:#1C1F24;--warn:#FBBF24;--warnbg:#3A2A0A;--ok:#5EEAD4}}
    body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
    main{max-width:1200px;margin:0 auto;padding:24px 16px}
    h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 8px}
    .muted{color:var(--muted)}.cards{display:flex;flex-wrap:wrap;gap:12px;margin:16px 0}
    .card{flex:1 1 160px;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
    .card b{display:block;font-size:22px;font-variant-numeric:tabular-nums}
    .scroll{overflow-x:auto;background:var(--card);border:1px solid var(--line);border-radius:10px}
    table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums}
    th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
    th{color:var(--muted);font-weight:500}.warn{color:var(--warn)}
    """
    warn_q = [q for q in qcs if q["status"] != "ok"]
    for q in warn_q:
        q["Fluss"] = _spark(q.get("_rate"))
    cards = [("Dateien", summary["files"]), ("gelesen", summary["loaded"]), ("Fehler", summary["errors"]),
             ("QC-Warnungen", summary["qc_warn"]), ("Beschriftung abweichend", summary["label_issues"])]
    parts = [f"<!doctype html><html lang='de'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>FlowCy Testbericht</title><style>{css}</style></head><body><main>",
             f"<h1>FlowCy Testbericht</h1><p class='muted'>Version {summary['version']} · {summary['finished']}</p>",
             "<div class='cards'>" + "".join(f"<div class='card'><span class='muted'>{k}</span><b>{v}</b></div>" for k, v in cards) + "</div>"]
    if conc:
        c = summary.get("concordance", {})
        parts.append("<h2>FlowJo-Konkordanz</h2><p>" + " · ".join(f"{k}: {v}" for k, v in c.items()) +
                     (f" · FlowCy = FlowKit: {summary.get('flowcy_eq_flowkit')}" if summary.get("flowcy_eq_flowkit") else "") + "</p>")
        bad = [r for r in conc if r.get("status") not in ("gleich",)]
        parts.append(_table(bad[:300], ["sample", "gate", "flowjo", "flowcy", "flowkit", "rel_diff", "status"],
                            {"rel_diff": lambda v: "" if v is None else f"{100 * v:.3f} %"}))
    parts.append("<h2>QC-Auffälligkeiten</h2>")
    parts.append(_table(warn_q, ["file", "status", "issues", "Fluss"]) if warn_q else "<p>Keine.</p>")
    if labels:
        parts.append("<h2>Abweichende Kanalbeschriftung</h2>")
        parts.append(_table(labels, ["file", "channel", "label", "majority_label", "majority_n"]))
    if autos:
        parts.append("<h2>Automatische Plausibilitäts-Strategie</h2><p class='muted'>Kein Ersatz für ein "
                     "gezeichnetes Gating – dient zum Erkennen auffälliger Proben.</p>")
        cols = [c for c in ["sample", "fsc_scale", "fsc_cut", "cells_pct", "singlets_pct_of_cells",
                            "live_label", "live_cut", "live_pct", "note"] if any(c in a for a in autos)]
        parts.append(_table(autos, cols))
    parts.append("<h2>Inventar</h2>")
    parts.append(_table(inv, ["file", "events", "channels", "cytometer", "date", "has_spill", "size_mb", "load_s"]))
    if errors:
        parts.append("<h2 class='warn'>Fehler</h2>" + _table(errors, ["file", "error"]))
    parts.append("</main></body></html>")
    pathlib.Path(path).write_text("".join(parts), encoding="utf-8")
