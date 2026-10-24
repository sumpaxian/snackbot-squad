"""Synthetische FCS-Dateien mit bekannter Wahrheit – für Tests und die Demo.

Jedes Ereignis trägt ein Etikett (Debris, Dublette, tot, CD4, CD8, DN, DP, CD3neg), so dass
Gating-Ergebnisse gegen die Wahrheit geprüft werden können. Die Fluoreszenz wird mit einer
bekannten Spillover-Matrix verschmiert und die Matrix als $SPILLOVER in die Datei geschrieben.
"""
from __future__ import annotations

import pathlib
import zlib
import numpy as np

CHANNELS = ["FSC-A", "FSC-H", "SSC-A", "CD3", "CD4", "CD8", "LD", "Time"]
LABELS = ["", "", "", "CD3 BV421", "CD4 FITC", "CD8 PE", "Zombie NIR", ""]
FLUOR = ["CD3", "CD4", "CD8", "LD"]
SPILL = np.array([
    [1.00, 0.05, 0.00, 0.00],
    [0.02, 1.00, 0.22, 0.00],
    [0.00, 0.03, 1.00, 0.04],
    [0.00, 0.00, 0.06, 1.00],
])


def _pos(rng, n, center, sd=0.25):
    return 10 ** rng.normal(center, sd, n)


def make_events(n=60_000, seed=0, cd4_frac=0.65, cd8_frac=0.28, dp_frac=0.02, flow_dip=None):
    rng = np.random.default_rng(seed)
    labels = rng.choice(["debris", "doublet", "cell"], size=n, p=[0.12, 0.06, 0.82]).astype(object)
    fsc = np.empty(n)
    ssc = np.empty(n)
    deb = labels == "debris"
    fsc[deb] = rng.gamma(2.0, 6000, deb.sum())
    ssc[deb] = rng.gamma(2.0, 5000, deb.sum())
    cells = ~deb
    fsc[cells] = rng.normal(80_000, 9_000, cells.sum())
    ssc[cells] = rng.normal(30_000, 6_000, cells.sum())
    ratio = np.where(labels == "doublet", rng.normal(0.55, 0.04, n), rng.normal(0.93, 0.02, n))
    fsch = np.clip(fsc * ratio, 0, None)
    fsc = np.where(labels == "doublet", fsc * 1.8, fsc)

    true = np.zeros((n, 4))
    idx = np.flatnonzero(labels == "cell")
    dead = rng.random(idx.size) < 0.10
    labels[idx[dead]] = "dead"
    live = idx[~dead]
    cd3pos = rng.random(live.size) < 0.70
    labels[live[~cd3pos]] = "CD3neg"
    t = live[cd3pos]
    u = rng.random(t.size)
    sub = np.where(u < cd4_frac, "CD4", np.where(u < cd4_frac + cd8_frac, "CD8",
                   np.where(u < cd4_frac + cd8_frac + dp_frac, "DP", "DN")))
    labels[t] = sub
    # Hintergrund: Autofluoreszenz um 0 mit Streuung
    true[:] = rng.normal(0, 120, (n, 4))
    true[t, 0] = _pos(rng, t.size, 3.9)
    for name, c4, c8 in (("CD4", True, False), ("CD8", False, True), ("DP", True, True)):
        m = labels == name
        if c4:
            true[m, 1] = _pos(rng, m.sum(), 3.7)
        if c8:
            true[m, 2] = _pos(rng, m.sum(), 4.0)
    m = labels == "dead"
    true[m, 3] = _pos(rng, m.sum(), 4.2)
    obs = true @ SPILL

    time = np.sort(rng.uniform(0, 300.0, n))
    if flow_dip:
        a, b = flow_dip
        keep = ~((time > a) & (time < b) & (rng.random(n) < 0.8))
        time, fsc, fsch, ssc, obs, labels = time[keep], fsc[keep], fsch[keep], ssc[keep], obs[keep], labels[keep]
    data = np.column_stack([fsc, fsch, ssc, obs, time])
    data[:, :3] = np.clip(data[:, :3], 0, 262143)
    return data, np.asarray(labels, dtype=str)


def spill_text():
    vals = ",".join(f"{v:g}" for v in SPILL.ravel())
    return f"{len(FLUOR)}," + ",".join(FLUOR) + "," + vals


def write_fcs(path, data, metadata=None):
    import flowio
    md = {"SPILLOVER": spill_text(), "CYT": "FlowCy Synthetik", "DATE": "14-OCT-2026"}
    md.update(metadata or {})
    for i in range(len(CHANNELS)):
        md[f"P{i + 1}R"] = "262144" if CHANNELS[i] != "Time" else "1024"
    path = pathlib.Path(path)
    with open(path, "wb") as fh:
        flowio.create_fcs(fh, data.astype(np.float32).ravel(), CHANNELS, opt_channel_names=LABELS, metadata_dict=md)
    return path


def make_study(folder, n_per_group=4, n=40_000):
    """Kontrolle vs. Behandlung (CD4↓, CD8↑), eine Probe mit Flusseinbruch."""
    folder = pathlib.Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    truth = {}
    for g, (c4, c8) in (("K", (0.66, 0.28)), ("B", (0.58, 0.35))):
        for i in range(1, n_per_group + 1):
            dip = (120, 170) if (g == "B" and i == 3) else None
            data, lab = make_events(n, seed=zlib.crc32(f"{g}{i}".encode()), cd4_frac=c4, cd8_frac=c8, flow_dip=dip)
            name = f"{g}{i:02d}_PBMC.fcs"
            write_fcs(folder / name, data)
            truth[name] = lab
    return truth
