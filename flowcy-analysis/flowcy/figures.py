"""Abbildungen für Publikationen: ein Plot je Probe, gleiche Achsen, nach Gruppen geordnet.

Stil nach Nature-Vorgaben: Helvetica/Arial 7 pt, Linien 0,5 pt, Gate-Linien schwarz,
Dichte in Pseudocolor wie FlowJo. Ausgabe als PDF (Vektor, Dichte als Rasterbild) und PNG.
"""
from __future__ import annotations

import pathlib

import numpy as np

PSEUDO = ["#1414C8", "#0A64F0", "#00B4F0", "#00D278", "#B4E600", "#FFD200", "#FF7800", "#E60000"]
VIRIDIS = ["#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725"]
GRAY = ["#BDBDBD", "#7A7A7A", "#3A3A3A", "#000000"]
MM = 1 / 25.4


def _style():
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2, "axes.labelsize": 7, "xtick.labelsize": 6,
        "ytick.labelsize": 6, "pdf.fonttype": 42, "svg.fonttype": "none",
    })


def _cmap(name):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, {"pseudo": PSEUDO, "viridis": VIRIDIS, "gray": GRAY}[name])


def gate_outline(g, tx, ty):
    """Gate-Umriss in Anzeigeskala der Plotachsen (Liste von Linienzügen)."""
    if g.kind == "polygon":
        v = np.asarray(g.vertices, float)
        dv = np.column_stack([g.xform.forward(v[:, 0]), g.yform.forward(v[:, 1])])
        pts = []
        for i in range(len(dv)):
            a, b = dv[i], dv[(i + 1) % len(dv)]
            pts.extend(a + (b - a) * t for t in np.linspace(0, 1, 12, endpoint=False))
        dv = np.asarray(pts + [pts[0]])
        raw = np.column_stack([g.xform.inverse(dv[:, 0]), g.yform.inverse(dv[:, 1])])
        return [np.column_stack([tx.forward(raw[:, 0]), ty.forward(raw[:, 1])])]
    if g.kind == "rectangle" and g.y:
        big = 1e12
        x0 = tx.forward(g.xmin) if g.xmin is not None else -big
        x1 = tx.forward(g.xmax) if g.xmax is not None else big
        y0 = ty.forward(g.ymin) if g.ymin is not None else -big
        y1 = ty.forward(g.ymax) if g.ymax is not None else big
        return [np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]], float)]
    if g.kind == "quadrant":
        xd, yd = float(tx.forward(g.xdiv)), float(ty.forward(g.ydiv))
        return [np.array([[xd, -1e12], [xd, 1e12]]), np.array([[-1e12, yd], [1e12, yd]])]
    return []


def sample_grid(project, node, x, y, samples=None, out=None, palette="pseudo", bins=200, panel_mm=32,
                title=None, show_pct=True):
    """Ein Panel je Probe; Zeilen = Gruppen. Gibt die Pfade der PDF- und PNG-Datei zurück."""
    _style()
    import matplotlib.pyplot as plt
    from scipy import ndimage as ndi

    names = samples or [n for n, e in project.samples.items() if not e.get("error")]
    groups = {}
    for n in names:
        groups.setdefault(project.samples[n].get("group") or "ohne Gruppe", []).append(n)
    ncol = max(len(v) for v in groups.values())
    nrow = len(groups)
    tx, xlo, xhi = project.axis(x)
    ty, ylo, yhi = project.axis(y)
    cmap = _cmap(palette)
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * panel_mm * MM, nrow * panel_mm * MM * 1.05),
                             squeeze=False, constrained_layout=True)
    labels = {c["name"]: c["display"] for c in project.channels()}
    xt = [(p, l) for _, p, l in tx.ticks(xlo, xhi) if xlo <= p <= xhi]
    yt = [(p, l) for _, p, l in ty.ticks(ylo, yhi) if ylo <= p <= yhi]
    for r, (grp, members) in enumerate(groups.items()):
        for c in range(ncol):
            ax = axes[r][c]
            if c >= len(members):
                ax.axis("off")
                continue
            name = members[c]
            s = project.sample(name)
            tree = project.tree_of(name, activate=False)
            masks = tree.masks(s)
            nid = node if node in tree.nodes else "root"
            rows = masks[nid]
            dx = tx.forward(s.column(x)[rows])
            dy = ty.forward(s.column(y)[rows])
            ok = np.isfinite(dx) & np.isfinite(dy)
            h, _, _ = np.histogram2d(dx[ok], dy[ok], bins=bins, range=[[xlo, xhi], [ylo, yhi]])
            h = ndi.gaussian_filter(h, 0.6)
            img = np.log1p(h).T
            img = np.ma.masked_where(h.T <= 0.05, img)
            ax.imshow(img, origin="lower", extent=[xlo, xhi, ylo, yhi], aspect="auto", cmap=cmap,
                      interpolation="nearest", rasterized=True)
            n_parent = int(rows.sum())
            for cid in tree.children(nid):
                g = tree.gate_for(cid, name)
                if g is None or g.kind == "boolean" or g.x != x or g.y != y:
                    continue
                for line in gate_outline(g, tx, ty):
                    ax.plot(line[:, 0], line[:, 1], color="black", lw=0.6)
                if show_pct and n_parent:
                    pct = 100 * masks[cid].sum() / n_parent
                    if g.kind == "quadrant":
                        pos = {"-+": (0.03, 0.97, "left", "top"), "++": (0.97, 0.97, "right", "top"),
                               "--": (0.03, 0.03, "left", "bottom"), "+-": (0.97, 0.03, "right", "bottom")}[g.quadrant]
                        ax.text(pos[0], pos[1], f"{pct:.1f}".replace(".", ","), transform=ax.transAxes,
                                ha=pos[2], va=pos[3], fontsize=6)
                    else:
                        ln = gate_outline(g, tx, ty)[0]
                        ln = ln[(ln[:, 0] > xlo) & (ln[:, 0] < xhi) & (ln[:, 1] > ylo) & (ln[:, 1] < yhi)]
                        if len(ln):
                            top = ln[np.argmax(ln[:, 1])]
                            ax.text(top[0], min(top[1] + 0.02 * (yhi - ylo), yhi), f"{pct:.1f}".replace(".", ","),
                                    ha="center", va="bottom", fontsize=6)
            ax.set_xlim(xlo, xhi)
            ax.set_ylim(ylo, yhi)
            ax.set_xticks([p for p, _ in xt], [l for _, l in xt])
            ax.set_yticks([p for p, _ in yt], [l for _, l in yt])
            if r < nrow - 1:
                ax.set_xticklabels([])
            if c > 0:
                ax.set_yticklabels([])
            ax.set_title(name.rsplit(".", 1)[0], fontsize=6, pad=2)
            if c == 0:
                ax.set_ylabel(f"{grp}\n{labels.get(y, y)}")
            if r == nrow - 1:
                ax.set_xlabel(labels.get(x, x))
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
    if title:
        fig.suptitle(title, fontsize=7, fontweight="bold")
    out = pathlib.Path(out or (pathlib.Path.home() / "FlowCy" / "Export" / f"{_slug(project.name)}_{_slug(x)}_{_slug(y)}"))
    out.parent.mkdir(parents=True, exist_ok=True)
    pdf, png = out.with_suffix(".pdf"), out.with_suffix(".png")
    fig.savefig(pdf, dpi=600)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    return {"pdf": str(pdf), "png": str(png), "rows": nrow, "cols": ncol}


def _slug(s):
    import re
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(s)).strip("_") or "x"
