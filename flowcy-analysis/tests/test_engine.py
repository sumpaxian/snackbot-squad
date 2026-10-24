import numpy as np
import pytest

from flowcy import GatingTree, PolygonGate, RectangleGate, EllipseGate, Transform, load_sample
from flowcy import synth
from flowcy.gates import points_in_polygon, quadrant_set, gate_from_dict
from flowcy.proposals import main_population, singlet_band, valley_threshold
from flowcy.qc import run_qc
from flowcy.transforms import default_transform
from flowcy.tree import robust_cv, ROOT


def test_fcs_roundtrip_and_spill(sample_and_truth):
    s, _ = sample_and_truth
    assert s.channels == synth.CHANNELS
    assert s.labels[3] == "CD3 BV421"
    chans, m = s.spill
    assert chans == synth.FLUOR
    np.testing.assert_allclose(m, synth.SPILL)
    idx = [s.index(c) for c in chans]
    # Kompensation ist die exakte Umkehr der Verschmierung
    np.testing.assert_allclose(s.events()[:, idx] @ m, s.raw[:, idx], rtol=1e-6, atol=1e-3)
    assert len(s.checksum) == 64


@pytest.mark.parametrize("t", [
    Transform("linear", {"min": 0, "max": 262144}),
    Transform("logicle", {"T": 262144, "W": 0.5, "M": 4.5, "A": 0}),
    Transform("asinh", {"T": 262144, "M": 4.5, "A": 0}),
    Transform("flowjo_biex", {}),
])
def test_transform_inverse(t):
    # FlowJos Biex mit negative=0 staucht stark negative Werte, daher dort ab −50
    lo = -50.0 if t.kind == "flowjo_biex" else -500.0
    x = np.array([lo, -10.0, 0.0, 10.0, 500.0, 5000.0, 100000.0])
    np.testing.assert_allclose(t.inverse(t.forward(x)), x, rtol=1e-4, atol=1e-2)
    y = t.forward(np.array([0.0, 262144.0]))
    assert y[0] < y[1]


def test_log_transform():
    t = Transform("log", {"T": 10000, "M": 4})
    np.testing.assert_allclose(t.forward([10000, 1000, 1]), [1.0, 0.75, 0.0])
    np.testing.assert_allclose(t.inverse(t.forward([3.0, 300.0])), [3.0, 300.0])


def test_logicle_ticks_are_ordered():
    t = Transform("logicle", {})
    pos = [p for _, p, _ in t.ticks()]
    assert pos == sorted(pos) and len(pos) >= 4


def test_point_in_polygon_concave():
    vx = np.array([0, 4, 4, 2, 0.0])
    vy = np.array([0, 0, 4, 2, 4.0])      # Einkerbung oben
    px = np.array([1, 3, 2, 2, 5.0])
    py = np.array([1, 1, 3, 1, 1.0])
    assert points_in_polygon(px, py, vx, vy).tolist() == [True, True, False, True, False]


def test_rectangle_and_ellipse(sample_and_truth):
    s, _ = sample_and_truth
    lin = Transform("linear", {"min": 0, "max": 262144})
    r = RectangleGate("FSC-A", "SSC-A", 40_000, 120_000, 10_000, 60_000, xform=lin, yform=lin)
    m = r.evaluate(s)
    x, y = s.column("FSC-A"), s.column("SSC-A")
    assert np.array_equal(m, (x >= 40_000) & (x < 120_000) & (y >= 10_000) & (y < 60_000))
    e = EllipseGate("FSC-A", "SSC-A", (80_000, 30_000), 0.1, 0.05, xform=lin, yform=lin)
    em = e.evaluate(s)
    dx, dy = x / 262144 - 80_000 / 262144, y / 262144 - 30_000 / 262144
    assert np.array_equal(em, (dx / 0.1) ** 2 + (dy / 0.05) ** 2 <= 1)


def build_strategy(s):
    """Strategie wie im Laboralltag: Zellen → Einzelzellen → lebend → CD3+ → CD4/CD8."""
    tree = GatingTree("Test")
    lin = default_transform("FSC-A", 262144)
    cells_g, _ = main_population(s, "FSC-A", "SSC-A", lin, lin, coverage=0.95)
    cells = tree.add("Zellen", cells_g)
    rows = np.flatnonzero(tree.masks(s)[cells])
    sing_g, _ = singlet_band(s, rows=rows)
    sing = tree.add("Einzelzellen", sing_g, cells)
    rows = np.flatnonzero(tree.masks(s)[sing])
    ld_t = default_transform("LD", 262144, s.column("LD"))
    live_g, _ = valley_threshold(s, "LD", ld_t, rows=rows, keep="below")
    live = tree.add("Lebend", live_g, sing)
    rows = np.flatnonzero(tree.masks(s)[live])
    cd3_t = default_transform("CD3", 262144, s.column("CD3"))
    cd3_g, _ = valley_threshold(s, "CD3", cd3_t, rows=rows, keep="above")
    cd3 = tree.add("CD3+", cd3_g, live)
    t4 = default_transform("CD4", 262144, s.column("CD4"))
    t8 = default_transform("CD8", 262144, s.column("CD8"))
    ids = {}
    for name, g in quadrant_set("CD4", "CD8", 1000, 1000, t4, t8,
                                names=["CD8+", "DP", "DN", "CD4+"]):
        ids[name] = tree.add(name, g, cd3)
    return tree, dict(cells=cells, sing=sing, live=live, cd3=cd3, **ids)


def test_gating_matches_truth(sample_and_truth):
    s, lab = sample_and_truth
    tree, ids = build_strategy(s)
    masks = tree.masks(s)

    def purity_recall(nid, truth_labels):
        """Reinheit im Gate; Ausbeute bezogen auf die Elternpopulation (Verluste weiter oben zählen nicht)."""
        m = masks[nid]
        parent = masks[tree.nodes[nid].parent]
        t = np.isin(lab, truth_labels)
        return (m & t).sum() / m.sum(), (m & t).sum() / (parent & t).sum()

    p, r = purity_recall(ids["cells"], ["cell", "dead", "CD3neg", "CD4", "CD8", "DN", "DP", "doublet"])
    assert p > 0.97 and r > 0.85
    p, r = purity_recall(ids["sing"], ["dead", "CD3neg", "CD4", "CD8", "DN", "DP"])
    assert p > 0.98 and r > 0.9
    p, r = purity_recall(ids["live"], ["CD3neg", "CD4", "CD8", "DN", "DP"])
    assert p > 0.98 and r > 0.95
    p, r = purity_recall(ids["cd3"], ["CD4", "CD8", "DN", "DP"])
    assert p > 0.98 and r > 0.95
    p, r = purity_recall(ids["CD4+"], ["CD4"])
    assert p > 0.97 and r > 0.95
    p, r = purity_recall(ids["CD8+"], ["CD8"])
    assert p > 0.97 and r > 0.95


def test_statistics_consistent(sample_and_truth):
    s, _ = sample_and_truth
    tree, ids = build_strategy(s)
    st = {r["id"]: r for r in tree.statistics(s, channels=["CD4"])}
    q = [ids[k] for k in ("CD4+", "CD8+", "DP", "DN")]
    # Quadranten zerlegen die Elternpopulation vollständig
    assert sum(st[i]["count"] for i in q) == st[ids["cd3"]]["count"]
    assert abs(sum(st[i]["pct_parent"] for i in q) - 100) < 1e-9
    assert st[ROOT]["count"] == s.n_events
    assert st[ids["CD4+"]]["median:CD4"] > st[ids["CD8+"]]["median:CD4"]


def test_tree_roundtrip_and_override(tmp_path, sample_and_truth):
    s, _ = sample_and_truth
    tree, ids = build_strategy(s)
    p = tmp_path / "s.flowcy.json"
    tree.save(p)
    t2 = GatingTree.load(p)
    a = [r["count"] for r in tree.statistics(s)]
    b = [r["count"] for r in t2.statistics(s)]
    assert a == b
    # Anpassung nur für diese Probe
    g = gate_from_dict(t2.nodes[ids["live"]].gate.to_dict())
    g.xmax = g.xmax / 4
    t2.set_gate(ids["live"], g, sample=s.name, user="test")
    st = {r["id"]: r for r in t2.statistics(s)}
    assert st[ids["live"]]["override"] and st[ids["live"]]["count"] < a[[i for i, _ in tree.order()].index(ids["live"])]
    t2.reset_override(ids["live"], s.name)
    assert [r["count"] for r in t2.statistics(s)] == a
    assert any("angepasst" in e["action"] for e in t2.log)


def test_polygon_matches_flowkit(sample_and_truth):
    """Unabhängige Referenz: dasselbe Polygon in FlowKit (Gating-ML-Semantik)."""
    import flowkit as fk
    s, _ = sample_and_truth
    T = Transform("logicle", {"T": 262144, "W": 0.5, "M": 4.5, "A": 0})
    verts = [(200, 300), (8000, 200), (60000, 40000), (3000, 90000), (-100, 5000)]
    g = PolygonGate("CD4", "CD8", verts, xform=T, yform=T)
    ours = int(g.evaluate(s, compensated=True).sum())

    gs = fk.GatingStrategy()
    gs.add_transform("lg", fk.transforms.LogicleTransform(param_t=262144, param_w=0.5, param_m=4.5, param_a=0))
    chans, m = s.spill
    gs.add_comp_matrix("spill", fk.Matrix(m, chans, chans))
    dims = [fk.Dimension("CD4", compensation_ref="spill", transformation_ref="lg"),
            fk.Dimension("CD8", compensation_ref="spill", transformation_ref="lg")]
    dv = [[float(T.forward(x)), float(T.forward(y))] for x, y in verts]
    gs.add_gate(fk.gates.PolygonGate("poly", dims, dv), gate_path=("root",))
    fks = fk.Sample(s.path, subsample=0)
    res = gs.gate_sample(fks)
    theirs = int(res.get_gate_count("poly"))
    assert abs(ours - theirs) <= 2, (ours, theirs)


def test_qc_flags_flow_dip(study):
    folder, _ = study
    bad = run_qc(load_sample(folder / "B03_PBMC.fcs"))
    good = run_qc(load_sample(folder / "K01_PBMC.fcs"))
    assert good["status"] == "ok", good["issues"]
    assert bad["status"] == "warn"
    seg = bad["flow"]["segments"]
    assert seg and seg[0]["start"] >= 110 and seg[-1]["end"] <= 180


def test_robust_cv_normal():
    rng = np.random.default_rng(0)
    v = rng.normal(1000, 100, 200_000)
    assert abs(robust_cv(v) - 10.0) < 0.2


def test_disk_cache_identical(study, tmp_path):
    folder, _ = study
    a = load_sample(folder / "K02_PBMC.fcs")
    b = load_sample(folder / "K02_PBMC.fcs", cache_dir=tmp_path)   # schreibt Cache
    c = load_sample(folder / "K02_PBMC.fcs", cache_dir=tmp_path)   # liest Cache
    assert isinstance(c.raw, np.memmap)
    for s in (b, c):
        assert s.channels == a.channels and s.labels == a.labels and s.ranges == a.ranges
        np.testing.assert_array_equal(np.asarray(s.raw), a.raw)
        np.testing.assert_array_equal(s.events(), a.events())
        np.testing.assert_array_equal(s.spill[1], a.spill[1])
    tree, _ = build_strategy(a)
    assert [r["count"] for r in tree.statistics(a)] == [r["count"] for r in tree.statistics(c)]
