import json
import os
import pathlib

import numpy as np
import pandas as pd
import pytest

from flowcy import GatingTree, RectangleGate, Transform, load_sample
from flowcy.gates import BooleanGate, points_in_polygon
from flowcy.importers import tree_from_spec
from flowcy.transforms import estimate_cofactor
from flowcy.validate import run

FK_DATA = os.environ.get("FLOWKIT_DATA")   # Ordner data/ aus github.com/whitews/FlowKit


def test_arcsinh_cofactor():
    t = Transform("arcsinh", {"cofactor": 8.45, "T": 262144})
    x = np.array([-300.0, -5, 0, 5, 184.8, 5000])
    np.testing.assert_allclose(t.inverse(t.forward(x)), x, rtol=1e-9, atol=1e-9)
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 7.0, 100_000)
    assert abs(estimate_cofactor(noise) - 7.0) < 0.2


def test_spec_import_log10_polygon_equals_direct(sample_and_truth):
    """Aus Bildern ausgelesene Polygone (log10-Raum) werden exakt wie dort ausgewertet."""
    s, _ = sample_and_truth
    verts = [[30000, 10000], [140000, 12000], [150000, 70000], [40000, 60000]]
    tree = tree_from_spec([{"name": "Zellen", "parent": None, "x": "FSC-A", "y": "SSC-A",
                            "scale": "log10", "vertices": verts}])
    ours = tree.masks(s)[tree.find("Zellen")]
    x, y = s.column("FSC-A"), s.column("SSC-A")
    ok = (x > 0) & (y > 0)
    v = np.log10(np.array(verts, float))
    direct = np.zeros(s.n_events, bool)
    direct[ok] = points_in_polygon(np.log10(x[ok]), np.log10(y[ok]), v[:, 0], v[:, 1])
    assert np.array_equal(ours, direct)


def test_spec_threshold_child(sample_and_truth):
    s, _ = sample_and_truth
    tree = tree_from_spec([
        {"name": "alle FSC", "parent": None, "x": "FSC-A", "min": 1},
        {"name": "lebend", "parent": "alle FSC", "x": "LD", "scale": "arcsinh", "cofactor": 100, "max": 2000},
    ])
    m = tree.masks(s)
    ld = s.column("LD")
    assert m[tree.find("lebend")].sum() == ((ld < 2000) & (s.column("FSC-A") >= 1)).sum()


def test_boolean_gate(sample_and_truth):
    s, _ = sample_and_truth
    lin = Transform("linear", {"min": 0, "max": 262144})
    tree = GatingTree()
    a = tree.add("A", RectangleGate("FSC-A", None, 50_000, None, xform=lin))
    b = tree.add("B", RectangleGate("SSC-A", None, 30_000, None, xform=lin))
    both = tree.add("A und B", BooleanGate("and", [a, b]))
    either = tree.add("A oder B", BooleanGate("or", [a, b]))
    nota = tree.add("nicht A", BooleanGate("not", [a]))
    m = tree.masks(s)
    assert np.array_equal(m[both], m[a] & m[b])
    assert np.array_equal(m[either], m[a] | m[b])
    assert np.array_equal(m[nota], ~m[a])
    t2 = GatingTree.from_dict(json.loads(json.dumps(tree.to_dict())))
    assert np.array_equal(t2.masks(s)[both], m[both])


def test_validator_is_read_only_and_reports(study, tmp_path):
    folder, _ = study
    before = sorted(p.name for p in pathlib.Path(folder).iterdir())
    out = tmp_path / "bericht"
    summary = run([folder], out, log=lambda *_: None)
    assert sorted(p.name for p in pathlib.Path(folder).iterdir()) == before
    assert summary["files"] == 8 and summary["errors"] == 0 and summary["qc_warn"] == 1
    for f in ("inventar.tsv", "qc.tsv", "auto_gating.tsv", "bericht.html", "protokoll.txt"):
        assert (out / f).exists(), f
    auto = pd.read_csv(out / "auto_gating.tsv", sep="\t")
    # Wahrheit der Synthetik: ~88 % Nicht-Debris, ~90 % lebend unter Einzelzellen
    assert auto["cells_pct"].between(85, 91).all()
    assert auto["live_pct"].between(87, 93).all()
    qc = pd.read_csv(out / "qc.tsv", sep="\t")
    assert qc.loc[qc.file == "B03_PBMC.fcs", "status"].item() == "warn"


def test_validator_refuses_output_inside_data(study):
    folder, _ = study
    with pytest.raises(SystemExit):
        run([folder], pathlib.Path(folder) / "ausgabe", log=lambda *_: None)


@pytest.mark.skipif(not FK_DATA, reason="FLOWKIT_DATA nicht gesetzt (siehe tests/README)")
@pytest.mark.parametrize("wsp,fcs", [
    ("simple_diamond_example/simple_diamond_example_quad_gate.wsp", "simple_diamond_example"),
    ("simple_diamond_example/test_data_diamond_biex_rect.wsp", "simple_diamond_example"),
    ("simple_diamond_example/test_data_diamond_asinh_rect.wsp", "simple_diamond_example"),
    ("simple_line_example/single_ellipse_51_events.wsp", "simple_line_example"),
    ("8_color_data_set/reused_quad_gate_with_child.wsp", "8_color_data_set/fcs_files"),
    ("8_color_data_set/8_color_ICS_boolean_gate_testing.wsp", "8_color_data_set/fcs_files"),
    ("8_color_data_set/8_color_ICS_with_ellipse.wsp", "8_color_data_set/fcs_files"),
])
def test_flowjo_concordance(wsp, fcs):
    from flowcy.flowjo import concordance
    d = pathlib.Path(FK_DATA)
    df = pd.DataFrame(concordance(d / wsp, d / fcs))
    assert not (df.status == "nicht übernommen").any(), df[df.status == "nicht übernommen"]
    # FlowCy rechnet exakt wie FlowKit (unabhängige Gating-ML-Implementierung)
    assert df.flowcy_eq_flowkit.dropna().all()
    # gegen FlowJos gespeicherte Werte: die große Mehrheit gleich oder innerhalb 0,5 %
    ref = df[df.status.isin(["gleich", "≈", "abweichend"])]
    if len(ref):
        assert (ref.status != "abweichend").mean() >= 0.8
