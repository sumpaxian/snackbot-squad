import numpy as np
import pytest
from fastapi.testclient import TestClient

from flowcy import server
from flowcy.project import Project


@pytest.fixture()
def client(study):
    folder, _ = study
    server.S.project = Project()
    c = TestClient(server.app)
    r = c.post("/api/open", json={"folder": str(folder)})
    assert r.status_code == 200, r.text
    return c


def test_open_and_samples(client):
    info = client.get("/api/info").json()
    assert info["samples"] == 8
    names = [c["name"] for c in info["channels"]]
    assert "FSC-A" in names and "CD4" in names
    lab = {c["name"]: c["label"] for c in info["channels"]}
    assert lab["CD4"] == "CD4 FITC"


def test_qc_and_exclusion(client):
    q = client.post("/api/qc/B03_PBMC.fcs").json()
    assert q["status"] == "warn"
    segs = [[s["start"], s["end"]] for s in q["flow"]["segments"]]
    before = client.get("/api/tree", params={"sample": "B03_PBMC.fcs"}).json()["nodes"][0]["count"]
    client.post("/api/exclude", json={"sample": "B03_PBMC.fcs", "segments": segs})
    after = client.get("/api/tree", params={"sample": "B03_PBMC.fcs"}).json()["nodes"][0]["count"]
    assert after < before


def test_groups_regex(client):
    r = client.post("/api/groups/regex", json={"regex": r"^([KB])\d", "mapping": {"K": "Kontrolle", "B": "Behandlung"}}).json()
    assert r["assigned"] == 8
    groups = {s["group"] for s in client.get("/api/samples").json()}
    assert groups == {"Kontrolle", "Behandlung"}


def test_gate_roundtrip_from_display(client):
    """Gate in Anzeigekoordinaten anlegen → Plot liefert dieselbe Form zurück."""
    plot = client.post("/api/plot", json={"sample": "K01_PBMC.fcs", "x": "FSC-A", "y": "SSC-A"}).json()
    X, Y = plot["x"], plot["y"]
    norm = [(0.2, 0.05), (0.5, 0.05), (0.5, 0.3), (0.2, 0.3)]
    disp = [[X["lo"] + a * (X["hi"] - X["lo"]), Y["lo"] + b * (Y["hi"] - Y["lo"])] for a, b in norm]
    nid = client.post("/api/gates", json={"parent": "root", "name": "Zellen", "sample": "K01_PBMC.fcs",
                                          "spec": {"kind": "polygon", "x": "FSC-A", "y": "SSC-A", "vertices": disp}}).json()["id"]
    plot = client.post("/api/plot", json={"sample": "K01_PBMC.fcs", "x": "FSC-A", "y": "SSC-A"}).json()
    g = next(g for g in plot["gates"] if g["id"] == nid)
    np.testing.assert_allclose(g["shape"]["points"], norm, atol=1e-9)
    assert g["shape"]["editable"] and 0 < g["pct"] < 100
    tree = client.get("/api/tree", params={"sample": "K01_PBMC.fcs"}).json()
    node = next(n for n in tree["nodes"] if n["id"] == nid)
    assert node["count"] == g["count"]


def test_quadrants_move_together_and_override(client):
    r = client.post("/api/gates", json={"parent": "root", "sample": "K01_PBMC.fcs",
                                        "spec": {"kind": "quadrant", "x": "CD4", "y": "CD8", "xdiv": 0.3, "ydiv": 0.3}}).json()
    assert len(r["ids"]) == 4
    tree = client.get("/api/tree", params={"sample": "K01_PBMC.fcs"}).json()
    qs = [n for n in tree["nodes"] if n["id"] in r["ids"]]
    assert sum(n["count"] for n in qs) == tree["nodes"][0]["count"]
    client.put(f"/api/gates/{r['ids'][0]}", json={"spec": {"kind": "quadrant", "x": "CD4", "y": "CD8", "xdiv": 0.5, "ydiv": 0.5},
                                                   "scope": "sample", "sample": "K01_PBMC.fcs"})
    t1 = client.get("/api/tree", params={"sample": "K01_PBMC.fcs"}).json()
    t2 = client.get("/api/tree", params={"sample": "K02_PBMC.fcs"}).json()
    assert all(n["override"] for n in t1["nodes"] if n["id"] in r["ids"])
    assert not any(n["override"] for n in t2["nodes"] if n["id"] in r["ids"])
    assert sum(n["count"] for n in t1["nodes"] if n["id"] in r["ids"]) == t1["nodes"][0]["count"]


def test_proposals_and_compare(client):
    client.post("/api/groups/regex", json={"regex": r"^([KB])\d", "mapping": {"K": "Kontrolle", "B": "Behandlung"}})
    p = client.post("/api/propose", json={"sample": "K01_PBMC.fcs", "kind": "population", "x": "FSC-A", "y": "SSC-A"}).json()
    nid = client.post("/api/gates", json={"parent": "root", "name": "Zellen", "spec": p["spec"]}).json()["id"]
    p2 = client.post("/api/propose", json={"sample": "K01_PBMC.fcs", "kind": "valley", "node": nid, "x": "CD3", "keep": "above"}).json()
    assert p2["shape"]["type"] == "range"
    cd3 = client.post("/api/gates", json={"parent": nid, "name": "CD3+", "spec": p2["spec"]}).json()["id"]
    q = client.post("/api/gates", json={"parent": cd3, "spec": {"kind": "quadrant", "x": "CD4", "y": "CD8",
                                                                 "xdiv": 0.45, "ydiv": 0.45}}).json()["ids"]
    tree = client.get("/api/tree").json()
    cd4 = next(n["id"] for n in tree["nodes"] if n["id"] in q and n["name"].startswith("CD4⁺ CD8⁻"))
    res = client.post("/api/compare", json={"node": cd4, "metric": "pct_parent", "reference": "Kontrolle"}).json()
    assert [g["n"] for g in res["groups"]] == [4, 4]
    t = res["tests"][0]
    # Synthetik: Behandlung hat weniger CD4 (58 % vs 66 %) → alle 4 unter allen 4, exaktes p = 2/70
    assert t["diff_median"] < -4
    assert abs(t["p"] - 2 / 70) < 1e-9
    csv = client.get("/api/export/stats.csv").text
    assert "pct_parent" in csv and "K01_PBMC.fcs" in csv


def test_axis_change_keeps_membership(client):
    plot = client.post("/api/plot", json={"sample": "K01_PBMC.fcs", "x": "CD4", "y": "CD8"}).json()
    X, Y = plot["x"], plot["y"]
    disp = [[X["lo"] + a * (X["hi"] - X["lo"]), Y["lo"] + b * (Y["hi"] - Y["lo"])] for a, b in [(0.5, 0.1), (0.9, 0.1), (0.9, 0.4), (0.5, 0.4)]]
    nid = client.post("/api/gates", json={"parent": "root", "name": "CD4", "sample": "K01_PBMC.fcs",
                                          "spec": {"kind": "polygon", "x": "CD4", "y": "CD8", "vertices": disp}}).json()["id"]
    n0 = next(n for n in client.get("/api/tree", params={"sample": "K01_PBMC.fcs"}).json()["nodes"] if n["id"] == nid)["count"]
    client.put("/api/axes/CD4", json={"kind": "arcsinh"})
    n1 = next(n for n in client.get("/api/tree", params={"sample": "K01_PBMC.fcs"}).json()["nodes"] if n["id"] == nid)["count"]
    assert n0 == n1
    plot = client.post("/api/plot", json={"sample": "K01_PBMC.fcs", "x": "CD4", "y": "CD8"}).json()
    g = next(g for g in plot["gates"] if g["id"] == nid)
    assert g["shape"]["editable"] is False and len(g["shape"]["points"]) > 4


def test_project_save_load(client, tmp_path):
    client.post("/api/gates", json={"parent": "root", "name": "FSC", "spec": {"kind": "rectangle", "x": "FSC-A", "y": None, "xmin": 0.1, "xmax": 0.9}})
    f = client.post("/api/project/save", json={"path": str(tmp_path / "p.flowcy-projekt.json")}).json()["file"]
    server.S.project = Project()
    client.post("/api/project/load", json={"path": f})
    names = [n["name"] for n in client.get("/api/tree").json()["nodes"]]
    assert "FSC" in names and client.get("/api/info").json()["samples"] == 8


def test_panels_separate_strategies(client):
    r = client.post("/api/panels/regex", json={"regex": r"^([KB])", "mapping": {"K": "Panel K", "B": "Panel B"}}).json()
    assert r["assigned"] == 8 and set(r["panels"]) >= {"Panel K", "Panel B"}
    client.post("/api/gates", json={"parent": "root", "name": "nur K", "sample": "K01_PBMC.fcs",
                                    "spec": {"kind": "rectangle", "x": "FSC-A", "y": None, "xmin": 0.1, "xmax": 0.9}})
    tk = [n["name"] for n in client.get("/api/tree", params={"sample": "K02_PBMC.fcs"}).json()["nodes"]]
    tb = [n["name"] for n in client.get("/api/tree", params={"sample": "B01_PBMC.fcs"}).json()["nodes"]]
    assert "nur K" in tk and "nur K" not in tb
    client.post("/api/panels/copy", json={"src": "Panel K", "dst": "Panel B"})
    tb = [n["name"] for n in client.get("/api/tree", params={"sample": "B01_PBMC.fcs"}).json()["nodes"]]
    assert "nur K" in tb
    rows = client.get("/api/results").json()
    assert {r["panel"] for r in rows} == {"Panel K", "Panel B"}


def test_figure_export(client, tmp_path):
    client.post("/api/groups/regex", json={"regex": r"^([KB])\d"})
    plot = client.post("/api/plot", json={"sample": "K01_PBMC.fcs", "x": "CD4", "y": "CD8"}).json()
    client.post("/api/gates", json={"parent": "root", "sample": "K01_PBMC.fcs",
                                    "spec": {"kind": "quadrant", "x": "CD4", "y": "CD8", "xdiv": 0.45, "ydiv": 0.45}})
    r = client.post("/api/export/figure", json={"x": "CD4", "y": "CD8", "out": str(tmp_path / "abb")}).json()
    assert r["rows"] == 2 and r["cols"] == 4
    assert (tmp_path / "abb.pdf").stat().st_size > 10_000 and (tmp_path / "abb.png").exists()
