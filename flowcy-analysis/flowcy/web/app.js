"use strict";
/* FlowCy Analysis – Oberfläche. Alle Rechnungen (Transformationen, Gates, Statistik)
   laufen im Python-Backend; hier nur Darstellung und Interaktion. */

const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else if (k === "class") el.className = v;
    else if (k === "html") el.innerHTML = v;
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const c of kids.flat()) if (c !== null && c !== undefined && c !== false) el.append(c.nodeType ? c : document.createTextNode(c));
  return el;
};
const SVGNS = "http://www.w3.org/2000/svg";
const s = (tag, attrs = {}) => {
  const el = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) el.setAttribute(k, v);
  return el;
};
const fmtInt = (n) => (n === null || n === undefined) ? "–" : Math.round(n).toLocaleString("de-DE");
const fmtPct = (p, d = 1) => (p === null || p === undefined || isNaN(p)) ? "–" : p.toLocaleString("de-DE", { minimumFractionDigits: d, maximumFractionDigits: d }) + " %";
const fmtNum = (v, d = 2) => (v === null || v === undefined || isNaN(v)) ? "–" : Number(v).toLocaleString("de-DE", { maximumFractionDigits: d });
const fmtP = (p) => p === null || p === undefined || isNaN(p) ? "–" : (p < 0.001 ? "< 0,001" : p.toLocaleString("de-DE", { maximumFractionDigits: 3 }));

async function api(method, url, body) {
  const r = await fetch(url, { method, headers: body ? { "Content-Type": "application/json" } : {}, body: body ? JSON.stringify(body) : undefined });
  if (!r.ok) {
    let msg = r.statusText;
    try { msg = (await r.json()).detail || msg; } catch (_) { /* Text */ }
    throw new Error(msg);
  }
  return r.json();
}
function toast(msg, ms = 3500) {
  const t = $("#toast"); t.textContent = msg; t.classList.remove("hidden");
  clearTimeout(toast._t); toast._t = setTimeout(() => t.classList.add("hidden"), ms);
}
function fail(e) { console.error(e); toast("Fehler: " + e.message, 6000); }
const store = {
  get(k, d) { try { const v = localStorage.getItem("flowcy:" + k); return v ? JSON.parse(v) : d; } catch (_) { return d; } },
  set(k, v) { try { localStorage.setItem("flowcy:" + k, JSON.stringify(v)); } catch (_) { /* privat */ } },
};

/* ---------- Farbskalen ---------- */
function lut(stops) {
  const rgb = stops.map((c) => [1, 3, 5].map((i) => parseInt(c.slice(i, i + 2), 16)));
  const out = new Uint8ClampedArray(256 * 4);
  for (let v = 1; v < 256; v++) {
    const t = (v - 1) / 254 * (rgb.length - 1);
    const i = Math.min(Math.floor(t), rgb.length - 2), f = t - i;
    for (let k = 0; k < 3; k++) out[v * 4 + k] = rgb[i][k] + (rgb[i + 1][k] - rgb[i][k]) * f;
    out[v * 4 + 3] = 255;
  }
  return out;
}
const PALETTES = {
  pseudo: ["#1414C8", "#0A64F0", "#00B4F0", "#00D278", "#B4E600", "#FFD200", "#FF7800", "#E60000"],
  viridis: ["#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725"],
  gray: ["#BDBDBD", "#7A7A7A", "#3A3A3A", "#000000"],
};
const LUTS = Object.fromEntries(Object.entries(PALETTES).map(([k, v]) => [k, lut(v)]));

/* ---------- Zustand ---------- */
const st = {
  view: "import", samples: [], channels: [], sample: null, tree: null, selNode: "root",
  tool: "edit", scope: "template", palette: store.get("palette", "pseudo"), bins: store.get("bins", 256),
  plots: [], selGate: null, projectName: "", detail: null, cmpNode: null,
};
let plotSeq = 0;

/* ---------- Navigation ---------- */
function show(view) {
  st.view = view;
  $$(".step").forEach((b) => b.classList.toggle("active", b.dataset.view === view));
  $$(".view").forEach((v) => v.classList.toggle("active", v.id === "view-" + view));
  if (view === "gating") renderGating();
  if (view === "compare") renderCompare();
}
$$(".step").forEach((b) => b.addEventListener("click", () => show(b.dataset.view)));

/* ---------- Projekt ---------- */
async function refreshInfo() {
  const info = await api("GET", "/api/info");
  st.channels = info.channels || [];
  st.projectName = info.project || "–";
  $("#projectName").textContent = info.project || "–";
  st.samples = await api("GET", "/api/samples");
  if (!st.sample || !st.samples.find((x) => x.name === st.sample)) st.sample = st.samples.find((x) => !x.error)?.name || null;
  const saved = store.get("plots:" + st.projectName, null);
  if (saved && saved.length) { st.plots = saved; plotSeq = Math.max(...saved.map((p) => p.id)) + 1; }
  else defaultPlots();
  renderImport();
  await loadTree();
}
function defaultPlots() {
  const names = st.channels.map((c) => c.name);
  const pick = (re) => names.find((n) => re.test(n));
  const fa = pick(/^FSC.*A$/i) || names[0], sa = pick(/^SSC.*A$/i) || names[1], fh = pick(/^FSC.*H$/i);
  st.plots = [{ id: plotSeq++, node: "root", x: fa, y: sa }];
  if (fh) st.plots.push({ id: plotSeq++, node: "root", x: fa, y: fh });
  savePlots();
}
function savePlots() { store.set("plots:" + st.projectName, st.plots); }

async function openFolder(path) {
  if (!path) return;
  try {
    toast("Lese Dateiköpfe…");
    await api("POST", "/api/open", { folder: path });
    store.set("lastFolder", path);
    await refreshInfo();
    toast(`${st.samples.length} Proben geladen`);
  } catch (e) { fail(e); }
}
$("#btnOpen").addEventListener("click", async () => {
  if (window.pywebview?.api?.choose_folder) {
    const p = await window.pywebview.api.choose_folder();
    if (p) { $("#folderPath").value = p; openFolder(p); }
  } else { $("#folderPath").focus(); toast("Pfad eingeben und „Öffnen“ wählen"); }
});
$("#btnOpenPath").addEventListener("click", () => openFolder($("#folderPath").value.trim()));
$("#folderPath").addEventListener("keydown", (e) => { if (e.key === "Enter") openFolder(e.target.value.trim()); });
$("#folderPath").value = store.get("lastFolder", "");
$("#btnDemo").addEventListener("click", async () => {
  try { toast("Erzeuge synthetische Beispielstudie…"); await api("POST", "/api/demo", {}); await refreshInfo(); } catch (e) { fail(e); }
});
$("#btnSave").addEventListener("click", async () => {
  try { const r = await api("POST", "/api/project/save", {}); toast("Gespeichert: " + r.file); } catch (e) { fail(e); }
});
$("#btnLoadProject").addEventListener("click", async () => {
  let p = null;
  if (window.pywebview?.api?.choose_file) p = await window.pywebview.api.choose_file("FlowCy-Projekt (*.json)");
  else p = await ask("Projekt laden", "Pfad zur .flowcy-projekt.json", "");
  if (!p) return;
  try { await api("POST", "/api/project/load", { path: p }); await refreshInfo(); toast("Projekt geladen"); } catch (e) { fail(e); }
});
$("#userName").value = store.get("user", "");
$("#userName").addEventListener("change", (e) => { store.set("user", e.target.value); api("POST", "/api/project/user", { user: e.target.value }).catch(fail); });
if ($("#userName").value) api("POST", "/api/project/user", { user: $("#userName").value }).catch(() => {});

/* ---------- Dialog ---------- */
function ask(title, label, value = "", extra = null) {
  return new Promise((resolve) => {
    $("#dlgTitle").textContent = title;
    $("#dlgInput").previousSibling.textContent = label + " ";
    $("#dlgInput").value = value;
    const ex = $("#dlgExtra"); ex.innerHTML = ""; if (extra) ex.append(extra);
    $("#dialog").classList.remove("hidden");
    setTimeout(() => { $("#dlgInput").focus(); $("#dlgInput").select(); }, 0);
    const done = (v) => { $("#dialog").classList.add("hidden"); $("#dlgForm").onsubmit = null; $("#dlgCancel").onclick = null; resolve(v); };
    $("#dlgForm").onsubmit = (e) => { e.preventDefault(); done($("#dlgInput").value.trim()); };
    $("#dlgCancel").onclick = () => done(null);
  });
}
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape" && !$("#dialog").classList.contains("hidden")) $("#dlgCancel").click();
});

/* =================== Import & QC =================== */
function spark(rate, w = 120, hgt = 24) {
  if (!rate || rate.length < 2) return null;
  const top = Math.max(...rate) || 1;
  const pts = rate.map((v, i) => `${(i * w / (rate.length - 1)).toFixed(1)},${(hgt - 2 - v / top * (hgt - 4)).toFixed(1)}`).join(" ");
  const svg = s("svg", { width: w, height: hgt, class: "spark", role: "img", "aria-label": "Flussrate über die Zeit" });
  svg.append(s("polyline", { points: pts, fill: "none", stroke: "currentColor", "stroke-width": 1.2 }));
  return svg;
}
function qcChip(qc) {
  if (!qc) return h("span", { class: "chip" }, "offen");
  if (qc.status === "warn") return h("span", { class: "chip warn", title: qc.issues.join("\n") }, "prüfen");
  if (qc.status === "info") return h("span", { class: "chip", title: qc.issues.join("\n") }, "eingeschränkt");
  return h("span", { class: "chip ok" }, "OK");
}
function renderImport() {
  const n = st.samples.length, ev = st.samples.reduce((a, b) => a + (b.events || 0), 0);
  const warn = st.samples.filter((x) => x.qc?.status === "warn").length;
  const groups = new Set(st.samples.map((x) => x.group).filter(Boolean));
  const err = st.samples.filter((x) => x.error).length;
  $("#importCards").replaceChildren(
    ...[["Proben", n], ["Ereignisse gesamt", fmtInt(ev)], ["Gruppen", groups.size], ["QC prüfen", warn], ["Lesefehler", err]]
      .map(([k, v]) => h("div", { class: "card" }, h("span", { class: "muted small" }, k), h("b", {}, String(v)))));
  $("#sampleCount").textContent = n ? `${n} Dateien` : "";
  const tb = $("#sampleTable tbody");
  tb.replaceChildren(...st.samples.map((x) => {
    const grp = h("input", { value: x.group || "", "aria-label": "Gruppe " + x.name, list: "groupNames" });
    grp.addEventListener("change", async () => {
      try { await api("PUT", `/api/samples/${encodeURIComponent(x.name)}/group`, { group: grp.value }); x.group = grp.value; renderImport(); } catch (e) { fail(e); }
    });
    const pan = h("input", { value: x.panel || "Standard", "aria-label": "Panel " + x.name, list: "panelNames", style: "width:120px" });
    pan.addEventListener("change", async () => {
      try { await api("PUT", `/api/samples/${encodeURIComponent(x.name)}/panel`, { panel: pan.value }); x.panel = pan.value; renderImport(); } catch (e) { fail(e); }
    });
    const qcBtn = h("button", { type: "button", onclick: () => runQc(x.name) }, "QC");
    return h("tr", { class: x.name === st.detail ? "sel" : "" },
      h("td", { class: "file" }, h("button", { type: "button", onclick: () => showDetail(x.name) }, x.name)),
      h("td", { class: "num" }, x.error ? "Fehler" : fmtInt(x.events)),
      h("td", {}, grp),
      h("td", {}, pan),
      h("td", {}, x.error ? h("span", { class: "chip warn", title: x.error }, "nicht lesbar") : qcChip(x.qc)),
      h("td", {}, spark(x.qc?.flow?.rate) || ""),
      h("td", {}, x.error ? "" : qcBtn));
  }));
  let dl = $("#groupNames"); if (!dl) { dl = h("datalist", { id: "groupNames" }); document.body.append(dl); }
  dl.replaceChildren(...[...groups].map((g) => h("option", { value: g })));
  let pl = $("#panelNames"); if (!pl) { pl = h("datalist", { id: "panelNames" }); document.body.append(pl); }
  pl.replaceChildren(...[...new Set(st.samples.map((x) => x.panel || "Standard"))].map((g) => h("option", { value: g })));
}
async function runQc(name) {
  try {
    const q = await api("POST", `/api/qc/${encodeURIComponent(name)}`);
    const x = st.samples.find((y) => y.name === name); x.qc = q;
    renderImport(); if (st.detail === name) showDetail(name);
  } catch (e) { fail(e); }
}
$("#btnQcAll").addEventListener("click", async () => {
  const todo = st.samples.filter((x) => !x.error);
  for (let i = 0; i < todo.length; i++) {
    toast(`QC ${i + 1}/${todo.length}: ${todo[i].name}`);
    try { todo[i].qc = await api("POST", `/api/qc/${encodeURIComponent(todo[i].name)}`); } catch (e) { fail(e); }
    renderImport();
  }
  toast("QC abgeschlossen");
});
async function showDetail(name) {
  st.detail = name; renderImport();
  const box = $("#sampleDetail");
  box.replaceChildren(h("p", { class: "muted" }, "Lade…"));
  try {
    const d = await api("GET", `/api/samples/${encodeURIComponent(name)}`);
    const x = st.samples.find((y) => y.name === name);
    const parts = [h("span", { class: "muted small" }, "Ausgewählte Datei"), h("h2", { style: "font-family:var(--mono)" }, name)];
    const fl = x.qc?.flow;
    if (fl?.available) {
      parts.push(h("h3", { class: "cap" }, "Flussrate über Zeit"));
      const W = 300, H = 100, top = Math.max(...fl.rate) || 1, t0 = fl.time[0], t1 = fl.time[fl.time.length - 1];
      const svg = s("svg", { viewBox: `0 0 ${W} ${H}`, preserveAspectRatio: "none", class: "flowchart", role: "img", "aria-label": "Flussrate" });
      for (const seg of fl.segments) {
        const a = (seg.start - t0) / (t1 - t0) * W, b = (seg.end - t0) / (t1 - t0) * W;
        svg.append(s("rect", { x: a, y: 0, width: Math.max(b - a, 1), height: H, fill: "#FDE7C2" }));
      }
      svg.append(s("polyline", { points: fl.rate.map((v, i) => `${i * W / (fl.rate.length - 1)},${H - 4 - v / top * (H - 10)}`).join(" "), fill: "none", stroke: "var(--accent)", "stroke-width": 1.5, "vector-effect": "non-scaling-stroke" }));
      parts.push(svg);
      if (fl.segments.length) {
        const excluded = (d.excluded_time || []).length > 0;
        parts.push(h("p", { class: "small" }, `${fl.segments.length} auffällige Abschnitte, ${fmtPct(100 * fl.affected_time_fraction, 0)} der Messzeit, ${fmtPct(100 * fl.affected_fraction)} der Ereignisse.`));
        parts.push(h("div", { class: "row" },
          h("button", { type: "button", class: excluded ? "" : "primary", onclick: async () => {
            await api("POST", "/api/exclude", { sample: name, segments: excluded ? [] : fl.segments.map((g) => [g.start, g.end]) });
            toast(excluded ? "Ausschluss aufgehoben" : "Abschnitte ausgeschlossen – gilt für alle Gates dieser Probe"); showDetail(name);
          } }, excluded ? "Ausschluss aufheben" : "Abschnitte ausschließen")));
      }
    } else if (x.qc) parts.push(h("p", { class: "small muted" }, x.qc.issues.join(" · ")));
    else parts.push(h("button", { type: "button", onclick: () => runQc(name) }, "QC ausführen"));
    if (x.qc?.issues?.length) parts.push(h("ul", { class: "small" }, ...x.qc.issues.map((i) => h("li", {}, i))));
    parts.push(h("h3", { class: "cap" }, "Metadaten"));
    const kv = h("dl", { class: "kv" });
    const md = { "$CYT": d.cytometer, "$DATE": d.date, "$BTIM": d.start, "$ETIM": d.end, "$TOT": fmtInt(d.events), "FCS": d.fcs_version, "$VOL": d.metadata.vol ? d.metadata.vol + " nL" : null, "Kanäle": d.channels.length };
    for (const [k, v] of Object.entries(md)) if (v) kv.append(h("dt", {}, k), h("dd", {}, String(v)));
    parts.push(kv);
    parts.push(h("h3", { class: "cap" }, "Kanäle"));
    parts.push(h("div", { class: "scroll", style: "max-height:220px" }, h("table", { class: "data" },
      h("thead", {}, h("tr", {}, h("th", {}, "$PnN"), h("th", {}, "$PnS"), h("th", { class: "num" }, "$PnR"))),
      h("tbody", {}, ...d.channels.map((c) => h("tr", {}, h("td", {}, c.name), h("td", {}, c.label || ""), h("td", { class: "num" }, fmtInt(c.range))))))));
    if (d.spill) {
      parts.push(h("h3", { class: "cap" }, "Spillover-Matrix (Datei)"));
      const t = h("table", { class: "matrix" });
      t.append(h("tr", {}, h("th", {}), ...d.spill.channels.map((c, i) => h("th", { title: d.spill.labels[i] }, c))));
      d.spill.matrix.forEach((row, i) => t.append(h("tr", {}, h("th", { title: d.spill.labels[i] }, d.spill.channels[i]),
        ...row.map((v, j) => {
          const pct = v * 100, bg = i === j ? "var(--soft)" : pct >= 10 ? "#F5B97A" : pct >= 3 ? "#FBDDB8" : pct > 0 ? "#FFF4E6" : "transparent";
          return h("td", { style: `background:${bg};color:${i !== j && pct >= 3 ? "#16181D" : "inherit"}` }, i === j ? "100" : pct.toFixed(1).replace(".", ","));
        }))));
      parts.push(h("div", { class: "scroll", style: "max-height:260px" }, t));
      parts.push(h("p", { class: "small muted" }, "In % · Zeile = Fluorochrom, Spalte = Detektor. Kompensation: Ereignisse × Inverse dieser Matrix."));
    }
    box.replaceChildren(...parts);
  } catch (e) { box.replaceChildren(h("p", { class: "warn" }, e.message)); }
}
$("#btnGrpRegex").addEventListener("click", async () => {
  const mapping = {};
  $("#grpMap").value.split(",").map((x) => x.trim()).filter(Boolean).forEach((pair) => { const [k, v] = pair.split("="); if (k && v) mapping[k.trim()] = v.trim(); });
  try { const r = await api("POST", "/api/groups/regex", { regex: $("#grpRegex").value, mapping }); toast(`${r.assigned} Proben zugeordnet`); await refreshInfo(); } catch (e) { fail(e); }
});
$("#btnPanelRegex").addEventListener("click", async () => {
  const mapping = {};
  $("#panelMap").value.split(",").map((x) => x.trim()).filter(Boolean).forEach((pair) => { const [k, v] = pair.split("="); if (k && v) mapping[k.trim()] = v.trim(); });
  try { const r = await api("POST", "/api/panels/regex", { regex: $("#panelRegex").value, mapping }); toast(`${r.assigned} Proben Panels zugeordnet`); await refreshInfo(); } catch (e) { fail(e); }
});
let grpCsv = "";
$("#grpFile").addEventListener("change", async (e) => {
  const f = e.target.files[0]; if (!f) return;
  grpCsv = await f.text();
  const first = grpCsv.split(/\r?\n/)[0];
  const sep = [";", "\t", ","].sort((a, b) => first.split(b).length - first.split(a).length)[0];
  const cols = first.split(sep).map((c) => c.replace(/^"|"$/g, ""));
  for (const id of ["#grpKey", "#grpCol"]) $(id).replaceChildren(...cols.map((c) => h("option", { value: c }, c)));
  if (cols.length > 1) $("#grpCol").value = cols[1];
});
$("#btnGrpTable").addEventListener("click", async () => {
  if (!grpCsv) { toast("Zuerst eine Tabelle wählen"); return; }
  try {
    const r = await api("POST", "/api/groups/table", { csv: grpCsv, key_col: $("#grpKey").value, group_col: $("#grpCol").value, regex: $("#grpIdRegex").value });
    toast(`${r.assigned} Proben zugeordnet`); await refreshInfo();
  } catch (e) { fail(e); }
});

/* =================== Gating =================== */
$$(".tool").forEach((b) => b.addEventListener("click", () => setTool(b.dataset.tool)));
const HINTS = {
  edit: "Auswählen: Gate anklicken, Eckpunkte oder Gate ziehen",
  polygon: "Polygon: Punkte setzen, Doppelklick schließt · Esc bricht ab",
  rect: "Rechteck: in 2D-Plot ziehen",
  quad: "Quadranten: Klick setzt den Mittelpunkt",
  range: "Bereich: im Histogramm waagerecht ziehen",
};
function setTool(t) {
  st.tool = t; $$(".tool").forEach((b) => b.classList.toggle("active", b.dataset.tool === t));
  $("#toolHint").textContent = HINTS[t]; cancelDraft();
}
$("#scope").addEventListener("change", (e) => { st.scope = e.target.value; });
$("#palette").value = st.palette;
$("#palette").addEventListener("change", (e) => { st.palette = e.target.value; store.set("palette", st.palette); legend(); st.plots.forEach((p) => p._data && paint(p)); });
$("#bins").value = String(st.bins);
$("#bins").addEventListener("change", (e) => { st.bins = +e.target.value; store.set("bins", st.bins); renderPlots(); });
function legend() { $("#legendBar").style.background = `linear-gradient(90deg, ${PALETTES[st.palette].join(",")})`; }
legend();
$("#btnAddPlot").addEventListener("click", () => addPlot(st.selNode));
$("#sampleFilter").addEventListener("input", renderSampleList);

function addPlot(node, x, y) {
  const used = new Set(st.plots.filter((p) => p.node === node).flatMap((p) => [p.x, p.y]));
  const fl = st.channels.filter((c) => c.kind === "fluor").map((c) => c.name);
  const free = fl.filter((c) => !used.has(c));
  st.plots.push({ id: plotSeq++, node: node || "root", x: x || free[0] || st.channels[0]?.name, y: y === undefined ? (free[1] || null) : y });
  savePlots(); renderPlots();
}
function renderSampleList() {
  const q = $("#sampleFilter").value.toLowerCase();
  const byGroup = {};
  st.samples.filter((x) => !x.error && x.name.toLowerCase().includes(q)).forEach((x) => (byGroup[x.group || "ohne Gruppe"] ||= []).push(x));
  const out = [];
  for (const [g, xs] of Object.entries(byGroup)) {
    out.push(h("div", { class: "grp" }, `${g} · ${xs.length}`));
    for (const x of xs) out.push(h("button", { type: "button", class: x.name === st.sample ? "sel" : "", "aria-pressed": x.name === st.sample ? "true" : "false", onclick: () => selectSample(x.name) },
      h("span", { class: "dot " + (x.qc?.status === "warn" ? "warn" : x.qc?.status === "ok" ? "ok" : ""), title: x.qc ? (x.qc.issues.join(" · ") || "QC bestanden") : "QC noch nicht ausgeführt" }),
      h("span", { class: "name" }, x.name)));
  }
  $("#sampleList").replaceChildren(...out);
}
async function selectSample(name) { st.sample = name; renderSampleList(); await loadTree(); renderPlots(); }
async function renderGating() {
  renderSampleList();
  const x = st.samples.find((y) => y.name === st.sample);
  $("#gatingTitle").textContent = st.sample || "Keine Probe";
  $("#gatingSub").textContent = x ? `${fmtInt(x.events)} Ereignisse${x.group ? " · " + x.group : ""}` : "Zuerst unter Import einen Ordner öffnen";
  $("#gatingPanel").textContent = "Panel: " + (x?.panel || "Standard");
  const w = $("#gatingWarn");
  if (x?.qc?.status === "warn") { w.textContent = "QC: " + x.qc.issues[0]; w.classList.remove("hidden"); } else w.classList.add("hidden");
  renderPlots();
}

/* ----- Hierarchie ----- */
async function loadTree() {
  try {
    st.tree = await api("GET", "/api/tree" + (st.sample ? "?sample=" + encodeURIComponent(st.sample) : ""));
    if (!st.tree.nodes.find((n) => n.id === st.selNode)) st.selNode = "root";
    renderTree();
  } catch (e) { fail(e); }
}
function nodeById(id) { return st.tree?.nodes.find((n) => n.id === id); }
function renderTree() {
  if (!st.tree) return;
  $("#tree").replaceChildren(...st.tree.nodes.map((n) => h("button", {
    type: "button", class: n.id === st.selNode ? "sel" : "", style: `padding-left:${8 + n.depth * 14}px`,
    "aria-pressed": n.id === st.selNode ? "true" : "false", onclick: () => selectNode(n.id),
    ondblclick: () => { if (n.kind && n.kind !== "boolean") addPlot(n.id); },
  }, h("span", {}, n.name, n.override ? h("span", { class: "ov", title: "für diese Probe angepasst" }, "angepasst") : null),
  h("span", { class: "pct" }, n.depth === 0 ? fmtInt(n.count) : fmtPct(n.pct_parent)))));
  renderNodeStats();
  const others = (st.tree.panels || []).filter((q) => q !== st.tree.panel);
  $("#copyFrom").replaceChildren(...others.map((q) => h("option", { value: q }, q)));
  $("#btnCopyPanel").disabled = !others.length;
  $("#log").replaceChildren(...st.tree.log.slice().reverse().slice(0, 30).map((e) => h("li", {}, h("time", {}, (e.time || "").slice(11, 16)), e.action + (e.user ? " · " + e.user : ""))));
  renderPopList();
}
function selectNode(id) { st.selNode = id; renderTree(); $$(".plot").forEach((el) => el.classList.toggle("focus", st.plots.find((p) => "p" + p.id === el.id)?.node === id)); st.plots.forEach((p) => p._data && drawGates(p)); }
function renderNodeStats() {
  const n = nodeById(st.selNode); const box = $("#nodeStats");
  if (!n) { box.replaceChildren(); return; }
  const rare = n.count !== null && n.count < 100 && n.depth > 0;
  const kv = h("dl", { class: "kv" },
    h("dt", {}, "Ereignisse"), h("dd", {}, fmtInt(n.count)),
    h("dt", {}, "% Eltern"), h("dd", {}, fmtPct(n.pct_parent, 2)),
    h("dt", {}, "% gesamt"), h("dd", {}, fmtPct(n.pct_total, 2)),
    h("dt", {}, "Zähl-CV"), h("dd", { class: rare ? "warn" : "" }, n.count_cv === null ? "–" : fmtPct(n.count_cv) + (rare ? " (selten)" : "")),
    n.per_ul ? h("dt", {}, "pro µL") : null, n.per_ul ? h("dd", {}, fmtNum(n.per_ul, 1)) : null,
    n.kind ? h("dt", {}, "Gate") : null, n.kind ? h("dd", {}, `${n.kind} · ${n.x}${n.y ? " × " + n.y : ""}`) : null);
  const btns = h("div", { class: "row" });
  if (n.id !== "root") {
    if (n.kind && n.kind !== "boolean") btns.append(h("button", { type: "button", onclick: () => addPlot(n.id) }, "Plot auf diese Population"));
    btns.append(h("button", { type: "button", onclick: async () => {
      const name = await ask("Gate umbenennen", "Name", n.name); if (!name) return;
      await api("POST", `/api/gates/${n.id}/rename`, { name }).catch(fail); await afterGateChange();
    } }, "Umbenennen"));
    if (n.override) btns.append(h("button", { type: "button", onclick: async () => {
      await api("POST", `/api/gates/${n.id}/reset`, { sample: st.sample }).catch(fail); await afterGateChange();
    } }, "Auf Vorlage zurücksetzen"));
    btns.append(h("button", { type: "button", onclick: async () => {
      if (!confirm(`Gate „${n.name}“ und alle Untergates löschen?`)) return;
      await api("DELETE", `/api/gates/${n.id}`).catch(fail); st.selNode = "root"; await afterGateChange();
    } }, "Löschen"));
    btns.append(h("button", { type: "button", onclick: () => booleanDialog(n) }, "Boolesch…"));
  } else btns.append(h("button", { type: "button", onclick: () => addPlot("root") }, "Plot auf alle Ereignisse"));
  box.replaceChildren(h("h3", {}, n.name), kv, btns);
}
async function booleanDialog(n) {
  const sibs = st.tree.nodes.filter((m) => m.parent === n.parent && m.kind && m.kind !== "boolean");
  const op = h("select", {}, h("option", { value: "and" }, "UND"), h("option", { value: "or" }, "ODER"), h("option", { value: "not" }, "NICHT (nur erstes)"));
  const list = h("div", { class: "stack" }, ...sibs.map((m) => h("label", { style: "flex-direction:row;align-items:center;gap:8px" }, h("input", { type: "checkbox", value: m.id, checked: m.id === n.id }), m.name)));
  const name = await ask("Boolesches Gate", "Name", "", h("div", { class: "stack" }, h("label", {}, "Verknüpfung", op), list));
  if (!name) return;
  const refs = $$("input:checked", list).map((i) => i.value);
  if (!refs.length) return;
  try { await api("POST", "/api/gates", { parent: n.parent, name, sample: st.sample, spec: { kind: "boolean", op: op.value, refs } }); await afterGateChange(); } catch (e) { fail(e); }
}
async function afterGateChange() { await loadTree(); renderPlots(); }

/* ----- Plots ----- */
function chOptions(sel, allowNone) {
  const o = st.channels.map((c) => h("option", { value: c.name, selected: c.name === sel }, c.display));
  if (allowNone) o.unshift(h("option", { value: "", selected: !sel }, "Histogramm"));
  return o;
}
function renderPlots() {
  const wrap = $("#plots");
  if (!st.sample || !st.tree) { wrap.replaceChildren(h("p", { class: "muted" }, "Keine Probe geladen.")); return; }
  wrap.replaceChildren(...st.plots.map((p) => plotCard(p)));
  st.plots.forEach((p) => loadPlot(p));
}
function plotCard(p) {
  const nodeSel = h("select", { "aria-label": "Elternpopulation" }, ...st.tree.nodes.filter((n) => n.kind !== "boolean" || n.id === "root").map((n) => h("option", { value: n.id, selected: n.id === p.node }, "  ".repeat(n.depth) + n.name)));
  nodeSel.addEventListener("change", () => { p.node = nodeSel.value; savePlots(); loadPlot(p); });
  const xs = h("select", { "aria-label": "x-Achse" }, ...chOptions(p.x));
  const ys = h("select", { "aria-label": "y-Achse" }, ...chOptions(p.y, true));
  xs.addEventListener("change", () => { p.x = xs.value; savePlots(); loadPlot(p); });
  ys.addEventListener("change", () => { p.y = ys.value || null; savePlots(); loadPlot(p); });
  const close = h("button", { type: "button", class: "x", "aria-label": "Plot entfernen", onclick: () => { st.plots = st.plots.filter((q) => q !== p); savePlots(); renderPlots(); } }, "×");
  const canvas = h("canvas");
  const svg = s("svg", { viewBox: "0 0 1000 1000", preserveAspectRatio: "none" });
  const wrapEl = h("div", { class: "canvas-wrap" }, canvas); wrapEl.append(svg);
  const yax = h("div", { class: "yaxis" }), xax = h("div", { class: "xaxis" });
  const foot = h("div", { class: "plot-foot" });
  const card = h("section", { class: "plot" + (p.node === st.selNode ? " focus" : ""), id: "p" + p.id, "aria-label": "Plot" },
    h("div", { class: "plot-head" }, nodeSel, xs, h("span", { class: "muted small" }, "×"), ys, close),
    h("div", { class: "plot-body" }, yax, wrapEl, h("span"), xax), foot);
  Object.assign(p, { _canvas: canvas, _svg: svg, _yax: yax, _xax: xax, _foot: foot, _wrap: wrapEl });
  bindPointer(p);
  return card;
}
async function loadPlot(p) {
  if (!st.sample) return;
  try {
    p._data = await api("POST", "/api/plot", { sample: st.sample, node: p.node, x: p.x, y: p.y, bins: st.bins });
    paint(p); drawAxes(p); drawGates(p); drawFoot(p);
  } catch (e) { p._foot && p._foot.replaceChildren(h("span", { class: "warn" }, e.message)); }
}
function paint(p) {
  const d = p._data, c = p._canvas, L = LUTS[st.palette];
  if (d.grid) {
    const n = d.bins, bytes = Uint8Array.from(atob(d.grid), (ch) => ch.charCodeAt(0));
    c.width = n; c.height = n;
    const ctx = c.getContext("2d"), img = ctx.createImageData(n, n);
    for (let i = 0; i < bytes.length; i++) { const v = bytes[i]; if (!v) continue; img.data.set(L.subarray(v * 4, v * 4 + 4), i * 4); }
    ctx.putImageData(img, 0, 0);
  } else {
    const W = 600, H = 600; c.width = W; c.height = H;
    const ctx = c.getContext("2d"); ctx.clearRect(0, 0, W, H);
    const hist = d.hist, n = hist.length;
    ctx.beginPath(); ctx.moveTo(0, H);
    hist.forEach((v, i) => ctx.lineTo(i / (n - 1) * W, H - v * (H - 20)));
    ctx.lineTo(W, H); ctx.closePath();
    ctx.fillStyle = PALETTES[st.palette][1] + "88"; ctx.fill();
    ctx.strokeStyle = PALETTES[st.palette][1]; ctx.lineWidth = 2; ctx.stroke();
  }
}
function drawAxes(p) {
  const d = p._data;
  const xl = h("span", { class: "lab", role: "button", tabindex: 0, title: "Achse einstellen", onclick: () => axisDialog(d.x.channel) }, d.x.label);
  p._xax.replaceChildren(...d.x.ticks.map((t) => h("span", { class: "t", style: `left:${t.pos * 100}%` }, t.label)), xl);
  if (d.y) {
    const yl = h("span", { class: "lab", role: "button", tabindex: 0, title: "Achse einstellen", onclick: () => axisDialog(d.y.channel) }, d.y.label);
    p._yax.replaceChildren(...d.y.ticks.map((t) => h("span", { class: "t", style: `top:${(1 - t.pos) * 100}%` }, t.label)), yl);
  } else p._yax.replaceChildren(h("span", { class: "lab" }, "Anteil"));
}
function drawFoot(p) {
  const d = p._data, f = p._foot;
  const btn = (label, fn) => h("button", { type: "button", onclick: fn }, label);
  const items = [h("span", {}, `${fmtInt(d.n)} Ereignisse${d.offscale ? ` · ${fmtInt(d.offscale)} außerhalb` : ""}`), h("span", { class: "grow" })];
  if (d.y) {
    items.push(btn("Alle Proben als Abbildung", async () => {
      try {
        toast("Erzeuge Abbildung aller Proben…", 4000);
        const r = await api("POST", "/api/export/figure", { node: p.node, x: p.x, y: p.y, palette: st.palette, title: nodeById(p.node)?.name });
        toast(`Gespeichert: ${r.pdf} (und PNG)`, 9000);
      } catch (e) { fail(e); }
    }));
    items.push(btn("Vorschlag: Population", () => propose(p, { kind: "population" })));
    if (/FSC.*A$/i.test(p.x) && /FSC.*H$/i.test(p.y || "")) items.push(btn("Vorschlag: Einzelzellen", () => propose(p, { kind: "singlet" })));
  } else {
    items.push(btn("Tal-Schwelle: darunter", () => propose(p, { kind: "valley", keep: "below" })));
    items.push(btn("darüber", () => propose(p, { kind: "valley", keep: "above" })));
  }
  f.replaceChildren(...items);
}
async function axisDialog(ch) {
  let ax; try { ax = await api("GET", `/api/axes/${encodeURIComponent(ch)}`); } catch (e) { return fail(e); }
  const kind = h("select", {}, ...[["linear", "linear"], ["log", "log10"], ["logicle", "biexponentiell (Logicle)"], ["arcsinh", "arcsinh mit Kofaktor"], ["asinh", "asinh (Gating-ML)"]].map(([v, l]) => h("option", { value: v, selected: v === ax.transform.kind }, l)));
  const cof = h("input", { type: "number", step: "any", value: ax.transform.params.cofactor ?? "", placeholder: "automatisch" });
  const W = h("input", { type: "number", step: "any", value: ax.transform.params.W ?? "", placeholder: "automatisch" });
  const lo = h("input", { type: "number", step: "any", placeholder: "Rohwert" });
  const hi = h("input", { type: "number", step: "any", placeholder: "Rohwert" });
  const extra = h("div", { class: "stack" }, h("label", {}, "Transformation", kind), h("label", {}, "Kofaktor (arcsinh)", cof), h("label", {}, "Breite W (Logicle)", W),
    h("div", { class: "row" }, h("label", {}, "Minimum", lo), h("label", {}, "Maximum", hi)),
    h("p", { class: "small muted" }, "Gilt projektweit für alle Proben. Gates behalten ihre eigene Transformation; die Zugehörigkeit ändert sich dadurch nicht."));
  const r = await ask("Achse " + ax.label, "Kanal", ch, extra);
  if (r === null) return;
  const body = {};
  if (kind.value !== ax.transform.kind || cof.value || W.value) {
    body.kind = kind.value; body.params = {};
    if (cof.value) body.params.cofactor = +cof.value;
    if (W.value) body.params.W = +W.value;
    if (kind.value === "logicle") body.params = { ...ax.transform.params, ...body.params };
  }
  if (lo.value !== "") body.raw_lo = +lo.value;
  if (hi.value !== "") body.raw_hi = +hi.value;
  try { await api("PUT", `/api/axes/${encodeURIComponent(ch)}`, body); renderPlots(); } catch (e) { fail(e); }
}

/* ----- Gates zeichnen und bearbeiten ----- */
const PX = (n) => n * 1000, PY = (n) => (1 - n) * 1000;
function drawGates(p) {
  const svg = p._svg, d = p._data; if (!d) return;
  svg.replaceChildren();
  const k = 1000 / (p._wrap.getBoundingClientRect().width || 400);
  const quadDone = new Set();
  for (const g of d.gates) {
    const sh = g.shape; if (!sh) continue;
    const sel = g.id === st.selNode;
    const cls = "gate-shape" + (sel ? " sel" : "") + (g.override ? " override" : "");
    let label = null;
    if (sh.type === "polygon") {
      const el = s("polygon", { class: cls, points: sh.points.map(([x, y]) => `${PX(x)},${PY(y)}`).join(" ") });
      el.dataset.gate = g.id; svg.append(el);
      const top = sh.points.reduce((a, b) => (b[1] > a[1] ? b : a));
      label = [PX(top[0]), PY(top[1]) - 10 * k];
      if (sel && sh.editable && st.tool === "edit") sh.points.forEach((pt, i) => { const hd = s("circle", { class: "handle", cx: PX(pt[0]), cy: PY(pt[1]), r: 5 * k, "stroke-width": 1.5 * k }); hd.dataset.gate = g.id; hd.dataset.vertex = i; svg.append(hd); });
    } else if (sh.type === "rect") {
      const x0 = PX(Math.max(sh.x0, -0.05)), x1 = PX(Math.min(sh.x1, 1.05)), y0 = PY(Math.min(sh.y1, 1.05)), y1 = PY(Math.max(sh.y0, -0.05));
      const el = s("rect", { class: cls, x: x0, y: y0, width: x1 - x0, height: y1 - y0 }); el.dataset.gate = g.id; svg.append(el);
      label = [x0 + 4 * k, y0 - 8 * k];
      if (sel && sh.editable && st.tool === "edit") [[x0, y0, "tl"], [x1, y1, "br"]].forEach(([cx, cy, c]) => { const hd = s("circle", { class: "handle", cx, cy, r: 5 * k, "stroke-width": 1.5 * k }); hd.dataset.gate = g.id; hd.dataset.corner = c; svg.append(hd); });
    } else if (sh.type === "range") {
      const x0 = PX(Math.max(sh.x0, -0.05)), x1 = PX(Math.min(sh.x1, 1.05));
      const el = s("rect", { class: cls, x: x0, y: 40, width: x1 - x0, height: 960 }); el.dataset.gate = g.id; svg.append(el);
      label = [x0 + 4 * k, 30];
      if (sel && sh.editable && st.tool === "edit") [[x0, "l"], [x1, "r"]].forEach(([cx, c]) => { const hd = s("circle", { class: "handle", cx, cy: 500, r: 6 * k, "stroke-width": 1.5 * k }); hd.dataset.gate = g.id; hd.dataset.corner = c; svg.append(hd); });
    } else if (sh.type === "quad") {
      const cx = PX(sh.x), cy = PY(sh.y);
      if (!quadDone.has(sh.group)) {
        quadDone.add(sh.group);
        const ln = s("path", { class: "gate-shape" + (d.gates.some((q) => q.shape?.group === sh.group && q.id === st.selNode) ? " sel" : ""), d: `M${cx} 0 V1000 M0 ${cy} H1000`, fill: "none" });
        ln.dataset.gate = g.id; svg.append(ln);
        if (st.tool === "edit" && sh.editable) { const hd = s("circle", { class: "handle", cx, cy, r: 6 * k, "stroke-width": 1.5 * k }); hd.dataset.gate = g.id; hd.dataset.corner = "c"; svg.append(hd); }
      }
      const pos = { "-+": [8 * k, 22 * k, "start"], "++": [1000 - 8 * k, 22 * k, "end"], "--": [8 * k, 1000 - 10 * k, "start"], "+-": [1000 - 8 * k, 1000 - 10 * k, "end"] }[sh.quadrant];
      const t = s("text", { class: "gate-label", x: pos[0], y: pos[1], "text-anchor": pos[2], "font-size": 11.5 * k, "stroke-width": 3 * k, "font-weight": sel ? 700 : 500 });
      t.textContent = `${g.name} · ${fmtPct(g.pct)}`; t.dataset.gate = g.id; svg.append(t);
      continue;
    }
    if (label) {
      const t = s("text", { class: "gate-label", x: Math.max(label[0], 4 * k), y: Math.max(label[1], 14 * k), "font-size": 11.5 * k, "stroke-width": 3 * k, "font-weight": sel ? 700 : 500 });
      t.textContent = `${g.name} · ${fmtPct(g.pct)}`; t.dataset.gate = g.id; svg.append(t);
    }
  }
  if (p._draft) drawDraft(p);
}
function norm(p, e) {
  const r = p._svg.getBoundingClientRect();
  return [Math.min(Math.max((e.clientX - r.left) / r.width, -0.05), 1.05), Math.min(Math.max(1 - (e.clientY - r.top) / r.height, -0.05), 1.05)];
}
function toDisp(ax, n) { return ax.lo + n * (ax.hi - ax.lo); }
function cancelDraft() { st.plots.forEach((p) => { p._draft = null; p._drag = null; if (p._data) drawGates(p); }); }
function drawDraft(p) {
  const dr = p._draft; const k = 1000 / (p._wrap.getBoundingClientRect().width || 400);
  let el;
  if (dr.type === "polygon") {
    const pts = dr.points.concat(dr.hover ? [dr.hover] : []);
    el = s("polyline", { class: "draft", points: pts.map(([x, y]) => `${PX(x)},${PY(y)}`).join(" ") });
    p._svg.append(el);
    dr.points.forEach(([x, y]) => p._svg.append(s("circle", { class: "handle", cx: PX(x), cy: PY(y), r: 4 * k, "stroke-width": 1.2 * k })));
  } else if (dr.type === "rect" || dr.type === "range") {
    const [a, b] = [dr.start, dr.end || dr.start];
    const x0 = PX(Math.min(a[0], b[0])), x1 = PX(Math.max(a[0], b[0]));
    const y0 = dr.type === "range" ? 40 : PY(Math.max(a[1], b[1])), y1 = dr.type === "range" ? 1000 : PY(Math.min(a[1], b[1]));
    p._svg.append(s("rect", { class: "draft", x: x0, y: y0, width: x1 - x0, height: y1 - y0 }));
  } else if (dr.type === "proposal") {
    const sh = dr.shape;
    if (sh.type === "polygon") p._svg.append(s("polygon", { class: "draft", points: sh.points.map(([x, y]) => `${PX(x)},${PY(y)}`).join(" ") }));
    else if (sh.type === "range") p._svg.append(s("rect", { class: "draft", x: PX(Math.max(sh.x0, -0.05)), y: 40, width: PX(Math.min(sh.x1, 1.05)) - PX(Math.max(sh.x0, -0.05)), height: 960 }));
  }
}
function bindPointer(p) {
  const svg = p._svg;
  svg.addEventListener("pointerdown", (e) => {
    if (!p._data) return;
    const pt = norm(p, e); const tool = st.tool; const is2d = !!p._data.y;
    const target = e.target; const gid = target.dataset?.gate;
    if (tool === "edit") {
      if (gid) { st.selNode = gid; renderTree(); }
      const g = gid && p._data.gates.find((x) => x.id === gid);
      if (g && g.shape.editable) {
        p._drag = { gate: g, start: pt, vertex: target.dataset.vertex !== undefined ? +target.dataset.vertex : null, corner: target.dataset.corner || null, shape: JSON.parse(JSON.stringify(g.shape)) };
        svg.setPointerCapture(e.pointerId);
      }
      drawGates(p); return;
    }
    if (tool === "polygon" && is2d) {
      if (!p._draft || p._draft.type !== "polygon") p._draft = { type: "polygon", points: [] };
      const pts = p._draft.points;
      // Doppelklick selbst erkennen: das SVG wird zwischen den Klicks neu gezeichnet
      const now = performance.now(), last = p._lastDown;
      p._lastDown = { t: now, pt };
      if (last && now - last.t < 400 && Math.hypot(last.pt[0] - pt[0], last.pt[1] - pt[1]) < 0.01 && pts.length >= 3) { finishPolygon(p); return; }
      if (pts.length >= 3) { const [fx, fy] = pts[0]; if (Math.hypot(fx - pt[0], fy - pt[1]) < 0.015) { finishPolygon(p); return; } }
      pts.push(pt); drawGates(p); return;
    }
    if (tool === "rect" && is2d) { p._draft = { type: "rect", start: pt }; svg.setPointerCapture(e.pointerId); return; }
    if (tool === "range" && !is2d) { p._draft = { type: "range", start: pt }; svg.setPointerCapture(e.pointerId); return; }
    if (tool === "quad" && is2d) { createGate(p, { kind: "quadrant", x: p.x, y: p.y, xdiv: toDisp(p._data.x, pt[0]), ydiv: toDisp(p._data.y, pt[1]) }, "Quadranten"); return; }
    if ((tool === "rect" || tool === "polygon" || tool === "quad") && !is2d) toast("Dieses Werkzeug braucht einen 2D-Plot");
    if (tool === "range" && is2d) toast("Bereich nur im Histogramm (y-Achse: Histogramm)");
  });
  svg.addEventListener("pointermove", (e) => {
    const pt = norm(p, e);
    if (p._drag) { moveDrag(p, pt); return; }
    if (!p._draft) return;
    if (p._draft.type === "polygon") p._draft.hover = pt; else p._draft.end = pt;
    drawGates(p);
  });
  svg.addEventListener("pointerup", async (e) => {
    if (p._drag) { const dr = p._drag; p._drag = null; if (dr.moved) await saveDrag(p, dr); return; }
    const dr = p._draft; if (!dr) return;
    if (dr.type === "rect" || dr.type === "range") {
      const a = dr.start, b = norm(p, e); p._draft = null;
      if (Math.abs(a[0] - b[0]) < 0.01 || (dr.type === "rect" && Math.abs(a[1] - b[1]) < 0.01)) { drawGates(p); return; }
      const X = p._data.x;
      const spec = { kind: "rectangle", x: p.x, y: dr.type === "rect" ? p.y : null, xmin: toDisp(X, Math.min(a[0], b[0])), xmax: toDisp(X, Math.max(a[0], b[0])) };
      if (dr.type === "rect") { spec.ymin = toDisp(p._data.y, Math.min(a[1], b[1])); spec.ymax = toDisp(p._data.y, Math.max(a[1], b[1])); }
      createGate(p, spec, dr.type === "rect" ? "Rechteck" : "Bereich");
    }
  });
  svg.addEventListener("dblclick", () => { if (p._draft?.type === "polygon" && p._draft.points.length >= 3) finishPolygon(p); });
}
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape") cancelDraft();
  const tag = (e.target.tagName || "").toLowerCase();
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === "z" && !["input", "select", "textarea"].includes(tag) && st.view === "gating") { e.preventDefault(); undo(); }
});
async function undo() {
  try { await api("POST", "/api/undo"); await afterGateChange(); toast("Rückgängig gemacht"); } catch (e) { toast(e.message); }
}
$("#btnUndo").addEventListener("click", undo);
function finishPolygon(p) {
  const pts = p._draft.points; p._draft = null;
  const uniq = pts.filter((q, i) => i === 0 || Math.hypot(q[0] - pts[i - 1][0], q[1] - pts[i - 1][1]) > 0.002);
  if (uniq.length < 3) { drawGates(p); toast("Mindestens drei Punkte"); return; }
  createGate(p, { kind: "polygon", x: p.x, y: p.y, vertices: uniq.map(([x, y]) => [toDisp(p._data.x, x), toDisp(p._data.y, y)]) }, "Population");
}
async function createGate(p, spec, defName) {
  // Quadranten erzeugen vier benannte Gates – kein Namensdialog nötig
  const name = spec.kind === "quadrant" ? "Quadranten" : await ask("Neues Gate", "Name", defName);
  drawGates(p);
  if (!name) return;
  try {
    await api("POST", "/api/gates", { parent: p.node, name, sample: st.sample, spec });
    await afterGateChange();
  } catch (e) { fail(e); }
}
function moveDrag(p, pt) {
  const dr = p._drag, sh = dr.gate.shape, o = dr.shape; dr.moved = true;
  const dx = pt[0] - dr.start[0], dy = pt[1] - dr.start[1];
  if (sh.type === "polygon") {
    if (dr.vertex !== null) sh.points[dr.vertex] = pt;
    else sh.points = o.points.map(([x, y]) => [x + dx, y + dy]);
  } else if (sh.type === "rect") {
    if (dr.corner === "tl") { sh.x0 = Math.min(pt[0], o.x1); sh.y1 = Math.max(pt[1], o.y0); }
    else if (dr.corner === "br") { sh.x1 = Math.max(pt[0], o.x0); sh.y0 = Math.min(pt[1], o.y1); }
    else { sh.x0 = o.x0 + dx; sh.x1 = o.x1 + dx; sh.y0 = o.y0 + dy; sh.y1 = o.y1 + dy; }
  } else if (sh.type === "range") {
    if (dr.corner === "l") sh.x0 = Math.min(pt[0], o.x1); else if (dr.corner === "r") sh.x1 = Math.max(pt[0], o.x0); else { sh.x0 = o.x0 + dx; sh.x1 = o.x1 + dx; }
  } else if (sh.type === "quad") {
    for (const g of p._data.gates) if (g.shape?.group === sh.group) { g.shape.x = pt[0]; g.shape.y = pt[1]; }
  }
  drawGates(p);
}
async function saveDrag(p, dr) {
  const sh = dr.gate.shape, X = p._data.x, Y = p._data.y; let spec;
  const open = (v, lo) => (v <= -0.049 || v >= 1.049) ? null : v;   // offene Grenzen bleiben offen
  if (sh.type === "polygon") spec = { kind: "polygon", x: p.x, y: p.y, vertices: sh.points.map(([x, y]) => [toDisp(X, x), toDisp(Y, y)]) };
  else if (sh.type === "rect") spec = { kind: "rectangle", x: p.x, y: p.y, xmin: open(sh.x0) === null ? null : toDisp(X, sh.x0), xmax: open(sh.x1) === null ? null : toDisp(X, sh.x1), ymin: open(sh.y0) === null ? null : toDisp(Y, sh.y0), ymax: open(sh.y1) === null ? null : toDisp(Y, sh.y1) };
  else if (sh.type === "range") spec = { kind: "rectangle", x: p.x, y: null, xmin: open(sh.x0) === null ? null : toDisp(X, sh.x0), xmax: open(sh.x1) === null ? null : toDisp(X, sh.x1) };
  else if (sh.type === "quad") spec = { kind: "quadrant", x: p.x, y: p.y, xdiv: toDisp(X, sh.x), ydiv: toDisp(Y, sh.y) };
  try {
    await api("PUT", `/api/gates/${dr.gate.id}`, { spec, scope: st.scope, sample: st.sample });
    await afterGateChange();
  } catch (e) { fail(e); }
}
async function propose(p, opts) {
  try {
    const r = await api("POST", "/api/propose", { sample: st.sample, node: p.node, x: p.x, y: p.y, ...opts });
    p._draft = { type: "proposal", shape: r.shape }; drawGates(p);
    const info = r.info, frac = info.fraction !== undefined ? fmtPct(100 * info.fraction) : "";
    const accept = h("button", { type: "button", class: "primary", onclick: async () => {
      p._draft = null; await createGate(p, r.spec, opts.kind === "singlet" ? "Einzelzellen" : opts.kind === "valley" ? (opts.keep === "below" ? "negativ" : "positiv") : "Population");
    } }, "Übernehmen");
    const reject = h("button", { type: "button", onclick: () => { p._draft = null; drawFoot(p); drawGates(p); } }, "Verwerfen");
    const desc = opts.kind === "valley" ? `Tal bei ${fmtNum(info.cut, 1)} (${info.modes} Gipfel) · ${frac}` : opts.kind === "singlet" ? `FSC-H/FSC-A ${fmtNum(info.lower, 3)}–${fmtNum(info.upper, 3)} · ${frac}` : `Dichtekontur · ${frac}`;
    p._foot.replaceChildren(h("span", {}, "Vorschlag: " + desc), h("span", { class: "grow" }), reject, accept);
  } catch (e) { fail(e); }
}
window.addEventListener("resize", () => st.plots.forEach((p) => p._data && drawGates(p)));

$("#btnCopyPanel").addEventListener("click", async () => {
  const src = $("#copyFrom").value; if (!src) return;
  if (!confirm(`Die Strategie von „${src}“ ersetzt die von „${st.tree.panel}“. Fortfahren?`)) return;
  try { await api("POST", "/api/panels/copy", { src, dst: st.tree.panel }); await afterGateChange(); toast("Strategie übernommen"); } catch (e) { fail(e); }
});
$("#btnImportWsp").addEventListener("click", async () => {
  let path = null;
  if (window.pywebview?.api?.choose_file) path = await window.pywebview.api.choose_file("FlowJo-Workspace (*.wsp)");
  else path = await ask("FlowJo-Workspace importieren", "Pfad zur .wsp-Datei", "");
  if (!path) return;
  if (!confirm("Die aktuelle Gating-Strategie wird durch die aus FlowJo ersetzt. Fortfahren?")) return;
  try {
    const r = await api("POST", "/api/import/wsp", { path });
    toast(r.skipped.length ? `Übernommen; nicht übernommen: ${r.skipped.join("; ")}` : "Workspace übernommen", 8000);
    await afterGateChange();
  } catch (e) { fail(e); }
});
$("#btnImportStrategy").addEventListener("click", () => $("#strategyFile").click());
$("#strategyFile").addEventListener("change", async (e) => {
  const f = e.target.files[0]; if (!f) return;
  try { await api("POST", "/api/import/strategy", { strategy: JSON.parse(await f.text()) }); await afterGateChange(); toast("Strategie übernommen"); } catch (err) { fail(err); }
});

/* =================== Vergleich & Export =================== */
function renderPopList() {
  if (!st.tree) return;
  const nodes = st.tree.nodes.filter((n) => n.id !== "root");
  if (!st.cmpNode || !nodes.find((n) => n.id === st.cmpNode)) st.cmpNode = nodes[0]?.id || null;
  $("#popList").replaceChildren(...nodes.map((n) => h("button", { type: "button", class: n.id === st.cmpNode ? "sel" : "", style: `padding-left:${8 + (n.depth - 1) * 12}px`, onclick: () => { st.cmpNode = n.id; renderPopList(); if (st.view === "compare") renderCompare(); } }, n.name)));
  const m = $("#metric"), cur = m.value;
  const base = [["pct_parent", "% der Elternpopulation"], ["pct_total", "% aller Ereignisse"], ["count", "Anzahl Ereignisse"], ["per_ul", "Ereignisse pro µL ($VOL)"]];
  const med = st.channels.filter((c) => c.kind === "fluor").map((c) => ["median:" + c.name, "Median " + c.display]);
  m.replaceChildren(...base.concat(med).map(([v, l]) => h("option", { value: v, selected: v === cur }, l)));
}
["#metric", "#refGroup", "#correction"].forEach((id) => $(id).addEventListener("change", () => renderCompare()));
async function renderCompare() {
  renderPopList();
  const groups = [...new Set(st.samples.map((x) => x.group || "ohne Gruppe"))];
  const rg = $("#refGroup");
  // Referenz: zuletzt gewählte, sonst eine Kontrollgruppe, sonst die erste
  const ctrl = groups.find((g) => /^(kontroll|control|ctrl|mock|ut\b|untreated|unbehandelt|vehikel|vehicle|pbs)/i.test(g));
  const curRef = groups.includes(rg.value) ? rg.value : (ctrl || groups[0]);
  rg.replaceChildren(...groups.map((g) => h("option", { value: g, selected: g === curRef }, g)));
  if (!st.cmpNode) { $("#cmpTitle").textContent = "Noch keine Gates"; return; }
  const n = nodeById(st.cmpNode);
  $("#cmpTitle").textContent = n?.name || "–";
  $("#cmpSub").textContent = (n?.path || "") + " · Panel " + st.tree.panel;
  try {
    toast("Berechne alle Proben…", 1500);
    const r = await api("POST", "/api/compare", { node: st.cmpNode, panel: st.tree.panel, metric: $("#metric").value, reference: rg.value, correction: $("#correction").value });
    drawStrip(r);
    drawTests(r);
    await drawResultTable();
  } catch (e) { fail(e); }
}
function drawStrip(r) {
  const W = 720, H = 340, L = 64, R = 16, T = 16, B = 46;
  const vals = r.samples.map((x) => x.value);
  const box = $("#strip");
  if (!vals.length) { box.replaceChildren(h("p", { class: "muted" }, "Keine Werte.")); return; }
  let lo = Math.min(...vals), hi = Math.max(...vals); const pad = (hi - lo) * 0.12 || Math.abs(hi) * 0.1 || 1; lo -= pad; hi += pad;
  const gs = r.groups.map((g) => g.group), gw = (W - L - R) / Math.max(gs.length, 1);
  const y = (v) => T + (H - T - B) * (1 - (v - lo) / (hi - lo));
  const svg = s("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Gruppenvergleich, ein Punkt je Probe" });
  svg.append(s("rect", { x: 0, y: 0, width: W, height: H, fill: "none" }));
  const nt = 5;
  for (let i = 0; i <= nt; i++) {
    const v = lo + (hi - lo) * i / nt, yy = y(v);
    svg.append(s("line", { x1: L, x2: W - R, y1: yy, y2: yy, stroke: "var(--line)", "stroke-dasharray": "3 3" }));
    const t = s("text", { x: L - 8, y: yy + 4, "text-anchor": "end", "font-size": 11, fill: "var(--muted)", "font-family": "var(--mono)" }); t.textContent = fmtNum(v, Math.abs(hi - lo) < 10 ? 2 : 0); svg.append(t);
  }
  const colors = ["#5B6270", "#D97706", "#1F5FAF", "#0F766E", "#B06AA8", "#CC3078", "#789CCC", "#8FBF7F"];
  r.groups.forEach((g, gi) => {
    const cx = L + gw * (gi + 0.5), col = colors[gi % colors.length];
    svg.append(s("line", { x1: cx - gw * 0.22, x2: cx + gw * 0.22, y1: y(g.median), y2: y(g.median), stroke: "var(--fg)", "stroke-width": 3 }));
    const pts = r.samples.filter((x) => x.group === g.group);
    pts.forEach((x, i) => {
      const jit = ((i % 5) - 2) * Math.min(gw * 0.06, 9);
      const warn = x.qc === "warn", rare = x.count_cv > 10;
      const c = s(rare ? "rect" : "circle", rare ? { x: cx + jit - 6, y: y(x.value) - 6, width: 12, height: 12 } : { cx: cx + jit, cy: y(x.value), r: 6 });
      c.setAttribute("fill", warn ? "var(--card)" : col); c.setAttribute("stroke", col); c.setAttribute("stroke-width", 2);
      const tt = s("title"); tt.textContent = `${x.sample}: ${fmtNum(x.value, 3)} (${fmtInt(x.count)} Ereignisse${warn ? ", QC-Warnung" : ""}${rare ? ", Zähl-CV > 10 %" : ""})`; c.append(tt);
      svg.append(c);
    });
    const t = s("text", { x: cx, y: H - B + 20, "text-anchor": "middle", "font-size": 12, fill: "var(--fg)", "font-weight": 500 }); t.textContent = g.group; svg.append(t);
    const t2 = s("text", { x: cx, y: H - B + 36, "text-anchor": "middle", "font-size": 11, fill: "var(--muted)" }); t2.textContent = `n = ${g.n}`; svg.append(t2);
  });
  svg.append(s("line", { x1: L, x2: L, y1: T, y2: H - B, stroke: "var(--muted)" }));
  box.replaceChildren(svg);
  $("#cmpNote").replaceChildren(h("p", {}, "Linie = Median. Offene Punkte: QC-Warnung. Quadrate: weniger als 100 Ereignisse (Zähl-CV > 10 %, Cossarizza et al. 2021)."));
}
function drawTests(r) {
  const rows = r.tests.map((t) => h("tr", {}, h("td", {}, `${t.b} vs. ${t.a}`), h("td", { class: "num" }, `${t.n_b} / ${t.n_a}`),
    h("td", { class: "num" }, fmtNum(t.diff_median, 3)), h("td", { class: "num" }, fmtNum(t.hodges_lehmann, 3)),
    h("td", { class: "num" }, fmtP(t.p)), h("td", { class: "num" }, fmtP(t.p_adj)), h("td", { class: "num" }, fmtP(t.min_p_possible))));
  const glob = r.global ? h("p", { class: "small" }, `${r.global.test}: H = ${fmtNum(r.global.statistic, 2)}, p = ${fmtP(r.global.p)}`) : null;
  $("#cmpTests").replaceChildren(h("div", { class: "tests" }, glob,
    r.tests.length ? h("div", { class: "scroll" }, h("table", { class: "data" }, h("thead", {}, h("tr", {}, ...["Vergleich", "n", "Δ Median", "Hodges-Lehmann", "p (Mann-Whitney)", "p korrigiert", "kleinstes mögliches p"].map((c) => h("th", {}, c)))), h("tbody", {}, ...rows))) : h("p", { class: "muted small" }, "Für Tests mindestens zwei Gruppen zuordnen (Import & QC)."),
    r.tests.length ? h("p", { class: "small muted" }, `Korrektur: ${r.tests[0].correction}. Exakter Test bei kleinen Gruppen; „kleinstes mögliches p“ zeigt die Grenze bei dieser Fallzahl.`) : null));
}
async function drawResultTable() {
  const rows = await api("GET", "/api/results");
  const nodes = st.tree.nodes.filter((n) => n.id !== "root");
  const by = {};
  rows.filter((r) => r.panel === st.tree.panel).forEach((r) => { (by[r.sample] ||= { group: r.group })[r.id] = r; });
  const t = $("#resultTable");
  t.replaceChildren(h("thead", {}, h("tr", {}, h("th", {}, "Probe"), h("th", {}, "Gruppe"), ...nodes.map((n) => h("th", { class: "num", title: n.path, style: n.id === st.cmpNode ? "color:var(--accent)" : "" }, n.name)))),
    h("tbody", {}, ...Object.entries(by).map(([name, r]) => h("tr", {}, h("td", { style: "font-family:var(--mono)" }, name), h("td", {}, r.group || ""),
      ...nodes.map((n) => h("td", { class: "num", style: n.id === st.cmpNode ? "background:var(--accent-soft)" : "" }, r[n.id] ? fmtNum(r[n.id].pct_parent, 2) : "–"))))));
}
$("#expSvg").addEventListener("click", () => {
  const svg = $("#strip svg"); if (!svg) return;
  const clone = svg.cloneNode(true);
  clone.setAttribute("xmlns", SVGNS);
  const css = getComputedStyle(document.documentElement);
  let txt = new XMLSerializer().serializeToString(clone);
  for (const v of ["--fg", "--muted", "--line", "--card", "--mono"]) txt = txt.split(`var(${v})`).join(css.getPropertyValue(v).trim());
  const a = h("a", { href: URL.createObjectURL(new Blob([txt], { type: "image/svg+xml" })), download: `flowcy_${(nodeById(st.cmpNode)?.name || "vergleich").replace(/\W+/g, "_")}.svg` });
  document.body.append(a); a.click(); a.remove();
});

/* ---------- Start ---------- */
refreshInfo().catch(fail);
