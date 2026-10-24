# FlowCy Analysis

Desktop-Anwendung zur Auswertung durchflusszytometrischer Daten (FCS 2.0–3.1).
Import und Qualitätskontrolle → Gating → Gruppenvergleich und Export.
Die Daten bleiben auf dem Rechner; Datenordner werden nie beschrieben.

## Starten

**macOS:** `FlowCy starten.command` doppelklicken. Beim ersten Start wird unter
`~/FlowCy/.venv` eine eigene Python-Umgebung angelegt (Python ≥ 3.10 nötig).

**Kommandozeile:**

    pip install -e ".[desktop]"
    flowcy app [Ordner]          # eigenes Fenster (pywebview), sonst Browser
    flowcy app --browser         # immer im Browser
    flowcy demo ~/FlowCy/Beispielstudie   # synthetische Studie mit bekannter Wahrheit

Projekte werden unter `~/FlowCy/Projekte/`, Abbildungen unter `~/FlowCy/Export/` gespeichert.

## Funktionen

| Schritt | Inhalt |
|---|---|
| Import & QC | Ordner einlesen (nur Dateiköpfe, schnell), Gruppen und Panels aus Dateinamen oder Tabelle, Fluss-QC mit Ausschluss auffälliger Zeitabschnitte, Sättigung, Metadaten, Spillover-Matrix |
| Gating | Pseudocolor-Dichte (FlowJo-artig), Viridis, Graustufen; Polygon, Rechteck, Quadranten, Bereich; Eckpunkte und Gates ziehen; boolesche Gates; Vorschläge (Dichtekontur, Einzelzell-Band, Tal-Schwelle) zum Annehmen; Anpassung nur für eine Probe; Achsen: linear, log10, Logicle, arcsinh mit Kofaktor, asinh; projektweit gleiche Achsen; eigene Strategie je Panel; FlowJo-Workspace und FlowCy-Strategie importieren; Änderungsprotokoll |
| Vergleich & Export | Punktdiagramm je Gruppe, Mann-Whitney-U (exakt bei kleinen Gruppen), Kruskal-Wallis, Holm/BH, Hodges-Lehmann, kleinstes mögliches p, Zähl-CV seltener Populationen; Ergebnistabelle; CSV, Strategie-JSON, SVG; Abbildung aller Proben mit identischen Achsen (PDF/PNG, 7 pt) |

## Test-Analyzer

Prüft FlowCy an vorhandenen Daten, ohne sie zu verändern:

    flowcy validate <Ordner> --out <eigener Ausgabeordner> [--gates strategie.flowcy.json] [--wsp workspace.wsp]

Ergebnis: `bericht.html`, Inventar, QC, Kanalbeschriftung, automatische Plausibilitäts-Strategie,
Statistik und – für jeden gefundenen FlowJo-Workspace – Konkordanz FlowCy vs. FlowJo vs. FlowKit.

## Rechenweg

- Einlesen über FlowKit/flowio, Kompensation `Ereignisse × Spillover⁻¹` (Matrix aus der Datei oder dem Workspace).
- Gates werden in Rohwerten gespeichert und in der Anzeigeskala ihrer Transformation geprüft
  (Gating-ML-Semantik; identisch mit FlowKit in allen geprüften Workspaces).
- Dichte in Klassen der transformierten Achse, Achsengrenzen datengetrieben (0,1.–99,9. Perzentil).
- Fluss-QC: Zeitfenster, robuste z-Werte (Median/MAD), zusammenhängende Abschnitte.

## Tests

    pip install -e ".[test]"
    pytest                                     # synthetische Daten mit bekannter Wahrheit
    FLOWKIT_DATA=/pfad/FlowKit/data pytest     # zusätzlich FlowJo-Konkordanz

Siehe `tests/README.md`.
