# FlowCy Analysis – Desktop-App-Design

Entwurf einer Desktop-Anwendung zur Auswertung durchflusszytometrischer Daten (FCS 3.0/3.1).
Interaktiver Canvas mit allen Screens: https://claude.ai/artifact/EjnxHX2Mx1KWmvBVhXq2YP

## Zielgruppe

Klinisch-wissenschaftliche Labore: Ärztinnen und Ärzte, Biologinnen und Biologen, Studienteams.
Anforderungen: reproduzierbares Gating, nachvollziehbare QC, publikations- und auditfähige Ergebnisse.

## Arbeitsablauf (3 Schritte, oben als Navigation)

| Schritt | Screen | Kernfunktionen |
|---|---|---|
| 1 | Import & QC | FCS-Import, Gruppenzuordnung, Flussrate-/Signal-QC je Datei, Metadaten (`$CYT`, `$TOT`, …), SHA-256-Prüfsumme, Spillover-Matrix |
| 2 | Gating | Probenliste mit QC-Status, 2D-Plots (Pseudocolor-Dichte), Gate-Werkzeuge, Gating-Hierarchie, Statistik je Gate, Änderungsprotokoll |
| 3 | Vergleich & Export | Strip-Plot je Population, Median-Differenz, exakter Mann-Whitney-U, Korrektur für multiples Testen, Ergebnistabelle, Export, MIFlowCyt-Checkliste |

## Layout

- Drei Spalten: Navigation/Kontext links, Arbeitsfläche Mitte, Details/Aktionen rechts.
- Kopfzeile: Projekt, Schritt-Navigation, primäre Aktion des Schritts.
- Werkzeugleiste nur im Gating: Auswahl, Rechteck, Polygon, Ellipse, Quadrant, Histogramm; Achsentransformation (Logicle, linear, log, ArcSinh).

## Gestaltung

- Schrift: IBM Plex Sans (UI), IBM Plex Mono (Zahlen, Dateinamen, Kanäle).
- Farben: neutrales Grau, Akzent `#1F5FAF`; Warnungen Bernstein `#92400E`/`#FEF3C7`, OK Petrol `#0F766E`.
- Dichteplots in Viridis (wahrnehmungslinear, farbfehlsichtig-sicher).
- Gruppenvergleich Grau vs. Orange (unterscheidbar auch in Helligkeit); QC-auffällige Proben als offene Punkte.
- Status nie nur über Farbe: Icon + Text.
- Zahlenformat deutsch (Tausenderpunkt, Dezimalkomma).

## Fachliche Leitplanken

- QC-Warnungen bleiben bis in Gating und Vergleich sichtbar (Beispiel B03, Flussrate-Einbruch).
- Gating-Vorlagen werden versioniert und auf Gruppen angewendet; jede Änderung landet im Protokoll.
- Statistik: Hinweis auf kleine Fallzahlen (bei n = 4/4 ist p = 0,029 das Minimum des exakten Tests).
- Export: CSV/XLSX, PDF-Bericht, SVG/PNG, Gating-ML 2.0, R-/Python-Skript (flowCore, FlowKit).

## Offen

- Gerätekonfiguration, Antikörper-Klone und Lots (Platzhalter `[GERÄT]`).
- Plattform (Electron, Tauri oder nativ), Lizenz- und Validierungsanforderungen (z. B. IVDR, GxP).
- Hochdimensionale Analysen (UMAP/t-SNE, Clustering) als möglicher vierter Schritt.

Alle Zahlen im Entwurf sind Beispielwerte.
