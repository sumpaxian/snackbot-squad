"""Kommandozeile: flowcy validate | qc | demo | app"""
from __future__ import annotations

import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="flowcy", description="FlowCy Analysis")
    sub = ap.add_subparsers(dest="cmd", required=True)

    v = sub.add_parser("validate", help="Test-Analyzer: vorhandene FCS-Dateien und FlowJo-Workspaces prüfen (nur lesend)")
    v.add_argument("folders", nargs="+", help="Ordner oder einzelne FCS-Dateien")
    v.add_argument("--out", required=True, help="eigener Ausgabeordner (nicht im Datenordner)")
    v.add_argument("--gates", help="Gating-Strategie (.flowcy.json) auf alle Dateien anwenden")
    v.add_argument("--wsp", nargs="*", default=[], help="FlowJo-Workspaces; ohne Angabe werden *.wsp in den Ordnern gesucht")
    v.add_argument("--pattern", default=".", help="regulärer Ausdruck für Dateinamen")
    v.add_argument("--max-files", type=int)
    v.add_argument("--no-checksum", action="store_true")
    v.add_argument("--channels", nargs="*", help="Kanäle für Median/rCV in der Statistik")

    d = sub.add_parser("demo", help="synthetische Beispielstudie erzeugen")
    d.add_argument("folder")

    a = sub.add_parser("app", help="Desktop-App starten")
    a.add_argument("folder", nargs="?", help="Ordner mit FCS-Dateien")
    a.add_argument("--port", type=int, default=0)
    a.add_argument("--browser", action="store_true", help="im Browser statt im eigenen Fenster öffnen")

    args = ap.parse_args(argv)
    if args.cmd == "validate":
        from .validate import run
        s = run(args.folders, args.out, gates=args.gates, wsp=args.wsp, pattern=args.pattern,
                max_files=args.max_files, checksum=not args.no_checksum, stats_channels=args.channels)
        print(s)
        return 0 if s["errors"] == 0 else 1
    if args.cmd == "demo":
        from .synth import make_study
        make_study(args.folder)
        print(f"Beispielstudie in {args.folder}")
        return 0
    if args.cmd == "app":
        from .desktop import launch
        return launch(args.folder, port=args.port, browser=args.browser)
    return 2


if __name__ == "__main__":
    sys.exit(main())
