#!/bin/bash
# FlowCy Analysis starten (macOS: Doppelklick im Finder).
# Legt beim ersten Start eine eigene Python-Umgebung unter ~/FlowCy/.venv an –
# bestehende Umgebungen (z. B. conda "bioinfo") werden nicht verändert.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
VENV="$HOME/FlowCy/.venv"
PY="$(command -v python3.12 || command -v python3.11 || command -v python3)"
if [ ! -x "$VENV/bin/python" ]; then
  echo "Erster Start: richte Python-Umgebung in $VENV ein (einige Minuten)…"
  mkdir -p "$HOME/FlowCy"
  "$PY" -m venv "$VENV"
  "$VENV/bin/pip" install --upgrade pip >/dev/null
fi
"$VENV/bin/pip" install -q -e "$HERE[desktop]"
exec "$VENV/bin/flowcy" app "$@"
