# Tests

    pip install -e ".[test]"
    pytest

Synthetische FCS-Dateien mit bekannter Wahrheit werden zur Laufzeit erzeugt (`flowcy/synth.py`).

FlowJo-Konkordanz: öffentliche Workspaces aus dem FlowKit-Repository (BSD-3, © Scott White).

    git clone --depth 1 https://github.com/whitews/FlowKit.git /tmp/FlowKit
    FLOWKIT_DATA=/tmp/FlowKit/data pytest

Ohne `FLOWKIT_DATA` werden diese Tests übersprungen.
