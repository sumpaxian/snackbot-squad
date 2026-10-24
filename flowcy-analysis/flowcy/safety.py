"""Schreibschutz: FlowCy schreibt nie in Cloud-Sync-Ordner (Dropbox u. ä.) oder Datenordner.

Rohdaten liegen oft in Dropbox, OneDrive, iCloud oder Google Drive; ein versehentlich dort
angelegtes Projekt, ein Export oder Cache würde sofort an alle Beteiligten synchronisiert.
"""
from __future__ import annotations

import os
import pathlib

PROTECTED_PARTS = ("dropbox", "onedrive", "icloud", "mobile documents", "google drive", "googledrive",
                   "cloudstorage", "box sync", "nextcloud", "owncloud")


class ProtectedPathError(PermissionError):
    pass


def _extra():
    return [pathlib.Path(p).expanduser().resolve() for p in os.environ.get("FLOWCY_READONLY", "").split(os.pathsep) if p]


def is_protected(path) -> bool:
    p = pathlib.Path(path).expanduser()
    try:
        p = p.resolve()
    except OSError:
        pass
    # Ordnernamen vergleichen, nicht Teilstrings des Pfads („Dropbox“, „Dropbox (Firma)“, „OneDrive - X“)
    if any(c.lower().startswith(PROTECTED_PARTS) for c in p.parts):
        return True
    return any(p == r or r in p.parents for r in _extra())


def check_writable(path):
    """Vor jedem Schreiben aufrufen. Wirft ProtectedPathError für geschützte Orte."""
    if is_protected(path):
        raise ProtectedPathError(
            f"FlowCy schreibt nicht in Cloud-Sync- oder Datenordner: {path}. "
            "Bitte einen lokalen Ordner wählen, z. B. ~/FlowCy/ oder ~/Downloads/.")
    return pathlib.Path(path).expanduser()
