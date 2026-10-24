"""Desktop-Start: lokaler Server plus eigenes Fenster (pywebview) oder Browser."""
from __future__ import annotations

import socket
import threading
import time
import webbrowser


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class Api:
    """Aus der Oberfläche aufrufbar: window.pywebview.api.*"""

    def choose_folder(self):
        import webview
        res = webview.windows[0].create_file_dialog(webview.FOLDER_DIALOG)
        return res[0] if res else None

    def choose_file(self, kinds="FlowJo-Workspace (*.wsp)"):
        import webview
        res = webview.windows[0].create_file_dialog(webview.OPEN_DIALOG, file_types=(kinds, "Alle Dateien (*.*)"))
        return res[0] if res else None

    def save_file(self, name):
        import webview
        res = webview.windows[0].create_file_dialog(webview.SAVE_DIALOG, save_filename=name)
        return res if isinstance(res, str) else (res[0] if res else None)


def serve(port):
    import uvicorn
    from .server import app
    uvicorn.run(app, host="127.0.0.1", port=port, log_level="warning")


def launch(folder=None, port=0, browser=False):
    port = port or _free_port()
    threading.Thread(target=serve, args=(port,), daemon=True).start()
    url = f"http://127.0.0.1:{port}/"
    for _ in range(100):
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.1):
                break
        except OSError:
            time.sleep(0.05)
    if folder:
        from .server import open_folder
        open_folder({"folder": folder})
    if not browser:
        try:
            import webview
            webview.create_window("FlowCy Analysis", url, width=1440, height=920, min_size=(1024, 680), js_api=Api())
            webview.start()
            return 0
        except ImportError:
            print("pywebview nicht installiert – öffne im Browser (pip install pywebview für ein eigenes Fenster)")
    webbrowser.open(url)
    print(f"FlowCy läuft auf {url} – beenden mit Strg+C")
    try:
        while True:
            time.sleep(3600)
    except KeyboardInterrupt:
        return 0
