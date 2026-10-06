#!/usr/bin/env python3
"""Serve the offline advisor at http://127.0.0.1:8000 (no third-party packages).

Endpoints: GET /api/config, GET /api/evidence, POST /api/recommend (JSON object).
Only loopback binding is accepted; the UI makes no external requests.
"""

from __future__ import annotations

import argparse
import json
import mimetypes
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit

import advisor

ROOT = Path(__file__).resolve().parent
WEB = ROOT / "web"
MAX_BODY = 32 * 1024


class Handler(BaseHTTPRequestHandler):
    server_version = "SequencingAdvisor/1.0"

    def _respond(self, status: int, data: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'")
        self.end_headers()
        self.wfile.write(data)

    def _json(self, status: int, value: dict) -> None:
        self._respond(status, json.dumps(value, ensure_ascii=False, allow_nan=False).encode("utf-8"), "application/json; charset=utf-8")

    def do_GET(self) -> None:
        path = unquote(urlsplit(self.path).path)
        try:
            if path == "/api/config":
                self._json(200, {"config": advisor.load_config(), "defaults": advisor.normalize_input({}),
                                 "scenarios": advisor.read_json(ROOT / "config" / "scenarios.json")})
                return
            if path == "/api/evidence":
                self._json(200, advisor.load_evidence())
                return
            if path == "/api/health":
                self._json(200, {"status": "ok", "model_version": advisor.MODEL_VERSION})
                return
            requested = (WEB / ("index.html" if path == "/" else path.lstrip("/"))).resolve()
            if not requested.is_relative_to(WEB.resolve()) or not requested.is_file():
                self._json(404, {"error": "Page not found."})
                return
            mime = mimetypes.guess_type(requested.name)[0] or "application/octet-stream"
            if mime.startswith("text/") or mime in ("application/javascript", "application/json"):
                mime += "; charset=utf-8"
            self._respond(200, requested.read_bytes(), mime)
        except (advisor.InputError, OSError, ValueError) as error:
            self._json(500, {"error": str(error)})

    def do_POST(self) -> None:
        if urlsplit(self.path).path != "/api/recommend":
            self._json(404, {"error": "Endpoint not found."})
            return
        # The local page uses JSON; require it and reject cross-origin requests.
        origin = self.headers.get("Origin")
        host = self.headers.get("Host", "")
        if origin and origin != f"http://{host}":
            self._json(403, {"error": "Cross-origin request rejected."})
            return
        if self.headers.get("Content-Type", "").split(";")[0].strip() != "application/json":
            self._json(415, {"error": "Use Content-Type: application/json."})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if not 0 < length <= MAX_BODY:
                self._json(413, {"error": "JSON request must contain 1–32768 bytes."})
                return
            scenario = json.loads(self.rfile.read(length).decode("utf-8"))
            self._json(200, advisor.recommend(scenario))
        except (advisor.InputError, ValueError, UnicodeError) as error:
            self._json(400, {"error": str(error)})
        except OSError as error:
            self._json(500, {"error": str(error)})

    def log_message(self, format: str, *args: object) -> None:
        print(f"{self.log_date_time_string()} {format % args}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--host", choices=("127.0.0.1", "localhost", "::1"), default="127.0.0.1")
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error("port must be from 1 to 65535")
    if args.host == "::1":
        import socket
        ThreadingHTTPServer.address_family = socket.AF_INET6
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"Sequencing Experiment Design Advisor: http://{args.host if args.host != '::1' else '[::1]'}:{args.port}")
    print("Press Ctrl+C to stop. No data leave this computer.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
