import time
import uuid
from threading import Thread
from typing import Any, Dict, Optional

import requests
from flask import Blueprint, Flask, request, session, url_for
from werkzeug.serving import make_server

api = Blueprint("api", __name__)


@api.route("/deposit/depositions", methods=["POST"])
def create() -> Dict[str, Any]:
    assert request.headers["Authorization"] == "Bearer test"
    return {
        "id": 1234,
        "links": {
            "bucket": url_for("api.bucket", _external=True),
        },
    }


@api.route("/bucket", methods=["PUT"], defaults={"filename": ""})
@api.route("/bucket/<filename>", methods=["PUT"])
def bucket(filename: str) -> Dict[str, Any]:
    assert request.headers["Authorization"] == "Bearer test"
    session["files"] = session.get("files", []) + [filename]
    return {}


@api.route("/deposit/depositions/<dep_id>/actions/publish", methods=["POST"])
def publish(dep_id: str) -> Dict[str, Any]:
    assert request.headers["Authorization"] == "Bearer test"
    return {
        "doi": f"10.5281/zenodo.{dep_id}",
        "files": [{"filename": f} for f in session.get("files", [])],
    }


class ZenodoMock:
    def __init__(self, port: int = 5050):
        self.port = port
        self.app = Flask(__name__)
        self.server = make_server("localhost", self.port, self.app)
        self.url = f"http://localhost:{self.port}"
        self.thread: Optional[Thread] = None

        self.app.secret_key = uuid.uuid4().hex

        @self.app.route("/alive", methods=["GET"])
        def _() -> str:
            return "True"

        self.app.register_blueprint(api, url_prefix="/api")

    def start(self) -> None:
        self.thread = Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

        alive = False
        tries = 0
        while not alive:
            if tries >= 50:  # noqa
                raise Exception(
                    "Failed to start and connect to mock server. "
                    f"Is port {self.port} in use by another application?"
                )
            tries += 1
            try:
                requests.get(f"{self.url}/alive", timeout=0.2)
                alive = True
            except (
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
            ):
                time.sleep(0.1)

    def stop(self) -> None:
        self.server.shutdown()
        if self.thread is not None:
            self.thread.join()
