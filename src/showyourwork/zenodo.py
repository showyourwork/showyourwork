import os
from functools import cached_property
from typing import Any, Dict, Optional

import requests

from showyourwork.version import __version__


class Zenodo:
    def __init__(
        self,
        deposition_id: Optional[str] = None,
        token: Optional[str] = None,
        sandbox: bool = False,
    ):
        self.deposition_id = deposition_id
        self._token = token
        self.sandbox = sandbox
        if sandbox:
            self.base_url = "https://sandbox.zenodo.org/api"
        else:
            self.base_url = "https://zenodo.org/api"
        self.session: Optional[requests.Session] = None

    @cached_property
    def token(self) -> Optional[str]:
        if self._token is None:
            if self.sandbox:
                token = os.environ.get("SANDBOX_TOKEN", None)
            else:
                token = os.environ.get("ZENODO_TOKEN", None)
            return token
        else:
            return self._token

    def __enter__(self) -> "Zenodo":
        self.session = requests.Session()
        self.session.headers["User-Agent"] = f"showyourwork/{__version__}"
        if self.token is not None:
            self.session.headers["Authorization"] = f"Bearer {self.token}"
        return self

    def __exit__(self, *_: Any) -> None:
        if self.session is not None:
            self.session.close()

    def request(
        self,
        method: str,
        path: Optional[str] = None,
        url: Optional[str] = None,
        check: bool = True,
        require_token: bool = True,
        **kwargs: Any,
    ) -> requests.Response:
        if require_token and self.token is None:
            raise ValueError(
                "A Zenodo access token is required but one was not provided"
            )
        if url is None:
            assert path is not None
            url = f"{self.base_url}{path}"
        if self.session is None:
            with self:
                response = self.session.request(method, url, **kwargs)  # type: ignore
        else:
            response = self.session.request(method, url, **kwargs)
        if check:
            response.raise_for_status()
        return response

    def create_new_draft(self, **metadata: Any) -> Dict[str, Any]:
        metadata["title"] = metadata.get("title", "Staged data from showyourwork")
        metadata["description"] = metadata.get(
            "description",
            f"Data automatically uploaded by version {__version__} of the "
            "<code>showyourwork</code> workflow. Please visit "
            "<a href='https://github.com/showyourwork'>github.com/showyourwork</a> "
            "for more information.",
        )
        metadata["upload_type"] = metadata.get("upload_type", "dataset")
        metadata["creators"] = metadata.get("creators", [{"name": "showyourwork"}])

        metadata = {"metadata": metadata}
        data = self.request("POST", "/deposit/depositions", json=metadata).json()
        self.deposition_id = data["id"]
        return data
