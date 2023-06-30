import hashlib
import json
import os
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from showyourwork.paths import PathLike, package_data, path_to_identifier
from showyourwork.plugins.staging.stages import Stage
from showyourwork.version import __version__


class ZenodoStage(Stage):
    def __init__(
        self,
        name: str,
        restore: bool,
        info_file: Optional[PathLike] = None,
        url: Optional[str] = None,
        sandbox: bool = False,
        token: Optional[str] = None,
        working_directory: Optional[PathLike] = None,
    ):
        super().__init__(name, restore, working_directory=working_directory)
        self._info_file = info_file
        if self.info_file.exists():
            with open(self.info_file, "r") as f:
                info = json.load(f)
            self.sandbox = info.get("doi", "").startswith("10.5072")
        else:
            self.sandbox = sandbox
        self._token = token
        if url is None:
            if sandbox:
                self.url = "https://sandbox.zenodo.org/api"
            else:
                self.url = "https://zenodo.org/api"
        else:
            self.url = url

    @property
    def info_file(self) -> Path:
        if self._info_file is None:
            return self.working_directory / f"{self.name}.zenodo.json"
        return Path(self._info_file)

    @property
    def draft_info_file(self) -> Path:
        return self.working_directory / f"{self.name}.draft.zenodo.json"

    def upload_info_file(self, file: PathLike) -> Path:
        ident = path_to_identifier(file)
        return self.working_directory / f"{self.name}.zenodo" / f"{ident}.upload.json"

    def snakefile(self) -> Path:
        return package_data("workflow", "rules", "zenodo.smk")

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

    @property
    def session(self) -> requests.Session:
        session = requests.Session()
        session.headers["User-Agent"] = f"showyourwork/v{__version__}"
        if self.token is not None:
            session.headers["Authorization"] = f"Bearer {self.token}"
        retry = Retry(
            backoff_factor=0.1,
            status_forcelist=[403],
            allowed_methods=["DELETE", "GET", "PUT", "POST"],
        )
        adapter = HTTPAdapter(max_retries=retry)
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        return session

    def request(
        self,
        method: str,
        path: Optional[str] = None,
        url: Optional[str] = None,
        check: bool = True,
        require_token: bool = True,
        session: Optional[requests.Session] = None,
        **kwargs: Any,
    ) -> requests.Response:
        if require_token and self.token is None:
            raise ValueError(
                "A Zenodo access token is required but one was not provided"
            )
        if url is None:
            assert path is not None
            url = f"{self.url}{path}"
        if session is None:
            with self.session as session:
                response = session.request(method, url, **kwargs)
        else:
            response = session.request(method, url, **kwargs)
        if check:
            try:
                response.raise_for_status()
            except requests.HTTPError:
                print(response.text)
                raise
        return response

    def create_draft(self, info_file: PathLike, **metadata: Any) -> None:
        metadata_proc: Dict[str, Any] = {
            "title": f"Staged showyourwork Workflow: {self.name}",
            "description": """
This is a a snapshot of the outputs of a showyourwork workflow
""",
            "creators": [{"name": f"showyourwork/v{__version__}"}],
            "upload_type": "dataset",
        }
        metadata_proc = dict(metadata_proc, **metadata)
        metadata_proc = {"metadata": metadata_proc}

        response = self.request(
            "POST",
            "/deposit/depositions",
            require_token=True,
            check=True,
            json=metadata_proc,
        )

        with open(info_file, "w") as f:
            json.dump(response.json(), f, indent=2)

    def upload_file(
        self, draft_info_file: PathLike, file: PathLike, upload_info_file: PathLike
    ) -> None:
        with open(draft_info_file, "r") as f:
            draft_info = json.load(f)

        bucket_url = draft_info["links"]["bucket"]
        ident = path_to_identifier(file)
        with open(file, "rb") as f:
            response = self.request(
                "PUT",
                url=f"{bucket_url}/{ident}",
                require_token=True,
                check=True,
                data=f,
            )

        with open(upload_info_file, "w") as f:
            json.dump(response.json(), f, indent=2)

    def publish_draft(self, draft_info_file: PathLike, info_file: PathLike) -> None:
        with open(draft_info_file, "r") as f:
            draft_info = json.load(f)
        dep_id = draft_info["id"]

        response = self.request(
            "POST",
            f"/deposit/depositions/{dep_id}/actions/publish",
            require_token=True,
            check=True,
        )

        # Save the draft data to the output file
        with open(info_file, "w") as f:
            json.dump(response.json(), f, indent=2)

    def new_record(self, info_file: PathLike, *files: PathLike, **metadata: Any) -> str:
        # Set default metadata for required fields
        metadata_proc: Dict[str, Any] = {
            "title": f"Staged showyourwork Workflow: {self.name}",
            "description": """
This is a a snapshot of the outputs of a showyourwork workflow
""",
            "creators": [{"name": f"showyourwork/v{__version__}"}],
            "upload_type": "dataset",
        }
        metadata_proc = dict(metadata_proc, **metadata)
        metadata_proc = {"metadata": metadata_proc}

        with self.session as session:
            # Create the deposition draft
            response = self.request(
                "POST",
                "/deposit/depositions",
                require_token=True,
                check=True,
                session=session,
                json=metadata_proc,
            )
            draft_data = response.json()
            dep_id = draft_data["id"]
            bucket_url = draft_data["links"]["bucket"]

            # Upload all of the target files
            for file in files:
                ident = path_to_identifier(file)
                with open(file, "rb") as f:
                    response = self.request(
                        "PUT",
                        url=f"{bucket_url}/{ident}",
                        require_token=True,
                        check=True,
                        session=session,
                        data=f,
                    )

            # Publish the record
            response = self.request(
                "POST",
                f"/deposit/depositions/{dep_id}/actions/publish",
                require_token=True,
                check=True,
                session=session,
            )

            # Save the draft data to the output file
            with open(info_file, "w") as f:
                json.dump(response.json(), f, indent=2)

            return dep_id

    def download_file(
        self, info_file: PathLike, file: PathLike, verify: bool = True
    ) -> None:
        with open(info_file, "r") as f:
            info = json.load(f)

        # Search the info file for the file we want to download
        ident = path_to_identifier(file)
        download_url = f"{info['links']['record_html']}/files/{ident}"
        files_info = list(
            filter(lambda f: f["filename"] == ident, info.get("files", []))
        )
        if len(files_info) != 1:
            raise RuntimeError(
                f"File {file} not found in record metadata file {info_file}"
            )
        file_info = files_info[0]

        # Stream from the download URL directly into the target file
        if verify:
            checksum = hashlib.md5()
        with self.request(
            "GET",
            url=download_url,
            require_token=False,
            check=True,
            params={"download": 1},
            stream=True,
        ) as response:
            with open(file, "wb") as output_file:
                for chunk in response.iter_content(chunk_size=4096):
                    if verify:
                        checksum.update(chunk)
                    output_file.write(chunk)

        if verify:
            if checksum.hexdigest() != file_info["checksum"]:
                raise RuntimeError(f"Checksum mismatch for downloaded file {file}")
