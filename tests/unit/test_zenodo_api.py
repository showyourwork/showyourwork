import json
from pathlib import Path

import pytest
import snakemake

from showyourwork.zenodo import Zenodo, get_dataset_dois, get_dataset_urls


# Mock snakemake workflow for tests
class DummyWorkflow:
    config = {"github_actions": False}


snakemake.workflow = DummyWorkflow()

# Test data file path
TEST_DATA_FILE = Path(__file__).parent / "data" / "TOI640b.json"


@pytest.fixture
def config():
    config = {}
    config["datasets"] = {
        "10.5281/zenodo.6515864": {
            "destination": "src/data/TOI640",
            "contents": {
                "README.md": None,
                "TOI640b.json": "src/data/TOI640/planet.json",
                "images.tar.gz": {
                    "README.md": None,
                    "S06": {"image.json": "src/data/TOI640/S06.json"},
                    "S07": {"image.json": "src/data/TOI640/S07.json"},
                },
                "lightcurves.zip": {
                    "lightcurves": {
                        "README.md": None,
                        "S06": {"lc.txt": "src/data/TOI640/S06.txt"},
                        "S07": {"lc.txt": "src/data/TOI640/S07.txt"},
                    }
                },
            },
        }
    }
    return config


def test_get_dataset_urls(config):
    # Test with a single good file
    urls = get_dataset_urls(["src/data/TOI640/planet.json"], config["datasets"])
    assert len(urls) == 1
    assert urls[0] == "https://zenodo.org/records/6515864"

    # Test with part of a zip file
    urls = get_dataset_urls(["src/data/TOI640/S06.json"], config["datasets"])
    assert urls[0] == "https://zenodo.org/records/6515864"

    # Test with a non-existing file
    urls = get_dataset_urls(["src/data/TOI640/S067.json"], config["datasets"])
    assert len(urls) == 0


def test_get_dataset_doi(config):
    # Test with a single good file
    dois = get_dataset_dois(["src/data/TOI640/planet.json"], config["datasets"])
    assert len(dois) == 1
    assert dois[0] == "10.5281/zenodo.6515864"

    # Test with part of a zip file
    dois = get_dataset_dois(["src/data/TOI640/S06.json"], config["datasets"])
    assert dois[0] == "10.5281/zenodo.6515864"

    # Test with a non-existing file
    dois = get_dataset_dois(["src/data/TOI640/S067.json"], config["datasets"])
    assert len(dois) == 0


def test_get_id_type(config):
    doi = next(iter(config["datasets"]))
    deposit = Zenodo(doi)
    assert deposit.get_id_type() == "version"


# TODO: Try and reduce uploads to zenodo in these tests, maybe using existing deposit
# TODO: Properly mark to test only when remote access?
# TODO: xfail when no access token?
def test_download_latest_draft():
    sandbox = Zenodo("sandbox")

    cache_dir = sandbox._download_latest_draft()
    assert cache_dir.is_dir()

    with open(cache_dir / ".metadata.json") as f:
        metadata = json.load(f)

    # Upload a file and publish to test automated creation
    draft = sandbox._get_draft()
    sandbox.upload_file_to_draft(draft, TEST_DATA_FILE, "testing")
    sandbox.publish()

    cache_dir_post = sandbox._download_latest_draft()
    assert cache_dir_post.is_dir()

    assert cache_dir == cache_dir_post

    with open(cache_dir_post / ".metadata.json") as f:
        metadata_post = json.load(f)

    # Ensure DOI is different, i.e. new draft was created
    assert (
        metadata_post["prereserve_doi"]["recid"] != metadata["prereserve_doi"]["recid"]
    )


def test_get_draft():
    sandbox = Zenodo("sandbox")

    # Test that we can recover the intiial draft
    draft = sandbox._get_draft()
    assert isinstance(draft, dict)
    assert "created" in draft
    assert not draft["submitted"]

    # Upload a file and publish to test automated creation
    sandbox.upload_file_to_draft(draft, TEST_DATA_FILE, "testing")
    sandbox.publish()

    # Make sure a new draft is created on request
    draft_post = sandbox._get_draft()
    assert not draft_post["submitted"]
    assert draft != draft_post
    assert draft["created"] < draft_post["created"]
