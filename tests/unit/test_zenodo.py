import pytest

from showyourwork import zenodo


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
    urls = zenodo.get_dataset_urls(["src/data/TOI640/planet.json"], config["datasets"])
    assert len(urls) == 1
    assert urls[0] == "https://zenodo.org/records/6515864"

    # Test with part of a zip file
    urls = zenodo.get_dataset_urls(["src/data/TOI640/S06.json"], config["datasets"])
    assert urls[0] == "https://zenodo.org/records/6515864"

    # Test with a non-existing file
    urls = zenodo.get_dataset_urls(["src/data/TOI640/S067.json"], config["datasets"])
    assert len(urls) == 0


def test_get_dataset_doi(config):
    # Test with a single good file
    dois = zenodo.get_dataset_dois(["src/data/TOI640/planet.json"], config["datasets"])
    assert len(dois) == 1
    assert dois[0] == "10.5281/zenodo.6515864"

    # Test with part of a zip file
    dois = zenodo.get_dataset_dois(["src/data/TOI640/S06.json"], config["datasets"])
    assert dois[0] == "10.5281/zenodo.6515864"

    # Test with a non-existing file
    dois = zenodo.get_dataset_dois(["src/data/TOI640/S067.json"], config["datasets"])
    assert len(dois) == 0
