"""
Downloads and installs ``tectonic`` from source.
This script is called from the ``tectonic`` rule.

"""
import json
import urllib
import urllib.request
import tarfile


# Params defined in `../rules/tectonic.smk`
TEMP = snakemake.params["TEMP"]
OS = snakemake.params["OS"]


def get_tectonic_link():
    """
    Get the download link for the latest Linux release of tectonic on GitHub.

    """
    link = None
    with urllib.request.urlopen(
        "https://api.github.com/repos/tectonic-typesetting/tectonic/releases"
    ) as url:
        data = json.loads(url.read().decode())
        for entry in data:
            if entry.get("tag_name", "") == "continuous":
                assets = entry.get("assets")
                for asset in assets:
                    if OS in asset.get("name", ""):
                        link = asset.get("browser_download_url")
                        return link


# Download the tarball
link = get_tectonic_link()
urllib.request.urlretrieve(link, "tectonic.tar.gz")


# Extract it
with tarfile.open("tectonic.tar.gz") as file:
    file.extractall(TEMP)