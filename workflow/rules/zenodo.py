"""
Collect information about all the datasets we'll need to upload to/download from Zenodo.

"""
from pathlib import Path
import copy
from sphinx_mock import *


class zenodo:
    """
    Class containing metadata for all Zenodo deposits.

    """

    #: Script for generating each dataset
    script = {}

    #: File name for each dataset
    file_name = {}

    #: Path to each dataset
    file_path = {}

    #: Use Zenodo sandbox for a particular dataset?
    sandbox = {}

    #: Zenodo token name (name of an env. var.) for each dataset
    token_name = {}

    #: Zenodo deposit title for each dataset
    deposit_title = {}

    #: Zenodo deposit description for each dataset
    deposit_description = {}

    #: Zenodo authors for each dataset
    deposit_creators = {}

    #: Id for each Zenodo download
    zenodo_id = {}

    #: Url for each Zenodo download
    zenodo_url = {}


# Get repo name for Zenodo metadata
repo = "/".join(get_repo_url().split("/")[-2:])


# Loop over datasets
for dataset in config["zenodo"]:

    # Get the dataset name & properties
    if type(dataset) is str:
        dependency = dataset
        dep_props = config["zenodo"][dependency]
    elif type(dataset) is OrderedDict:
        dependency = list(dataset)[0]
        dep_props = dict(dataset[dependency])
    else:
        raise ValueError("Cannot parse `zenodo` entry in `showyourwork.yml`.")

    # Sandbox?
    zenodo.sandbox[dependency] = dep_props.get("sandbox", False)
    if zenodo.sandbox[dependency]:
        zenodo.zenodo_url[dependency] = "sandbox.zenodo.org"
    else:
        zenodo.zenodo_url[dependency] = "zenodo.org"

    # Generate & upload settings
    zenodo.script[dependency] = dep_props.get("script", None)
    zenodo.file_name[dependency] = str(Path(dependency).name)
    zenodo.file_path[dependency] = str(Path(dependency).parent)
    zenodo.token_name[dependency] = dep_props.get("token_name", "ZENODO_TOKEN")
    zenodo.deposit_title[dependency] = dep_props.get(
        "title", f"{repo}:{zenodo.file_name[dependency]}"
    )
    zenodo.deposit_description[dependency] = dep_props.get(
        "description", f"File uploaded from {repo}."
    )
    zenodo.deposit_creators[dependency] = dep_props.get(
        "creators", get_repo_url().split("/")[-2]
    )

    # Download settings
    zenodo.zenodo_id[dependency] = dep_props.get("id", None)

    # Checks
    if zenodo.script[dependency] is None and zenodo.zenodo_id[dependency] is None:
        raise ValueError(
            f"Please provide a `script` or a Zenodo `id` for dependency "
            f"{dependency}."
        )
    elif zenodo.script[dependency] and zenodo.zenodo_id[dependency]:
        raise ValueError(
            f"Please provide either a `script` *or* a Zenodo `id` for "
            f"dependency {dependency} (but not both)."
        )

    # Name of the dummy file we'll use to track the upload/download
    dot_zenodo_file = posix(f"{dependency}.zenodo")

    # Make the deposit a dependency of the PDF
    # so we can add Zenodo links to the figure caption
    if not dot_zenodo_file in files.dot_zenodo:
        files.dot_zenodo.append(dot_zenodo_file)

    # Track download-only files
    if zenodo.zenodo_id[dependency] and not dependency in files.zenodo_files_manual:
        files.zenodo_files_manual.append(dependency)

    # Track upload-download files
    if zenodo.script[dependency] and not dependency in files.zenodo_files_auto:
        files.zenodo_files_auto.append(dependency)


# Loop over files w/ dependencies
dependencies = copy.deepcopy(config["dependencies"])
for file in dependencies:

    # Loop over dependencies for this file
    for dependency in dependencies[file]:

        # Name of the dummy file we'll use to track the upload/download
        dot_zenodo_file = posix(f"{dependency}.zenodo")

        # If this is a Zenodo dataset, make it a dependency of the figure
        if dot_zenodo_file in files.dot_zenodo:
            if not f"{dependency}.zenodo" in config["dependencies"][file]:
                config["dependencies"][file].append(f"{dependency}.zenodo")