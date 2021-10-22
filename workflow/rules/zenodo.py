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

    #: List of dependencies of each dataset
    generate_deps = {}

    #: Shell command for generating each dataset
    generate_shell = {}

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
        dep_name = dataset
        dep_props = config["zenodo"][dep_name]
    elif type(dataset) is OrderedDict:
        dep_name = list(dataset)[0]
        dep_props = dict(dataset[dep_name])
    else:
        raise ValueError("Cannot parse `zenodo` entry in `showyourwork.yml`.")

    # User settings for generating / uploading / downloading the dependency
    generate = dep_props.get("generate", {})
    download = dep_props.get("download", {})
    if generate and download:
        raise ValueError(
            "Cannot specify both `generate` and `download` for a Zenodo dataset."
        )
    elif not generate and not download:
        raise ValueError(
            f"Must specify either `generate` or `download` for Zenodo dataset {dep_name}."
        )

    # Generate & upload settings
    zenodo.generate_deps[dep_name] = [
        posix(Path("src/figures") / gd) for gd in generate.get("dependencies", [])
    ]
    zenodo.generate_shell[dep_name] = generate.get("shell", None)
    if generate and zenodo.generate_shell[dep_name] is None:
        raise ValueError(f"Please provide a `shell` command for dependency {dep_name}.")
    zenodo.file_name[dep_name] = str(Path(dep_name).name)
    zenodo.file_path[dep_name] = str(Path("src/figures") / Path(dep_name).parent)
    zenodo.sandbox[dep_name] = generate.get("sandbox", False)
    zenodo.token_name[dep_name] = generate.get("token_name", "ZENODO_TOKEN")
    zenodo.deposit_title[dep_name] = generate.get("title", f"{repo}:{dep_name}")
    zenodo.deposit_description[dep_name] = generate.get(
        "description", f"File uploaded from {repo}."
    )
    zenodo.deposit_creators[dep_name] = generate.get(
        "creators", get_repo_url().split("/")[-2]
    )

    # Download settings
    zenodo.zenodo_id[dep_name] = download.get("id", None)
    download_sandbox = download.get("sandbox", False)
    if download_sandbox:
        zenodo.zenodo_url[dep_name] = "sandbox.zenodo.org"
    else:
        zenodo.zenodo_url[dep_name] = "zenodo.org"
    if download and zenodo.zenodo_id[dep_name] is None:
        raise ValueError(f"Please provide a Zenodo `id` for dependency {dep_name}.")

    # Name of the dummy file we'll use to track the upload/download
    dot_zenodo_file = posix(relpaths.figures / f"{dep_name}.zenodo")

    # Make the deposit a dependency of the PDF
    # so we can add Zenodo links to the figure caption
    if not dot_zenodo_file in files.dot_zenodo:
        files.dot_zenodo.append(dot_zenodo_file)

    # Track download-only files
    if download and not dep_name in files.zenodo_files_manual:
        files.zenodo_files_manual.append(dep_name)

    # Track upload-download files
    if generate and not dep_name in files.zenodo_files_auto:
        files.zenodo_files_auto.append(dep_name)


# Loop over figures
figure_dependencies = copy.deepcopy(config["figure_dependencies"])
for fig in figure_dependencies:

    # Loop over dependencies for this figure
    for dep_name in figure_dependencies[fig]:

        # Name of the dummy file we'll use to track the upload/download
        dot_zenodo_file = posix(relpaths.figures / f"{dep_name}.zenodo")

        # If this is a Zenodo dataset, make it a dependency of the figure
        if dot_zenodo_file in files.dot_zenodo:
            if not f"{dep_name}.zenodo" in config["figure_dependencies"][fig]:
                config["figure_dependencies"][fig].append(f"{dep_name}.zenodo")