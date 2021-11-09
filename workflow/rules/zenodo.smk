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

    #: Contents of each dataset (if it's a tarball)
    deposit_contents = {}

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
    zenodo.script[dependency] = dep_props.get("script", files.unknown)
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
    zenodo.deposit_contents[dependency] = dep_props.get("contents", [])

    # Download settings
    zenodo.zenodo_id[dependency] = dep_props.get("id", None)

    # Checks
    if zenodo.script[dependency] != files.unknown and zenodo.zenodo_id[dependency]:
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
    if zenodo.zenodo_id[dependency]:
        if not dependency in files.zenodo_files_manual:
            files.zenodo_files_manual.append(dependency)
    # Track upload-download files
    else:
        if not dependency in files.zenodo_files_auto:
            files.zenodo_files_auto.append(dependency)

    # Set up rules to compress and extract the tarballs
    tf = zenodo.file_name[dependency]
    if not tf.endswith(".tar.gz"):
        # Not a tarball
        if len(zenodo.deposit_contents[dependency]) > 0:
            raise ValueError(
                f"File `{tf}` is not a tarball, so `contents` cannot be specified in `showyourwork.yml`."
            )
    else:
        # This is a tarball
        if len(zenodo.deposit_contents[dependency]) == 0:
            raise ValueError(
                f"Must specify `contents` for tar file `{tf}` in `showyourwork.yml`."
            )

        # Make the tarball an explicit dependency of any figure that
        # depends on its contents; this ensures we get the Zenodo link
        # next to the figure caption when building the PDF
        for file in config["dependencies"]:
            for dep in config["dependencies"][file]:
                if dep in zenodo.deposit_contents[dependency]:
                    if dependency not in config["dependencies"][file]:
                        config["dependencies"][file].append(dependency)

        if config["CI"]:
            # Dynamically create a rule to unpack the tarball
            rule:
                input:
                    dependency
                output:
                    zenodo.deposit_contents[dependency]
                shell:
                    "tar -xzvf {input}"
        else:
            # Dynamically create a rule to generate the tarball
            rule:
                input:
                    zenodo.deposit_contents[dependency]
                output:
                    dependency
                shell:
                    "tar -czvf {output} {input}"
            

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

# On CI, we delete all datasets at the end so we don't cache them;
# this ensures we're always generating the figures based on the
# latest version of the deposit.
if config["CI"]:
    rule remove_zenodo_datasets:
        run:
            for file in files.zenodo_files_auto:
                if Path(file).exists():
                    Path(file).unlink()