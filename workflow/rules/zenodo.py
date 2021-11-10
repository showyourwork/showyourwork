"""
Collect information about all the datasets we'll need to upload to/download from Zenodo.

"""
from pathlib import Path
import copy
from sphinx_mock import *
import jinja2
import re


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


# List of dynamic rules, to be included in the `Snakefile`
dynamic_rules = []


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

        # Jinja env
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(
                abspaths.workflow / "resources" / "templates"
            ),
        )

        if config["CI"] or zenodo.zenodo_id[dependency]:

            # Dynamically create a rule to unpack the tarball
            rulename = re.sub("[^0-9a-zA-Z]+", "_", f"extract_{dependency}")
            with open(abspaths.temp / f"{rulename}.smk", "w") as f:
                smk = env.get_template("extract.smk").render(
                    rulename=rulename,
                    input=dependency,
                    contents=zenodo.deposit_contents[dependency],
                )
                print(smk, file=f)
            dynamic_rules.append(rulename)

        else:

            if zenodo.script[dependency] != files.unknown:

                # Dynamically create a rule to generate the datasets if
                # the user provided a `script`
                rulename = re.sub(
                    "[^0-9a-zA-Z]+", "_", f"run_{zenodo.script[dependency]}"
                )
                zenodo_script_name = Path(zenodo.script[dependency]).name
                zenodo_script_path = Path(zenodo.script[dependency]).parents[0]
                shell_cmd = f"cd {zenodo_script_path} && python {zenodo_script_name}"
                with open(abspaths.temp / f"{rulename}.smk", "w") as f:
                    smk = env.get_template("run.smk").render(
                        rulename=rulename,
                        input=zenodo.script[dependency],
                        contents=zenodo.deposit_contents[dependency],
                        shell_cmd=shell_cmd,
                    )
                    print(smk, file=f)
                dynamic_rules.append(rulename)

            # Dynamically create a rule to generate the tarball
            rulename = re.sub("[^0-9a-zA-Z]+", "_", f"compress_{dependency}")
            with open(abspaths.temp / f"{rulename}.smk", "w") as f:
                smk = env.get_template("compress.smk").render(
                    rulename=rulename,
                    output=dependency,
                    contents=zenodo.deposit_contents[dependency],
                )
                print(smk, file=f)
            dynamic_rules.append(rulename)


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