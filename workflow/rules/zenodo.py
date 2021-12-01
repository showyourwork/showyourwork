"""
Collect information about all the datasets we'll need to upload to/download from Zenodo.

"""
from pathlib import Path
import copy
from showyourwork.workflow.helpers.exceptions import ShowyourworkException
from sphinx_mock import *
from collections import OrderedDict
import jinja2
import re


class zenodo:
    """
    Class containing metadata for all Zenodo deposits.

    """

    # Metadata for this class
    name = "zenodo"
    default_token_name = "ZENODO_TOKEN"
    url = "zenodo.org"

    #: Script for generating each dataset
    script = {}

    #: File name for each dataset
    file_name = {}

    #: Path to each dataset
    file_path = {}

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

    #: Zenodo id input by the user
    deposit_input_id = {}

    #: Zenodo deposit concept id
    deposit_concept_id = {}

    #: Zenodo deposit version id
    deposit_version_id = {}

    #: Zenodo deposit draft id
    deposit_draft_id = {}


class zenodo_sandbox(zenodo):
    """
    Class containing metadata for all Zenodo Sandbox deposits.

    """

    # Metadata for this class
    name = "zenodo_sandbox"
    default_token_name = "ZENODO_SANDBOX_TOKEN"
    url = "sandbox.zenodo.org"


# Get repo name for Zenodo metadata
repo = "/".join(get_repo_url().split("/")[-2:])


# List of dynamic rules, to be included in the `Snakefile`
dynamic_rules = []


# Loop through `zenodo` and `zenodo_sandbox` entries in the yaml file
for zrepo in [zenodo, zenodo_sandbox]:

    # Loop over datasets
    for dataset in config[zrepo.name]:

        # Get the dataset name & properties
        if type(dataset) is str:
            dependency = dataset
            dep_props = config[zrepo.name][dependency]
        elif type(dataset) is OrderedDict:
            dependency = list(dataset)[0]
            dep_props = dict(dataset[dependency])
        else:
            raise ValueError(
                f"Cannot parse `{zrepo.name}` entry in `showyourwork.yml`."
            )

        # Settings
        zrepo.script[dependency] = dep_props.get("script", files.unknown)
        zrepo.file_name[dependency] = str(Path(dependency).name)
        zrepo.file_path[dependency] = str(Path(dependency).parent)
        zrepo.token_name[dependency] = dep_props.get(
            "token_name", zrepo.default_token_name
        )
        zrepo.deposit_title[dependency] = dep_props.get(
            "title", f"{repo}:{zrepo.file_name[dependency]}"
        )
        zrepo.deposit_description[dependency] = dep_props.get(
            "description", f"File uploaded from {repo}."
        )
        zrepo.deposit_creators[dependency] = dep_props.get(
            "creators", get_repo_url().split("/")[-2]
        )
        zrepo.deposit_contents[dependency] = dep_props.get("contents", [])

        # Name of the dummy file we'll use to track the upload/download
        dot_zenodo_file = posix(f"{dependency}.zenodo")

        # Make the deposit a dependency of the PDF
        # so we can add Zenodo links to the figure caption
        if not dot_zenodo_file in files.dot_zenodo:
            files.dot_zenodo.append(dot_zenodo_file)

        # Determine if `id` is a concept or version Zenodo id.
        # This step requires an internet connection the first time it is run
        # for a given ID; subsequent runs read the cached info (which should
        # never change!)
        id = dep_props.get("id", None)
        if id is None:
            raise ShowyourworkException(
                f"Zenodo dependency {dependency} does not have an id.",
                context=f"Please provide an id for the dependency {dependency} "
                "in the `showyourwork.yml` config file. Refer to the docs for "
                "more information.",
                brief=f"Zenodo dependency {dependency} does not have an id.",
            )

        cache_file = abspaths.temp / f"{id}.{zrepo.url}"
        if not cache_file.exists():
            # Function defined in `../helpers/zenodo.py`
            id_type = get_id_type(
                dep_props.get("id", None), zrepo.url, zrepo.token_name[dependency]
            )
            with open(cache_file, "w") as f:
                print(id_type, file=f)
        else:
            with open(cache_file, "r") as f:
                id_type = f.readline().replace("\n", "")

        # Require that the user provides either a
        # concept id (for upload/download) or a
        # version id (download only)
        if id_type == "concept":

            # Showyourwork manages this dependency
            if not dependency in files.zenodo_files_auto:
                files.zenodo_files_auto.append(dependency)

        elif id_type == "version":

            # This dependency was manually uploaded to Zenodo
            # and is static (download only)
            if not dependency in files.zenodo_files_manual:
                files.zenodo_files_manual.append(dependency)

            # Disallow metadata specification in this case,
            # since it will either be wrong or redundant; the
            # id is all we need!
            if (
                "title" in dep_props
                or "description" in dep_props
                or "creators" in dep_props
            ):
                raise ShowyourworkException(
                    "The `title`, `description`, and `creator` fields cannot be "
                    "provided for a Zenodo version id.",
                    brief="The `title`, `description`, and `creator` fields cannot be "
                    "provided for a Zenodo version id.",
                    context=f"The Zenodo id {id_dict['input']} specified in the config "
                    f"file for dependency {dependency} is an id tied to a specific version "
                    "of a deposit. As such, there is no need to specify the `title`, "
                    "`description`, or `creators` fields for this deposit, since they "
                    "are already populated on Zenodo. Please remove them from the "
                    "config file. Alternatively, if you would like showyourwork to "
                    "manage the deposit, and to upload a new version when the dependency "
                    "changes, you should provide a *concept* id for this record. "
                    "Please refer to the docs for more information on this, or "
                    "read more about Zenodo versioning here: https://help.zenodo.org/#versioning.",
                )

        else:

            raise ShowyourworkException(
                f"The Zenodo id {id} specified in the config "
                f"file for dependency {dependency} "
                "is not a valid concept or version id.",
                brief=f"The Zenodo id {id} is not a valid concept or version id.",
                context=f"The Zenodo id {id} specified in the config "
                f"file for dependency {dependency} "
                "does not seem to be a valid concept or version id. Zenodo dependencies "
                "specified in the `showyourwork.yml` config file must include an `id` that "
                "corresponds to a specific published version of a record (a `version id`) "
                "or an `id` that corresponds to *all* versions of a record (known as a "
                "`concept id`). In the former case, the dependency is considered to be "
                "static: showyourwork will only ever download it (and not attempt to "
                "regenerate it if it is stale). In the latter case, the "
                "dependency is considered to be dynamic; if any of the inputs have changed, "
                "and the workflow is being run locally (as opposed to on GitHub Actions), "
                "showyourwork will attempt to regenerate it and upload a new version to "
                "Zenodo (if the user has the right authentication). "
                "The easiest way to obtain a new concept id is to run `make reserve` on "
                "the command line, or to find the DOI to cite all versions on the webpage "
                "for the Zenodo deposit. You can read more about versioning here: "
                "https://help.zenodo.org/#versioning.",
            )

        # Set up rules to compress and extract the tarballs
        tf = zrepo.file_name[dependency]
        if not tf.endswith(".tar.gz"):

            # Not a tarball
            if len(zrepo.deposit_contents[dependency]) > 0:
                raise ValueError(
                    f"File `{tf}` is not a tarball, so `contents` cannot be specified in `showyourwork.yml`."
                )

        else:

            # This is a tarball
            if len(zrepo.deposit_contents[dependency]) == 0:
                raise ValueError(
                    f"Must specify `contents` for tar file `{tf}` in `showyourwork.yml`."
                )

            # Jinja env
            env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(
                    abspaths.workflow / "resources" / "templates"
                ),
            )

            # Dynamically create a rule to extract the tarball
            extract_rulename = re.sub("[^0-9a-zA-Z]+", "_", f"extract_{dependency}")
            with open(abspaths.temp / f"{extract_rulename}.smk", "w") as f:
                smk = env.get_template("extract.smk").render(
                    rulename=extract_rulename,
                    input=dependency,
                    contents=zrepo.deposit_contents[dependency],
                )
                print(smk, file=f)
            dynamic_rules.append(extract_rulename)

            # Dynamically create a rule to generate the datasets,
            # but only if the user provided a `script`
            if zrepo.script[dependency] != files.unknown:
                generate_rulename = re.sub(
                    "[^0-9a-zA-Z]+", "_", f"run_{zrepo.script[dependency]}"
                )
                zenodo_script_name = Path(zrepo.script[dependency]).name
                zenodo_script_path = Path(zrepo.script[dependency]).parents[0]

                # TODO: Allow non-python scripts here!
                shell_cmd = f"cd {zenodo_script_path} && python {zenodo_script_name}"

                with open(abspaths.temp / f"{generate_rulename}.smk", "w") as f:
                    smk = env.get_template("run.smk").render(
                        rulename=generate_rulename,
                        input=zrepo.script[dependency],
                        contents=zrepo.deposit_contents[dependency],
                        shell_cmd=shell_cmd,
                    )
                    print(smk, file=f)
                dynamic_rules.append(generate_rulename)
            else:
                generate_rulename = None

            # Dynamically create a rule to compress the tarball
            compress_rulename = re.sub("[^0-9a-zA-Z]+", "_", f"compress_{dependency}")
            with open(abspaths.temp / f"{compress_rulename}.smk", "w") as f:
                smk = env.get_template("compress.smk").render(
                    rulename=compress_rulename,
                    output=dependency,
                    contents=zrepo.deposit_contents[dependency],
                )
                print(smk, file=f)
            dynamic_rules.append(compress_rulename)

            # Specify a `ruleorder` for each rule to resolve the ambiguity
            if generate_rulename is not None:
                ruleorder_rulename = re.sub(
                    "[^0-9a-zA-Z]+", "_", f"ruleorder_{dependency}"
                )
                with open(abspaths.temp / f"{ruleorder_rulename}.smk", "w") as f:

                    if config["CI"]:

                        # TODO: More rule order entries?
                        smk = "\n".join(
                            [
                                f"ruleorder: {extract_rulename} > {generate_rulename}",
                            ]
                        )

                    else:

                        # TODO: More rule order entries?
                        smk = "\n".join(
                            [
                                f"ruleorder: {generate_rulename} > {extract_rulename}",
                            ]
                        )

                    print(smk, file=f)
                dynamic_rules.append(ruleorder_rulename)

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