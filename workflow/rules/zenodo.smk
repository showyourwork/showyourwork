import os
from pathlib import Path


# Get repo name for Zenodo metadata
repo = "/".join(get_repo_url().split("/")[-2:])


# Are we on GitHub Actions?
ON_GITHUB_ACTIONS = os.getenv("CI", "false") == "true"


# Loop over figures
for fig in figure_dependencies:

    # Loop over dependencies for this figure
    for dep in figure_dependencies[fig]:


        # Get the dependency name and any instructions on how to generate it
        if type(dep) is OrderedDict:
            dep_name = list(dep)[0]
            if dep[dep_name] is None:
                # This is a static dependency, with no rules on how to generate
                # it or upload it to Zenodo. Let's move on.
                continue
            dep_props = dict(dep[dep_name])
        elif type(figure_dependencies[fig]) is dict:
            dep_name = dep
            dep_props = figure_dependencies[fig][dep]
        else:
            # This is a static dependency, with no rules on how to generate
            # it or upload it to Zenodo. Let's move on.
            continue


        # User settings for generating / uploading / downloading the dependency
        generate = dep_props.get("generate", {})
        download = dep_props.get("download", {})
        if generate and download:
            raise ValueError("Cannot specify both `generate` and `download` for a figure dependency.")

        # Generate & upload settings
        generate_deps = [POSIX(Path("src/figures") / gd) for gd in generate.get("dependencies", [])]
        generate_shell = generate.get("shell", None)
        if generate and generate_shell is None:
            raise ValueError(f"Please provide a `shell` command for dependency {dep_name}.")
        file_name = str(Path(dep_name).name)
        file_path = str(Path("src/figures") / Path(dep_name).parent)
        sandbox = generate.get("sandbox", False)
        token_name = generate.get("token_name", "ZENODO_TOKEN")
        deposit_title = generate.get("title", f"{repo}:{dep_name}")
        deposit_description = generate.get("description", f"File uploaded from {repo}.")
        deposit_creators = generate.get("creators", get_repo_url().split("/")[-2])
        
        # Download settings
        zenodo_id = download.get("id", None)
        if download and zenodo_id is None:
            raise ValueError(f"Please provide either a Zenodo `id` for dependency {dep_name}.")

        # Now either download the dataset (if an `id` was provided or if 
        # we are on GH Actions) or generate & upload it (locally).
        # Note that we do NOT cache the dataset on GH Actions, since we don't
        # want to risk running into the cache limit (5 GB I think).
        if download:

            rule:
                message:
                    f"Downloading dependency file {dep_name} from Zenodo..."
                output:
                    temp(f"src/figures/{dep_name}") if ON_GITHUB_ACTIONS else f"src/figures/{dep_name}",
                    POSIX(f"{file_path}/{file_name}.zenodo")
                shell:
                    " && ".join(
                        [
                            f"curl https://zenodo.org/record/{zenodo_id}/files/{dep_name} --output {{output[0]}}", 
                            f"echo 'https://zenodo.org/record/{zenodo_id}' > {file_path}/{file_name}.zenodo"
                        ]
                    )

        elif ON_GITHUB_ACTIONS:

            rule:
                message:
                    f"Downloading dependency file {dep_name} from Zenodo..."
                output:
                    temp(f"src/figures/{dep_name}"),
                    POSIX(f"{file_path}/{file_name}.zenodo")
                conda:
                    POSIX(USER / "environment.yml")
                params:
                    action="download",
                    file_name=file_name,
                    file_path=file_path,
                    deposit_title=deposit_title,
                    sandbox=sandbox,
                    token_name=token_name
                script:
                    "../scripts/zenodo.py"

        elif generate:

            rule:
                message:
                    f"Generating figure dependency file {dep_name}..."
                input:
                    generate_deps
                output:
                    f"src/figures/{dep_name}"
                conda:
                    POSIX(USER / "environment.yml")
                shell:
                    f"cd src/figures && {generate_shell}"

            rule:
                message:
                    f"Uploading dependency file {dep_name} to Zenodo..."
                input:
                    f"src/figures/{dep_name}"
                output:
                    POSIX(f"{file_path}/{file_name}.zenodo")
                conda:
                    POSIX(USER / "environment.yml")
                params:
                    action="upload",
                    file_name=file_name,
                    file_path=file_path,
                    deposit_title=deposit_title,
                    deposit_description=deposit_description,
                    deposit_creators=deposit_creators,
                    sandbox=sandbox,
                    token_name=token_name,
                    generate=generate,
                    repo_url="{}/tree/{}".format(get_repo_url(), get_repo_sha())
                script:
                    "../scripts/zenodo.py"

        # Make the deposit a dependency of the figure
        figure_dependencies[fig].append(f"{dep_name}.zenodo")

        # Make it a dependency of the PDF as well
        # so we can add Zenodo links to the figure caption
        zenodo_files.append(POSIX(FIGURES / f"{dep_name}.zenodo"))