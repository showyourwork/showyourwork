"""
Extracts files from a `.tar.gz` tarball.

"""
from pathlib import Path
import tarfile
import sys


# Hack to import our custom exception
sys.path.insert(0, str(Path(__file__).parents[1]))
from helpers.exceptions import ShowyourworkException


# Grab user info
tarball = Path(snakemake.input[0])
contents = [Path(file) for file in snakemake.output]

# Open the tarball
with tarfile.open(tarball) as f:

    # Get the actual contents of the tarball
    actual_contents = [Path(file) for file in f.getnames()]

    # Try to extract the user-specified files, one at a time
    for file in contents:

        if file in actual_contents:

            # The file exists and the path is correct
            actual_file = str(file)
            print(f"Extracting {tarball}/{actual_file} -> {file}...")
            f.extract(actual_file)

        elif Path(file.name) in actual_contents:

            # The file exists at the root
            actual_file = file.name
            print(f"Extracting {tarball}/{actual_file} -> {file}...")
            f.extract(actual_file, file.parents[0])

        elif file.relative_to(Path("src") / "data") in actual_contents:

            # The file is nested under a folder at the root
            actual_file = str(file.relative_to(Path("src") / "data"))
            print(f"Extracting {tarball}/{actual_file} -> {file}...")
            f.extract(actual_file, Path("src") / "data")

        else:

            raise ShowyourworkException(
                f"Unable to extract file `{file}` from tarball `{tarball}`.",
                script="extract.py",
                rule_name="extract",
                brief=f"An error occurred while attempting to unpack `{tarball}`.",
                context=f"The file `{file}` specified under the contents of the "
                f"tarball `{tarball}` in the `showyourwork.yml` config file could "
                "not be located within the tarball. Please check that the full path "
                "to the file is correct. See the docs page on the `showyourwork.yml` "
                "file for more information.",
            )