import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--repo-path")
parser.add_argument("--output")
parser.add_argument("files", nargs="*")
args = parser.parse_args()

repo_path = Path(args.repo_path)
output_directory = Path(args.output)

for file in args.files:
    output_filename = output_directory / file
    output_filename = output_filename.with_name(output_filename.name + ".png")
    output_filename.parent.mkdir(parents=True, exist_ok=True)

    # Determine the original aspect ratio of the image
    r = subprocess.run(
        ["convert", file, "-format", "%[fx:w/h]", "info:"],
        check=True,
        capture_output=True,
    )
    aspect = float(r.stdout.decode())
    width = aspect * 500

    # Resize the image and save it as a png in the output directory
    subprocess.run(
        [
            "convert",
            "-resize",
            f"{width}x500",
            "-background",
            "white",
            "-alpha",
            "remove",
            "-bordercolor",
            "black",
            "-border",
            "15",
            file,
            str(output_filename),
        ],
        check=True,
    )
