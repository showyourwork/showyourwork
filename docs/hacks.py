import subprocess
import sys
from pathlib import Path


# Add our module to the path
sys.path.insert(0, str(Path(__file__).absolute().parents[1]))


class RemoveSubmodulesHeading(StopIteration):
    pass


class RemoveContentsHeading(StopIteration):
    pass


# Create our API autodoc pages w/ a little customization
for file in Path("api").rglob("*.rst"):
    file.unlink()
subprocess.run("sphinx-apidoc -feMT -d 5 -o api ../showyourwork", shell=True)
with open("api/showyourwork.rst", "r") as f:
    lines = f.readlines()
lines[0] = "Developer API\n"
with open("api/showyourwork.rst", "w") as f:
    f.writelines(lines)
for file in Path("api").rglob("*.rst"):
    with open(file, "r") as f:
        lines = f.readlines()
    lines[0] = (
        lines[0]
        .replace(" package", "")
        .replace(" module", "")
        .replace("showyourwork.", "")
    )
    while 1:
        try:
            subpackages = False
            for l, line in enumerate(lines):
                if line.startswith("Subpackages"):
                    subpackages = True
                    raise RemoveContentsHeading
                elif line.startswith("Submodules"):
                    if subpackages:
                        raise RemoveSubmodulesHeading
                    else:
                        raise RemoveContentsHeading
                elif line.startswith("Contents"):
                    raise RemoveContentsHeading
        except RemoveSubmodulesHeading:
            lines = lines[: l - 1] + lines[l + 6 :]
        except RemoveContentsHeading:
            lines = lines[:l] + lines[l + 3 :]
        else:
            break
    with open(file, "w") as f:
        f.writelines(lines)