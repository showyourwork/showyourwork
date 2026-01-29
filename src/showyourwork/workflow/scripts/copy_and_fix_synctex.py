import gzip
from pathlib import Path

from showyourwork import paths

if __name__ == "__main__":
    # Read the generated SyncTeX
    with gzip.open(snakemake.input[0], "rb") as f:
        data = f.read()

    # Rewrite the input paths to be relative to the user's tex directory if
    # that file exists
    with gzip.open(snakemake.output[0], "wb") as f:
        for line in data.split(b"\n"):
            if line.decode().startswith("Input:"):
                parts = line.decode().split(":", 2)
                path = parts[-1].strip()
                if len(path):
                    path = Path(path)
                    try:
                        path = path.relative_to(paths.user().compile)
                    except ValueError:
                        pass
                    else:
                        path = paths.user().tex / path
                        if path.exists():
                            line = ":".join(  # noqa
                                parts[:-1] + [path.as_posix()]
                            ).encode()
            line += b"\n"  # noqa
            f.write(line)
