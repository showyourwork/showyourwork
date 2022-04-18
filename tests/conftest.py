import shutil
from pathlib import Path


def pytest_sessionstart(session):
    # Clean the sandbox
    for folder in (Path(__file__).parents[0] / "sandbox").glob("*"):
        if folder.is_dir():
            shutil.rmtree(folder)