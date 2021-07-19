"""Install script for `showyourwork`."""
from setuptools import find_packages, setup
from pathlib import Path

setup(
    name="showyourwork",
    author="Rodrigo Luger",
    author_email="rodluger@gmail.com",
    url="https://github.com/rodluger/showyourwork.py",
    description="a workflow for reproducible scientific articles",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(),
    use_scm_version={
        "write_to": Path("showyourwork") / "showyourwork_version.py",
        "write_to_template": '__version__ = "{version}"\n',
    },
    install_requires=["setuptools_scm", "questionary"],
    entry_points={
        "console_scripts": ["showyourwork=showyourwork.entry_points:main"]
    },
    setup_requires=["setuptools_scm"],
    include_package_data=True,
    zip_safe=False,
)
