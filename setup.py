"""Install script for `showyourwork`."""
from pathlib import Path
from setuptools import find_packages, setup

setup(
    name="showyourwork",
    author="Rodrigo Luger",
    author_email="rodluger@gmail.com",
    url="https://github.com/showyourwork/showyourwork",
    description="A workflow for open-source scientific articles",
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(),
    use_scm_version={
        "write_to": Path("showyourwork") / "_showyourwork_version.py",
        "write_to_template": '__version__ = "{version}"\n',
    },
    install_requires=[
        "setuptools_scm",
        "jinja2>=2.11.1",
        "click",
        "pyyaml",
        "requests",
        "cookiecutter",
        "packaging",
        "lastversion",
    ],
    extras_require={
        "tests": ["pytest", "pytest-asyncio-cooperative"],
    },
    entry_points={
        "console_scripts": ["showyourwork=showyourwork.cli:entry_point"]
    },
    setup_requires=["setuptools_scm"],
    include_package_data=True,
    zip_safe=False,
)