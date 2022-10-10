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
        "setuptools_scm>=7.0.1",
        "jinja2>=2.11.1",
        "click>=8.1.3",
        "pyyaml>=6.0",
        "requests>=2.28.0",
        "cookiecutter>=2.1.0",
        "packaging>=21.3",
    ],
    extras_require={
        "tests": ["pytest>=7.0.0", "pytest-asyncio-cooperative>=0.28.0"]
    },
    entry_points={
        "console_scripts": ["showyourwork=showyourwork.cli:entry_point"]
    },
    setup_requires=["setuptools_scm>=7.0.1"],
    include_package_data=True,
    zip_safe=False,
)
