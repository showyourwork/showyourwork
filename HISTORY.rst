.. :changelog:

0.3.0 (June 2022)
+++++++++++++++++

- Brand new release, featuring a complete re-write of the workflow. Below is 
  a list of the *major* changes:
- ``showyourwork!`` is now a pip-installable Python package. It is no longer 
  a git submodule. The ``showyourwork.yml`` file must now specify which version
  of ``showyourwork`` to use.
- There is no longer a ``Makefile``. Articles should be built using the 
  ``showyourwork`` command, which creates and activates a clean conda environment
  containing ``Snakemake`` and all dependencies needed to run the pipeline for
  each article. Users therefore no longer need to install ``Snakemake`` or ``mamba``
  in the base environment.
- New articles can now be created using the ``showyourwork setup`` command rather than
  via a GitHub repository template.
- The syntax for many of the settings in the ``showyourwork.yml`` config file has
  changed, particularly for specifying Zenodo datasets.
- The directory structure for article repositories has changed slightly. Figure scripts
  should now be placed in the ``src/scripts`` directory (used to be ``src/figures``).
  The TeX files should now be placed in the ``src/tex`` directory (used to be ``src``).
  Figure output files should now be generated in ``src/tex/figures`` (used to be ``src/figures``).
- To help with the transition to the new directory structure, new repositories include a
  file ``src/scripts/paths.py`` that specifies absolute Pathlib paths to the main directories
  in the repository.
- Figure scripts are no longer inferred from the ``\label`` command in LaTeX. Instead,
  users should specify the script associated with a given figure using the new ``\script``
  command.
- Added support for programmatically-generated files that can be included in the TeX
  manuscript via the new ``\variable`` command.
- Users must now manually add ``\usepackage{showyourwork}`` to their LaTeX manuscript.
- Overhauled the way ``showyourwork`` integrates with Zenodo. Static datasets should now
  be specified using their full DOI. Dynamic datasets are deprecated in favor of "cached"
  datasets. These are intermediate results that get cached on Zenodo Sandbox alongside
  a hash of the rules and all upstream dependencies used to generate them, making it
  possible to automatically restore results from the cloud in a way that preserves the
  full reproducibility of the workflow.
- Implemented (experimental) integration with Overleaf projects, allowing users to pull
  changes to the manuscript and push changes to the figures.
- Drastically improved the command-line interface, suppressing most of the noise generated
  by Snakemake in favor of succinct informational messages describing the build process.
  All messages now get logged to files in ``.showyourwork/logs`` for easier debugging.
  Similarly, improved error catching and added informational error messages for most of the failure
  modes of the workflow.
- Several other tweaks, bugfixes, and improvements. Lots of changes to the back end to make
  ``showyourwork`` easier to develop, maintain, and extend!

0.2.3 (2022-02-21)
++++++++++++++++++

- Bump ``jinja2`` version to fix issue with ``markupsafe``
- Changelog: `v0.2.2...v0.2.3 <https://github.com/showyourwork/showyourwork/compare/v0.2.2...v0.2.3>`_

0.2.2 (2022-01-05)
++++++++++++++++++

- Tweaks to the logo (now an actual font!)
- Added option to exclude files from article cache on CI
- Fixed behavior of figures labeled with an asterisk (e.g., ``\label{fig*:...}``)
- Changelog: `v0.2.1...v0.2.2 <https://github.com/showyourwork/showyourwork/compare/v0.2.1...v0.2.2>`_

0.2.1 (2021-12-18)
++++++++++++++++++

- Implemented custom DAG generation (cleaner, prettier)
- Added several entries to the Projects page on the docs
- Streamlined arXiv tarball generation step
- Added a basic `make lint` command to check for repo issues
- Users can now disable caching on CI by setting the cache number to ``null``
- Switch to installing ``graphviz`` with ``conda``
- Pinned all ``conda`` dependencies to specific versions
- Changelog: `v0.2.0...v0.2.1 <https://github.com/showyourwork/showyourwork/compare/v0.2.0...v0.2.1>`_

0.2.0 (2021-12-07)
++++++++++++++++++

- Major changes to the Zenodo interface! Please read the Zenodo section of the documentation on the
  `showyourwork.yml file <https://showyourwork.readthedocs.io/en/v0.2.0/config/>`_
  for details on what changed. The tl;dr is that all deposits now require either a **concept** or
  a **version** id (read more about that here: https://help.zenodo.org/#versioning); this id
  now uniquely identifies the deposit (previously, we relied on the uniqueness of the deposit
  title & creators).
- Added support for non-Python scripts to generate datasets and other dependencies
- Fixed issue with large datasets lingering in the arxiv tarball
- Added a ``make fast`` option to reproduce the results locally w/out running expensive steps
- Signficiant improvements to the documentation, now with detailed info on LaTeX features
- Added a ``marginicon`` command for custom margin icons next to figures
- Added a progress bar to Zenodo uploads
- Now re-downloading datasets on GitHub Actions if a newer version is available
- Better path resolution when extracting tarballs
- Changelog: `v0.1.35...v0.2.0 <https://github.com/showyourwork/showyourwork/compare/v0.1.35...v0.2.0>`_

0.1.35 (2021-11-22)
+++++++++++++++++++

- Fixed issue with unnecessary downloads of Zenodo datasets on CI.
- Fixed issue that prevented downloading the Zenodo datasets if the Zenodo API key belonged to someone other than the creator of the deposit.
- Fixed issue that caused the conda env creation to fail. We are now installing snakemake-minimal and pinning the mamba version; this is only a temporary solution.
- Changelog: `v0.1.34...v0.1.35 <https://github.com/showyourwork/showyourwork/compare/v0.1.34...v0.1.35>`_

0.1.34 (2021-11-18)
+++++++++++++++++++

- Now uploads a build artifact whenever the workflow fails on CI for easier debugging.
- Changelog: `v0.1.33...v0.1.34 <https://github.com/showyourwork/showyourwork/compare/v0.1.33...v0.1.34>`_

0.1.33 (2021-11-17)
+++++++++++++++++++

- Fixed issue with unnecessary reruns of figure scripts
- Changelog: `v0.1.32...v0.1.33 <https://github.com/showyourwork/showyourwork/compare/v0.1.32...v0.1.33>`_

0.1.32 (2021-11-17)
+++++++++++++++++++

- Fixed `issue #57 <https://github.com/showyourwork/showyourwork/issues/57>`_.
- Better documentation for the config file and the ``Snakefile``.
- Changelog: `v0.1.31...v0.1.32 <https://github.com/showyourwork/showyourwork/compare/v0.1.31...v0.1.32>`_

0.1.31 (2021-11-15)
+++++++++++++++++++

- Migrated to the new Zenodo API; previously the limit for uploading files was 100 MB (now 50 GB).
- Changelog: `v0.1.30...v0.1.31 <https://github.com/showyourwork/showyourwork/compare/v0.1.30...v0.1.31>`_

0.1.30 (2021-11-12)
+++++++++++++++++++

- Added an example on how to use jinja templating to simplify the ``showyourwork.yml`` config file.
- **Developers:** Undo the reset build cache operation from the previous patch, since this causes race conditions when
  accessing the cache during the unit tests (since we are concurrently running dozens of actions on a single repo!)
- Changelog: `v0.1.29...v0.1.30 <https://github.com/showyourwork/showyourwork/compare/v0.1.29...v0.1.30>`_

0.1.29 (2021-11-10)
+++++++++++++++++++

- Support for creation/download of Zenodo tarballs.
- Implements the idea in `#48 <https://github.com/showyourwork/showyourwork/issues/48>`_ for specifying custom manuscript dependencies.
- Bugfix for rules that subclass the main showyourwork figure rule.
- Implements the idea in `#47 <https://github.com/showyourwork/showyourwork/issues/47>`_ for custom Zenodo dataset generation.
- **Developers:** Now resetting the build cache before each unit test on ``showyourwork-example`` and then re-running the cached build.
- Changelog: `v0.1.28...v0.1.29 <https://github.com/showyourwork/showyourwork/compare/v0.1.28...v0.1.29>`_

0.1.28 (2021-11-09)
+++++++++++++++++++

- Added support for non-Python scripts; users can now define instructions in the YAML config file to execute other kinds of scripts.
- Implemented better error messages when figure scripts fail.
- Allow users to specify a ``graphicspath`` for all figures in the document.
- Allow users to customize the name of the manuscript (it no longer needs to be called ``ms.tex``).
- Changelog: `v0.1.27...v0.1.28 <https://github.com/showyourwork/showyourwork/compare/v0.1.27...v0.1.28>`_

0.1.27 (2021-11-03)
+++++++++++++++++++

- Added support for installing a minimal TeX distribution so that TeX can be rendered in matplotlib; see Custom workflows.
- Changelog: `v0.1.26...v0.1.27 <https://github.com/showyourwork/showyourwork/compare/v0.1.26...v0.1.27>`_

0.1.26 (2021-11-02)
+++++++++++++++++++

- Fixed issue causing documentation builds to fail
- Changelog: `v0.1.25...v0.1.26 <https://github.com/showyourwork/showyourwork/compare/v0.1.25...v0.1.26>`_

0.1.25 (2021-11-02)
+++++++++++++++++++

- Fixed issue that prevented ORCID badges from showing up when building the PDF on GitHub Actions
- Changelog: `v0.1.24...v0.1.25 <https://github.com/showyourwork/showyourwork/compare/v0.1.24...v0.1.25>`_

0.1.24 (2021-11-02)
+++++++++++++++++++

- Fixed issue with ``os.get_terminal_size`` breaking CI builds when displaying error messages
- Changelog: `v0.1.23...v0.1.24 <https://github.com/showyourwork/showyourwork/compare/v0.1.23...v0.1.24>`_

0.1.23 (2021-11-02)
+++++++++++++++++++

- Added explicit support for MNRAS and A&A LaTeX document classes
- Improved support for new Apple M1 chips
- Fixed options clash for package ``hyperref``
- Changelog: `v0.1.22...v0.1.23 <https://github.com/showyourwork/showyourwork/compare/v0.1.22...v0.1.23>`_

0.1.22 (2021-11-02)
+++++++++++++++++++

- Updated LaTeX package ``fontawesome`` to ``fontawesome5``
- **Developers:** Can now run tests on PR branches to generate `showyourwork-example-dev` branches
- Changelog: `v0.1.21...v0.1.22 <https://github.com/showyourwork/showyourwork/compare/v0.1.21...v0.1.22>`_

0.1.21 (2021-11-01)
+++++++++++++++++++

- Fixed minor issue with error messages for custom figures
- Improved documentation page on projects that use ``showyourwork``
- Changelog: `v0.1.20...v0.1.21 <https://github.com/showyourwork/showyourwork/compare/v0.1.20...v0.1.21>`_

0.1.20 (2021-10-28)
+++++++++++++++++++

- Fixed issue with figure link formatting when enabling linenumbers in AASTeX
- Made `arxiv_tarball_exclude` paths relative to the repository root
- Added a `make update` option to update ``showyourwork`` to the latest release.
- Changelog: `v0.1.19...v0.1.20 <https://github.com/showyourwork/showyourwork/compare/v0.1.19...v0.1.20>`_

0.1.19 (2021-10-25)
+++++++++++++++++++

- Fixed typo that causes Zenodo integration to fail.
- Changelog: `v0.1.18...v0.1.19 <https://github.com/showyourwork/showyourwork/compare/v0.1.18...v0.1.19>`_

0.1.18 (2021-10-25)
+++++++++++++++++++

- Added more informative error messages that are displayed at the very *end* of the build logs.
  Still more work to be done on this front, but error logs should now be much easier to parse.
- Implemented the new Zenodo config structure in the ``showyourwork.yml`` file, as per
  `#31 <https://github.com/showyourwork/showyourwork/issues/31>`_.
- Changelog: `v0.1.17...v0.1.18 <https://github.com/showyourwork/showyourwork/compare/v0.1.17...v0.1.18>`_

0.1.17 (2021-10-22)
+++++++++++++++++++

- Changed the way Zenodo dependencies are provided in the ``showyourwork.yml`` file. Dependencies like
  datasets should still be listed as entries under the corresponding figure scripts in ``figure_dependencies``,
  but all information on how to ``generate`` or ``download`` them should now go in a separate top-level
  ``zenodo:`` key. This makes it much easier to, e.g., specify datasets used by multiple figures.
  Please see the ``Custom workflows`` section of the docs for more information.
- Improved the API documentation.
- Changelog: `v0.1.16...v0.1.17 <https://github.com/showyourwork/showyourwork/compare/v0.1.16...v0.1.17>`_

0.1.16 (2021-10-22)
+++++++++++++++++++

- **Template repo update:** Pared down the ``Makefile`` in the template repository. This now calls
  a ``Makefile`` in the ``showyourwork`` submodule (this repo), which contains all the directives.
  This makes it easier to improve/update the workflow, since we can just update ``showyourwork``.
- Changelog: `v0.1.15...v0.1.16 <https://github.com/showyourwork/showyourwork/compare/v0.1.15...v0.1.16>`_

0.1.15 (2021-10-21)
+++++++++++++++++++

- **Template repo update:** Added options to the ``Makefile`` to generate a report and a DAG.
  Added a submodule setup check; if the user didn't init the showyourwork submodule, does it
  automatically before building.
- Changelog: `v0.1.14...v0.1.15 <https://github.com/showyourwork/showyourwork/compare/v0.1.14...v0.1.15>`_

0.1.14 (2021-10-21)
+++++++++++++++++++

- Remove duplicated Zenodo links from figure captions
- Changelog: `v0.1.13...v0.1.14 <https://github.com/showyourwork/showyourwork/compare/v0.1.13...v0.1.14>`_

0.1.13 (2021-10-21)
+++++++++++++++++++

- Fixed API documentation
- Fixed error with `arxiv_tarball_exclude` and arxiv tarball issue (`#21 <https://github.com/showyourwork/showyourwork/issues/21>`_)
- Changelog: `v0.1.12...v0.1.13 <https://github.com/showyourwork/showyourwork/compare/v0.1.12...v0.1.13>`_

0.1.12 (2021-10-20)
+++++++++++++++++++

- Revert code that prevents the Snakefile from being loaded more than once. Turns out that is
  expected behavior, and is required in order for the module import syntax to work!
- Switched to adding checks within the ``zenodo.py`` script to prevent dependencies from getting
  ingested multiple times.
- Changelog: `v0.1.11...v0.1.12 <https://github.com/showyourwork/showyourwork/compare/v0.1.11...v0.1.12>`_

0.1.11 (2021-10-20)
+++++++++++++++++++

- Fix bug preventing figures from being cached properly when one script generates multiple figures
- Fixed issues due to Snakefile being loaded multiple times
- Auto-populate the ``projects`` page on the docs via a GitHub API search on every release
- Changelog: `v0.1.10...v0.1.11 <https://github.com/showyourwork/showyourwork/compare/v0.1.10...v0.1.11>`_

0.1.10 (2021-10-20)
+++++++++++++++++++

- Cleaned up the workflow, separating rules into their own files with better documentation.
- Added a fix for nested figures (figures under subdirectories in the ``src/figures`` folder).
- Fixed issue with multiple Zenodo datasets causing the build to fail.
- Added support for figures in figure* environments.
- Fixed issue with occasional missing </HTML> closing tags in the showyourwork XML tree.
- Added some API documentation; more coming soon.
- Changelog: `v0.1.9...v0.1.10 <https://github.com/showyourwork/showyourwork/compare/v0.1.9...v0.1.10>`_

0.1.9 (2021-10-18)
++++++++++++++++++

- **Template repo update:** Added a ``Makefile`` for quick article generation; added docs on how to use it.
- Changelog: `v0.1.8...v0.1.9 <https://github.com/showyourwork/showyourwork/compare/v0.1.8...v0.1.9>`_

0.1.8 (2021-10-18)
++++++++++++++++++

- Added "One script, multiple figures" example
- Improved the documentation for script dependencies and datasets
- Fixed a bug when downloading deposits from Zenodo
- Added release testing
- Changelog: `v0.1.7...v0.1.8 <https://github.com/showyourwork/showyourwork/compare/v0.1.7...v0.1.8>`_

0.1.7 (2021-10-18)
++++++++++++++++++

- Added explicit support for Zenodo-hosted datasets.
- **Template repo update:** Added the environment variable ``ZENODO_TOKEN`` to ``.github/workflows/showyourwork.yml``.
- Changelog: `v0.1.6...v0.1.7 <https://github.com/showyourwork/showyourwork/compare/v0.1.6...v0.1.7>`_

0.1.6 (2021-10-14)
++++++++++++++++++

- Added documentation for the ``expensive-figure`` example.
- Changelog: `v0.1.5...v0.1.6 <https://github.com/showyourwork/showyourwork/compare/v0.1.5...v0.1.6>`_

0.1.5 (2021-10-14)
++++++++++++++++++

- Added the ``expensive-figure`` example for computationally expensive figure generation.
- Changelog: `v0.1.4...v0.1.5 <https://github.com/showyourwork/showyourwork/compare/v0.1.4...v0.1.5>`_

0.1.4 (2021-10-13)
++++++++++++++++++

- Initial release of the workflow.