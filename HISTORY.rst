.. :changelog:

0.1.19 (2021-10-25)
+++++++++++++++++++

- Fixed typo that causes Zenodo integration to fail.
- Changelog: `v0.1.18...v0.1.19 <https://github.com/rodluger/showyourwork/compare/v0.1.18...v0.1.19>`_

0.1.18 (2021-10-25)
+++++++++++++++++++

- Added more informative error messages that are displayed at the very *end* of the build logs.
  Still more work to be done on this front, but error logs should now be much easier to parse.
- Implemented the new Zenodo config structure in the ``showyourwork.yml`` file, as per
  `#31 <https://github.com/rodluger/showyourwork/issues/31>`_.
- Changelog: `v0.1.17...v0.1.18 <https://github.com/rodluger/showyourwork/compare/v0.1.17...v0.1.18>`_

0.1.17 (2021-10-22)
+++++++++++++++++++

- Changed the way Zenodo dependencies are provided in the ``showyourwork.yml`` file. Dependencies like
  datasets should still be listed as entries under the corresponding figure scripts in ``figure_dependencies``,
  but all information on how to ``generate`` or ``download`` them should now go in a separate top-level
  ``zenodo:`` key. This makes it much easier to, e.g., specify datasets used by multiple figures.
  Please see the ``Custom workflows`` section of the docs for more information.
- Improved the API documentation.
- Changelog: `v0.1.16...v0.1.17 <https://github.com/rodluger/showyourwork/compare/v0.1.16...v0.1.17>`_

0.1.16 (2021-10-22)
+++++++++++++++++++

- **Template repo update:** Pared down the ``Makefile`` in the template repository. This now calls
  a ``Makefile`` in the ``showyourwork`` submodule (this repo), which contains all the directives.
  This makes it easier to improve/update the workflow, since we can just update ``showyourwork``.
- Changelog: `v0.1.15...v0.1.16 <https://github.com/rodluger/showyourwork/compare/v0.1.15...v0.1.16>`_

0.1.15 (2021-10-21)
+++++++++++++++++++

- **Template repo update:** Added options to the ``Makefile`` to generate a report and a DAG.
  Added a submodule setup check; if the user didn't init the showyourwork submodule, does it
  automatically before building.
- Changelog: `v0.1.14...v0.1.15 <https://github.com/rodluger/showyourwork/compare/v0.1.14...v0.1.15>`_

0.1.14 (2021-10-21)
+++++++++++++++++++

- Remove duplicated Zenodo links from figure captions
- Changelog: `v0.1.13...v0.1.14 <https://github.com/rodluger/showyourwork/compare/v0.1.13...v0.1.14>`_

0.1.13 (2021-10-21)
+++++++++++++++++++

- Fixed API documentation
- Fixed error with `arxiv_tarball_exclude` and arxiv tarball issue (`#21 <https://github.com/rodluger/showyourwork/issues/21>`_)
- Changelog: `v0.1.12...v0.1.13 <https://github.com/rodluger/showyourwork/compare/v0.1.12...v0.1.13>`_

0.1.12 (2021-10-20)
+++++++++++++++++++

- Revert code that prevents the Snakefile from being loaded more than once. Turns out that is
  expected behavior, and is required in order for the module import syntax to work!
- Switched to adding checks within the ``zenodo.py`` script to prevent dependencies from getting
  ingested multiple times.
- Changelog: `v0.1.11...v0.1.12 <https://github.com/rodluger/showyourwork/compare/v0.1.11...v0.1.12>`_

0.1.11 (2021-10-20)
+++++++++++++++++++

- Fix bug preventing figures from being cached properly when one script generates multiple figures
- Fixed issues due to Snakefile being loaded multiple times
- Auto-populate the `projects` page on the docs via a GitHub API search on every release
- Changelog: `v0.1.10...v0.1.11 <https://github.com/rodluger/showyourwork/compare/v0.1.10...v0.1.11>`_

0.1.10 (2021-10-20)
+++++++++++++++++++

- Cleaned up the workflow, separating rules into their own files with better documentation.
- Added a fix for nested figures (figures under subdirectories in the ``src/figures`` folder).
- Fixed issue with multiple Zenodo datasets causing the build to fail.
- Added support for figures in figure* environments.
- Fixed issue with occasional missing </HTML> closing tags in the showyourwork XML tree.
- Added some API documentation; more coming soon.
- Changelog: `v0.1.9...v0.1.10 <https://github.com/rodluger/showyourwork/compare/v0.1.9...v0.1.10>`_

0.1.9 (2021-10-18)
++++++++++++++++++

- **Template repo update:** Added a ``Makefile`` for quick article generation; added docs on how to use it.
- Changelog: `v0.1.8...v0.1.9 <https://github.com/rodluger/showyourwork/compare/v0.1.8...v0.1.9>`_

0.1.8 (2021-10-18)
++++++++++++++++++

- Added "One script, multiple figures" example
- Improved the documentation for script dependencies and datasets
- Fixed a bug when downloading deposits from Zenodo
- Added release testing
- Changelog: `v0.1.7...v0.1.8 <https://github.com/rodluger/showyourwork/compare/v0.1.7...v0.1.8>`_

0.1.7 (2021-10-18)
++++++++++++++++++

- Added explicit support for Zenodo-hosted datasets.
- **Template repo update:** Added the environment variable ``ZENODO_TOKEN`` to ``.github/workflows/showyourwork.yml``.
- Changelog: `v0.1.6...v0.1.7 <https://github.com/rodluger/showyourwork/compare/v0.1.6...v0.1.7>`_

0.1.6 (2021-10-14)
++++++++++++++++++

- Added documentation for the ``expensive-figure`` example.
- Changelog: `v0.1.5...v0.1.6 <https://github.com/rodluger/showyourwork/compare/v0.1.5...v0.1.6>`_

0.1.5 (2021-10-14)
++++++++++++++++++

- Added the ``expensive-figure`` example for computationally expensive figure generation.
- Changelog: `v0.1.4...v0.1.5 <https://github.com/rodluger/showyourwork/compare/v0.1.4...v0.1.5>`_

0.1.4 (2021-10-13)
++++++++++++++++++

- Initial release of the workflow.
