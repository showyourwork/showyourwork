.. :changelog:

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
