Introduction
============

The ``showyourwork`` workflow is intended to help authors publish the code that generated
the figures and results in a scientific article. It ensures that the compiled article
PDF is always in sync with all of the code used to generate it. It does this
automatically—and seamlessly—with the help of the Snakemake workflow management system,
the tectonic typesetting engine, and Github Actions CI. The basic philosophy behind
showyourwork is this: scientific papers should exist as GitHub repositories comprised of
LaTeX files, figure scripts, rules to access datasets, a platform/environment specification,
and nothing else. Anyone should be able to re-generate the article PDF from scratch at
the click of a button.
