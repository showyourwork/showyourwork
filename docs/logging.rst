Logging
=======

Log files live in ``.showyourwork/logs``. There is typically a log file
for |showyourwork| messages and errors (``showyourwork.log``), one for Snakemake
messages and errors, including the full verbose Snakemake workflow output 
(``snakemake.log``), and one for LaTeX messages and errors (``tectonic.log``).

If you run into any errors in your build, these log files might help you debug!
Note that you can increase the verbosity of Snakemake by passing the ``--verbose``
command line option to ``showyourwork build``.

Finally, you can also set ``verbose: true`` in ``showyourwork.yml``, which
causes everything that usually gets printed to ``snakemake.log`` to be
printed to the terminal in real time.