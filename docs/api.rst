Developer API
=============

.. note::

    The API documentation is still under development.
    Please check back soon for more details.


Modules
-------

These modules are located in ``showyourwork/worfklow/rules`` and their contents
are imported into the global namespace when executing the workflow. 
The functions and variables defined in these modules are used in several of 
the rules.

config.py
^^^^^^^^^
.. automodule:: rules.config

errors.py
^^^^^^^^^
.. automodule:: rules.errors

exceptions.py
^^^^^^^^^^^^^
.. automodule:: rules.exceptions

files.py
^^^^^^^^
.. automodule:: rules.files

functions.py
^^^^^^^^^^^^
.. automodule:: rules.functions

paths.py
^^^^^^^^
.. automodule:: rules.paths

zenodo.py
^^^^^^^^^
.. automodule:: rules.zenodo


Rules
-----

These rules are located in files of the same name in the 
``showyourwork/workflow/rules`` directory with the extension ``.smk``.
They are the building blocks of the ``Snakemake`` workflow, and can even be
called externally as ``make <rule_name>``.

.. include:: rules.rst


Scripts
-------

These scripts are located in ``showyourwork/workflow/scripts`` and are
called from some of the rules defined above.

.. include:: scripts.rst