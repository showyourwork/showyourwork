Zenodo integration
==================

.. note:: Coming soon!

The ``showyourwork`` workflow automatically caches the results of intermediate
steps in your pipeline on Zenodo (see :doc:`zenodo` for details). In order to
do this, it needs access to a Zenodo API token, which you can 
`generate here <https://zenodo.org/account/settings/applications/tokens/new>`_
(you'll need to set up a Zenodo account first if you don't already have one).
Name the token something informative (like ``showyourwork api token``) and make
sure to give it ``deposit:actions`` and ``deposit:write`` permissions. Copy the
token and store it somewhere secure. 

.. warning::

    Never commit your Zenodo API token (or any API token) directly to your
    repository!
