Zenodo integration
==================

|showyourwork| integrates with the `Zenodo <https://zenodo.org>`_ service to
easily allow users to host datasets and simulation results used in the
workflow. There are two main ways in which this integration occurs: with
static datasets and with dynamic datasets.

Static datasets
---------------

.. note:: Page under construction. More information coming soon!


Dynamic datasets
----------------

The |showyourwork| workflow has the ability to cache the results of intermediate
steps in your pipeline on Zenodo. In order to do this, it needs access to a Zenodo API token, which you can 
`generate here <https://zenodo.org/account/settings/applications/tokens/new>`_
(you'll need to set up a Zenodo account first if you don't already have one).
Name the token something informative (like ``showyourwork api token``) and make
sure to give it ``deposit:actions`` and ``deposit:write`` permissions. Copy the
token and store it somewhere secure. 

.. warning::

    Never commit your Zenodo API token (or any API token) directly to your
    repository!

.. note:: Page under construction. More information coming soon!