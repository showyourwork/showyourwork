Zenodo integration
==================

|showyourwork| integrates with the `Zenodo <https://zenodo.org>`_ and
`Zenodo Sandbox <https://sandbox.zenodo.org>`_ services to
easily allow users to host datasets and simulation results used in the
workflow. There are two main ways in which this integration occurs: with
static datasets and with dynamic datasets.


Static datasets
---------------

If your workflow depends on data that cannot be programmatically generated,
such as data collected from a telescope, a mass spectrometer, a DNA
sequencer, etc., that data should be made available to anyone trying to
reproduce your results. Instead of committing the dataset directly to the
repository, we recommend you archive it on an online open-access 
file-hosting service. You can use whichever service you want, but if you use
Zenodo, |showyourwork| does all the communicating back-and-forth for you. All
you need to do is specify the ID of the (public) archive and some information about the
files your workflow needs in the ``showyourwork.yml`` config file.
For detailed information about how to work with static Zenodo datasets, please
check out :ref:`config.datasets`.


Dynamic datasets
----------------

The |showyourwork| workflow has the ability to cache the results of intermediate
steps in your pipeline on either Zenodo or Zenodo Sandbox. 
This is useful for workflows that entail
running lengthy computations, simulations, etc., that third-party users may
not want to run on their own. It's also useful for builds on GitHub Actions,
which has limited compute resources and a timeout of a few hours. 

The way |showyourwork| deals with these cases is to cache these lengthy
computations on Zenodo or Zenodo Sandbox alongside a record of all the
inputs that went into generating the cached output. If, on subsequent runs
of the workflow, the inputs remain unchanged, |showyourwork| will simply
download the cached results from Zenodo or Zenodo Sandbox, *maintaining the guarantee that
the output you get follows determinstically from the given inputs*.

Before we get into how this works and how to take advantage of it, let's
discuss how to set up the integration by generating Zenodo and Zenodo Sandbox API tokens.


Setting up Zenodo integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, caching takes place on Zenodo Sandbox, which is equal to Zenodo in all
respects except that records on that service are *temporary*. This makes it
perfect for caching intermediate results during the development cycle, where things
are likely to change a lot and you may not want to assign a static, permanent
DOI with any particular dataset.

In order for |showyourwork| to communicate with Zenodo Sandbox, it needs access to an API token.
You can generate it by clicking
`here <https://sandbox.zenodo.org/account/settings/applications/tokens/new>`__.
If you don't yet have a Zenodo Sandbox account, you'll be asked to create one.

Name the token something informative and make
sure to give it ``deposit:actions`` and ``deposit:write`` permissions. Copy the
token and store it somewhere secure. 

.. warning::

    Never commit your Zenodo API token (or any API token) directly to your
    repository!

Then, on your local computer, create an environment variable called
``$SANDBOX_TOKEN`` with value equal to the API token you just generated. 
You can either do this manually in the terminal, e.g.,

.. code-block:: bash

    export SANDBOX_TOKEN=YYYYYY

or by adding that line to your shell config file (``.bashrc``, ``.zshrc``, etc.)
and re-starting your session.
In order for |showyourwork| to have access to Zenodo Sandbox when running on GitHub
Actions, you'll also have to provide this value as a secret (name it ``SANDBOX_TOKEN``); 
read more about those 
`here <https://docs.github.com/en/actions/security-guides/encrypted-secrets>`_.

If you've done all that, the next time you create a new article repository
using :ref:`syw_setup`, pass the ``--cache`` option and |showyourwork| will automatically create a Zenodo
Sandbox draft deposit which it will use to cache your intermediate results. Note that
you can also manually create a draft deposit by running ``showyourwork cache create``
(see :ref:`syw_cache`, and below, for details).


Intermediate results
^^^^^^^^^^^^^^^^^^^^

Earlier, we mentioned that the Zenodo integration allows users to cache intermediate
results in their workflow.
But hang on--what's an "intermediate result"? The standard procedure for generating
figures using |showyourwork| is to define a figure script that generates the
figure output (PDF, PNG, etc.) and to specify that script using the ``\script``
command in your TeX file (see :ref:`latex_script`). There's no intermediate step
there--we simply go from figure script to figure output (and then to the article
PDF).

Suppose, however, that our figure script involves some lengthy computation,
integration, or simulation. Every time we change anything in that script,
|showyourwork| will attempt to re-run the entire computation when asked to
build your article PDF. This is good for reproducibility---i.e., to always
ensure the output is up to date with the inputs---but it is extremely
wasteful in cases where, e.g., we wish to tweak some aspect of the plot,
like the color of a line in a graph. 

To avoid this, we can split our script into two: one that runs the simulation
and saves the results, and one that loads the results and plots them. Let's
discuss how to do this by considering an example.


An expensive workflow
^^^^^^^^^^^^^^^^^^^^^

Consider a workflow containing the following figure script and TeX file:


.. code-block:: python
    :caption: **File:** ``src/scripts/figure.py``

    import simulation
    import matplotlib.pyplot as plt
    import paths

    # Run the simulation for some inputs
    simulation.run(x=10, y=25)
    data = simulation.get_results()

    # Plot the results
    fig, ax = plt.subplots(1)
    ax.plot(data, color="k")
    fig.savefig(paths.figures / "figure.pdf")


.. code-block:: TeX
    :caption: **File:** ``src/tex/ms.tex``
    
    ...

    \begin{figure}[ht!]
        \script{figure.py}
        \begin{centering}
            \includegraphics{figures/figure.pdf}
            \caption{Simulation results.}
            \label{fig:figure}
        \end{centering}
    \end{figure}

    ...


where ``simulation`` is some custom package we're using to run
an expensive simulation. As we mentioned above, changing anything in the
file ``src/scripts/figure.py``, including something as trivial as the plot
line color, will result in a re-run of the entire simulation the next time
we build the article.


The streamlined version
^^^^^^^^^^^^^^^^^^^^^^^

We would like to streamline our workflow by decoupling the plotting step
from the simulation step. We can do this by introducing a new script, which
we'll call ``simulation.py``, that runs and saves the result of the simulation.
Then, in ``figure.py``, we load the result and plot our figure:

.. code-block:: python
    :caption: **File:** ``src/scripts/simulation.py``

    import simulation
    import numpy as np
    import paths

    # Run the simulation for some inputs
    simulation.run(x=10, y=25)
    data = simulation.get_results()

    # Save the results
    np.savetxt(paths.data / "simulation.dat", data)


.. code-block:: python
    :caption: **File:** ``src/scripts/figure.py``

    import numpy as np
    import matplotlib.pyplot as plt
    import paths

    # Load the data
    data = np.loadtxt(paths.data / "simulation.dat")

    # Plot the results
    fig, ax = plt.subplots(1)
    ax.plot(data, color="k")
    fig.savefig(paths.figures / "figure.pdf")


Our workflow is now separable: changes to ``figure.py`` will not result
in the re-execution of the simulation, as they are merely plotting changes.
The simulation will only be re-executed if we change something in ``simulation.py``,
like the input arguments to our ``simulation.run()`` function.

In order to get this all to work, we need to tell |showyourwork| two things:
(1) the script ``figure.py`` has a dependency called ``simulation.dat`` and
(2) the dependency ``simulation.dat`` can be generated by running the script
``simulation.py``. We accomplish this by (1) editing the config file:

.. code-block:: yaml
    :caption: **File:** ``showyourwork.yml``

    dependencies:
        src/scripts/figure.py:
            - src/data/simulation.dat

(see :ref:`config.dependencies` for details) and (2) adding a custom
rule to our Snakefile:

.. code-block:: python
    :caption: **File:** ``Snakefile``

    rule simulation:
        output:
            "src/data/simulation.dat"
        script:
            "src/scripts/simulation.py"

(see :doc:`snakefile` for details).


Caching the intermediate result
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The workflow above is now separable, but we're still not caching anything.
If we commit and push it to GitHub, the runner will still have to execute
``simulation.py`` in order to generate ``simulation.dat``; the same goes for
third-party users who have cloned your repository. Adding caching functionality
can be done by adding a single line to the ``Snakefile``:

.. code-block:: python
    :caption: **File:** ``Snakefile``

    rule simulation:
        output:
            "src/data/simulation.dat"
        cache:
            True
        script:
            "src/scripts/simulation.py"


which tells |showyourwork| to cache the output of that rule (``simulation.dat``).
Normally, if we were just running this in a regular Snakemake pipeline, this
would result in the data file getting cached in some local hidden folder. The
next time you run your workflow, Snakemake will check to see if any of the inputs
to the ``simulation`` rule changed and, if not, it will restore ``simulation.dat``
from the cache (if it's needed).

|showyourwork| builds on this functionality by also caching the file ``simulation.dat``
on Zenodo Sandbox, allowing the results to be restored on *any* computer running your
workflow (as long as they have the correct ``SANDBOX_TOKEN``; but more on this
in a moment). This means that, provided you have run your workflow locally first, 
the runner on GitHub Actions will never have to execute ``simulation.py``, as
it can just download the result from Zenodo Sandbox. Recall that this procedure still
guarantees that you'll get the *same result* as if you had run your entire
simulation (provided your workflow is deterministic), since a cache is only
restored if *none* of the upstream inputs to a rule have changed.

The cached files (and the hashes of the rule inputs)
are stored in a Zenodo Sandbox deposit draft with concept ID specified
in your ``zenodo.yml`` config file. If you navigate to Zenodo Sandbox in your
browser and log in, you should see a draft with a title like 
``Data for user/repo [main]``, where ``user/repo`` is your repository slug
and ``main`` is the current branch. At any given time, you can only have
one draft per deposit, so if you change any of the inputs to your rule (e.g., if
you change the file ``simulation.py``), the draft will get overwritten with
a new version of the cache. Note, also, that drafts are *private*: only
users with access to your account can see their files.


.. note::

    If you switch branches, or if you set up a repository without caching
    functionality and would like to add it, you can create a new Zenodo Sandbox deposit 
    for the current branch by running

    .. code-block:: bash

        showyourwork cache create

    .. raw:: html

        <br/>


Freezing the cache
^^^^^^^^^^^^^^^^^^

Only one draft can exist at a given time on Zenodo Sandbox. This is inconvenient
if, during development, you find the need to frequently switch back and forth
between different sets of inputs (such as a few different randomizer seeds for
a simulation). |showyourwork| can help with this! At any time, you can "freeze"
the latest draft on Zenodo Sandbox by running

.. code-block:: bash

    showyourwork cache freeze

which publishes your draft into an actual deposit on Zenodo Sandbox. The next
time you edit your workflow, a new draft will be created. Future builds can
then draw from either the published version *or* the draft. This process can
be repeated as often as you need, so that at any given time the workflow is
able to restore intermediate results from any number of deposits on Zenodo
Sandbox. 

One useful side effect of freezing a deposit is that an API token is not
required to access it. This is especially convenient for working with
outside collaborators, who will not in general have access to your
``SANDBOX_TOKEN``. If the cache has been frozen, they will be able to 
restore output from the cache both locally and when building pull requests
on GitHub Actions.

.. note::

    It is important to keep in mind that Zenodo Sandbox deposits cannot be deleted.
    While the service is temporary---and can be wiped at any time to free up
    space on the Zenodo servers---users should take care to not abuse it!


Publishing the cache
^^^^^^^^^^^^^^^^^^^^

When you're ready to publish or distribute your article to the outside world
--and you're confident the inputs to your cached rules won't change again--
you should consider migrating your cached outputs from Zenodo Sandbox to Zenodo. 
You can easily do this by running

.. code-block:: bash

    showyourwork cache publish

in the top level of your repo. This will transfer the latest draft of your
Zenodo Sandbox deposit over to Zenodo and publish it, giving it a 
permanent DOI (digital object identifier) and making it visible to unauthenticated users.
Once you do this, anyone can take advantage of the caching functionality, for
all of eternity (well, for as long as Zenodo continues to exist).

Note that in order for this step to work, you'll need proper authentication
in the form of an API token stored as the environment variable ``ZENODO_TOKEN``.
This is different from the Zenodo Sandbox token we generated above. If you
don't have one, you may generate one by clicking
`here <https://sandbox.zenodo.org/account/settings/applications/tokens/new>`__.
If you don't yet have a Zenodo account, you'll be asked to create one.

.. note::

    Once you publish your deposit, further changes to a cached rule's inputs
    will result in a new draft being created on Zenodo Sandbox (much like what
    happens when you ``freeze`` your cache). Future runs of your workflow
    will be able to restore the cache from any of the frozen versions on Zenodo Sandbox, 
    published versions on Zenodo, or from the latest draft on Zenodo Sandbox.

.. warning::

    Published Zenodo deposits are permanent! There is no way to delete a Zenodo deposit once
    it's published, as it now has a perennial DOI associated with it. Therefore,
    it is important that users be responsible in their use of this service!
    If you simply want a few different versions of your cache while you develop
    your project, please use the ``freeze`` feature instead (see above).


Deleting the cache
^^^^^^^^^^^^^^^^^^

You can delete the latest Zenodo Sandbox cache draft for the current branch by running

.. code-block:: bash

    showyourwork cache delete

Note that, as we mentioned above, you can't delete deposits once they have
been published!