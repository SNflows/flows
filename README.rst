==============
Flows Pipeline
==============
.. image:: https://zenodo.org/badge/241705955.svg
   :target: https://zenodo.org/badge/latestdoi/241705955
.. image:: https://github.com/SNflows/flows/actions/workflows/tests.yml/badge.svg?branch=devel
    :target: https://github.com/SNflows/flows/actions/workflows/tests.yml
.. image:: https://codecov.io/github/SNflows/flows/branch/devel/graph/badge.svg?token=H8CQGPG0U6
    :target: https://codecov.io/github/SNflows/flows
.. image:: https://hitsofcode.com/github/SNflows/flows?branch=devel
    :alt: Hits-of-Code
    :target: https://hitsofcode.com/view/github/SNflows/flows?branch=devel
.. image:: https://img.shields.io/github/license/SNflows/flows.svg
    :alt: license
    :target: https://github.com/SNflows/flows/blob/devel/LICENSE

Installation instructions
=========================
* Go to the directory where you want the Python code to be installed and simply download it or clone it via *git* as::

  >>> git clone https://github.com/SNflows/flows.git .

* Required dependencies can be installed using the following command. It is recommended to do this in a dedicated `virtualenv <https://virtualenv.pypa.io/en/stable/>`_ or similar:

  >>> pip install -r requirements.txt
  >>> pip install -r requirements_dev.txt  # for tests/development

* Last step is to create a config-file. Create a file named "config.ini" and place it in the "flows" directory. Make sure that the file can only be read by you (``chmod 0600 config.ini``)!
  This file can contain all the settings for running the pipeline. A minimal file for working with the pipeline is

  .. code-block:: ini

      [api]
      token = <my api token>

      [TNS]
      api_key = <AUFLOWS_BOT API key>

  Where your API token can be found on the Flows webpage.


How to run tests
================
You can test your installation by going to the root directory where you cloned the repository and run the command::

>>> pytest

Full configuration file
=======================
Text coming soon...

.. code-block:: ini

    [api]
    token =
    photometry_cache =

    [photometry]
    archive_local =
    output =

    [casjobs]
    wsid =
    password =

    [TNS]
    api_key =
    bot_id =
    bot_name =
    user_id =
    user_name =

    [ztf]
    output_photometry =

    [database]
    username =
    password =

# Making a release
Bump sem-version when Devel is ready to merge.
Merge Devel into Master, and ensure tests are passing.
Create tag on Master corresponding to right semversion.
Merge Master into devel to propagate tag.
Create release on GH releases tab if all tests passing.
