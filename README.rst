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

.. code-block:: bash
    git clone https://github.com/SNflows/flows.git
    cd flows

* Required dependencies can be installed using the following command. It is recommended to do this in a dedicated `virtualenv <https://virtualenv.pypa.io/en/stable/>`_ or similar:

.. code-block:: bash
    python3 -m venv env
    source env/bin/activate
    pip install --upgrade pip
    pip install -r requirements.txt

* In addition, in order to run tests and do development, do

.. code-block:: bash
    pip install -r dev_requirements.txt  # in same virtual environment as above
    _pyver=$(find env/lib/ -type d -mindepth 1 -maxdepth 1 | cut -d '/' -f 3)
    ln -s ../env/lib/${_pyver}/site-packages/tendrils/utils/config.ini flows/config.ini
    wget -O tests/input/2020aatc/SN2020aatc_K_20201213_495s.fits.gz https://anon.erda.au.dk/share_redirect/FJGx69KFvg
    wget -O tests/input/2020lao/59000.96584_h_e_20200531_33_1_1_1_2020lao_LT_gp.fits.gz https://anon.erda.au.dk/share_redirect/E98lmqOVWf
    wget -O tests/input/2020lao/subtracted/59000.96584_h_e_20200531_33_1_1_1_2020lao_LT_gpdiff.fits.gz https://anon.erda.au.dk/share_redirect/bIxyzrRXbg
    wget -O tests/input/2021wyw/ADP.2021-10-15T11_40_06.553.fits.gz https://anon.erda.au.dk/share_redirect/Gr8p2K7ph5


TODO: Reformulate following bullet point and check/decide on which config paths are/should be available...:

* **Changed with ``tendrils`` API.** If using ``tendrils``, follow the steps below, but then let ``tendrils`` know of the config file location. Alternatively, individual config file elements can be set programatically using `tendrils` and will be saved to a config file automatically. Last step is to create a config-file. Create a file named "config.ini" and place it in the "flows" directory. Make sure that the file can only be read by you (``chmod 0600 config.ini``)!
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

.. code-block:: bash

    pytest

Full configuration file
=======================
Text coming soon...

.. code-block:: ini

    ##########################################################
    ### --Tendrils configurations used at FLOWS run-time-- ###
    ### Configurations without a leading '#' are required  ###
    ### and must be specified by the user.                 ###
    ### Configurations with a leading '#' are optional,    ###
    ### and their assigned values are documentation of     ###
    ### their default values in the FLOWS pipeline.        ###
    ### Default values of optional configurations with a   ###
    ### leading '$' signify environment variables resolved ###
    ### at run-time; their fallbacks are given following   ###
    ### a '/'.                                             ###
    ##########################################################

    [api]
    # photometry_cache = None
    # pipeline = False
    token = None

    # casjobs:
    #   wsid and password required for run_catalogs.py,
    #   user registration at
    #   https://galex.stsci.edu/casjobs/CreateAccount.aspx
    #   wsid can be found at
    #   https://galex.stsci.edu/casjobs/changedetails.aspx
    #   after login
    [casjobs]
    # wsid = $CASJOBS_WSID/None
    # password = $CASJOBS_PASSWORD/None

    # database:
    #   username and password required for run_catalogs.py,
    #   the user is a registered user in the flows database
    #   with access to the 'adastra' schema
    [database]
    # username = $AUDBUsername/None
    # password = $AUDBPassword/None

    [photometry]
    archive_local = None
    # output = .

    # TNS:
    #   api_key required for run_querytns.py,
    #   user registration at
    #   https://www.wis-tns.org/user
    #   api_key is that of a TNS bot; ask a flows group
    #   member for one
    #   if user_id and user_name are not given, fallback
    #   to a TNS bot's bot_id and bot_name, which must
    #   match with api_key
    [TNS]
    # api_key = None
    # bot_id = 191396
    # bot_name = AUFLOWS_BOT2
    # user_id = None
    # user_name = None

    [URL]
    # base_url = https://flows.phys.au.dk/api/
    # catalogs_url = reference_stars.php
    # catalogs_missing_url = catalog_missing.php
    # cleanup_photometry_status_url = cleanup_photometry_status.php
    # datafiles_url = datafiles.php
    # filters_url = filters.php
    # lightcurves_url = lightcurve.php
    # photometry_upload_url = upload_photometry.php
    # photometry_url = download_photometry.php
    # set_photometry_status_url = set_photometry_status.php
    # sites_url = sites.php
    # targets_post_url = targets_add.php
    # targets_url = targets.php
    # verify_ssl = True

    [ztf]
    # output_photometry = .

Making a release
================

 - Bump sem-version when Devel is ready to merge in file = VERSION (v1.0.0). Checkout devel. Edit Version. Push devel.
 - Merge Devel into Master (Create PR from Devel -> Master), wait until tests are passing. Create issues if not. Then Merge.
 - Create tag on Master corresponding to right semversion. This means, checkout master. Pull master locally. Create tag using git tag called "v1.0.0" or whatever the sem-version. Push local tag to GitHub.
 - Merge Master into devel to propagate tag (Create PR on GitHub).
 - Create release on GH releases tab if all tests passing.
