==============
Flows Pipeline
==============

.. image:: https://img.shields.io/github/license/rhandberg/flows.svg
    :alt: license
    :target: https://github.com/rhandberg/flows/blob/devel/LICENSE

Installation instructions
=========================
* Go to the directory where you want the Python code to be installed and simply download it or clone it via *git* as::

  >>> git clone https://github.com/rhandberg/flows.git .

* All dependencies can be installed using the following command. It is recommended to do this in a dedicated `virtualenv <https://virtualenv.pypa.io/en/stable/>`_ or similar:

  >>> pip install -r requirements.txt
  
* Last step is to create a config-file. Create a file named "config.ini" and place it in the "flows" directory. Make sure that the file can only be read by you (chmod 0600 config.ini)!
  This file can contain all the settings for running the pipeline. A minimal file for working with the pipeline is 
  
.. code-block:: text
    [api]
    token = <my api token>


How to run tests
================
You can test your installation by going to the root directory where you cloned the repository and run the command::

>>> pytest
