.. _installation-label:

Installation
============

System requirements
-------------------

*PAIREF* is a multiplatformal Python module and works on GNU/Linux, Windows and macOS. *PAIREF* depends on the `CCP4 Software Suite <http://www.ccp4.ac.uk/>`_ or `PHENIX <https://www.phenix-online.org/documentation/reference/refinement.html>`_. Both contain the `Computational Crystallography Toolbox <https://cci.lbl.gov/cctbx_docs/index.html>`_ with Python 2.7 and `pip <https://pip.pypa.io/en/stable/>`_).

Installation
------------

Installation can be done easily using *pip*. Open terminal (GNU/Linux, macOS) or CCP4 console (Windows) and execute as the system administrator a following command:

.. code ::

   cctbx.python -m pip install pairef --no-deps

.. note::
   Be sure that you have used the command `cctbx.python` instead of `python`.

If you do not have the administrator permissions, it is still possible to install the module in the `"user site" location <https://www.python.org/dev/peps/pep-0370/>`_. Just add an option :code:`--user`:

.. code ::

   cctbx.python -m pip install pairef --user --no-deps

.. note::
   If you have installed *PAIREF* as root or you added the user site-packages directory in your shell *PATH*, you will not need to write the whole :code:`cctbx.python -m pairef ARGUMENTS` but only :code:`pairef ARGUMENTS`.

Packages are also available at `PyPI <https://pypi.org/project/pairef/>`_.

Reinstallation
--------------

Reinstallation and upgrade to a new version can be done also using *pip*. If you have the administrator permissions, you can use the following command:

.. code ::

   cctbx.python -m pip install pairef --upgrade --force-reinstall --no-deps

If not, add the option :code:`--user`.

Uninstallation
--------------

Run command :code:`cctbx.python -m pip uninstall pairef`. If you have installed the package as system administrator, you must run this command as administrator, too.
