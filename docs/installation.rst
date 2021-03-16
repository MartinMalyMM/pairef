.. _installation-label:

Installation
============

System requirements
-------------------

*PAIREF* is a multiplatform Python module and works on GNU/Linux, Windows and macOS. *PAIREF* depends on the `CCP4 Software Suite <http://www.ccp4.ac.uk>`_ or `PHENIX <https://www.phenix-online.org>`_. Both contain the `Computational Crystallography Toolbox <https://cci.lbl.gov/cctbx_docs/index.html>`_ with Python 2.7.

Installation
------------

Installation can be done easily using *pip* that is usually distributed in CCTBX with Python 2.7. If not, see the :ref:`trouble-label` section below. In GNU/Linux or macOS, just open the `terminal <https://en.wikipedia.org/wiki/Terminal_emulator>`_. In Windows, find the CCP4 console or Phenix Command Prompt in the Start menu and open it (see the screenshots below).

.. image:: _static/ccp4console.gif
.. image:: _static/phenix_command_prompt.gif

Now execute the following command that installs *PAIREF* on your computer:

.. code ::

   cctbx.python -m pip install pairef --no-deps --user

.. note::
   Be sure that you have used the command `cctbx.python` instead of `python`.

If you have administrator permissions, it is possible to install *PAIREF* for all the users. Just skip the option :code:`--user`:

.. code ::

   cctbx.python -m pip install pairef --no-deps

.. note::
   If you have installed *PAIREF* as root or you added the user site-packages directory in your shell *PATH*, you will not need to write the whole :code:`cctbx.python -m pairef ARGUMENTS` but only :code:`pairef ARGUMENTS`.

Packages are also available at the `PyPI <https://pypi.org/project/pairef/>`_ repository.

.. _trouble-label:

Troubleshooting
---------------

Error: No module *pip* available
::::::::::::::::::::::::::::::::

This problem appears often while using Python from CCP4 on Ubuntu 20.04 or newer that does not contain *pip* for Python 2.7 in repositories. To install *pip* in this case, we follow the `instructions in pip documentation <https://pip.pypa.io/en/stable/installing/#installing-with-get-pip-py>`_: download the installation script *get-pip.py* and run it. 

.. code ::

   wget https://bootstrap.pypa.io/2.7/get-pip.py
   ccp4-python get-pip.py

Then it should be possible to install *PAIREF* using the command on the top of this page.

Command prompt freezes on Windows
:::::::::::::::::::::::::::::::::

"*PAIREF or its installation just hangs in the command prompt. Nothing is going on.*"

This can suddenly happen on MS Windows. It is caused by the `Quick Edit Mode of the Command Prompt <https://social.msdn.microsoft.com/Forums/en-US/bf9f97a1-ebbb-4f35-bbb6-6af740a71c76/how-to-disable-command-window-quick-edit-mode-once-and-for-all?forum=vcgeneral>`_. Press Enter or click in the prompt window, the process will continue running.

cctbx.python: Command not found
:::::::::::::::::::::::::::::::

*PAIREF* depends on the `CCP4 Software Suite <http://www.ccp4.ac.uk/>`_ or `PHENIX <https://www.phenix-online.org>`_. The appropriate paths to the executables must be set well in your working shell (typically bash or tcsh).

If you use CCP4, according to the `documentation <http://legacy.ccp4.ac.uk/docs.php#commandline>`_, run :code:`source /path/to/ccp4-<version>/bin/ccp4.setup-sh` (in bash/dash/zsh shells) or :code:`source /path/to/ccp4-<version>/bin/ccp4.setup-csh` (in csh/tcsh shells).

If you use Phenix, according to the `documentation <https://www.phenix-online.org/documentation/install-setup-run.html#setting-up-the-command-line-environment>`_, run :code:`. /usr/local/phenix-<version>/phenix_env.sh` (in bash shell) or :code:`source /usr/local/phenix-<version>/phenix_env` (in csh/tcsh shells).

Then it should be possible to install *PAIREF* using the command on the top of this page.

Reinstallation
--------------

Reinstallation and upgrade to a new version can be done also using *pip*, you can use the following command:

.. code ::

   cctbx.python -m pip install pairef --no-deps --user --upgrade --force-reinstall


If you have administrator permissions, skip the option :code:`--user`.

Uninstallation
--------------

Run command :code:`cctbx.python -m pip uninstall pairef`. If you have installed the package as system administrator, you must run this command as an administrator, too.
