*PAIREF*
========

*Automatic PAIRed REFinement protocol*

*PAIREF* is a tool for macromolecular crystallographers that performs the *PAIRed REFinement protocol* [1]_ automatically to estimate the optimal high-resolution cutoff. It is developed in Python and can be installed as a module into the `Computational Crystallography Toolbox <https://cci.lbl.gov/cctbx_docs/index.html>`_. It provides graphical and command-line interface that executes all the needed calculations. Parameters of refinement can be specified in detail to put all the calculations under full control of the user. Obtained results are presented as plots and tables in HTML log file. *PAIREF* supports `REFMAC5 <http://www.ccp4.ac.uk/html/refmac5.html>`_ (part of the `CCP4 Software Suite <http://www.ccp4.ac.uk/>`_) and `phenix.refine <https://www.phenix-online.org/documentation/reference/refinement.html>`_ (part of the `PHENIX <https://www.phenix-online.org/documentation/reference/refinement.html>`_) for structure model refinement.

.. image:: README_images/example_gui.gif

.. image:: README_images/example_head.gif

.. image:: README_images/example_free_work.gif

.. [1] `P.A. Karplus, K. Diederichs: "Linking crystallographic model and data quality." (2012) Science, 336(6084):1030-3. <https://science.sciencemag.org/content/336/6084/1030>`_

Installation and system requirements
------------------------------------

*PAIREF* depends on the `CCP4 Software Suite <http://www.ccp4.ac.uk/>`_ or `PHENIX <https://www.phenix-online.org/documentation/reference/refinement.html>`_. Both contain the `Computational Crystallography Toolbox <https://cci.lbl.gov/cctbx_docs/index.html>`_ with Python and `pip <https://pip.pypa.io/en/stable/>`_). *PAIREF* works with both Python 2 and Python 3.

*PAIREF* can be easily installed running command :code:`cctbx.python -m pip install pairef --user --no-deps` in terminal (GNU/Linux, macOS) or CCP4Console (Windows). More information are available in `documentation <https://pairef.fjfi.cvut.cz/docs/installation.html>`_. Check also the *PAIREF* homepage at `<https://pairef.fjfi.cvut.cz/>`_ and `PyPI repository <https://pypi.org/project/pairef/>`_.

Example
-------

To run paired refinement of a model (previously refined at 1.81 Å) for a series of cutoffs (1.7, 1.6, and 1.5 Å), execute a following command:

.. code ::

   cctbx.python -m pairef --XYZIN model_1-81A.pdb --HKLIN data_1-5A.mtz --HKLIN_UNMERGED data_1-5A_unmerged.mtz -i 1.81 -r 1.7,1.6,1.5

For detailed information about other program parameters, read the documentation available at `<http://pairef.fjfi.cvut.cz/docs>`_.

Credits and contact
-------------------

Please refer: M. Maly, K. Diederichs, J. Dohnalek, P. Kolenko: `Paired refinement under the control of PAIREF <https://journals.iucr.org/m/issues/2020/04/00/mf5044/index.html>`_ (2020) *IUCrJ* **7**

*PAIREF* is developed by Martin Malý in collaboration of Czech Technical University, Czech Academy of Sciences, and University of Konstanz. In case of any questions or problems, please do not hesitate and write us: `martin.maly@fjfi.cvut.cz <mailto:martin.maly@fjfi.cvut.cz>`_.
