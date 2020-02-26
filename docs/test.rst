.. _test-label:

Tests
=====

*PAIREF* package includes tests to check potential errors in the module code and setting. Framework `pytest <https://docs.pytest.org/en/latest/>`_ is used and has to be installed (this can be done running a following command: :code:`cctbx.python -m pip install pytest --user`. Unfortunately, tests are not supported for Windows now.

To run the test, download the code from `PyPI repository <https://pypi.org/project/pairef/>`_ and open terminal in a folder `test`. There, all the tests can be executed using *e.g.* a following command:

.. code::

   cctbx.python -m pytest -vv -s
