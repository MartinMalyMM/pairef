.. _gui-label:

Graphical interface
===================

*PAIREF* provides a graphical interface to make an execution of calculations easier. It can be launched running command:

.. code ::

   cctbx.python -m pairef --gui

.. note::
   If you have installed *PAIREF* as root or you added the user site-packages directory in your shell *PATH*, simply execute only :code:`pairef-gui`.

The user interface depens on `PyQt4 <https://wiki.python.org/moin/PyQt>`_. This framework is available in the *CCP4* Software Suite. If you use this suite and you got error messages while launching the graphical interface, try to execute :code:`ccp4-python -m pairef --gui`. If you have installed only *PHENIX*, you will have to install `PyQt4 <https://wiki.python.org/moin/PyQt>`_ manually.
