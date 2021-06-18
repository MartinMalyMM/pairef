.. _gui-label:

Graphical interface
===================

*PAIREF* provides a graphical user interface (GUI) to make an execution of calculations easier. It can be launched in the terminal (GNU/Linux, macOS) or CCP4 console (Windows) running command:

.. code ::

   ccp4-python -m pairef --gui

.. note::
   If you got an error message :code:`No module named pairef`, check whether *PAIREF* is properly installed in *CCP4* environment - read the section :ref:`troubleccp4-label`.

.. note::
   If you have installed *PAIREF* as root or you added the user site-packages directory in your shell *PATH*, you should be able to launch GIU simply executing :code:`pairef-gui`.

The graphical user interface depends on `PyQt4 <https://wiki.python.org/moin/PyQt>`_. This framework is available in the *CCP4* Software Suite - not in *PHENIX*. If you wanted to use GUI and *phenix.refine*, it is possible with the following setup: Firstly, set the paths to the executables for both *CCP4* and *PHENIX* - the instructions are written on the following page: :ref:`trouble-label`. In Windows, use the CCP4 console and the paths to *PHENIX* executables can be configured via running the script *phenix_env.bat*, so just execute :code:`"C:\Program Files\Phenix\phenix-installer-1.19.2-4158-intel-windows-x86_64\phenix_env.bat"` in the CCP4 console. With the correct path setting, you should be able to execute :code:`ccp4-python -m pairef --gui` and run *PAIREF* while using *phenix.refine*.
