.. _using-label:

Using *PAIREF*
==============

Getting started
---------------

Open terminal (GNU/Linux, macOS) or CCP4 console (Windows) and `go to <https://en.wikipedia.org/wiki/Cd_(command)>`_ the folder where your structure model and diffraction data are saved. For example, let's assume that your structure model *nuclease_model.pdb* has been refined using data file *data_2A.mtz* at 2 Å resolution and all the files are located in a folder :code:`/home/test/project/nuclease/dataset7/files/`. To change the current working directory, run a command:

.. code ::

   cd /home/test/project/nuclease/dataset7/files/

However, you have also prepared full-resolution merged diffraction data file *data_full_resolution.mtz* and unmerged diffraction data file *XDS_ASCII_full.HKL* (both at 1.5 Å).

.. note::
   The files *data_full_resolution.mtz* and *data_2A.mtz* should contain consistent free reflection sets.
   
   The usage of unmerged data is not obligatory, however, it is strongly recommended as they are required for the *CC** calculation. Various file formats are supported (*.HKL* from *XDS*, *.mtz*, *.sca*).

Now, you would like to perform the paired refinement protocol and use step-by-step following high resolution limits: 1.9 Å, 1.8 Å, 1.7 Å, 1.6 Å, and 1.5 Å. To execute these calculations, run a command:

.. code ::

   cctbx.python -m pairef --HKLIN data_full_resolution.mtz --XYZIN nuclease_model.pdb -u XDS_ASCII_full.HKL -i 2 -r 1.9,1.8,1.7,1.6,1.5 -p nuclease

Then a new folder *pairef_nuclease* is created in the folder where the command has been executed and all the log files, new structure models, *etc.*, will be saved there. Open a file *PAIREF_nuclease.html* in a web browser to see the current progress, results, plots, and statistics.

*PAIREF* will refine the input structure model in *REFMAC5* (default 10 cycles) against data up to 1.9 Å. Then it will calculate statistics relating to the refined model (also using *REFMAC5*) and plot graphs. After that, the refined model will be further refined against data up to 1.8 Å and its relating statistics will be computed. This will be also performed using the remaining high resolution diffraction limits 1.7 Å, 1.6 Å, and 1.5 Å. In the end, merging statisting will be calculated.

Detailed specification of refinement parameters
-----------------------------------------------

To obtain meaningful results, the refinement setting during the paired refinement protocol should be very similar to the setting that has been used in previous refinement steps. `PAIREF` provides many option to run all the calculations under full control of the user.

Using external *CIF file* (*LIBIN*)
+++++++++++++++++++++++++++++++++++

If a *CIF file* with external restrains has been used in previous refinement steps, it should be specified to be used also in the paired refinement. This can be specified with an option :code:`--LIBIN some_restrains.cif` (assuming that the file is saved in the folder where `PAIREF` is executed.

Weighting term
++++++++++++++

The weight of the X-ray term can be specified using option :code:`-w value`, *e.g.* :code:`-w 0.5`.

Number of refinement cycles
+++++++++++++++++++++++++++

The number of refinement cycles that is be performed in every resolution step can be controlled using an option :code:`--ncyc value`, *e.g.* :code:`--ncyc 20`. The default setting is 10 cycles.

TLS refinement
++++++++++++++

It is possible to perform a TLS refinement in `PAIREF` before a restrained refinement. An input TLS file is speciffied by an option :code:`--TLSIN`. A new TLS output file generated during a refinement is then used in the next refinement run (using data up to a higher resolution) as the TLS input file. To avoid this default behaviour and use the same  in all the refinement runs, use an option :code:`--TLSIN-keep`. A number of TLS refinement cycles can be set (*e.g.* for 5 cycles: :code:`--TLS-ncyc 5`), 10 cycles are performed by default.

*REFMAC5* parameters  (*Com file*)
++++++++++++++++++++++++++++++++++

A *Com file* (*command file*) is a text file describing refinement parameters for *REFMAC5*. Important parameters are *e.g.* a `weight matrix <http://www.ccp4.ac.uk/html/refmac5/keywords/keywords_5_5.html#Weight>`_ or a number of refinement cycles `ncyc <http://www.ccp4.ac.uk/html/refmac5/keywords/xray-principal.html#ncyc>`_.

To obtain a *command file* of particular refinement job in CCP4, select the last refinement job and press *ReRun Job..* . In a newly opened dialog, do not press *Run Now* but select *Run & View Com File* (details are described in the `CCP4 documentation <http://www.ccp4.ac.uk/dist/checkout/ccp4i/help/general/runjob.html>`_). Then a new dialog is opened -- the text at the bottom is the content of the *command file*. Select it, press *Ctrl+C* to copy it, paste it in a text editor and save it as *e. g.* *setting.com* in the folder where the diffraction data and model are placed.

The *command file* is specified with an option :code:`-c setting.com` where *setting.com* is the file containing parameters for *REFMAC5*.

.. note::
   Even thought the weight of the X-ray term or the number of refinement cycles are set in the *Com file* (*REFMAC5* keywords :code:`WEIGht MATRix` and :code:`NCYC`), the values specified by the options :code:`--weight` and :code:`--ncyc` have the higher priority.

Modification of input structure model
-------------------------------------

The input structure model can be modified and refined at the starting resolution before the paired refinement. These options should be used if the structure has been refined in another software (*e.g.* *phenix.refine*, or another version of *REFMAC5* than it is currently used, or the bias of previous free reflection selection is present. The number of refinement cycles at the starting resolution is be controlled by the option :code:`--prerefinement-ncyc` (20 cycles by default).

Possible modifications of the structure model:

* reset ADPs their mean value: :code:`--prerefinement-reset-bfactor`,
* add a value to the ADPs: :code:`--prerefinement-add-to-bfactor ADD_TO_BFACTOR`,
* set ADPs to a value: :code:`--prerefinement-set-bfactor`,
* perturb the atomic coordinates by an average of a value (0.25 Å by default): :code:`--prerefinement-shake-sites [SHAKE_SITES]`,
* no modification :code:`--prerefinement-no-modification`.

Summary of program options
--------------------------

.. code ::

   $ cctbx.python -m pairef -h
   usage: cctbx.python -m pairef --XYZIN XYZIN --HKLIN HKLIN [-u HKLIN_UNMERGED]
                                 [--LIBIN LIBIN] [-c COMIN] [-p PROJECT]
                                 [-r RES_SHELLS] [-n N_SHELLS] [-s STEP]
                                 [-i RES_INIT] [-f FLAG] [-w WEIGHT]
                                 [--ncyc NCYC] [--constant-grid] [--complete]
                                 [-h] [--prerefinement-ncyc PREREFINEMENT_NCYC]
                                 [--prerefinement-reset-bfactor]
                                 [--prerefinement-add-to-bfactor ADD_TO_BFACTOR]
                                 [--prerefinement-set-bfactor SET_BFACTOR]
                                 [--prerefinement-shake-sites [SHAKE_SITES]]
                                 [--prerefinement-no-modification]
   
   Automatic PAIRed REFinement protocol
   
   optional arguments specifying input files:
     --XYZIN XYZIN, --xyzin XYZIN
                           PDB file with current structure model
     --HKLIN HKLIN, --hklin HKLIN
                           MTZ file with processed diffraction data
     -u HKLIN_UNMERGED, --unmerged HKLIN_UNMERGED
                           unmerged processed diffraction data file
                           (XDS_ASCII.HKL)
     --LIBIN LIBIN, --libin LIBIN
                           CIF file geometric restraints
     --TLSIN TLSIN, --tlsin TLSIN
                           input TLS file
     -c COMIN, --comfile COMIN
                           configuration Com file for REFMAC5
   
   other optional arguments:
     -p PROJECT, --project PROJECT
                           project name
     -r RES_SHELLS         explicit definition of high resolution shells - values
                           must be divided using commas without any spaces and
                           written in decreasing order, e.g. 2.1,2.0,1.9
     -n N_SHELLS           number of high resolution shells to be added step by
                           step. Using this argument, setting of argument -s is
                           required.
     -s STEP, --step STEP  width of the added high resolution shells (in
                           angstrom). Using this argument, setting of argument -n
                           is required.
     -i RES_INIT           initial high resolution diffraction limit (in
                           angstrom) - if it is not necessary, do not use this
                           option, the script should find resolution
                           automatically in PDB file
     -f FLAG, --flag FLAG  definition which FreeRflag set will be excluded during
                           refinement (set 0 default)
     -w WEIGHT, --weight WEIGHT
                           manual definition of weighting term for REFMAC5
     --ncyc NCYC           number of refinement cycles that will be performed in
                           every resolution step
     --constant-grid       keep the same FFT grid through the whole paired
                           refinement.
     --complete            perform complete cross-validation (use all available
                           free reflection sets)
     --TLS-ncyc TLS_NCYC   number of cycles of TLS refinement (10 cycles by
                           default)
     --TLSIN-keep          keep using the same TLS input file in all the
                           refinement runs
     -h, --help            show this help message and exit
   
   optional arguments specifying structure model modification:
     --prerefinement-ncyc PREREFINEMENT_NCYC
                           number of refinement cycles to be performed as pre-
                           refinement of the input structure model before paired
                           refinement (the initial high resolution limit is
                           used). Pre-refinement is performed by default in case
                           of the complete cross-validation protocol. Other
                           related options are --prerefinement-reset-bfactor,
                           --prerefinement-add-to-bfactor, --prerefinement-set-
                           bfactor, --prerefinement-shake-sites, and
                           --prerefinement-no-modification. These options can be
                           useful when the structure has been refined in another
                           version of REFMAC5 than it is currently used or when
                           you want to reset the impact of used free reflections.
     --prerefinement-reset-bfactor
                           reset atomic B-factors of the input structure model to
                           the mean value. This is done by default in the case of
                           the completecross-validation protocol.
     --prerefinement-add-to-bfactor ADD_TO_BFACTOR
                           add the given value to B-factors of the input
                           structure model
     --prerefinement-set-bfactor SET_BFACTOR
                           set atomic B-factors of the input structure model to
                           the given value.
     --prerefinement-shake-sites [SHAKE_SITES]
                           randomize coordinates of the input structure model
                           with the given mean error value. This is done by
                           default in the case of the complete cross-validation
                           protocol - mean error 0.25.
     --prerefinement-no-modification
                           do not modify the input structure model before the
                           complete cross-validation protocol
   
   Dependencies: CCP4 package containing CCTBX with Python 2.7

Example: 

 * Structure model: *nuclease_model.pdb* (has been previously refined at 2.0 Å),
 * Diffraction data -- merged: *data_full_resolution.mtz* (data up to 1.5 Å),
 * Diffraction data -- unmerged: *XDS_ASCII_full.HKL* (data up to 1.5 Å),
 * High resolution limits: 1.9 Å, 1.8 Å, 1.7 Å, 1.6 Å, and 1.5 Å;
 * External restrains: *ligands.cif*,
 * Command file including external harmonics (*REFMAC5* parameters): *setting.com*.
 * X-ray weight: 0.04
 * Number of refinement cycles to be performed during every resolution step: 15
 * Project name: *nuclease*,

.. code ::

   cctbx.python -m pairef --HKLIN data_full_resolution.mtz --XYZIN nuclease_model.pdb -u XDS_ASCII_full.HKL --LIBIN ligands.cif -c setting.com -i 2 -r 1.9,1.8,1.7,1.6,1.5 -w 0.04 --ncyc 15 -p nuclease

The command file *setting.com* is the following text file:

.. code ::

   make -
       check NONE
   refi -
       resi MLKF -
       meth CGMAT -
       bref MIXED
   scal -
       type SIMP -
       LSSC -
       ANISO -
       EXPE
   solvent YES
   external harmonic residues from 3 B to 4 B sigma 0.03
   exte dist first chain A resi 777 atom CD second chain A resi 777 atom OE1 value 1.20 sigma 0.01
   PNAME nuclease
   DNAME nuclease_42

Advanced options
----------------

Complete cross-validation
+++++++++++++++++++++++++

To run the paired refinement protocol for each individual free reflections set (*e.i.* to perform the complete cross-validation), use an option :code:`--complete`. The input structure model is modified to remove the bias of previous free reflection selection. The default setting is: 

* the atomic coordinates are perturbed by an average of 0.25 Å,
* ADPs are set to their average value. 

The modified model is then refined at the starting resolution, the number of refinement cycles is controlled by an option :code:`--prerefinement-ncyc` (20 cycles by default). To disable the automatic modification, use an option :code:`--prerefinement-no-modification`. For further information about the input model modification, see the section `Modification of input structure model`_.

Constant FFT-grid
+++++++++++++++++

To keep the highest resolution FFT-grid in all the calculations, run `PAIREF` with an option :code:`--constant-grid`. The grid is then controlled by a *REFMAC5* keyword `SHANnon_factor <http://www.ccp4.ac.uk/html/refmac5/keywords/xray-general.html#shan>`_.

Problems
--------

Something is not working? Are you worried that you did not understand well? Is an important feature missing? Do you like our project? Do not hesitate -- please write us: `martin.maly@fjfi.cvut.cz <mailto:martin.maly@fjfi.cvut.cz>`_.
