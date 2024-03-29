# coding: utf-8
from __future__ import print_function
import sys
import os
import subprocess
import datetime
from math import sqrt, pow
from collections import OrderedDict  # Python 2.7
from .settings import warning_dict, date_time, settings
from .commons import twodec, twodecname, fourdec, extract_from_file
from .commons import warning_my, Popen_my, pick_work_free_from_csv_line


BINS_LOW = 10


def which(program):
    """Checks if `program` exists and finds its location. Analogy of the
    `which` GNU/Linux command.

    Args:
        program (str): Name of an executable

    Returns:
        str: Path of an executable location
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program) or is_exe(program + ".exe"):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
            exe_file = os.path.join(path, program + ".exe")
            if is_exe(exe_file):
                return exe_file
    return None


def welcome(args, pairef_version):
    '''Print introduction information about the module and
    input parameters.

    Args:
        args (parser): Input arguments (including e. g. name of the project) \
                       parsed by `argparse` via function process_arguments()
        pairef_version (str)

    Returns:
        bool: True
    '''
    import getpass
    import socket


    print("""
 __   _  ___ __  ___ ___
 )_) /_)  )  )_) )_  )_
/   / / _(_ / \ (__ (
""")
    print("automatic PAIRed REFinement protocol")
    print("version: " + pairef_version)
    print("run date and time: " + date_time)
    print("user@host: " + getpass.getuser() + "@" + socket.gethostname())
    print("")
    print('Please cite: "Paired refinement under the control of PAIREF"')
    print("M. Maly, K. Diederichs, J. Dohnalek, P. Kolenko (2020) IUCrJ 7")
    print("")
    print("Command line arguments: " + " ".join(sys.argv[1:]))
    print("")
    print("Program has been executed with following input parameters:")
    if args.refmac:
        print(" * Refinement software: REFMAC5")
    if args.phenix:
        print(" * Refinement software: phenix.refine")
    print(" * XYZIN: " + args.xyzin)
    print(" * HKLIN: " + args.hklin)
    if args.hklin_unmerged:
        print(" * HKLIN unmerged: " + args.hklin_unmerged)
    if args.libin:
        print(" * LIBIN: " + args.libin)
    if args.tlsin:
        print(" * TLSIN: " + args.tlsin)
    print(" * Project name: " + str(args.project))
    # if args.step:
    #     print(" * Resolution step in angstroem: " + str(args.step))
    # if args.n_shells:
    #     print(" * Number of res. shells: " + str(args.n_shells))
    if args.res_shells:
        print(" * Resolution shells: " + str(args.res_shells))
    if args.weight:
        print(" * Weight matrix: " + str(args.weight))
    if args.tls_ncyc:
        print(" * Number of number of cycles of TLS refinement: "
              "" + str(args.tls_ncyc))
    if args.ncyc:
        print(" * Number of refinement cycles that will be performed in "
              "every resolution step: " + str(args.ncyc))
    if args.prerefinement_ncyc:
        print(" * Number of pre-refinement cycles that will be performed "
              "before the paired refinement protocol: "
              "" + str(args.prerefinement_ncyc))
    if args.complete_cross_validation:
        print(" * Complete cross-validation will be performed.")

    if (args.complete_cross_validation or args.no_modification or
            args.reset_bfactor or args.add_to_bfactor or args.set_bfactor or
            args.shake_sites):
        print(" * Modification of the input structure model:")
    if args.reset_bfactor:
        print("   - Reset B-factors to the mean value")
    if args.add_to_bfactor:
        print("   - Add value to B-factors: " +
              twodec(args.add_to_bfactor))
    if args.set_bfactor:
        print("   - Set B-factors to the value: " + twodec(args.set_bfactor))
    if args.shake_sites:
        print("   - Randomize coordinates with the "
              "given mean error value: " + twodec(args.shake_sites))

    if args.constant_grid:
        print(" * The same FFT grid will be kept through the whole paired "
              "refinement.")
    if args.comin:
        print(" * Com file for REFMAC5: " + args.comin)
    if args.defin:
        print(" * Keyword file for phenix.refine: " + args.defin)
    if args.test:
        print(" * Light-testing mode (REFMAC5 will not be executed).")
    print("")
    return True


def create_workdir(project):
    '''Create directory `pairef_`project``, if exists, create directory
    `pairef_`project`_new` (recursively).

    Args:
        project (str): Name of the project

    Returns:
        str: Name of the created working directory
    '''
    workdir = "pairef_" + project
    # Make directory for the project
    try:
        os.mkdir(workdir)
    except OSError:
        warning_my("workdir",
                   "This directory contains directory "
                   "`pairef_" + project + "` "
                   "already, choosing another directory name.")
        workdir = create_workdir(project + "_new")
    return workdir


class output_log:
    "Set to write `sys.stdout` to screen and also in file `PAIREF_out.log`."
    # Not working for STDERR! TODO! write on internet...
    def __init__(self, stdout, filename):
        self.stdout = stdout
        self.logfile = open(filename, 'a')

    def write(self, text):
        self.stdout.write(text)
        self.logfile.write(text)

    def close(self):
        self.stdout.close()
        self.logfile.close()

    def flush(self):
        self.stdout.close()
        self.logfile.close()


def def_res_shells(args, refinement, res_high_mtz, res_low=999):
    """Determine high resolution shells and number of low resolution bins.

    The first value of list should be the resolution of data which
    were used for the refinement of the input structure model.
    If explicit definition (args.res_shells) is set, test its correctness.
    If it is valid, use, if not, define it automatically (shell step 0.05 A)

    Args:
        args (parser): Input arguments (including e. g. name of the project) \
                       parsed by `argparse` via function process_arguments()
        refinement (str): "refmac" or "phenix"
        res_high_mtz (float): High resolution diffraction limit of args.hklin
        res_low (float): Low resolution diffraction limit of args.hklin

    Returns:
        (tuple):
            * shells (*list*)
            * n_bins_low (*int*)
            * n_flag_sets (*int*)
            * default_shells_definition (*bool*)
    """
    # If the resolution of input model == resolution of diffr. data, abort
    if twodec(args.res_init) <= twodec(res_high_mtz):
        sys.stderr.write(
            "ERROR: Given input MTZ file " + args.hklin + " contain data "
            "only up to resolution " + twodec(res_high_mtz) + " A that is "
            "the same or lower than the initial high-resolution diffraction "
            "limit " + twodec(args.res_init) + " A. Nothing to do.\n"
            "Aborting.\n")
        sys.exit(1)
    # Estimation of
    #   1. a number of low resolution bins and
    #   2. a number of free reflection sets
    # If not successful, just divide into 12 resolution bins
    # and consider 20 free reflection sets.
    n_i_obs = 0
    n_i_obs_low = 0
    n_flag_sets = 0

    if refinement == "refmac":
        # using `mtzdump` from CCP4 to get both n_i_obs_low and n_flag_sets
        tool = "mtzdump"
        p = Popen_my(["mtzdump", "HKLIN", args.hklin],
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                     shell=settings["sh"])
        com = "STATS NBIN 1 RESO " + twodec(res_low) + " " \
            "" + twodec(args.res_init) + "\n end\n"
        output, err = p.communicate(com)
        # rc = p.returncode
        try:
            for line in output.splitlines():
                if " * Number of Reflections = " in line:
                    n_i_obs = int(line.split()[-1])
                elif len(line.split()) == 12:
                    if "free" in line.split()[-1].lower():
                        n_flag_sets = int(float(line.split()[3])) - \
                            int(float(line.split()[2])) + 1
                elif " No. of reflections used in FILE STATISTICS" in line:
                    n_i_obs_low = int(line.split()[-1])
                    break
        except ValueError:
            n_bins_low = 12
            warning_my("mtzdump",
                       "Definition of binning "
                       "and search for number of sets of free reflection "
                       "using mtzdump were not successful. "
                       "Data will be divided into " + str(n_bins_low) + " "
                       "resolution bins. "
                       "Is the input MTZ file " + args.hklin + " OK?")
    elif refinement == "phenix":
        from iotbx.reflection_file_reader import any_reflection_file
        from iotbx import mtz
        tool = "CCTBX"
        # 1. get n_i_obs_low using CCTBX
        n_i_obs_low_list = []
        n_i_obs_list = []
        hkl_in = any_reflection_file(file_name=args.hklin)
        miller_arrays = hkl_in.as_miller_arrays()
        for column in miller_arrays:
            if "xray" in str(column.observation_type()):
                n_i_obs_list.append(column.size())
                column_res_init = column.resolution_filter(
                    d_max=float(res_low), d_min=float(args.res_init))
                n_i_obs_low_list.append(column_res_init.size())
        n_i_obs = min(n_i_obs_list)
        n_i_obs_low = min(n_i_obs_low_list)
        # 2. get n_flag_sets using CCTBX
        mtz_object = mtz.object(file_name=args.hklin)
        try:  # Python 3
            import io
            out = io.StringIO()
            mtz_object.show_summary(out=out)
        except TypeError:  # Python 2.7 (CCP4 7, Phenix 1.20)
            import cStringIO
            out = cStringIO.StringIO()
            mtz_object.show_summary(out=out)
        out = out.getvalue()
        out_lines = out.splitlines()
        for i in range(len(out_lines)):
            if "free" in out_lines[i].lower():
                j = i
        if "j" in vars():
            n_flag_sets = int(float(out_lines[j].split()[4])) - \
                int(float(out_lines[j].split()[3])) + 1

    if n_flag_sets == 0:
        n_flag_sets = 20
        if args.complete_cross_validation:
            warning_my("flags",
                       "Search for number of sets of free reflection "
                       "using " + tool + " was not successful. "
                       "It will be assumed that 20 sets of free reflection "
                       "are present in the input MTZ file " + args.hklin + ".")
    else:
        if args.complete_cross_validation:
            print(str(n_flag_sets) + " sets of free reflection were found and "
                  "will be used in complete cross-validation.")
    if n_i_obs == 0 or n_i_obs_low == 0:
        n_bins_low = 12
        warning_my("binning",
                   "Definition of binning using " + tool + " was not success"
                   "ful. Data will be divided into " + str(n_bins_low) + " "
                   "resolution bins. "
                   "Is the input MTZ file " + args.hklin + " OK?")
    else:  # If run of mtzdump was succesful
        n_bins_low_start = 12
        n_i_obs_thinner_shell_start = 0
        n_bins_low = n_bins_low_start
        n_i_obs_thinner_shell = n_i_obs_thinner_shell_start
        while n_i_obs_thinner_shell < 2000:
            n_bins_low = n_bins_low - 1
            n_i_obs_thinner_shell = \
                sqrt(2) * n_i_obs_low / pow(n_bins_low, 1.5)
    if n_bins_low <= 1:
        n_bins_low = 2
    # Now `n_bins_low` is ready

    if args.res_shells:
        # Explicit definition of high resolution shells
        default_shells_definition = False
        shells_high = args.res_shells.split(",")
        for i in range(len(shells_high)):
            try:
                shells_high[i] = round(float(shells_high[i]), 2)
            except ValueError:
                sys.stderr.write(
                    "ERROR: Explicit definition of high "
                    "resolution shells (option -r) is not correct. ")
                sys.stderr.write(
                    "Values must be divided using commas without "
                    "any spaces (e.g. 2.1,2.0,1.9).")
                sys.stderr.write("\nAborting.\n")
                sys.exit(1)
            if i != 0:
                if float(shells_high[i - 1]) <= float(shells_high[i]):
                    sys.stderr.write(
                        "ERROR: Explicit definition of high "
                        "resolution shells (option -r) is not correct. ")
                    sys.stderr.write(
                        "Values must be set in the decreasing order "
                        "(e.g. 2.1,2.0,1.9).")
                    sys.stderr.write("\nAborting.\n")
                    sys.exit(1)
        for i in range(len(shells_high)):
            if twodec(shells_high[i]) <= twodec(res_high_mtz):
                shells_high = shells_high[:i] + [round(float(res_high_mtz), 2)]
                break
        if shells_high[0] >= args.res_init:
            sys.stderr.write(
                "ERROR: Explicit definition of high "
                "resolution shells (option -r) is not correct. ")
            sys.stderr.write(
                "Resolution of the first shell (" + twodec(shells_high[0]) + ""
                " A), "
                "which was explicitely defined, is lower than or the same as "
                "the resolution of data which were used for refinement of the "
                "input structure model (" + twodec(args.res_init) + " A).")
            sys.stderr.write("\nAborting.\n")
            sys.exit(1)

    elif args.step and args.n_shells:
        # Defined step (in A) and number of shells to be added
        default_shells_definition = False
        res_fin = args.res_init - args.step * args.n_shells
        if res_fin < 0:
            sys.stderr.write("ERROR: Current setup of resolution shells is "
                             "not valid. High resolution limit of the last "
                             "shell would be negative.\n")
            sys.stderr.write("\nRESOLUTION_INITIAL - (RESOLUTION_STEP * NUM"
                             "BER_OF_SHELLS) = " + twodec(args.res_init) + ""
                             " - ( " + str(args.step) + " * "
                             "" + str(args.n_shells) + ""
                             " ) = " + twodec(res_fin) + " < 0")
            sys.stderr.write("\nAborting.\n")
            sys.exit(1)
        shells_high = []
        for i in range(args.n_shells):
            new_shell = args.res_init - (i + 1) * args.step
            if twodec(new_shell) >= twodec(res_high_mtz):
                shells_high.append(args.res_init - (i + 1) * args.step)
            elif i == 0:  # and new_shell < res_high_mtz:
                shells_high = round(float(twodec(res_high_mtz)), 2)
                break
            else:  # i >= 1 and new_shell < res_high_mtz:
                if twodec(res_high_mtz) != twodec(shells_high[i - 1]):
                    shells_high.append(round(float(twodec(res_high_mtz)), 2))
                break
    else:
        # Default setting - 0.05A wide high resolution shells
        shell_width = 0.05
        default_shells_definition = True
        shells_high = []
        res_cur = args.res_init - shell_width
        while round(res_cur, 3) >= round(res_high_mtz, 3):
            shells_high = shells_high + [res_cur]
            res_cur = res_cur - shell_width
        shells_high = shells_high[:-1]
        shells_high = shells_high + [res_high_mtz]

        # Automatic estimation of high res. shell (using `mtzdump` from CCP4)
        # (shells with increasing number of reflections)
        # default_shells_definition = True
        # n_i_obs_high = n_i_obs - n_i_obs_low
        # n_i_obs_remaining = n_i_obs_high
        # n_bins_high = 0
        # i = 1
        # n_i_obs_shells_high = n_i_obs_thinner_shell * sqrt(n_bins_low + i)
        # while n_i_obs_remaining > n_i_obs_shells_high:
            # i += 1
            # n_i_obs_shell_high = n_i_obs_thinner_shell * sqrt(n_bins_low + i)
            # n_i_obs_remaining = n_i_obs_remaining - n_i_obs_shell_high
            # n_bins_high += 1
        # p = Popen_my(["mtzdump", "HKLIN", args.hklin],
                             # stdin=subprocess.PIPE,
                             # stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # com = "STATS NBIN " + str(n_bins_high) + " RESO " \
            # "" + twodec(args.res_init) + "" \
            # " " + twodec(res_high_mtz) + "\n end\n"
        # output, err = p.communicate(com)
        # shells_high = []
        # try:
            # for line in output.splitlines():
                # if " PARTIAL FILE STATISTICS for resolution range" in line:
                    # shells_high.append(
                        # round(1/sqrt(float(line.split()[-1])), 2))
                    # #  there is no break
        # except:
            # sys.stderr.write(
                # "ERROR: Definition of high resolution shells  using mtzdump "
                # "was not successful. Is the MTZ file OK?")
            # sys.stderr.write(
                # "Use manual definition of the high resolution shells to "
                # "overcome this problem (option -r).")
            # sys.stderr.write("\nAborting.\n")
            # sys.exit(1)
    if shells_high == []:
        sys.stderr.write(
            "ERROR: The given definition of high-resolution shells requires "
            "more data in high resolution than are available in the given "
            "input MTZ file " + args.hklin + " that contain data "
            "only up to resolution " + twodec(res_high_mtz) + " A. Nothing to do.\n"
            "Aborting.\n")
        sys.exit(1)
    shells = [args.res_init] + shells_high
    return shells, n_bins_low, n_flag_sets, default_shells_definition


def res_high_from_xyzin(xyzin, format=".pdb"):
    """Finds a line containing `RESOLUTION RANGE HIGH` in the file `xyzin`,
    picks the last word of the line (that should be the high resolution)
    and rounds it to two decimals. If it is not successful, returns -1.

    Args:
        xyzin (str): Input structure model (PDB or mmCIF format)

    Returns:
        float: initial high resolution limit rounded to two decimals or -1 if \
               it was not found
    """
    if "cif" in format:
        res_high = extract_from_file(
            filename=xyzin, searched="_refine.ls_d_res_high", skip_lines=0,
            n_lines=1, nth_word=-1, not_found='N/A', get_first=True)
    else:  # pdb
        res_high = extract_from_file(
            filename=xyzin, searched="RESOLUTION RANGE HIGH", skip_lines=0,
            n_lines=1, nth_word=-1, not_found='N/A', get_first=True)
    if res_high == 'N/A':
        return -1
    try:
        res_high = round(float(res_high), 2)
    except ValueError:
        return -1
    return res_high


def res_from_mtz(hklin):
    """Finds the low resolution diffraction limit of data `hklin` using CCTBX.

    Args:
        hklin (str): Name of diffraction data MTZ file

    Returns:
        (tuple):
            * res_low (*float*): Low resolution diffraction limit
            * res_high (*float*): High resolution diffraction limit
    """
    from iotbx.reflection_file_reader import any_reflection_file
    hkl_in = any_reflection_file(file_name=hklin)
    miller_arrays = hkl_in.as_miller_arrays()
    res_low = None
    res_high = None
    for column in miller_arrays:
        if "xray" in str(column.observation_type()):
            res_low = column.d_max_min()[0]
            res_high = column.d_max_min()[1]
            break
    return res_low, res_high


def res_from_hklin_unmerged(hklin_unmerged):
    """
    Finds a resolution range for given unmerged diffr. data file 
    `hklin_unmerged`.
    In the case of an ASCII file from XDS, find it by
    searching the option `INCLUDE_RESOLUTION_RANGE` in the file.
    In the case of an MTZ file, find it using `xia2`.
    In the case of a SCA file, do not check.

    Args:
        hklin_unmerged (str): Name of the unmerged diffr. data file

    Returns:
        (tuple):
            * res_low (*float*): Low resolution diffraction limit
            * res_high (*float*): High resolution diffraction limit
    """
    res_high_from_hklin_unmerged = None
    res_low_from_hklin_unmerged = None
    # may be an MTZ file
    if  "mtz" in hklin_unmerged.split(".")[-1].lower():
        ### Solution 1 - CCTBX
        ### res_low_from_hklin_unmerged, res_high_from_hklin_unmerged = \
        ###     res_from_mtz(hklin_unmerged)  # This needs to much RAM...

        ## Solution 2 - mtzdump
        ## p = Popen_my(["mtzdump", "HKLIN", hklin_unmerged],
        ##                      stdin=subprocess.PIPE,
        ##                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ## com = "end\n"
        ## output, err = p.communicate(com)
        ## for i in range(len(output.splitlines())):
        ##     if " *  Resolution Range :" in output.splitlines()[i]:
        ##         res_low_from_hklin_unmerged = \
        ##             float(output.splitlines()[i + 2].split()[-5])
        ##         res_high_from_hklin_unmerged = \
        ##             float(output.splitlines()[i + 2].split()[-3])
        ##         break

        # Solution 3 - xia2 - works with both CCP4 and phenix
        from xia2.Wrappers.CCP4.Mtzdump import Mtzdump
        m = Mtzdump()
        m.set_hklin(hklin_unmerged)
        m.dump()
        res_low_from_hklin_unmerged, res_high_from_hklin_unmerged = \
            m.get_resolution_range()
        
    # may be an XDS ASCII file
    elif "hkl" in hklin_unmerged.split(".")[-1].lower():
        hklin_unmerged_head = []
        try:  # Python 2 and 3 support
            xrange
        except NameError:
            xrange = range
        with open(hklin_unmerged, "r") as hklin_unmerged_file:
            hklin_unmerged_head = \
                [next(hklin_unmerged_file) for x in xrange(1000)]
        for line in hklin_unmerged_head:
            if "INCLUDE_RESOLUTION_RANGE" in line:
                res_high_from_hklin_unmerged = float(line.split()[-1])
                res_low_from_hklin_unmerged = float(line.split()[-2])
                break
    if not res_high_from_hklin_unmerged or not res_low_from_hklin_unmerged:
        warning_my("merging_stats",
                   "The check of a resolution range of the unmerged data "
                   "file " + hklin_unmerged + " was not successful. ")
        res_low_from_hklin_unmerged = float("inf")
        res_high_from_hklin_unmerged = 0
    return res_low_from_hklin_unmerged, res_high_from_hklin_unmerged


def check_refinement_software(args, versions_dict, refinement="refmac"):
    """
    Check if the input structure model `args.xyzin` was refined in REFMAC5 
    or phenix.refine and in which version (PDB and mmCIF format is accepted).
    If it was not refined or if it was refined in another
    version of REFMAC5 or phenix.refine than is now installed, write warning.

    Args:
        args.xyzin (str): Input structure model (PDB or mmCIF format)
        versions_dict (dict): Dictionary containing a key `refmac_version` \
                              or `phenix_version`.
        refmac_version_installed (str): Version of REFMAC5 or phenix.refine \
                                        that is installed

    Returns:
        str: version of REFMAC5 or phenix.refine  that was used for the \
             refinement of `args.xyzin`, if it was not found, "N/A" is returned
    """
    if refinement == "refmac":
        refinement_name = "REFMAC5"
        version_installed = versions_dict["refmac_version"]
        if "cif" in settings["pdbORmmcif"]:
            refmac_confirm = extract_from_file(
                filename=args.xyzin, searched="_software.name", skip_lines=0,
                n_lines=1, nth_word=-1, not_found="N/A", get_first=True)
            if refmac_confirm == "refmac":
                version_xyzin = extract_from_file(
                    filename=args.xyzin, searched="_software.version",
                    skip_lines=0, n_lines=1, nth_word=-1, not_found="N/A",
                    get_first=True)
                version_xyzin = version_xyzin.replace("'", "")
            else:
                version_xyzin = "N/A"
        else:  # pdb
            version_xyzin = extract_from_file(
                filename=args.xyzin, searched="REFMAC", skip_lines=0,
                n_lines=1, nth_word=-1, not_found="N/A", get_first=True)
    elif refinement == "phenix":
        refinement_name = "phenix.refine"
        version_installed = versions_dict["phenix_version"]
        if "cif" in settings["pdbORmmcif"]:
            phenix_confirm = extract_from_file(
                filename=args.xyzin, searched="phenix.refine", skip_lines=0,
                n_lines=1, nth_word=1, not_found="N/A", get_first=True)
            if phenix_confirm == "phenix.refine":  # TO DO check
                version_xyzin = extract_from_file(
                    filename=args.xyzin, searched="phenix.refine",
                    skip_lines=0, n_lines=1, nth_word=2, not_found="N/A",
                    get_first=True)
            else:
                version_xyzin = "N/A"          
        else:  # pdb
            version_xyzin = extract_from_file(
                filename=args.xyzin,
                searched="REMARK   3   PROGRAM     : PHENIX",
                skip_lines=0, n_lines=1, not_found="N/A", get_first=True)
            version_xyzin = version_xyzin[0]
            if version_xyzin != "N/A":
                version_xyzin = version_xyzin.split()[5]
                version_xyzin = version_xyzin.replace("(", "")
                version_xyzin = version_xyzin.replace(")", "")
        if version_xyzin != "N/A":
            version_xyzin = version_xyzin.split("_")[0]
    if args.prerefinement_ncyc:  # TO DO is it correct?
        return version_xyzin
    if version_xyzin == "N/A":
        warning_my("refinement_not_before",
                   "The input structure model `" + args.xyzin + "` seems not "
                   "to be refined in " + refinement_name + ". "
                   "The obtained results could be misleading. "
                   "Consider refinement of the structure model "
                   "" + args.xyzin + " in " + refinement_name + ". This could be "
                   "performed using the arguments --" + refinement + " "
                   "--prerefinement-ncyc")
    else:
        if version_xyzin != version_installed:
            warning_my("refinement_version_mismatch",
                       "The version of " + refinement_name + " used for "
                       "refinement of the  input structure model "
                       "`" + args.xyzin + "` ("
                       "" + version_xyzin + ") is not the same as the "
                       "version of " + refinement_name + " that is now "
                       "installed and used during paired refinement ("
                       "" + version_installed + "). "
                       "The results from PAIREF could be misleading. "
                       "Consider refinement of the structure model "
                       "`" + args.xyzin + "` using the installed version of"
                       " " + refinement_name + ". This could be performed "
                       "using the arguments --" + refinement + " "
                       "--prerefinement-ncyc")
    return version_xyzin


def res_opt(shell, args, refinement="refmac"):
    """
    Finds optical resolution running command
    :code:`sfcheck -f project_Rflag_shellA.mtz -m project_Rflag_shellA.pdb`
    and writes the gained value into project_Optical_resolution.csv
    PDB format is required, does not work with mmCIF.

    Args:
        shell (float): Current resolution shell
        args (parser): Input arguments (including e. g. name of the project) \
                       parsed by `argparse` via function process_arguments()

    Returns
        float: optical resolution
    """
    prefix = args.project + "_R" + str(args.flag).zfill(2) + "_" + \
        twodecname(shell) + "A"
    if refinement == "phenix":
        prefix += "_001"
    hklin = prefix + ".mtz"
    xyzin = prefix + ".pdb"
    command = ["sfcheck", "-f", hklin, "-m", xyzin]
    logfilename = prefix + "_sfcheck.out"
    with open(logfilename, "w") as logfile:
        p = Popen_my(
            command, stdout=logfile, stderr=logfile, shell=settings["sh"])
        out, err = p.communicate()

    res_opt = 0
    with open(logfilename, "r") as logfile:
        for line in logfile.readlines():
            if "Optical Resolution" in line:
                res_opt = float(line.split()[-1])
                break
    if os.path.isfile("sfcheck.xml"):
        os.rename("sfcheck.xml", prefix + "_sfcheck.xml")
    if os.path.isfile("sfcheck_xxxx.ps"):
        os.rename("sfcheck_xxxx.ps", prefix + "_sfcheck.ps")
    if os.path.isfile("sfcheck.log"):
        os.rename("sfcheck.log", prefix + "_sfcheck.log")
    # TODO: warning
    # TODO help
    csvfilename = args.project + "_Optical_resolution.csv"
    if not os.path.isfile(csvfilename):
        with open(csvfilename, "w") as csvfile:
            csvfile.write("# Nominal resolution          Optical resolution\n")
    with open(csvfilename, "a") as csvfile:
        csvfile.write(
            twodec(shell) + 26 * " " + twodec(res_opt) + "\n")
    return float(twodec(res_opt))


def calculate_merging_stats(hklin_unmerged, shells, project, bins_low,
                            res_low_from_hklin_unmerged=float("inf"),
                            res_high_from_hklin_unmerged=0):
    """
    For given file `hklin_unmerged`, calculate the merging statistics
    using CCTBX.

    Args:
        hklin_unmerged (str): Name of the unmerged diffraction data file
        shells (list)
        project (str): Name of the project
        bins_low (list)
        costant_reflections_in_bins (bool)

    Returns:
        str: Name of a CSV file where the calculated statistics has been saved
    """
    print("\n     * Calculating merging statistics...", end="")

    ## Calculate statistics
    def calculate_merging_stats_run_cctbx(project, hklin, res_high,
                                          res_low=None, n_bins=None,
                                          data_labels=None):
        import iotbx.merging_statistics
        import gc
        gc.collect()

        # Determine data labels (can be ambiguous in MTZ file)
        if "mtz" in hklin.split(".")[-1] and not data_labels:
            try:
                from iotbx.reflection_file_reader import any_reflection_file
                miller_arrays = any_reflection_file(hklin).as_miller_arrays()
                labels_i = []
                for label in miller_arrays:
                    if "xray.intensity" in str(label.observation_type()):
                        label_here = str(label.info()).split(":")[-1]
                        if ",merged" in label_here:
                            label_here = label_here.replace(",merged", "")
                        labels_i.append(label_here)
                if len(labels_i) == 0:
                    print("No intensity arrays were found in the file " + hklin)
                    print(".\nMerging statistics could not be calculated.")
                    return False
                elif len(labels_i) == 1:  # Default behaviour is OK
                    # label_imean = label_imean[0]
                    pass  # data_labels = None
                else:  # Multiple intensity arrays - specifying one:
                    for label in labels_i:
                        if "I,SIGI" in label:
                            data_labels = label
                            break
                if not data_labels:
                    data_labels = labels_imean[0]
            except:
                pass  # Try to continue...
        i_obs = iotbx.merging_statistics.select_data(file_name=hklin,
                                                     data_labels=data_labels)
        result = iotbx.merging_statistics.dataset_statistics(
            i_obs=i_obs,
            # crystal_symmetry=symm,
            d_min=res_high,   # My change
            d_max=res_low,    # My change
            n_bins=n_bins,
            log=None
            )
        logfilename = project + "_merging_stats_" + twodecname(res_high) + "" \
            "A.log"
        with open(logfilename, "w") as logfile:
            result.show(out=logfile, header=False)
        return data_labels

    bins_total_proposed = bins_low + shells
    bins_total = []
    labels = None
    for i in range(len(bins_total_proposed)-1):
        if res_low_from_hklin_unmerged < bins_total_proposed[i + 1] \
                or res_high_from_hklin_unmerged > bins_total_proposed[i]:
            warning_my("merging_stats",
                       "Merging statistics of the data in "
                       "bin " + twodec(bins_total_proposed[i]) + "-"
                       "" + twodec(bins_total_proposed[i + 1]) + " A "
                       "could not be calculated as the unmerged data "
                       "(file " + hklin_unmerged + ") do not contain "
                       "data in this resolution range.")
        else:
            bins_total.append(bins_total_proposed[i])
            bins_total.append(bins_total_proposed[i + 1])
            # Remove duplicates from bins_total and keep order
            bins_total = list(OrderedDict.fromkeys(bins_total))
            labels = calculate_merging_stats_run_cctbx(
                project, hklin_unmerged, res_high=bins_total[i + 1],
                res_low=bins_total[i], n_bins=1, data_labels=labels)
        print(" .", end="")
    if labels:
        print("\n       Using labels=" + str(labels))
    print("")

    # Prepare a csv file header
    csvfilename = project + "_merging_stats.csv"
    with open(csvfilename, "w") as csvfile:
        csvfile.writelines("#shell d_max  d_min   #obs  #uniq   mult.  %comp"
                           "       <I>  <I/sI>    r_mrg   r_meas    r_pim   "
                           "r_anom   cc1/2   cc_ano     cc* \n")

    ## Collect statistics, calculate CC*-values and save them to CSV file

    def calculate_CCstar(lines, shell=1):
        """For each string of lines array, find 11th word which should be
        value of CChalf and calculate CC*. Put the result into lines
        array. Add number of bin to lines array.

        Args:
            lines (list): list of strings
            shell (int)

        Returns:
            list
        """
        for i in range(len(lines)):
            # Compute CC*
            bin_CChalf = lines[i].split()[-2]
            if float(bin_CChalf) < 0:
                bin_CCstar = "N/A"
                warning_my("CC*", "A CC*-value for a particular shell "
                           "could not be calculated as it is undefined "
                           "for a negative value of CC1/2. "
                           "(CC1/2 = " + str(bin_CChalf) + ")")
            else:
                try:
                    bin_CCstar = \
                        sqrt(2 * float(bin_CChalf) / (1 + float(bin_CChalf)))
                except ValueError:
                    bin_CCstar = "N/A"
                    warning_my("CC*", "A CC*-value for a particular "
                               "shell could not be calculated. "
                               "(CC1/2 = " + str(bin_CChalf) + ")")
            lines[i] = lines[i].replace("\n", "")
            lines[i] = lines[i] + "   " + fourdec(bin_CCstar) + "\n"
            # Insert number of the particular shell
            lines[i] = str(shell).zfill(2) + "    " + lines[i]
            shell = shell + 1
        return lines

    # Insert statistics relating up to resolution res_init
    for i in range(len(bins_total) - 1):
        logfilename = project + "_merging_stats_" \
            "" + twodecname(bins_total[i + 1]) + "A.log"
        line = extract_from_file(logfilename,
                                 searched="Statistics by resolution bin:",
                                 skip_lines=2, n_lines=1)
        os.remove(logfilename)  # Little clean-up
        line = calculate_CCstar(line, shell=i + 1)
        with open(csvfilename, "a") as csvfile:
            csvfile.writelines(line)
    return csvfilename


def run_baverage(project, xyzin, res_init):
    """Finds average value of B-factors of all the atoms in the structure
    model `xyzin` using `baverage` from the CCP4 package.
    This creates files with a prefix `project_twodecname(res_init)A_baverage`.
    Then returns mean B-factor for all the atoms of the structure model to
    the obtained value.

    Args:
        project (str): Name of the project
        xyzin (str): Filename of structure model in PDB or mmCIF format
        res_init (float): Resolution of the input structure model

    Returns:
        float: Mean B-factor for all the atoms
    """
    prefix = project + "_" + twodecname(res_init) + "A_baverage"
    logout = prefix + ".log"
    if "cif" in settings["pdbORmmcif"]:
        xyzout = prefix + settings["pdbORmmcif"]
    else:  # pdb
        xyzout = prefix + settings["pdbORmmcif"]
    command = ["baverage", "XYZIN", xyzin, "RMSTAB",
               prefix + ".tab", "XYZOUT", xyzout]
    com = "end\n"
    with open(logout, "w") as logfile:
        p = Popen_my(command, stdin=subprocess.PIPE,
                     stdout=logfile, shell=settings["sh"])
        p.communicate(com)

    baverage = extract_from_file(filename=logout,
                                 searched="AVERAGE B VALUE FOR ALL ATOMS",
                                 skip_lines=0,
                                 n_lines=1,
                                 nth_word=-1)
    baverage = float(baverage)
    print("Average B-factor for all atoms: " + twodec(baverage))
    # if add_to_bfactor:
        # baverage = str(float(baverage) + add_to_bfactor)
    #     bfac_set = baverage + add_to_bfactor
    #     print("B-factor after application of the parameter "
    #           "--add-to-bfactor: " + twodec(bfac_set))
    # else:
    #     bfac_set = baverage
    # bfac_set = baverage

    # Seems that it is not working well for anisotropic B-factors
    # http://www.ccp4.ac.uk/html/baverage.html
    # BLIMIT <bminmc> <bmaxmc> <bminsc> <bmaxsc>
    # prefix = project + "_" + twodecname(res_init) + "A_baveraged"
    # logout = prefix + ".log"
    # xyzout = prefix + ".pdb"
    # com = "BLIMIT " + baverage + " " + baverage + " " + baverage + " " \
    #     "" + baverage + "\nend\n"
    # command = ["baverage", "XYZIN", xyzin, "RMSTAB",
    #            prefix + ".tab", "XYZOUT", xyzout]
    # with open(logout, "w") as logfile:
    #     p = Popen_my(command, stdin=subprocess.PIPE,
    #                          stdout=logfile)
    #    p.communicate(com)
    # return xyzout
    # return bfac_set
    return baverage


def run_bmean_iotbx(project, xyzin):
    """Finds average value of B-factors of all the atoms in the structure
    model `xyzin` using `pdb` from `iotbx`.
    Then returns mean B-factor for all the atoms of the structure model to
    the obtained value.

    Args:
        project (str): Name of the project
        xyzin (str): Filename of structure model in PDB or mmCIF format

    Returns:
        float: Mean B-factor for all the atoms
    """
    from iotbx import pdb
    pdb_inp = pdb.input(file_name=xyzin)
    atoms = pdb_inp.atoms()
    bfactors = atoms.extract_b()
    baverage = float(bfactors.format_mean("%5.2f"))
    return baverage


def run_pdbtools(args, baverage=0):
    """Modify the input structure model `args.xyzin` by `mmtbx.pdbtools`.
    The procces is controlled by `args.reset_bfactor`, `args.add_to_bfactor`,
    `args.set_bfactor`, and `args.shake_sites`. Reset of B-factors requires
    set of `baverage`. The modified structure model:
    `project_twodecname(res_init)A_modified`.

    Args:
        args (parser): Input arguments (including e. g. name of the project) \
                       parsed by `argparse` via function process_arguments()
        baverage (float): Mean B-factor for all the atoms

    Returns:
        str: Filename of the modified structure model (PDB or mmCIF format)
    """
    import shutil
    if not (args.reset_bfactor or args.add_to_bfactor or args.set_bfactor or
            args.shake_sites):
        # Nothing to do
        return args.xyzin

    prefix_phil = args.project + "_" + twodecname(args.res_init) + "A"
    prefix = args.project + "_" + twodecname(args.res_init) + "A_modified"
    logout = prefix + ".log"
    if "cif" in settings["pdbORmmcif"]:
        xyzout = prefix + ".cif"
    else:
        xyzout = prefix + ".pdb"

    phil_import_successful = None
    try:  # new versions of CCP4 or PHENIX
        from libtbx import phil
        from mmtbx.programs import pdbtools
        from iotbx.data_manager import DataManager
        phil_import_successful = True
    except ImportError:
        pass
    if phil_import_successful:
        dm = DataManager()
        dm.set_overwrite(True)
        model = dm.get_model(args.xyzin)
        phil_master = phil.parse(pdbtools.Program.master_phil_str,
                                 process_includes=True)
        obj_work = phil_master.extract()
        adp_obj = phil_master.extract().modify.adp[0]
        # prepare arguments
        if args.reset_bfactor and baverage:
            bfactor = baverage
            if args.add_to_bfactor:
                bfactor += float(args.add_to_bfactor)
                print("B-factor after application of the parameter "
                      "--prerefinement-add-to-bfactor: " + twodec(bfactor))
            adp_obj.set_b_iso = float(bfactor)
            # pdbtools_args.append("set_b_iso=" + twodec(bfactor))
        elif args.set_bfactor:
            adp_obj.set_b_iso = float(args.set_bfactor)
            # pdbtools_args.append("set_b_iso=" + twodec(args.set_bfactor))
        if args.add_to_bfactor and not args.reset_bfactor:
            adp_obj.shift_b_iso = float(args.add_to_bfactor)
            # pdbtools_args.append("shift_b_iso=" + twodec(args.add_to_bfactor))
        if adp_obj.set_b_iso or adp_obj.shift_b_iso:
            obj_work.modify.adp.append(adp_obj)
        if args.shake_sites:
            sites_obj = phil_master.extract().modify.sites[0]
            sites_obj.shake = float(args.shake_sites)
            obj_work.modify.sites.append(sites_obj)
            # pdbtools_args.append("shake=" + twodec(args.shake_sites))
        obj_work.output.prefix = prefix_phil
        obj_work.output.suffix = "_modified"
        return_phil = phil_master.format(obj_work)
        print("Modification of the input structure model using pdbtools...")
        Pdbtools = pdbtools.Program(data_manager=dm, params=obj_work)
        Pdbtools.run()
        return xyzout

    else:  #   old versions of CCP4 or PHENIX (no phil)
        pdbtools_args = []
        if args.reset_bfactor and baverage:
            bfactor = baverage
            if args.add_to_bfactor:
                bfactor += args.add_to_bfactor
                print("B-factor after application of the parameter "
                      "--prerefinement-add-to-bfactor: " + twodec(bfactor))
            pdbtools_args.append("set_b_iso=" + twodec(bfactor))
        elif args.set_bfactor:
            pdbtools_args.append("set_b_iso=" + twodec(args.set_bfactor))
        if args.add_to_bfactor and not args.reset_bfactor:
            pdbtools_args.append("shift_b_iso=" + twodec(args.add_to_bfactor))
        if args.shake_sites:
            pdbtools_args.append("shake=" + twodec(args.shake_sites))

        pdbtools_args.append("file_name=" + xyzout)
        print("Modification of the input structure model - pdbtools arguments: " + 
            " ".join(pdbtools_args))
        if which("phenix.pdbtools"):
            with open(logout, "w") as logfile:
                p = Popen_my(["phenix.pdbtools", args.xyzin] + pdbtools_args,
                             stdout=logfile, shell=settings["sh"])
                p.communicate()
        else:
            import mmtbx.command_line.pdbtools
            pdbtools_args.append("model_file_name=" + args.xyzin)
            with open(logout, "w") as logfile:
                mmtbx.command_line.pdbtools.run(pdbtools_args, out=logfile,
                                                replace_stderr=False)
    
    if os.path.isfile(xyzout):
        if "cif" in settings["pdbORmmcif"]:
            # refmac required .mmcif (.cif does not work)
            shutil.copy2(xyzout, xyzout[:-4] + ".mmcif")
            xyzout = xyzout[:-4] + ".mmcif"
        return xyzout
    else:
        sys.stderr.write("ERROR: File " + xyzout + " has not been created "
                         "by pdbtools. Check the log file " + logout + "."
                         "\nAborting.\n")
        sys.exit(1)


def suggest_cutoff(args, shells, n_bins_low, flag):
    shells_high = shells[1:]

    # Pick overall R-values
    csvfilename = args.project + "_R-values.csv"
    Rwork_overall_list = []
    Rfree_overall_list = []
    errors_Rwork_overall_list = []
    errors_Rfree_overall_list = []
    errors = False
    with open(csvfilename, "r") as csvfile:
        for line in csvfile.readlines():
            Rwork_overall_list, Rfree_overall_list, errors_Rwork_overall_list, \
                errors_Rfree_overall_list, continue_sign = \
                pick_work_free_from_csv_line(
                    line, Rwork_overall_list, Rfree_overall_list,
                    errors_Rwork_overall_list, errors_Rfree_overall_list, errors)
            if continue_sign:
                continue

    # Pick CCstar (and CC1/2)
    csvfilename_merg = args.project + "_merging_stats.csv"
    if os.path.isfile(csvfilename_merg):
        CChalf_list = []
        CCstar_list = []
        with open(csvfilename_merg, "r") as csvfile:
            for line in csvfile.readlines()[1 + n_bins_low:]:  #  high-res only
                if line.lstrip()[0] == "#":  # If it is a comment
                    continue                 # do not load data
                try:
                    CCstar_list.append(float(line.split()[-1]))
                except ValueError:
                    CCstar_list.append((float("nan")))
                try:
                    CChalf_list.append(float(line.split()[-3]))
                except ValueError:
                    CChalf_list.append((float("nan")))

    # Pick R-values and CCwork for the highest resolution shells
    if not args.complete_cross_validation:
        Nfree_shell_list = []
        Rwork_shell_list = []
        Rfree_shell_list = []
        CCwork_shell_list = []
        # CCfree_shell_list = []
        for shell in shells_high:
            csvfilename = args.project + "_R" + str(flag).zfill(2) + "_" \
                "" + twodecname(shell) + "A.csv"
            with open(csvfilename, "r") as csvfile:
                # for line in csvfile.readlines()[7 + n_bins_low:]:  #  high-res only
                line = csvfile.readlines()[-1]  # highest-res only
            try:
                Nfree_shell_list.append(float(line.split()[5]))
            except ValueError:
                Nfree_shell_list.append((float("nan")))
            try:
                Rwork_shell_list.append(float(line.split()[6]))
            except ValueError:
                Rwork_shell_list.append((float("nan")))
            try:
                Rfree_shell_list.append(float(line.split()[7]))
            except ValueError:
                Rfree_shell_list.append((float("nan")))
            try:
                CCwork_shell_list.append(float(line.split()[8]))
            except ValueError:
                CCwork_shell_list.append((float("nan")))

    # Rate shells
    # 1  overall Rfree decreased
    # 2  overall Rfree did not increase more than 0.0002 (including) and
    #    Rwork increased
    # 3  0.40 >= Rfree > 0.45 in the highest resolution shell but Nfree < 50
    # 4  Rfree > 0.45 in the highest resolution shell but Nfree < 50
    # ---
    # 5  0.40 >= Rfree > 0.45 in the highest resolution shell
    # 6  0.40 >= Rwork > 0.45 in the highest resolution shell
    # 7  overall Rfree increased > 0.0002
    # ---
    # 8  Rfree >= 0.45 in the highest resolution shell
    # 9  Rwork >= 0.45 in the highest resolution shell
    # 10 CCwork > CC*
    # 11 CC1/2 <= 0, CC* undefined
    # 12 overall Rwork increased > 0.0100 (increased pretty much)
    # cutoff = shells[0]
    # rating = [[]] * len(shells_high)
    # reason = [[]] * len(shells_high)
    # accepted = [None] * len(shells_high)
    rating = []
    reason = []
    for i, shell in enumerate(shells_high):
        rating.append([])
        reason.append([])
        reason_phrase =  " while using data in the shell " + \
            twodec(shells[i]) + "-" + twodec(shells[i + 1]) + " A"
        # If CC* is undefined or smaller than CCwork
        if os.path.isfile(csvfilename_merg):
            if CChalf_list[i] <= 0 or CChalf_list[i] == float("nan"):
                rating[i].append(11)
                reason[i].append("CC1/2 in high resolution is negative or "
                                 "undefined" + reason_phrase)
            else:
                if not args.complete_cross_validation:
                    if CCstar_list[i] < CCwork_shell_list[i] or \
                            CCstar_list[i] == float("nan"):
                        rating[i].append(10)
                        reason[i].append("CC* in high resolution is lower "
                                         "than CCwork" + reason_phrase)
        # If an R-value >= 0.45 or >= 40
        if not args.complete_cross_validation:
            # Rwork
            if Rwork_shell_list[i] >= 0.45 or \
                    Rwork_shell_list[i] == float("nan"):
                rating[i].append(9)
                reason[i].append("Rwork in high resolution is higher than "
                                 "0.45" + reason_phrase)
            elif Rwork_shell_list[i] >= 0.40:
                rating[i].append(6)
                reason[i].append("Rwork in high resolution is higher than "
                                 "0.40" + reason_phrase)
            # Rfree
            if Nfree_shell_list[i] < 50:
                warning_phrase = \
                    "here are only " + str(int(Nfree_shell_list[i])) + "" + \
                    " < 50 free reflections in the resolution shell " + \
                    "" + twodec(shells[i]) + "-" + \
                    "" + twodec(shells[i + 1]) + " A. Values of " + \
                    "statistics Rfree and CCfree in this shell " + \
                    "could be misleading. Consider setting " + \
                    "thicker resolution shells."
            if Rfree_shell_list[i] >= 0.45 or \
                    Rfree_shell_list[i] == float("nan"):
                if Nfree_shell_list[i] >= 50:
                    rating[i].append(8)
                    reason[i].append("Rfree in high resolution is higher than "
                                     "0.45" + reason_phrase)
                else:
                    rating[i].append(4)
                    reason[i].append(
                        "Rfree in high resolution is higher than "
                        "0.45" + reason_phrase + ". But t" + warning_phrase)
                    warning_my(
                        "lowNfree" + twodec(shell), "T" + warning_phrase)
            elif Rfree_shell_list[i] >= 0.40:
                if Nfree_shell_list[i] >= 50:
                    rating[i].append(5)
                    reason[i].append("Rfree in high resolution is higher than "
                                     "0.40" + reason_phrase)
                else:
                    rating[i].append(3)
                    reason[i].append(
                        "Rfree in high resolution is higher than "
                        "0.40" + reason_phrase + ". But t" + warning_phrase)
                    warning_my(
                        "lowNfree" + twodec(shell), "T" + warning_phrase)
        # Differences in overall R-values
        if Rwork_overall_list[i] > 0.01:
            rating[i].append(12)
            reason[i].append("Overall Rwork increased pretty much" + \
                             reason_phrase)
        elif Rfree_overall_list[i] <= 0.000009:
            rating[i].append(1)
            reason[i].append("Overall Rfree decreased" + reason_phrase)
        elif Rfree_overall_list[i] <= 0.000209 and Rwork_overall_list[i] > 0:
            rating[i].append(2)
            reason[i].append("Overall Rwork increased and Rfree remained " + \
                             "constant" + reason_phrase)
        else:
            rating[i].append(7)
            reason[i].append("Overall Rfree increased" + reason_phrase)
    # Suggest cutoff(s)
    cutoff = [shells[0], shells[0]]
    accepted = []  # order in list: [strict, benevolent] algorithm
    reason_phrase_bad_previous = "But statistics deteriorate in a previous resolution shell."
    for i, shell in enumerate(shells_high):
        accepted.append([None, None])
        # strict
        if i != 0 and not accepted[i - 1][0]:
            accepted[i][0] = False
            reason[i].insert(0, reason_phrase_bad_previous)
        elif max(rating[i]) < 5:
            accepted[i][0] = True
            cutoff[0] = shell
        elif max(rating[i]) >= 5:
            accepted[i][0] = False
        # benevolent
        if max(rating[i]) >= 7:
            accepted[i][1] = False
        elif max(rating[i]) < 7:
            if i >= 2 and not accepted[i - 1][1] and not accepted[i - 2][1]:
                accepted[i][1] = False
                reason[i].insert(0, reason_phrase_bad_previous)
            elif i != 0 and not accepted[i - 1][1] and max(rating[i - 1]) > 7:
                accepted[i][1] = False
                reason[i].insert(0, reason_phrase_bad_previous)
            elif i != 0 and not accepted[i - 1][1] and max(rating[i - 1]) == 7:
                # Analyse R-values at initial resolution
                # Pick them
                csvfilename = args.project + "_Rgap.csv"
                Rwork_overall_alt_list = []
                Rfree_overall_alt_list = []
                with open(csvfilename, "r") as csvfile:
                    for line in csvfile.readlines():
                        if line.lstrip()[0] == "#":  # If it is a comment
                            continue                 # do not load data
                        try:
                            Rwork_overall_alt_list.append(float(line.split()[1]))
                        except ValueError:
                            Rwork_overall_alt_list.append((float("nan")))
                        try:
                            Rfree_overall_alt_list.append(float(line.split()[2]))
                        except ValueError:
                            Rfree_overall_alt_list.append((float("nan")))
                # Be aware: Rwork_overall_alt_list and Rfree_overall_alt_list
                # has an extra element in the beggining - for the initial resol.!
                # (in comparison with other lists of statistics in this function)
                Rfree_overall_alt_diff = \
                    Rfree_overall_alt_list[i + 1] - Rfree_overall_alt_list[i - 1]
                Rwork_overall_alt_diff = \
                    Rwork_overall_alt_list[i + 1] - Rwork_overall_alt_list[i - 1]
                # Decide about the current and previous shell
                if Rfree_overall_alt_diff < 0:
                    accepted[i - 1][1] = True
                    reason[i - 1][-1] += ". But the next shell compensates it."
                    accepted[i][1] = True
                    cutoff[1] = shell
                elif Rfree_overall_alt_diff <= 0.000209 and Rwork_overall_alt_diff > 0:
                    accepted[i - 1][1] = True
                    reason[i - 1][-1] += ". But the next shell compensates it."
                    accepted[i][1] = True
                    cutoff[1] = shell
                else:
                    accepted[i][1] = False
                    reason[i].insert(0, reason_phrase_bad_previous)
            else:
                accepted[i][1] = True
                cutoff[1] = shell
    with open("PAIREF_cutoff.txt", "w") as f:
        f.write(twodec(cutoff[0]))
    # print(twodec(cutoff[0]), twodec(cutoff[1]))
    # for i, shell in enumerate(shells_high):
    #     print(twodec(shell) + "     " + str(accepted[i][0]) + "     " + \
    #           str(accepted[i][1]) + "     " + str(reason[i]))
    return(cutoff, accepted, reason)
