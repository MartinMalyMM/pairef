# coding: utf-8
from __future__ import print_function
import argparse
import sys
import os
import platform
import shutil
from .settings import warning_dict, settings
from .preparation import welcome, create_workdir, output_log, def_res_shells
from .preparation import which, res_high_from_xyzin, res_from_mtz, res_opt
from .preparation import calculate_merging_stats, run_baverage, run_pdbtools
from .preparation import res_from_hklin_unmerged, check_refinement_software
from .commons import twodec, twodecname, warning_my, try_symlink
from .refinement import collect_stat_OVERALL
from .refinement import collect_stat_OVERALL_AVG
from .refinement import collect_stat_BINNED
from .graphs import matplotlib_bar, matplotlib_line, write_log_html


RES_LOW = 50


class MyArgumentParser(argparse.ArgumentParser):
    """ Helper class for `argparse`

    It adds an attribute `add_argument_with_check` that check if file
    given in argument exist but not open it.
    Inspired by `https://codereview.stackexchange.com/questions/28608/
    checking-if-cli-arguments-are-valid-files-directories-in-python`

    Moreover, help message is displayed in the case of parser.error().
    """
    def __is_valid_file(self, parser, arg):
        """Checks if file
        given in argument `arg` exists but does not open it.
        If not, abort.

        Args:
            self
            parser: parser of `argparse`
            arg (str): argument of `argparse`

        Returns:
            str: Name of the checked file
        """
        if not os.path.isfile(arg):
            parser.error('The file {} does not exist!'.format(arg))
        else:
            # File exists so return the filename
            return arg

    def add_argument_with_check(self, *args, **kwargs):
        """New attribute for `argparse` that checks if file
        given in argument exist but does not open it.
        """
        # Look for your FILE settings
        # type = lambda x: self.__is_valid_file(self, x) # PEP8 E731
        def type(x):
            return self.__is_valid_file(self, x)
        kwargs['type'] = type
        self.add_argument(*args, **kwargs)

    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


# https://stackoverflow.com/questions/14117415/
# in-python-using-argparse-allow-only-positive-integers
def check_positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive "
                                         "int value" % value)
    return ivalue


def check_non_negative_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid non-negative "
                                         "int value" % value)
    return ivalue


def check_positive_float(value):
    ivalue = float(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive "
                                         "float value" % value)
    return ivalue


def process_arguments(input_args):
    '''Processes input arguments using `argparse`.

    Args:
        input_args (list): Input arguments

    Returns:
        list: Processed arguments
    '''
    # Input processing
    # parser = argparse.ArgumentParser(
    parser = MyArgumentParser(
        description="Automatic PAIRed REFinement protocol",
        epilog='Dependencies: CCP4 Software Suite or PHENIX containing CCTBX '
        'with Python 2.7',
        prog='ccp4-python -m pairef',
        add_help=False)

    # Arguments
    parser._optionals.title = 'optional arguments specifying input files'
    # This is a bit dirty hack but  the input files has to be added to `parser`
    # (not to a group) in order to add_argument_with_check() could work
    parser.add_argument(
        '--GUI', '--gui', dest='gui', action='store_true',
        help='Start graphical user interface (usually requires '
        "to be executed as ccp4-python, not as cctbx.python)")
    ###########################   GUI   ######################################
    if "--GUI" in input_args or "--gui" in input_args:   
        args, args_unknown = parser.parse_known_args(input_args)
        return args
    ###########################   GUI   ######################################
        
    parser.add_argument_with_check(
        '--XYZIN', '--xyzin', dest='xyzin',
        help='PDB or mmCIF file with current structure model', required=True)
    parser.add_argument_with_check(
        '--HKLIN', '--hklin', dest='hklin',
        help='MTZ file with processed diffraction data', required=True)
    parser.add_argument_with_check(
        '-u', "--unmerged", dest='hklin_unmerged',
        help='unmerged processed diffraction data file (e.g. XDS_ASCII.HKL '
        'or data_unmerged.mtz)')
    parser.add_argument_with_check(
        '--LIBIN', '--libin', dest='libin',
        help='CIF file geometric restraints')
    parser.add_argument_with_check(
        '--TLSIN', '--tlsin', dest='tlsin',
        help='input TLS file (only for REFMAC5)')
    parser.add_argument_with_check(
        '-c', "--comfile", dest='comin',  # TODO - better...
        help='configuration Com file with keywords for REFMAC5')
    parser.add_argument_with_check(
        '-d', "--def", dest='defin',
        help='configuration def file with keywords for phenix.refine')

    group1 = parser.add_mutually_exclusive_group()
    #    'optional arguments specifying refinement software')
    group1.add_argument(
        '-R', '--refmac', dest='refmac',
        help='Use REFMAC5 (default)', action='store_true')
    group1.add_argument(
        '-P', '--phenix', dest='phenix',
        help='Use phenix.refine', action='store_true')

    group2 = parser.add_argument_group('other optional arguments')
    group2.add_argument(
        '-p', "--project", dest='project', help='project name')
    group2.add_argument(
        '-r', dest='res_shells',
        help='explicit definition of high resolution shells - '
        'values must be divided using commas without any '
        'spaces and written in decreasing order, e.g. '
        '2.1,2.0,1.9')
    group2.add_argument(
        '-n', dest='n_shells',
        help='number of high resolution shells to be added step by step. '
        'Using this argument, setting of argument -s is required.',
        type=check_positive_int)
    group2.add_argument(
        '-s', "--step", dest='step',
        help='width of the added high resolution shells (in angstrom). '
        'Using this argument, setting of argument -n is required.',
        type=float)
    group2.add_argument(
        '-i', dest='res_init',
        help='initial high-resolution diffraction limit (in angstrom) '
        '- if it is not necessary, '
        'do not use this option, the script should find resolution '
        'automatically in PDB or mmCIF file',
        type=float)
    group2.add_argument(
        '-f', "--flag", dest='flag',
        help="definition which FreeRflag set will be excluded during "
        "refinement (set 0 default)", type=check_non_negative_int)
    group2.add_argument(
        '-w', "--weight", dest='weight',
        help="manual definition of weighting term (only for REFMAC5)",
        type=float)
    group2.add_argument(
        "--ncyc", dest='ncyc',
        help="number of refinement cycles that will be performed in every "
        "resolution step", type=check_positive_int)
    group2.add_argument(
        '--constant-grid', dest='constant_grid',
        help="keep the same FFT grid through the whole paired refinement "
        "(only for REFMAC5)", action='store_true')
    group2.add_argument(
        '--complete', dest='complete_cross_validation',
        help="perform complete cross-validation (use all available free "
        "reflection sets)", action='store_true')
    group2.add_argument(
        "--TLS-ncyc", "--tls-ncyc", dest='tls_ncyc',
        help="number of cycles of TLS refinement (10 cycles by default, "
        "only for REFMAC5)",
        type=check_positive_int)
    group2.add_argument(
        '--TLSIN-keep', "--tlsin-keep", dest='tlsin_keep',
        help="keep using the same TLS input file in all the refinement runs "
        "(only for REFMAC5)", action='store_true')
    group2.add_argument(
        "--open-browser", action="store_true", dest='open_browser',
        help="open web browser to show results "
        "(requires to be executed as ccp4-python, not as cctbx.python)")
    group2.add_argument(
        "-h", "--help", action="help", help="show this help message and exit")

    group3 = parser.add_argument_group(
        'optional arguments specifying structure model modification')
    group3.add_argument(
        "--prerefinement-ncyc", dest='prerefinement_ncyc',
        help="number of refinement cycles to be performed as pre-refinement "
        "of the input structure model before paired "
        "refinement (the initial high resolution limit is used). Pre-refine"
        "ment is performed by default in case of the complete cross-"
        "validation protocol. Other related options are "
        "--prerefinement-reset-bfactor, --prerefinement-add-to-bfactor, "
        "--prerefinement-set-bfactor, --prerefinement-shake-sites, and "
        "--prerefinement-no-modification. "
        "These options can be useful when the structure has "
        "been refined in another version of REFMAC5 or phenix.refine "
        "than it is currently used or when you want to reset the impact of "
        "used free reflections.",
        type=check_positive_int)
    group3.add_argument(
        "--prerefinement-reset-bfactor", dest='reset_bfactor',
        help="reset atomic B-factors of the input structure model to "
        "the mean value. This is done by default in the case of the complete"
        "cross-validation protocol.",
        # "This argument requires the argument --prerefinement_ncyc.",
        action='store_true')
    group3.add_argument(
        '--prerefinement-add-to-bfactor', dest='add_to_bfactor',
        help="add the given value to B-factors"
        " of the input structure model", type=float)
    group3.add_argument(
        "--prerefinement-set-bfactor", dest='set_bfactor',
        help="set atomic B-factors of the input structure model "
        "to the given value.",
        # "This argument requires the argument --prerefinement_ncyc.",
        type=check_positive_float)
    group3.add_argument(
        "--prerefinement-shake-sites", dest='shake_sites',
        help="randomize coordinates of the input structure model with the "
        "given mean error value. This is done by default in the case of the "
        "complete cross-validation protocol - mean error 0.25.",
        # "This argument requires the argument --prerefinement_ncyc.",
        type=check_positive_float,
        nargs='?', const=0.25)
    # nargs='?' means 0-or-1 arguments
    # const=0.25 sets the default when there are 0 arguments
    group3.add_argument(
        '--prerefinement-no-modification', dest='no_modification',
        help="do not modify the input structure model before the complete "
        "cross-validation protocol", action='store_true')

    parser.add_argument(
        '-t', dest='test', help=argparse.SUPPRESS, action='store_true')
    parser.add_argument(
        '-q', dest='quick', help=argparse.SUPPRESS, action='store_true')

    # If the arguments are not set correctly, show help and exit
    if (input_args == [__file__]) or \
       (input_args == [os.path.basename(__file__)]) or \
       (len(input_args) == 0):
        parser.print_help()  # help on empty input
        sys.exit(1)
    # Passing list of params to parse_args without script name
    if (__file__ in input_args) or \
       (os.path.basename(__file__) in input_args) or \
       (input_args[0].split("/")[-1] == "__main__.py") or \
       (input_args[0].split("/")[-1] == "pairef"):
        input_args = input_args[1:]
    args, args_unknown = parser.parse_known_args(input_args)
    # CCP4Console on Windows adds an unknown argument (abs path)\__main__.py
    for arg_unknown in args_unknown:
        if "__main__.py" not in arg_unknown:
            print("WARNING: Unknown argument: " + arg_unknown)
    # Management of conflicts
    #
    # if ((args.add_to_bfactor and
    #         not (args.complete_cross_validation or
    #              args.reset_bfactor))):
    #     parser.error("The argument --prerefinement-add-to-bfactor requires the"
    #                  " argument --complete or --prerefinement-reset-bfactor.")
    if args.complete_cross_validation and isinstance(args.flag, (int, long)):
        parser.error("It is a non-sense to use the option -f with the option "
                     "--complete.")
    if (args.no_modification and
            (args.reset_bfactor or args.add_to_bfactor or
             args.set_bfactor or args.shake_sites)):
        parser.error("It is a non-sense to use the option "
                     "--prerefinement-no-modification "
                     "with the options --prerefinement-reset-bfactor, "
                     "--prerefinement-add-to-bfactor, "
                     "--prerefinement-set-bfactor, or "
                     "--prerefinement-shake-sites.")
    if ((not args.prerefinement_ncyc and
            not args.complete_cross_validation) and
            (args.reset_bfactor or args.add_to_bfactor or
             args.set_bfactor or args.shake_sites)):
        parser.error("It is not supported to modify the input structure model "
                     "(arguments --prerefinement-reset-bfactor, "
                     "--prerefinement-add-to-bfactor, "
                     "--prerefinement-set-bfactor, or "
                     "--prerefinement-shake-sites) and not perform "
                     "pre-refinement at the initial resolution - use the "
                     "required argument --prerefinement-ncyc")
    if (args.no_modification and not args.complete_cross_validation):
        print("NOTE: The option --prerefinement-no-modification is "
              "redundant when the option --complete is not used.")
    if (args.set_bfactor and (args.reset_bfactor or args.add_to_bfactor)):
        parser.error("The option --prerefinement-set-bfactor cannot be "
                     "combined with the options --prerefinement-reset-bfactor "
                     "and --prerefinement-add-to-bfactor.")
    if args.complete_cross_validation and args.reset_bfactor:
        print("NOTE: The argument --prerefinement-reset-bfactor is "
              "redundant when the argument --complete is used - the B-factors "
              "are reseted to their average value by default.")
    if args.complete_cross_validation and args.no_modification:
        warning_my("no_modification", "Modification of the input structure "
                   "model is turned off - be sure that the input structure "
                   "model has been pertubed. If not, the results will be "
                   "biased.")
    # if ((args.reset_bfactor and
    #      not args.prerefinement_ncyc)):
    #     parser.error("The argument --prerefinement-reset-bfactor requires "
    #                  "the argument --prerefinement-ncyc.")
    if (args.step and not args.n_shells) or \
            (not args.step and args.n_shells):  # TODO
        parser.error("Both - the number and the width - of high resolution "
                     "shells (options -n and -s) "
                     "must be set in one time (not only one of "
                     "these options).")
        args.reset_bfactor = False
    # Input structure model modification - default behaviour
    if args.complete_cross_validation and not args.no_modification:
        if not (args.reset_bfactor or args.add_to_bfactor or
                args.set_bfactor or args.shake_sites):
            args.reset_bfactor = True
            args.shake_sites = 0.25
    # REFMAC5   X   phenix.refine
    if args.phenix and args.constant_grid:
        parser.error("Option --constant-grid is currently supported only for"
                     "REFMAC5, not for phenix.refine.")
    if args.phenix and args.weight:
        parser.error("Weight can be directly specifyied only if "
                     "phenix.refine is set as refinement software. Specify "
                     "keywords in a def file and use an option --def to set "
                     "parameters of phenix.refine.")
    if args.phenix and args.comin:
        parser.error("It is not possible to provide Com file for REFMAC5 when "
                     "phenix.refine is set as refinement software. Use an "
                     "option --def instead to specify keywords for "
                     "phenix.refine  or set refinement in REFMAC5 using an "
                     "option --refmac.")
    if not args.phenix and args.defin:
        parser.error("It is not possible to provide def file for phenix."
                     "refine when REFMAC5 is set as refinement software. Use "
                     "an option --comin instead to specify keywords for "
                     "REFMAC5 or set refinement in phenix.refine using an "
                     "option --phenix.")
    if args.phenix and (args.tlsin or args.tls_ncyc or args.tlsin_keep):
        parser.error("Specific options for TLS refinement are valid only for "
                     "REFMAC5. For phenix.refine, specify a refinement "
                     "strategy and TLS groups (keywords "
                     "refinement.refine.strategy and refinement.refine.adp) "
                     "in a configuration file (option --def).")
    # TLS
    if (args.tlsin_keep or args.tls_ncyc) and not args.tlsin:
        parser.error("Input TLS file must be specified (option --TLSIN) while "
                     "using the options --TLSIN-keep or --TLS-ncyc.")
    return(args)


def main(args):
    """The main function of the `pairef` module.

    Args:
        args: Input arguments processed by `argparse`
    """
    # Check software versions (matplotlib should be checked later)
    # if int(platform.python_version_tuple()[0]) != 2:
    #     sys.stderr.write("ERROR: This version of pairef module requires "
    #                      "Python 2.7 from CCTBX.\n")
    #     sys.exit(1)
    try:
        import cctbx.miller
    except ImportError:
        sys.stderr.write("ERROR: This version of pairef module requires "
                         "Python 2.7 from CCTBX.\n"
                         "It has to be executed using command:"
                         "ccp4-python -m pairef ARGUMENTS\n"
                         "or\n"
                         "cctbx.python -m pairef ARGUMENTS\n"
                         "Aborting.\n")
        sys.exit(1)

    # Decide which refinement software will be used
    if args.phenix:
        refinement = "phenix"
        refinement_name = "phenix.refine"
        from .refinement import refinement_phenix
        from .refinement import collect_stat_binned_phenix_low
    else:
        refinement = "refmac"  # REFMAC5 as default
        refinement_name = "REFMAC5"
        from .refinement import refinement_refmac
        from .refinement import collect_stat_binned_refmac_low

    # Check of the needed executables - works only on Linux
    if platform.system() == 'Linux':
        if refinement == "refmac":  # CCP4 & REFMAC5
            cryst_package = "CCP4 Software Suite"
            required_executables = ["refmac5", "baverage", "mtzdump", "sfcheck"]
            ## if args.update_waters:
            ##     required_executables.append("findwaters")
        elif refinement == "phenix":
            cryst_package = "PHENIX software suite"
            required_executables = ["phenix.refine"]
            # modules: from xia2.Wrappers.CCP4.Mtzdump import Mtzdump
            #          iotbx
        for required_executable in required_executables:
            if not which(required_executable) and not args.test:
                sys.stderr.write("ERROR: PAIREF requires installed `"
                                 "" + required_executable + "` (a part "
                                 "of the " + cryst_package + ") but it is not "
                                 "executable."
                                 "\nAborting.\n")
                sys.exit(1)

    # If a project name is not set, assign something
    if not args.project:
        args.project = "project"

    # Decide suffix for structure model (PDB X mmCIF) based on args.xyzin
    xyzin_suffix = args.xyzin.split(".")[-1]
    if xyzin_suffix == "mmcif" or xyzin_suffix == "cif":
        if refinement == "phenix":
            settings["pdbORmmcif"] = ".cif"
        else:  # refmac
            settings["pdbORmmcif"] = ".mmcif"
    else:
        settings["pdbORmmcif"] = ".pdb"

    # Create new working directory (name related to the project)
    workdir = create_workdir(args.project)

    # Set to write STDOUT to screen and file
    writer = output_log(sys.stdout, workdir + '/PAIREF_out.log')
    sys.stdout = writer

    versions_dict = {"refmac_version": "N/A",  # It will be found later
                     "phenix_version": "N/A",  # It will be found later
                     "pairef_version": "1.2.1"}

    # Show information about the module and input parameters
    welcome(args, versions_dict["pairef_version"])

    # Find resolution range of merged data
    res_low, res_high_mtz = res_from_mtz(args.hklin)  # uses CCTBX
    if not res_high_mtz:
        sys.stderr.write("ERROR: High resolution limit of data "
                         "" + args.hklin + " could not be found."
                         "\nAborting.\n")
    if not res_low:
        warning_my("low_res", "Low resolution limit could not be found."
                   "Setting it to a value " + RES_LOW + " A.")
        res_low = RES_LOW
    elif res_low > 999:
        warning_my("low_res", "Founded low resolution limit is too large "
                   "(" + str(res_low) + " A). Setting it to a value "
                   "" + str(RES_LOW) + "A.")
        res_low = RES_LOW
    else:
        print("Resolution of the merged diffraction data "
              "" + args.hklin + ": " + twodec(res_low) + "-"
              "" + twodec(res_high_mtz) + " A")

    if args.hklin_unmerged:
        # Find resolution range of unmerged data
        res_low_from_hklin_unmerged, res_high_from_hklin_unmerged = \
            res_from_hklin_unmerged(args.hklin_unmerged)
        if res_high_from_hklin_unmerged != 0 and \
                res_low_from_hklin_unmerged != float("inf"):
            print("Resolution of the unmerged diffraction data "
                  "" + args.hklin_unmerged + ": "
                  "" + twodec(res_low_from_hklin_unmerged) + "-"
                  "" + twodec(res_high_from_hklin_unmerged) + " A")

    # Set initial high resolution limit
    if args.res_init:
        print("Manual setting of initial high resolution limit will be "
              "used: " + twodec(args.res_init) + " A.")
    else:
        args.res_init = res_high_from_xyzin(
            args.xyzin, format=settings["pdbORmmcif"])
        if args.res_init < 0:
            sys.stderr.write(
                "ERROR: An attempt to determine a resolution of data which "
                "were used for refinement of the structure model "
                "" + args.xyzin + " was not successful. Please specify the "
                "initial high resolution limit manually using -i option."
                "\nAborting.\n")
            sys.exit(1)
        print("Initial high resolution limit found in the structure model "
              "" + args.xyzin + ": " + twodec(args.res_init) + " A.")

    # Check that resolution shell setting has sence
    # and determine resolution shells
    shells, n_bins_low, n_flag_sets, default_shells_definition = \
        def_res_shells(args, refinement, res_high_mtz, res_low)
    print("High resolution diffraction limits:", end=" ")
    for shell in shells[1:-1]:  # Skip the initial high resolution limit
        print(twodec(shell) + " A", end=", ")
    print(twodec(shells[-1]) + " A")  # Formatting issue

    # Set FreeRflag sets
    if args.complete_cross_validation:
        if n_flag_sets <= 2:
            sys.stderr.write(
                "Given input MTZ file " + args.hklin + " has too low number "
                "of free reflection sets (" + str(n_flag_sets) + "). k-fold "
                "cross-validation cannot be performed.\n"
                "Aborting.\n")
            sys.exit(1)
        flag_sets = range(n_flag_sets)
        if args.quick:  # Faster testing
            flag_sets = range(3)
    else:
        # Python 2.x requires the bracket (int, long)
        # https://stackoverflow.com/questions/3501382/
        # checking-whether-a-variable-is-an-integer-or-not
        if not isinstance(args.flag, (int, long)):
            args.flag = 0
        flag_sets = [args.flag]
        print(" * Data with FreeRflag set " + str(args.flag) + " will be "
              "excluded during refinement.")

    # Check the input files?

    # Copy input files to the working directory
    # and take only basename of the filenames
    in_files = ["hklin", "xyzin"]
    in_files_optional = ["libin", "comin", "defin", "tlsin"]
    for f in in_files_optional:
        if vars(args)[f]:
            in_files.append(f)
    for f in in_files:
        shutil.copy2(vars(args)[f], workdir)
        vars(args)[f] = os.path.basename(vars(args)[f])
    # Symlink HKLIN_unmerged
    if args.hklin_unmerged:
        try_symlink(os.path.join(os.getcwd(), args.hklin_unmerged),
                    os.path.join(os.getcwd(), workdir,
                                 os.path.basename(args.hklin_unmerged)))
        args.hklin_unmerged = os.path.basename(args.hklin_unmerged)
    print("")
    # Change the working directory
    os.chdir(workdir)
    print("Current working directory: " + os.getcwd())

    write_log_html(shells, [], args, versions_dict, flag_sets)
    htmlfilepath = os.path.abspath("PAIREF_" + args.project + ".html")
    print("------> RESULTS AND THE CURRENT STATUS OF CALCULATIONS ARE LISTED "
          "IN A HTML LOG FILE "
          "" + htmlfilepath)
    
    if args.open_browser and "ccp4" in sys.executable:  # cctbx.python fails
        import webbrowser
        print("Opening web browser...")
        webbrowser.open(htmlfilepath)
    print("")

    # Modification of the input structure model - Define starting XYZIN
    if args.no_modification:
        xyzin_start = args.xyzin
    else:
        if args.complete_cross_validation or args.reset_bfactor:
            baverage = run_baverage(args.project, args.xyzin, args.res_init)
            xyzin_start = run_pdbtools(args, baverage)
        else:
            xyzin_start = run_pdbtools(args)
    shutil.copy2(args.xyzin,
                 args.project + "_" + twodecname(shells[0]) + "A" + \
                 settings["pdbORmmcif"])
    if args.test:
        sys.exit(0)

    print("\nRefinement using " + refinement_name + ":\n")
    res_cur = shells[0]
    if args.complete_cross_validation or args.prerefinement_ncyc:
        print("   * Performing pre-refinement at "
              "" + twodec(res_cur) + " A resolution...")
    else:
        print("   * Calculating initial statistics at "
              "" + twodec(res_cur) + " A resolution...")
    for flag in flag_sets:
        if refinement == "refmac":
            results = refinement_refmac(res_cur=res_cur,
                                        res_prev=args.xyzin,
                                        res_high=shells[0],
                                        args=args,
                                        n_bins_low=n_bins_low,
                                        mode="first",
                                        res_low=res_low,
                                        res_highest=shells[-1],
                                        flag=flag,
                                        xyzin_start=xyzin_start)
            # bfac_set=bfac_set)
            versions_dict["refmac_version"] = results["version"]
        elif refinement == "phenix":
            n_bins = n_bins_low
            results = refinement_phenix(res_cur=res_cur,
                                        res_prev=args.xyzin,
                                        res_high=shells[0],
                                        args=args,
                                        n_bins=n_bins_low,
                                        mode="first",
                                        res_low=res_low,
                                        res_highest=shells[-1],
                                        flag=flag,
                                        xyzin_start=xyzin_start)
            versions_dict["phenix_version"] = results["version"]
        collect_stat_OVERALL([res_cur], args, flag, refinement)
        if args.complete_cross_validation or args.prerefinement_ncyc:
            matplotlib_line(
                shells=[shells[0]],
                project=args.project,
                statistics=["Rwork_cyc", "Rfree_cyc"],
                n_bins_low=n_bins_low,
                title=r"$\mathrm{" + twodec(shells[0]) + r"\ \AA\ -" +
                "\ flag\ " + str(flag) + "}$",
                filename_suffix="R" + str(flag).zfill(2) + "_" +
                twodecname(shells[0]) + "A_stats_vs_cycle", flag=flag,
                refinement=refinement)
            write_log_html(shells, [], args,
                           versions_dict, flag_sets, shells[0])

    if not args.complete_cross_validation:
        src = args.project + "_R" + str(flag).zfill(2) + "_Rgap.csv"
        dst = args.project + "_Rgap.csv"
        try_symlink(src, dst)

    matplotlib_line(shells=[shells[0]],
                    project=args.project,
                    statistics=["Rgap"],
                    n_bins_low=n_bins_low,
                    title=r"$\it{R}_{\mathrm{free}}-"
                    r"\it{R}_{\mathrm{work}}$",
                    filename_suffix="Rgap", flag=flag)
    shells_ready_with_res_init = [shells[0]]
    write_log_html(shells, shells_ready_with_res_init, args,
                   versions_dict, flag_sets)

    # Check which software (and which version) has been used
    # for refinement of the file XYZIN
    if not (args.complete_cross_validation or \
            args.prerefinement_ncyc):
        check_refinement_software(args, versions_dict, refinement)

    if refinement == "refmac":
        logfilename = args.project + "_R" + str(flag_sets[0]).zfill(2) + "_" \
            "" + twodecname(shells[0]) + "A_comparison" \
            "_at_" + twodecname(shells[0]) + "A.log"
        mtzfilename = \
            args.project + "_R" + str(flag_sets[0]).zfill(2) + "_" \
            "" + twodecname(shells[0]) + "A.mtz"
        bins_low = collect_stat_binned_refmac_low(
            logfilename, mtzfilename, args.hklin, n_bins_low, res_low, flag)[1]
    elif refinement == "phenix":
        pdbfilename = args.project + "_R" + str(flag_sets[0]).zfill(2) + "_" \
            "" + twodecname(shells[0]) + "A_comparison" \
            "_at_" + twodecname(shells[0]) + "A_001.pdb"
        bins_low = collect_stat_binned_phenix_low(pdbfilename, n_bins_low)[0]
    bins_low = [float(bin) for bin in bins_low]
    if not args.complete_cross_validation:
        collect_stat_BINNED([res_cur], args.project, args.hklin,
                                   n_bins_low, flag, res_low, refinement)
        if which("sfcheck"):
            res_opt(shells[0], args, refinement)
        matplotlib_line(shells=[shells[0]],
                        project=args.project,
                        statistics=["res_opt"],
                        n_bins_low=n_bins_low,
                        title="Optical resolution",
                        filename_suffix="Optical_resolution", flag=flag)
        matplotlib_line(shells=[shells[0]],
                        project=args.project,
                        statistics=["Rwork"],
                        n_bins_low=n_bins_low,
                        title=r"$\it{R}_{\mathrm{work}}$",
                        filename_suffix="Rwork", flag=flag)
        matplotlib_line(shells=[shells[0]],
                        project=args.project,
                        statistics=["Rfree"],
                        n_bins_low=n_bins_low,
                        title=r"$\it{R}_{\mathrm{free}}$",
                        filename_suffix="Rfree", flag=flag)
        matplotlib_line(shells=[shells[0]],
                        project=args.project,
                        statistics=["CCwork", "CC*"],
                        n_bins_low=n_bins_low,
                        title=r"CC$_\mathrm{work}$",
                        filename_suffix="CCwork", flag=flag)
        matplotlib_line(shells=[shells[0]],
                        project=args.project,
                        statistics=["CCfree", "CC*"],
                        n_bins_low=n_bins_low,
                        title=r"CC$_\mathrm{free}$",
                        filename_suffix="CCfree", flag=flag)
        matplotlib_line(shells=[shells[0]],
                        project=args.project,
                        statistics=["n_work", "n_free"],
                        n_bins_low=n_bins_low,
                        title="Number of reflections in resol. bins",
                        filename_suffix="No_work_free_reflections",
                        flag=flag, multiscale=True)
    write_log_html(shells, shells_ready_with_res_init, args,
                   versions_dict, flag_sets)

    for i in range(len(shells) - 1):
        # TODO: check files
        res_cur = shells[i + 1]
        res_prev = shells[i]

        # Real refinement
        print("\n   * Refining using data up to "
              "" + twodec(shells[i + 1]) + " A resolution...")
        for flag in flag_sets:
            if refinement == "refmac":
                results = refinement_refmac(res_cur=res_cur,
                                            res_prev=res_prev,
                                            res_high=shells[i + 1],
                                            args=args,
                                            n_bins_low=n_bins_low,
                                            mode="refine",
                                            res_low=res_low,
                                            res_highest=shells[-1],
                                            flag=flag)
            elif refinement == "phenix":
                n_bins += 1
                results = refinement_phenix(res_cur=res_cur,
                                            res_prev=res_prev,
                                            res_high=shells[i + 1],
                                            args=args,
                                            n_bins=n_bins,
                                            mode="refine",
                                            res_low=res_low,
                                            res_highest=shells[-1],
                                            flag=flag)
            matplotlib_line(
                shells=[res_cur],
                project=args.project,
                statistics=["Rwork_cyc", "Rfree_cyc"],
                n_bins_low=n_bins_low,
                title=r"$\mathrm{" + twodec(res_cur) + r"\ \AA\ -" +
                "\ flag\ " + str(flag) + "}$",
                filename_suffix="R" + str(flag).zfill(2) + "_" +
                twodecname(res_cur) +
                "A_stats_vs_cycle", flag=flag,
                refinement=refinement)
            write_log_html(shells, shells_ready_with_res_init, args,
                           versions_dict, flag_sets, res_cur)
            print("       Calculating statistics of the refined structure "
                  "model...", end="")
            if refinement == "refmac":
                # Statistics up to prev. res. limit
                results = refinement_refmac(res_cur=res_cur,
                                            res_prev=res_prev,
                                            res_high=shells[i],
                                            args=args,
                                            n_bins_low=n_bins_low,
                                            mode="prev_pair",
                                            res_low=res_low,
                                            res_highest=shells[-1],
                                            flag=flag)
                # Statistics for `n_bins_low` shells up to init. res. limit
                results = refinement_refmac(res_cur=res_cur,
                                            res_prev=res_prev,
                                            res_high=shells[0],
                                            args=args,
                                            n_bins_low=n_bins_low,
                                            mode="comp",
                                            res_low=res_low,
                                            res_highest=shells[-1],
                                            flag=flag)
            elif refinement == "phenix":
                # Statistics up to prev. res. limit
                results = refinement_phenix(res_cur=res_cur,
                                            res_prev=res_prev,
                                            res_high=shells[i],
                                            args=args,
                                            n_bins=n_bins-1,
                                            mode="prev_pair",
                                            res_low=res_low,
                                            res_highest=shells[-1],
                                            flag=flag)
                # Statistics for `n_bins_low` shells up to init. res. limit
                results = refinement_phenix(res_cur=res_cur,
                                            res_prev=res_prev,
                                            res_high=shells[0],
                                            args=args,
                                            n_bins=n_bins_low,
                                            mode="comp",
                                            res_low=res_low,
                                            res_highest=shells[-1],
                                            flag=flag)
            collect_stat_OVERALL(shells[:i + 2], args, flag, refinement)
            if not args.complete_cross_validation:
                # Update csv files
                symlinks_src = [
                    args.project + "_R" + str(flag).zfill(2) +
                    "_R-values.csv",
                    args.project + "_R" + str(flag).zfill(2) + "_Rgap.csv"
                    ]
                symlinks_dst = [
                    args.project + "_R-values.csv",
                    args.project + "_Rgap.csv"
                    ]
                for src, dst in zip(symlinks_src, symlinks_dst):
                    try_symlink(src, dst)
                # Optical resolution
                if which("sfcheck"):
                    res_opt(res_cur, args, refinement)
                # Statistics for high resolution shells
                n_high_resolution_shells_ready = i + 1
                for j in range(n_high_resolution_shells_ready):
                    if refinement == "refmac":
                        results = refinement_refmac(res_cur=shells[i + 1],
                                                    res_prev=shells[i],
                                                    res_high=shells[j + 1],
                                                    args=args,
                                                    n_bins_low=n_bins_low,
                                                    mode="comp",
                                                    res_low=shells[j],
                                                    res_highest=shells[-1],
                                                    flag=flag)
                    elif refinement == "phenix":
                        results = refinement_phenix(res_cur=shells[i + 1],
                                                    res_prev=shells[i],
                                                    res_high=shells[j + 1],
                                                    args=args,
                                                    n_bins=1,
                                                    mode="comp",
                                                    res_low=shells[j],
                                                    res_highest=shells[-1],
                                                    flag=flag)
        shells_ready_with_res_init = shells[:i + 2]
        print("")
        if args.complete_cross_validation:
            collect_stat_OVERALL_AVG(shells_ready_with_res_init,
                                            args.project, flag_sets)
        else:
            collect_stat_BINNED(
                shells_ready_with_res_init, args.project, args.hklin,
                n_bins_low, flag, res_low, refinement)

        print("       Updating graphs...")
        matplotlib_bar(args)
        if args.complete_cross_validation:
            matplotlib_bar(args=args, flag_sets=flag_sets,
                           ready_shells=shells_ready_with_res_init)
        else:
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["res_opt"],
                            n_bins_low=n_bins_low,
                            title="Optical resolution",
                            filename_suffix="Optical_resolution",
                            flag=flag)
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["Rwork"],
                            n_bins_low=n_bins_low,
                            title=r"$\it{R}_{\mathrm{work}}$",
                            filename_suffix="Rwork", flag=flag)
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["Rfree"],
                            n_bins_low=n_bins_low,
                            title=r"$\it{R}_{\mathrm{free}}$",
                            filename_suffix="Rfree", flag=flag)
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["CCwork", "CC*"],
                            n_bins_low=n_bins_low,
                            title=r"CC$_\mathrm{work}$",
                            filename_suffix="CCwork", flag=flag)
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["CCfree", "CC*"],
                            n_bins_low=n_bins_low,
                            title=r"CC$_\mathrm{free}$",
                            filename_suffix="CCfree", flag=flag)
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["n_work", "n_free"],
                            n_bins_low=n_bins_low,
                            title="Number of reflections in resol. bins",
                            filename_suffix="No_work_free_reflections",
                            flag=flag, multiscale=True)
        matplotlib_line(shells=[shells[0]],  # ???
                        project=args.project,
                        statistics=["Rgap"],
                        n_bins_low=n_bins_low,
                        title=r"$\it{R}_{\mathrm{free}}-"
                        r"\it{R}_{\mathrm{work}}$",
                        filename_suffix="Rgap", flag=flag)
        write_log_html(shells, shells_ready_with_res_init, args,
                       versions_dict, flag_sets)

    # If unmerged data are in disposal, calculate CC1/2 and CC*
    # for future graphs of CCwork, CCfree
    if args.hklin_unmerged:
        write_log_html(shells, shells, args, versions_dict, flag_sets)
        calculate_merging_stats(args.hklin_unmerged, shells, args.project,
                                bins_low, res_low_from_hklin_unmerged,
                                res_high_from_hklin_unmerged)
        matplotlib_line(shells=[shells[0]], project=args.project,
                        statistics=["Rmerge", "Rmeas", "Rpim"],
                        n_bins_low=n_bins_low, title="$\it{R}$-values",
                        filename_suffix="Rmerge_Rmeas_Rpim")
        matplotlib_line(shells=[shells[0]], project=args.project,
                        statistics=["<I/sI>", "<I>"],
                        n_bins_low=n_bins_low, title="Average intensities",
                        filename_suffix="Intensities", multiscale=True)
        matplotlib_line(shells=[shells[0]], project=args.project,
                        statistics=["Completeness", "Multiplicity"],
                        n_bins_low=n_bins_low,
                        title="Completeness and multiplicity",
                        filename_suffix="Comp_Mult", multiscale=True)
        matplotlib_line(shells=[shells[0]], project=args.project,
                        statistics=["CChalf", "CC*"], n_bins_low=n_bins_low,
                        title="Correlation coefficent", filename_suffix="CC")
        matplotlib_line(shells=[shells[0]], project=args.project,
                        statistics=["n_unique", "n_obs"],
                        n_bins_low=n_bins_low,
                        title="Number of reflections in resol. bins",
                        filename_suffix="No_reflections",
                        multiscale=True)

        if not args.complete_cross_validation:
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["CCwork", "CC*"],
                            n_bins_low=n_bins_low,
                            title=r"CC$_\mathrm{work}$",
                            filename_suffix="CCwork")
            matplotlib_line(shells=shells_ready_with_res_init,
                            project=args.project,
                            statistics=["CCfree", "CC*"],
                            n_bins_low=n_bins_low,
                            title=r"CC$_\mathrm{free}$",
                            filename_suffix="CCfree")
        write_log_html(shells, shells, args, versions_dict, flag_sets,
                       ready_merging_statistics=True, done=True)
    else:
        write_log_html(shells, shells, args, versions_dict, flag_sets,
                       done=True)

    if warning_dict:
        print("\nCalculation ended.")
        print("These warning messages appeared during calculation:")
        for key in warning_dict:
            print(warning_dict[key])
    else:
        print("\nCalculation ended successfully.")
    print("\nResults are listed "
          "in logfile " + os.getcwd() + "/PAIREF_" + args.project + ".html\n")
    return


def run_pairef(input_args=None):
    """**THE LAUNCHING FUNCTION** - to process input arguments
    using :func:`launcher.process_arguments` and launch
    the :func:`launcher.main` function.

    Args:
        input_args (list): List of input arguments - parameters
    """
    # Load init setting
    # init()
    # Proccess arguments
    # (Option for specifying arguments via `input_args` is due to testing)
    if not input_args:
        input_args = sys.argv
    args = process_arguments(input_args)

    if args.gui:
        # Run GUI in PyQt
        from .gui import gui
        gui()
    else:
        # Run the protocol
        main(args)
    return
 
