# coding: utf-8
from __future__ import print_function
from __future__ import division
import os
import sys
import re
import subprocess
import shutil
from math import sqrt
from .settings import warning_dict, settings
from .commons import twodec, twodecname, fourdec, extract_from_file, warning_my


def refinement_refmac(res_cur,
                      res_prev,
                      res_high,
                      args,
                      n_bins_low,
                      mode="refine",
                      res_low=0,
                      res_highest=0,
                      flag=0,
                      xyzin_start="",
                      bfac_set=0):  # TODO: hklfree
    """Refine using REFMAC5.

    Resolution limit is controlled by parameters `res_high` and `res_low`.
    Names of files are controlled by parameters `res_cur` and `res_high`.

    There are 4 different modes:

    * `mode="refine"`

          refine structure `args.project_twodecname(res_prev)A.pdb`;
          results with prefix `args.project_twodecname(res_cur)A`;
          `ncyc` given by `comfile` (default: `ncyc 10`)

    * `mode="comp"`

          use structure `args.project_twodecname(res_cur)A.pdb`;
          results with prefix `args.project_twodecname(res_cur)
          A_comparison_at_twodecname(res_high)A`;
          `ncyc 0`

    * `mode="prev_pair"`

          use structure `args.project_twodecname(res_cur)A.pdb`;
          results with prefix `args.project_twodecname(res_cur)
          "A_comparison_at_twodecname(res_high)A_prev_pair`;
          `ncyc 0`

    * `mode="first"`

          use structure `xyzin_start`);
          results with prefix `structure args.project_twodecname(res_cur)A.pdb`
          plus extra copy of log `structure args.project_
          twodecname(res_cur)A_comparison_at_twodecname(res_high)A.log`;
          `ncyc` controlled by an option `--prerefinement-ncyc` (0 cycles
          by default, 20 cycles by default for the complete
          cross-validation);
          usually `res_cur = res_init`, `res_high = res_init`

    If the input structure model was in mmCIF format, files with suffix .mmcif
    (not .pdb) are used.

    Args:
        res_cur (float)
        res_prev (float or str)
        res_high (float)
        args (parser): Input arguments (including e. g. name of the project) \
                       parsed by `argparse` via function process_arguments()
        n_bins_low (int)
        mode (str): [`"refine"`, `"comp"`, `"prev_pair"`, `"first"`]
        res_low (float)
        res_highest (float)
        flag (int)
        xyzin_start (str): Filename of the model to be refined (valid only \
                           for `mode="first"`)
        bfac_set (float): Value of B-factor that will be set to all atoms
                          before refinement (not used now)

    Returns:
        (dict):
            Dictionary containing names of files that have been created
            by REFMAC5 and a version of REFMAC5, *e. i.*
            `HKLOUT`, `XYZOUT`, `LOGOUT`, and `version` (all `str`)
    """
    if res_low:
        reso = twodec(res_low) + " " + twodec(res_high)
    else:
        reso = twodec(res_high)
        res_low = "Dmax"

    if not args.comin:
        com = """
make -
    check NONE
refi -
    type REST reso """ + reso + """ -
    resi MLKF -
    meth CGMAT -
    bref MIXED
scal -
    type SIMP -
    LSSC -
    ANISO -
    EXPE
solvent YES
"""
    else:                                 # args.comin is str, name of the file
        with open(args.comin, "r") as comfile:  # comfile is the opened file
            com = comfile.read()   # com is str (modified content of the file)
        re_end = re.compile(re.escape("end"), re.IGNORECASE)
        com = re_end.sub("", com)
        # com = re.sub("(?i)end","", com) # case-insensitive
        com += "\n refi reso " + reso + " \n"
        if not any(["ncyc" in line.lower() for line in com.splitlines()]) \
           and not args.ncyc:
            args.ncyc = 10
    prefix = args.project + "_R" + str(flag).zfill(2) + "_" \
        "" + twodecname(res_cur) + "A"
    if mode == "comp":
        print(" .", end="")
        # print("       Calculating statistics of the refined structure model."
        ncyc_not_zero = False
        prefix += "_comparison" \
            "_at_" + twodecname(res_high) + "A"
        com += "\n ncyc 0"
        com += "\n bins " + str(n_bins_low)
    elif mode == "prev_pair":
        if args.complete_cross_validation:
            print(" .")
        else:
            print(" .", end="")
        ncyc_not_zero = False
        prefix += "_comparison_at_" + twodecname(res_high) + "A_prev_pair"
        com += "\n ncyc 0"
        com += "\n bins " + str(n_bins_low)
    elif mode == "first":
        prefix += ""
        if args.complete_cross_validation:
            ncyc_not_zero = True
            if args.prerefinement_ncyc:
                com += "\n ncyc " + str(args.prerefinement_ncyc)
            elif args.ncyc:
                com += "\n ncyc " + str(args.ncyc)
            elif not args.comin:
                com += "\n ncyc 20"
        else:
            if args.prerefinement_ncyc:
                ncyc_not_zero = True
                com += "\n ncyc " + str(args.prerefinement_ncyc)
            else:
                ncyc_not_zero = False
                com += "\n ncyc 0"
        if args.quick:
            ncyc_not_zero = True
            com += "\n ncyc 1"
        com += "\n bins " + str(n_bins_low)
        xyzin = xyzin_start
    elif mode == "refine":
        prefix += ""
        com += "\n bins " + str(n_bins_low)
        ncyc_not_zero = True
        if args.ncyc:
            com += "\n ncyc " + str(args.ncyc)
        elif not args.comin:
            com += "\n ncyc 10"
        if args.quick:
            com += "\n ncyc 1"
    if bfac_set:
        com += "\n BFACtor SET " + twodec(bfac_set)
    if args.weight:
        com += "\n weight matrix " + fourdec(args.weight)
    # Shannon factor
    if args.constant_grid and res_highest:
        shannon_factor = 1.5 * res_high / res_highest
        com += "\n SHANnon_factor " + fourdec(shannon_factor)
    # TLS - TLSC line
    if args.tlsin and ncyc_not_zero:
        if args.tls_ncyc:
            com += "\n refi tlsc " + str(args.tls_ncyc)
        else:
            com += "\n refi tlsc 10"
    com += "\n monitor MEDIum"
    com += "\n free " + str(flag)
    com += "\n END"
    hklout = prefix + ".mtz"
    if "cif" in settings["pdbORmmcif"]:
        xyzout = prefix + settings["pdbORmmcif"]
        xyzout_secondary = prefix + ".pdb"
    else:  # pdb
        xyzout = prefix + settings["pdbORmmcif"]
        xyzout_secondary = prefix + ".mmcif"
    libout = prefix + ".cif"
    logout = prefix + ".log"
    if mode == "refine":
        xyzin = args.project + "_R" + str(flag).zfill(2) + "_" + \
            twodecname(res_prev) + "A" + settings["pdbORmmcif"]
    elif mode == "comp" or mode == "prev_pair":
        xyzin = args.project + "_R" + str(flag).zfill(2) + "_" + \
            twodecname(res_cur) + "A" + settings["pdbORmmcif"]
    command = ["refmac5", "HKLIN", args.hklin, "XYZIN", xyzin, "HKLOUT",
               hklout, "XYZOUT", xyzout, "LIBOUT", libout]
    # Optional LIBIN
    if args.libin:
        command.append("LIBIN")
        command.append(args.libin)
    # Optional TLSIN (and TLSOUT)
    if args.tlsin:
        if ncyc_not_zero:
            command.append("TLSIN")
            if args.tlsin_keep:
                command.append(args.tlsin)
            else:
                if mode == "first":
                    tlsin = args.tlsin
                elif mode == "refine":
                    tlsin = args.project + "_R" + str(flag).zfill(2) + "_" + \
                        twodecname(res_prev) + "A.tlsout"
                command.append(tlsin)
            command.append("TLSOUT")
            tlsout = prefix + ".tlsout"
            command.append(tlsout)
        elif mode == "first":  # and ncyc 0
            shutil.copy2(args.tlsin, prefix + ".tlsout")

    if (mode == "refine" or
            (mode == "first" and args.complete_cross_validation)):
        if args.complete_cross_validation:
            print("     – FreeRflag set " + str(flag))
        print("       Running command:")
        print("       " + " ".join(command))
    with open(logout, "w") as logfile:
        p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=logfile)
        #                    encoding='utf8')  # Probably required in Python 3
        p.communicate(com)
    for fileout in [logout, hklout, xyzout]:
        if not os.path.isfile(fileout):
            sys.stderr.write("ERROR: File " + fileout + " has not been created"
                             " by REFMAC5. Check a file " + logout + " for "
                             "the details.\nAborting.\n")
            sys.exit(1)
    # Copy log and tls while running REFMAC5 at the starting resolution
    if mode == "first":
        prefix_copy = prefix + "_comparison_at_" + twodecname(res_high) + "A"
        shutil.copy2(logout, prefix_copy + ".log")
        shutil.copy2(hklout, prefix_copy + ".mtz")
    version = extract_from_file(logout, "version", 0, 1, nth_word=5,
                                get_first=True)
    results = {"HKLOUT": hklout, "XYZOUT": xyzout, "LOGOUT": logout,
               "version": version}
    if mode == "comp" or mode == "prev_pair":
        files_to_be_removed = [xyzout, xyzout_secondary, libout]
        if "tlsout" in vars():
            to_be_removed += [tlsout]
        for filename in files_to_be_removed:
            if os.path.exists(filename):
                os.remove(filename)
    return results


def refinement_phenix(res_cur,
                      res_prev,
                      res_high,
                      args,
                      n_bins=1,
                      mode="refine",
                      res_low=0,
                      res_highest=0,
                      flag=0,
                      xyzin_start="",
                      bfac_set=0,
                      label=None):  # TODO: hklfree
    """Refine using phenix.refine.

    Resolution limit is controlled by parameters `res_high` and `res_low`.
    Names of files are controlled by parameters `res_cur` and `res_high`.

    There are 4 different modes:

    * `mode="refine"`

          refine structure `args.project_twodecname(res_prev)A_001.pdb`;
          results with prefix `args.project_twodecname(res_cur)A_001`;
          `number_of_macro_cycles` given by `deffile` (default: 3)

    * `mode="comp"`

          use structure `args.project_twodecname(res_cur)A_001.pdb`;
          results with prefix `args.project_twodecname(res_cur)
          A_comparison_at_twodecname(res_high)A_001`;
          `number_of_macro_cycles=1`, `refine.strategy=None`,
          ordered_solvent=False

    * `mode="prev_pair"`

          use structure `args.project_twodecname(res_cur)A_001.pdb`;
          results with prefix `args.project_twodecname(res_cur)
          "A_comparison_at_twodecname(res_high)A_prev_pair_001`;
          `number_of_macro_cycles=1`, `refine.strategy=None`,
          ordered_solvent=False

    * `mode="first"`

          use structure `xyzin_start`;
          results with prefix 
          `structure args.project_twodecname(res_cur)A_001.pdb`
          plus extra copy of log `structure args.project_
          twodecname(res_cur)A_comparison_at_twodecname(res_high)A_001.log`;
          `number_of_macro_cycles` controlled by an option 
          `--prerefinement-ncyc` (1 cycles with no strategy and no ordered
          solvent by default, 6 cycles by default for the complete
          cross-validation);
          usually `res_cur = res_init`, `res_high = res_init`

    If the input structure model was in mmCIF format, files with suffix .cif
    (not .pdb) are used.

    Args:
        res_cur (float)
        res_prev (float or str)
        res_high (float)
        args (parser): Input arguments (including e. g. name of the project) \
                       parsed by `argparse` via function process_arguments()
        n_bins (int)
        mode (str): [`"refine"`, `"comp"`, `"prev_pair"`, `"first"`]
        res_low (float)
        res_highest (float)
        flag (int)
        xyzin_start (str): Filename of the model to be refined (valid only \
                           for `mode="first"`)
        bfac_set (float): Value of B-factor that will be set to all atoms
                          before refinement (not used now)

    Returns:
        (dict):
            Dictionary containing names of files that have been created
            by phenix.refine and a version of phenix.refine, *e. i.*
            `HKLOUT`, `XYZOUT`, `LOGOUT`, and `version` (all `str`)
    """
    prefix = args.project + "_R" + str(flag).zfill(2) + "_" \
        "" + twodecname(res_cur) + "A"
    com = "refinement.input.xray_data.high_resolution=" + twodec(res_high)
    if res_low:
        com += "\nrefinement.input.xray_data.low_resolution=" + twodec(res_low)
    if args.defin:
        with open(args.defin, "r") as deffile:
            deff = deffile.read()
            if not any(["number_of_macro_cycles" in line for line in deff.splitlines()]) \
               and not args.ncyc:
                args.ncyc = 3

    if mode == "comp":
        print(" .", end="")
        # print("       Calculating statistics of the refined structure model."
        ncyc_not_zero = False
        prefix += "_comparison" \
            "_at_" + twodecname(res_high) + "A"
        com += "\nrefinement.main.number_of_macro_cycles=1"
        com += "\nrefinement.refine.strategy=None"
        com += "\nrefinement.main.ordered_solvent=False"
    elif mode == "prev_pair":
        if args.complete_cross_validation:
            print(" .")
        else:
            print(" .", end="")
        ncyc_not_zero = False
        prefix += "_comparison_at_" + twodecname(res_high) + "A_prev_pair"
        com += "\nrefinement.main.number_of_macro_cycles=1"
        com += "\nrefinement.refine.strategy=None"
        com += "\nrefinement.main.ordered_solvent=False"
    elif mode == "first":
        prefix += ""
        if args.complete_cross_validation:
            ncyc_not_zero = True
            if args.prerefinement_ncyc:
                com += "\nrefinement.main.number_of_macro_cycles=" + \
                    str(args.prerefinement_ncyc)
            elif args.ncyc:
                com += "\nrefinement.main.number_of_macro_cycles=" + str(args.ncyc)
            elif not args.comin:
                com += "\nrefinement.main.number_of_macro_cycles=6"
        else:
            if args.prerefinement_ncyc:
                ncyc_not_zero = True
                com += "\nrefinement.main.number_of_macro_cycles=" + \
                    str(args.prerefinement_ncyc)
            else:
                ncyc_not_zero = False
                com += "\nrefinement.main.number_of_macro_cycles=1"
                com += "\nrefinement.refine.strategy=None"
                com += "\nrefinement.main.ordered_solvent=False"
        if args.quick:
            ncyc_not_zero = True
            com += "\nrefinement.main.number_of_macro_cycles=1"
        xyzin = xyzin_start
    elif mode == "refine":
        prefix += ""
        ncyc_not_zero = True
        if args.ncyc:
            com += "\nrefinement.main.number_of_macro_cycles=" + str(args.ncyc)
        elif not args.comin:
            com += "\nrefinement.main.number_of_macro_cycles=3"
        if args.quick:
            com += "\nrefinement.main.number_of_macro_cycles=1"
    com += "\nrefinement.output.n_resolution_bins=" + str(n_bins)
    com += "\nrefinement.output.prefix=" + prefix
    com += "\nrefinement.input.xray_data.r_free_flags.test_flag_value=" + \
        str(flag)
    if hasattr(args, "label"):  # always false for the very first mode="refine"
        com += "\nrefinement.input.xray_data.labels=" + args.label
    # Output filenames
    if "cif" in settings["pdbORmmcif"]:
        xyzout = prefix + "_001" + settings["pdbORmmcif"]
        xyzout_secondary = prefix + "_001" + ".pdb"
    else:  # pdb
        xyzout = prefix + "_001" + settings["pdbORmmcif"]
        xyzout_secondary = prefix + "_001" + ".cif"
    hklout = prefix + "_001.mtz"
    logout = prefix + "_001.log"
    geoout = prefix + "_001.geo"
    outout = prefix + "_001.out"
    params = prefix + "_001.params"
    if mode == "refine":
        xyzin = args.project + "_R" + str(flag).zfill(2) + "_" + \
            twodecname(res_prev) + "A_001" + settings["pdbORmmcif"]
    elif mode == "comp" or mode == "prev_pair":
        xyzin = args.project + "_R" + str(flag).zfill(2) + "_" + \
            twodecname(res_cur) + "A_001" + settings["pdbORmmcif"]
    command = ["phenix.refine", args.hklin, xyzin]
    if args.libin:
        command.append(args.libin)
    if args.defin:
        command.append(args.defin)
    with open(params, "w") as pars:
        pars.write(com)
    command.append(params)

    if (mode == "refine" or
            (mode == "first" and args.complete_cross_validation)):
        if args.complete_cross_validation:
            print("     – FreeRflag set " + str(flag))
        print("       Running command:")
        print("       " + " ".join(command))
    with open(outout, 'w') as out:
        p = subprocess.Popen(command, stdout=out, stderr=out)
        #                    encoding='utf8')  # Probably required in Python 3
        p.communicate()
    print("xyzout", xyzout)
    for fileout in [logout, hklout, xyzout, outout]:
        if not os.path.isfile(fileout):
            error = True
            if mode == "first" and not hasattr(args, "label"):
                label, labels_all = refinement_phenix_get_label(outout)
                if label:
                    warning_my(
                        "amb_labels", "Multiple equally suitable arrays of "
                        "observed xray data found. Possible choices: " + \
                        labels_all + " . Automatically choosing "
                        "refinement.input.xray_data.labels=" + label)
                    args.label = label  # a bit dirty hack but it works
                    results = refinement_phenix(res_cur=res_cur,
                                                res_prev=res_prev,
                                                res_high=res_high,
                                                args=args,
                                                n_bins=n_bins,
                                                mode="first",
                                                res_low=res_low,
                                                res_highest=res_highest,
                                                flag=flag,
                                                xyzin_start=xyzin_start)
                    error = False
            if error:
                sys.stderr.write("ERROR: File " + fileout + " has not been "
                                 "created by phenix.refine. Check a file " + \
                                 outout + " for the details."
                                 "\nAborting.\n")
                sys.exit(1)
    # Copy log and tls while running phenix.refine at the starting resolution
    if mode == "first":
        pdbout = prefix + "_001" + ".pdb"
        prefix_copy = prefix + "_comparison_at_" + twodecname(res_high) + "A"
        shutil.copy2(pdbout, prefix_copy + "_001.pdb")
    version = extract_from_file(logout, "Version", 0, 1, nth_word=1,
                                get_first=True)
    results = {"HKLOUT": hklout, "XYZOUT": xyzout, "LOGOUT": logout,
               "version": version}
    if mode == "comp" or mode == "prev_pair":
        files_to_be_removed = [hklout, geoout]
        if "tlsout" in vars():
            to_be_removed += [tlsout]
        for filename in files_to_be_removed:
            if os.path.exists(filename):
                os.remove(filename)
    return results


def refinement_phenix_get_label(outout):
    """Get possible choices of refinement.input.xray_data.labels from 
    standard output from phenix.refine (saved in file `outout`) if
    multiple equally suitable arrays of observed xray data were found.

    Args:
        out (str): Filename

    Returns:
        (tuple): tuple containing the chosen label and all the found labels - \
                 both are (`str`)"""
    signal = extract_from_file(
        filename=outout,
        searched="Multiple equally suitable arrays of observed xray data foun",
        skip_lines=0, n_lines=0, not_found="N/A")
    if signal == ["N/A"]:
        return None, None
    with open(outout, "r") as o:
        lines = o.read().splitlines()
    for i in range(len(lines)):
        if "Possible choices:" in lines[i]:
            j = i
            choices = []
    if "j" not in vars():
        return None, None
    while lines[j] != "":
        j += 1
        choices.append(lines[j].strip())
    labels_all = " ".join(choices)
    for choice in reversed(choices):
        if "+" in choice and len(choices) > 1:
            choices.remove(choice)
    label = choices[0]
    return label, labels_all


def collect_stat_BINNED(shells, project, hklin, n_bins_low, flag,
                        res_low, refinement="refmac"):
    """Collects statistics of a particular structure model depending on
    resolution (*e. i.* values are binned) and saves them in a CSV file.
    Statistics are picked from a REFMAC5 logfile relating to the model
    `project_RXX_twodecname(shells[-1])A.pdb`, where `XX` is a number
    of a flag.

    This function is called by the function `main()` in file `launcher.py`.

    Args:
        shells (list) - Contains high resolution diffraction limits
                        from the initial to the last which were used
                        for that model (float).
        project (str): Name of the project
        hklin (str): Name of diffraction data MTZ file
        n_bins_low (int)
        flag (int)
        res_low (float): Low resolution of diffraction data `hklin`
                         (required for calling the function
                         `collect_stat_refmac_log_low()`)
        refinement (str): "refmac" or "phenix"

    Returns:
        str: Name of the created CSV file
    """
    print("       Collecting statistics from logfiles...")
    prefix = project + "_R" + str(flag).zfill(2) + "_" \
        "" + twodecname(shells[-1]) + "A"
    if refinement == "refmac":
        # Pick statistics for `n_bins_low` res. shells up to initial diffr. limit
        logfilename = prefix + "_comparison" \
            "_at_" + twodecname(shells[0]) + "A.log"
        mtzfilename = prefix + ".mtz"
        bin_res_mean, bin_res_low, bin_res_high, bin_Nwork, bin_Nfree, \
            bin_Rwork, bin_Rfree, bin_CCwork, bin_CCfree = \
            collect_stat_binned_refmac_low(
                logfilename, mtzfilename, hklin, n_bins_low, res_low, flag)
        # Pick overall values for data up to previous diffraction limit
        # (to be comparable pairwisely)
        filename_prefix = prefix
        overall_Rwork, overall_Rfree = \
            collect_stat_overall_refmac(prefix + ".log", flag)
        overall_CCwork, overall_CCfree = calculate_correlation(mtzfilename,
                                                               hklin, flag)
        overall_CCwork = fourdec(overall_CCwork)
        overall_CCfree = fourdec(overall_CCfree)
    elif refinement == "phenix":
        # Pick statistics for `n_bins_low` res. shells up to initial diffr. limit
        pdbfilename = prefix + "_comparison" \
            "_at_" + twodecname(shells[0]) + "A_001.pdb"
        bin_res_low, bin_res_high, bin_Nwork, bin_Nfree, \
            bin_Rwork, bin_Rfree, bin_CCwork, bin_CCfree = \
            collect_stat_binned_phenix_low(pdbfilename, n_bins_low, res_low)
        # Pick overall values for data up to previous diffraction limit
        # (to be comparable pairwisely)
        pdbfilename = prefix + "_001.pdb"
        overall_Rwork, overall_Rfree = \
            collect_stat_overall_phenix(pdbfilename)
        overall_CCwork = "N/A"
        overall_CCfree = "N/A"
        overall_CCavg = "N/A"
    # Write the statistics to a csv file
    csvfilename = prefix + ".csv"
    with open(csvfilename, "w") as csvfile:
        csvfile.write("# Statistics of refined structure model calculated "
                      "at resolution range " + twodec(res_low) + "-"
                      "" + twodec(shells[-1]) + " A.\n")
        csvfile.write("# Rwork = " + overall_Rwork + " \n")
        csvfile.write("# Rfree = " + overall_Rfree + " \n")
        csvfile.write("# CCwork = " + overall_CCwork + " \n")
        csvfile.write("# CCfree = " + overall_CCfree + " \n")
        csvfile.write("# CCavg = N/A \n")
        csvfile.write("#\n")
        csvfile.write("#Shell Res_low - Res_hig      Nwork    Nfree  Rwork    "
                      "Rfree   CCwork  CCfree \n")
    collect_stat_write(csvfilename, bin_res_low, bin_res_high,
                              bin_Nwork, bin_Nfree,
                              bin_Rwork, bin_Rfree,
                              bin_CCwork, bin_CCfree)
    # Pick stats. for the high resolution shells and write them to the CSV file
    for i in range(len(shells) - 1):
        if refinement == "refmac":
            bin_res_low = [twodec(shells[i])]
            bin_res_high = [twodec(shells[i + 1])]
            logfilename = prefix + "_comparison" \
                "_at_" + twodecname(shells[i + 1]) + "A.log"
            mtzfilename = prefix + ".mtz"
            bin_Nwork, bin_Nfree, bin_Rwork, bin_Rfree, bin_CCwork, \
                bin_CCfree = collect_stat_binned_refmac_high(
                    logfilename, mtzfilename, hklin, n_bins_low,
                    float(twodec(shells[i])), float(twodec(shells[i + 1])), flag)
        elif refinement == "phenix":
            pdbfilename = prefix + "_comparison" \
                "_at_" + twodecname(shells[i + 1]) + "A_001.pdb"
            bin_res_low, bin_res_high, bin_Nwork, bin_Nfree, bin_Rwork, bin_Rfree, \
                bin_CCwork, bin_CCfree = collect_stat_binned_phenix_high(pdbfilename, n_bins_low)
        collect_stat_write(
            csvfilename, bin_res_low, bin_res_high, bin_Nwork, bin_Nfree,
            bin_Rwork, bin_Rfree, bin_CCwork, bin_CCfree,
            shell_number=n_bins_low + 1 + i)
    return csvfilename


def collect_stat_binned_refmac_low(logfilename, mtzfilename, hklin,
                                   n_bins_low, res_low=999, flag=0):
    """Picks and returns statistics values in the given `REFMAC5` logfile.
    Logfile supposed to contain information from `n_bins_low` shells.

    This function is called by the function `collect_stat_BINNED()`
    from this file and `main()` from file `launcher.py`.

    Args:
        logfilename (str): Name of a `REFMAC5` logfile
        n_bins_low (int): Number of low resolution bins

    Returns:
        (tuple): tuple containing statistics
        values `bin_Nwork`, `bin_Nfree`, `bin_Rwork`,
        `bin_Rfree`, `bin_CCwork`, and `bin_CCfree` (all are `str`)
    """
    bin_4SSQLL_mean = []
    bin_res_mean = []
    bin_res_low = []
    bin_res_high = []
    bin_Nwork = []
    bin_Nfree = []
    bin_Rwork = []
    bin_Rfree = []
    bin_CCwork = []
    bin_CCfree = []

    # # Definition of bins, numbers of reflections, R-values
    logfile_lines = extract_from_file(
        logfilename, "Things for loggraph, R factor and others  ",
        11, n_bins_low)
    # Refmac gives only mean 4SSQLL, so res_high and res_low (in angtroem)
    # calculations are neccessary
    shell_step_4SSQLL = \
        float(logfile_lines[1].split()[0]) - float(logfile_lines[0].split()[0])
    for i in range(n_bins_low):
        bin_4SSQLL_mean.append(logfile_lines[i].split()[0])
        bin_res_mean.append(str(twodec(sqrt(1 / float(bin_4SSQLL_mean[i])))))
        if i == 0:
            # This is more accurate and there is no danger of zero division
            res_low_cur = res_low
        else:
            res_low_cur = sqrt(
                1 / (abs(float(bin_4SSQLL_mean[i]) - shell_step_4SSQLL / 2)))
        bin_res_low.append(twodec(res_low_cur))
        res_high_cur = sqrt(
            1 / (float(bin_4SSQLL_mean[i]) + shell_step_4SSQLL / 2))            
        bin_res_high.append(twodec(res_high_cur))
        bin_Nwork.append(logfile_lines[i].split()[1])
        try:
            bin_Nfree.append(logfile_lines[i].split()[7])
            bin_Rwork.append(logfile_lines[i].split()[5])
            bin_Rfree.append(logfile_lines[i].split()[10])
        except IndexError:
            bin_Nfree.append("N/A")
            bin_Rwork.append("N/A")
            bin_Rfree.append("N/A")
            warning_my("low_R", "All R-values relating to the resolution"
                       " bin " + twodec(bin_res_low[-1]) + "-"
                       "" + twodec(bin_res_high[-1]) + " were not calculated "
                       "successfully. "
                       "For further details, see file " + logfilename + ".")
        CCwork, CCfree = calculate_correlation(
            mtzfilename, hklin, flag, res_low_cur, res_high_cur)
        bin_CCwork.append(fourdec(CCwork))
        bin_CCfree.append(fourdec(CCfree))
    return(bin_res_mean, bin_res_low, bin_res_high,
           bin_Nwork, bin_Nfree,
           bin_Rwork, bin_Rfree,
           bin_CCwork, bin_CCfree)


# TODO - what to do if CC-values are not available, strange log format
def collect_stat_binned_refmac_high(
    logfilename, mtzfilename, hklin, n_bins_low, res_low, res_high, flag=0):
    """Picks and returns statistics values in the given `REFMAC5` logfile.
    Logfile supposed to contain information from 1 shells.

    This function is called by the function `collect_stat_BINNED()`.

    Args:
        logfilename (str): Name of a `REFMAC5` logfile
        n_bins_low (int): Number of low resolution bins

    Returns:
        (tuple): tuple containing statistics
        values `bin_Nwork`, `bin_Nfree`, `bin_Rwork`,
        `bin_Rfree`, `bin_CCwork`, and `bin_CCfree` (all are `str`)
    """
    # R-values and CC-values
    bin_Rwork = [extract_from_file(logfilename, "Overall R factor", 0, 1, -1)]
    bin_Rfree = [extract_from_file(logfilename, "Free R factor", 0, 1, -1)]
    if "*" in bin_Rwork:
        bin_Rwork = "N/A"  # TODO: Seems it is not working now...
        bin_Rfree = "N/A"
        warning_my("high_R", "R-values of a particular shell "
                   "were not calculated. For "
                   "further details, see file " + logfilename + ".")
    CCwork, CCfree = calculate_correlation(
        mtzfilename, hklin, flag, res_low, res_high)
    bin_CCwork = [fourdec(CCwork)]
    bin_CCfree = [fourdec(CCfree)]
    bin_Nwork = [extract_from_file(logfilename, "Number of used reflections",
                 0, 1, -1)]

    # A number of reflections in bin
    bin_Nfree = []
    logfile_lines = extract_from_file(
        logfilename, "Things for loggraph, R factor and others  ",
        11, n_bins_low)
    try:
        for i in range(n_bins_low):
            bin_Nfree.append(int(logfile_lines[i][49:57]))
        bin_Nfree = [str(sum(bin_Nfree))]
    except IndexError:
        bin_Nfree = ["N/A"]
    return(bin_Nwork, bin_Nfree, bin_Rwork, bin_Rfree,
           bin_CCwork, bin_CCfree)


def collect_stat_binned_phenix_low(pdbfilename, n_bins_low, res_low=999):
    """Picks and returns statistics values in the given `REFMAC5` logfile.
    Logfile supposed to contain information from `n_bins_low` shells.

    This function is called by the function `collect_stat_BINNED()`
    from this file and `main()` from file `launcher.py`.

    Args:
        pdbfilename (str): Name of a `REFMAC5` logfile
        n_bins_low (int): Number of low resolution bins

    Returns:
        (tuple): tuple containing statistics
        values `bin_Nwork`, `bin_Nfree`, `bin_Rwork`,
        `bin_Rfree`, `bin_CCwork`, and `bin_CCfree` (all are `str`)
    """
    bin_res_low = []
    bin_res_high = []
    bin_Nwork = []
    bin_Nfree = []
    bin_Rwork = []
    bin_Rfree = []
    bin_CCwork = []
    bin_CCfree = []
    # # Definition of bins, numbers of reflections, R-values
    pdbfile_lines = extract_from_file(
        pdbfilename, "REMARK   3   BIN  RESOLUTION RANGE  COMPL.    "
        "NWORK NFREE   RWORK  RFREE  CCWORK CCFREE", 1, n_bins_low)
    for i in range(n_bins_low):
        bin_res_low.append(pdbfile_lines[i].split()[3])
        bin_res_high.append(pdbfile_lines[i].split()[5])
        bin_Nwork.append(pdbfile_lines[i].split()[7])
        bin_Nfree.append(pdbfile_lines[i].split()[8])
        bin_Rwork.append(pdbfile_lines[i].split()[9])
        bin_Rfree.append(pdbfile_lines[i].split()[10])
        bin_CCwork.append(pdbfile_lines[i].split()[11])
        bin_CCfree.append(pdbfile_lines[i].split()[12])
    return(bin_res_low, bin_res_high, bin_Nwork, bin_Nfree,
           bin_Rwork, bin_Rfree, bin_CCwork, bin_CCfree)


def collect_stat_binned_phenix_high(pdbfilename, n_bins_low):
    """Picks and returns statistics values in the given `REFMAC5` logfile.
    Logfile supposed to contain information from 1 shells.

    This function is called by the function `collect_stat_BINNED()`.

    Args:
        logfilename (str): Name of a `REFMAC5` logfile
        n_bins_low (int): Number of low resolution bins

    Returns:
        (tuple): tuple containing statistics
        values `bin_Nwork`, `bin_Nfree`, `bin_Rwork`,
        `bin_Rfree`, `bin_CCwork`, and `bin_CCfree` (all are `str`)
    """
    line = extract_from_file(
        pdbfilename, "REMARK   3   BIN  RESOLUTION RANGE  COMPL.    "
       "NWORK NFREE   RWORK  RFREE  CCWORK CCFREE", 1, 1, not_found="N/A")
    line = line[0]
    bin_res_low = [line.split()[3]]
    bin_res_high = [line.split()[5]]
    bin_Nwork = [line.split()[7]]
    bin_Nfree = [line.split()[8]]
    bin_Rwork = [line.split()[9]]
    bin_Rfree = [line.split()[10]]
    bin_CCwork = [line.split()[11]]
    bin_CCfree = [line.split()[12]]
    return(bin_res_low, bin_res_high, bin_Nwork, bin_Nfree,
           bin_Rwork, bin_Rfree, bin_CCwork, bin_CCfree)


def collect_stat_write(csvfilename, bin_res_low, bin_res_high,
                              bin_Nwork, bin_Nfree,
                              bin_Rwork, bin_Rfree,
                              bin_CCwork, bin_CCfree,
                              shell_number=1):
    """Saves the given statistics values to a CSV file.

    This function is called by the function `collect_stat_BINNED()`.

    Args:
        csvfilename (str): Filename of the CSV file to be created.
        bin_* (list): Lists with statistics values
        shell_number (int): Number of the shell corresponding to the
                            statistics values

    Returns:
        str: Filename of the created CSV file.
    """
    with open(csvfilename, "a") as csvfile:
        for i in range(len(bin_Rwork)):
            csvfile.write(
                str(shell_number).zfill(2) + "        " + bin_res_low[i] + ""
                " - " + bin_res_high[i] + "         " + bin_Nwork[i] + "     "
                "" + bin_Nfree[i] + "     " + bin_Rwork[i] + "     "
                "" + bin_Rfree[i] + "     " + bin_CCwork[i] + "     "
                "" + bin_CCfree[i] + "\n")
            if float(bin_Nfree[i]) < 50:
                warning_my("Nfree", "There are only " + str(bin_Nfree[i]) + ""
                           " < 50 free reflections in the resolution shell "
                           "" + twodec(bin_res_low[i]) + "-"
                           "" + twodec(bin_res_high[i]) + " A. Values of "
                           "statistics Rfree and CCfree in this shell "
                           "could be misleading. Consider setting "
                           "thicker resolution shells.")
            shell_number = shell_number + 1
    return csvfilename


def collect_stat_OVERALL(shells, args, flag, refinement="refmac"):
    """Collects overall statistics of a particular structure model
    from a REFMAC5 logfile or from PDB file (phenix.refine).
    Model which is dealing with: `project_RXX_twodecname(shells[-1])A.pdb`
    (REFMAC5) or `project_RXX_twodecname(shells[-1])A_001.pdb` (phenix.refine),
    where `XX` is a number of a flag.

    This function is called by the function `main()` in file `launcher.py`.

    Args:
        shells (list):  Contains high resolution diffraction limits
                        from the initial to the last which were used
                        for that model (float).
        project (str): Name of the project
        flag (int):
        refinement (str): "refmac" or "phenix"

    Returns:
        (tuple): tuple containing names of created CSV files
    """
    csvfilenames = []
    prefix = args.project + "_R" + str(flag).zfill(2)
    # Pick overall values - values calculated on the same data - "paired"
    if len(shells) >= 2:
        if refinement == "refmac":
            # get Rwork, Rfree, CCwork, CCfree, CCavg
            filename_prefix = prefix + "_" + twodecname(shells[-2]) + "A"
            Rwork_before, Rfree_before = \
                collect_stat_overall_refmac(filename_prefix + ".log", flag)
            CCwork_before, CCfree_before = calculate_correlation(
                filename_prefix + ".mtz", args.hklin, flag) # floats
            filename_prefix = prefix + "_" + twodecname(shells[-1]) + "A_" \
                "comparison_at_" + twodecname(shells[-2]) + "A_prev_pair"
            Rwork_after, Rfree_after = \
                collect_stat_overall_refmac(filename_prefix + ".log", flag)
            mtzfilename = prefix + "_" + twodecname(shells[-1]) + "A.mtz"
            CCwork_after, CCfree_after = calculate_correlation(
                mtzfilename, args.hklin, flag, res_high=float(shells[-2]))  # floats
        elif refinement == "phenix":
            # get Rwork, Rfree                (not CCwork, CCfree, CCavg)
            pdbfilename = prefix + "_" + twodecname(shells[-2]) + "A_001.pdb"
            Rwork_before, Rfree_before = \
                collect_stat_overall_phenix(pdbfilename)
            pdbfilename = prefix + "_" + twodecname(shells[-1]) + "A_comparison" \
                "_at_" + twodecname(shells[-2]) + "A_prev_pair_001.pdb"
            Rwork_after, Rfree_after = \
                collect_stat_overall_phenix(pdbfilename)

        # Rwork, Rfree
        Rwork_change = fourdec(float(Rwork_after) - float(Rwork_before))
        Rfree_change = fourdec(float(Rfree_after) - float(Rfree_before))
        #
        csvfilename = prefix + "_R-values.csv"
        csvfilenames.append(csvfilename)
        # Write the begining of the csv file if it does not exist yet
        if not os.path.isfile(csvfilename):
            with open(csvfilename, "w") as csvfile:
                csvfile.write("# Shell      Rwork(init) Rwork(fin) Rwork(diff)"
                              "   Rfree(init) Rfree(fin) Rfree(diff)\n")
        # Formating issues
        space_Rwork = ""
        space_Rfree = ""
        if float(Rwork_change) >= 0:
            space_Rwork = " "
        if float(Rfree_change) >= 0:
            space_Rfree = " "
        # Write values in the csv file
        with open(csvfilename, "a") as csvfile:
            csvfile.write(twodec(shells[-2]) + "A->" + twodec(shells[-1]) + "A"
                          "      " + Rwork_before + "     " + Rwork_after + ""
                          "     " + space_Rwork + Rwork_change + "        "
                          "" + Rfree_before + "     " + Rfree_after + "     "
                          "" + space_Rfree + Rfree_change + "\n")

    # Pick overall values of the current structure model
    # and find Rfree-Rwork gap (at the initial resolution)
    if refinement == "refmac":     # get Rwork, Rfree
        logfilename = prefix + "_" + twodecname(shells[-1]) + "A" \
            "_comparison_at_" + twodecname(shells[0]) + "A.log"
        Rwork, Rfree = collect_stat_overall_refmac(logfilename, flag)
    elif refinement == "phenix":   # get Rwork, Rfree
        pdbfilename = prefix + "_" + twodecname(shells[-1]) + "A" \
            "_comparison_at_" + twodecname(shells[0]) + "A_001.pdb"
        Rwork, Rfree = collect_stat_overall_phenix(pdbfilename)
    Rgap = fourdec(float(Rfree) - float(Rwork))
    csvfilename_gap = prefix + "_Rgap.csv"
    # Write the begining of the csv file if it does not exist yet
    if not os.path.isfile(csvfilename_gap):
        with open(csvfilename_gap, "w") as csvfile:
            csvfile.write("# Resolution   Rwork   Rfree   Rfree-Rwork\n")
    with open(csvfilename_gap, "a") as csvfile:
        csvfile.write(twodec(shells[-1]) + "          " + Rwork + "   "
                      "" + Rfree + "   " + Rgap + "\n")
    csvfilenames.append(csvfilename_gap)
    return tuple(csvfilenames)


def collect_stat_overall_refmac(logfilename, flag=0):
    """Picks and returns overall Rwork, Rfree from a given `REFMAC5` logfile.

    This function is called by the functions `collect_stat_OVERALL()`
    and `collect_stat_BINNED()`.

    Args:
        logfilename (str): Filename of a `REFMAC5` logfile

    Returns:
        (tuple): tuple containing statistics values `Rwork` and `Rfree` \
                 (both are `str`)
    """
    Rwork = extract_from_file(logfilename, "R factor", 0, 1, -1)
    Rfree = extract_from_file(logfilename, "R free", 0, 1, -1)
    if "*" in Rwork:
        Rwork = "N/A"
        Rfree = "N/A"
        warning_my("overall_R", "Some overall R-values were not calculated. "
                   "For further details, see file " + logfilename + ".")
    return Rwork, Rfree


def calculate_correlation(hkl_calc, hklin,
                          flag=0, res_low=None, res_high=None):
    """Calculates CCwork and CCfree using `sftools`.

    Args:
        hkl_calc (str): Name of the MTZ file from a REFMAC5 run
        hklin (str): Name of the MTZ file with diffraction data
        flag (int): free reflection flag set
        res_low (float): low-resolution cutoff
        res_high (float): high-resolution cutoff

    Returns:
        (tuple): tuple containing `CCwork` and `CCfree` \
                 (both are `float` or `str`: "N/A")
    """
    from iotbx.reflection_file_reader import any_reflection_file
    from scitbx.array_family import flex

    # i_obs from HKLIN
    i_obs = None
    miller_arrays = any_reflection_file(file_name=hklin).as_miller_arrays()
    for i, column in enumerate(miller_arrays):
        if column.is_xray_intensity_array() and \
                column.anomalous_flag() == False:
            i_obs_label = str(column.info()).split(".mtz:")[1].split(",")[0]
            i_obs = column
            break
    if not i_obs:
        print("Error - intensities were not found in " + hklin)  # TO DO
        return "N/A", "N/A"

    # prepare sftools command and run
    if res_low and res_high:  # else: full resolution of MTZ file
        com_resolution = " resolution " + fourdec(res_low) + " " + fourdec(res_high)
    elif res_high:
        com_resolution = " resolution 999 " + fourdec(res_high)
    else:
        com_resolution = ""
    p = subprocess.Popen(["sftools"], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    com = "read " + hkl_calc + " col 1 FC_ALL\n" \
        "calc col FC_ALLsq = col FC_ALL col FC_ALL *\n" \
        "read " + hklin + " col " + i_obs_label + "\n" \
        "select only col 1 = " + str(flag) + "\n" \
        "correl col " + i_obs_label + " FC_ALLsq shells 1" + com_resolution + "\n" \
        "select only not col 1 = " + str(flag) + "\n" \
        "correl col " + i_obs_label + " FC_ALLsq shells 1" + com_resolution + "\n" \
        "stop\n" \
        "y\n"
    output, err = p.communicate(com)

    # pick values
    CCwork = None
    CCfree = None
    output_lines = output.splitlines()
    for i in range(len(output.splitlines())):
        if "WARNING" in output.splitlines()[i]:
            print(output.splitlines()[i])
            return "N/A", "N/A"
        elif " $TABLE: Correlation vs. resolution:" in output.splitlines()[i] and \
                "$$" in output.splitlines()[i + 2] and \
                "$$" in output.splitlines()[i + 4] and \
                "$$" in output.splitlines()[i + 6] and \
                "NSHELL  DMIN    DMAX    RFACT   CORREL    <F1>    <F2>    NREF $$" in output.splitlines()[i + 3]:
            if not CCfree:
                try:
                    CCfree = float(output.splitlines()[i + 5].split()[4])
                    Nfree = int(output.splitlines()[i + 5].split()[-1])
                except ValueError:
                    CCfree = "N/A"
                    Nfree = "N/A"
            else:  # CCfree already picked
                try:
                    CCwork = float(output.splitlines()[i + 5].split()[4])
                    Nwork = int(output.splitlines()[i + 5].split()[-1])
                except ValueError:
                    CCwork = "N/A"
                    Nwork = "N/A"
    if not CCfree:
        CCfree = "N/A"
    if not CCwork:
        CCwork = "N/A"
    return CCwork, CCfree

    
def collect_stat_overall_phenix(pdbfilename):
    """Picks and returns overall Rwork, Rfree values from a given PDB file 
    refined in phenix.refine.

    This function is called by the functions `collect_stat_OVERALL()`
    and `collect_stat_BINNED()`.

    Args:
        pdbfilename (str): Filename of a PDB structure model refined in \
                           phenix.refine

    Returns:
        (tuple): tuple containing statistics values `Rwork`, `Rfree` \
                 (all are `str`)
    """
    Rwork = extract_from_file(
        pdbfilename,
        "REMARK   3   R VALUE            (WORKING SET) : ", 0, 1, -1)
    Rfree = extract_from_file(
        pdbfilename,
        "REMARK   3   FREE R VALUE                     : ", 0, 1, -1)
    return Rwork, Rfree


def collect_stat_OVERALL_AVG(shells, project, flag_sets):
    """Calculates and saves average overall values from CSV files prepared by
    the function `collect_stat_OVERALL()`.

    This function is called by the function `main()` in file `launcher.py`.

    Args:
        shells (list): Contains high resolution diffraction limits
                       from the initial to the last which were used
                       for that model (float).
        project (str): Name of the project
        flag_sets (list): List of free reflection flag sets (int)

    Returns:
        (tuple): tuple containing names of created CSV files
    """
    import numpy as np
    # === Overall Rwork, Rfree ===
    Rwork_init = []
    Rwork_fin = []
    Rwork_diff = []
    Rfree_init = []
    Rfree_fin = []
    Rfree_diff = []
    for flag in flag_sets:
        csvfilename = project + "_R" + str(flag).zfill(2) + "_R-values.csv"
        with open(csvfilename, "r") as csvfile:
            csvfile_last_line = csvfile.readlines()[-1]
            Rwork_init.append(float(csvfile_last_line.split()[1]))
            Rwork_fin.append(float(csvfile_last_line.split()[2]))
            Rwork_diff.append(float(csvfile_last_line.split()[3]))
            Rfree_init.append(float(csvfile_last_line.split()[4]))
            Rfree_fin.append(float(csvfile_last_line.split()[5]))
            Rfree_diff.append(float(csvfile_last_line.split()[6]))
    Rwork_init_avg = fourdec(np.mean(Rwork_init))
    Rwork_fin_avg = fourdec(np.mean(Rwork_fin))
    Rwork_diff_avg = fourdec(np.mean(Rwork_diff))
    Rfree_init_avg = fourdec(np.mean(Rfree_init))
    Rfree_fin_avg = fourdec(np.mean(Rfree_fin))
    Rfree_diff_avg = fourdec(np.mean(Rfree_diff))
    Rwork_stdev = fourdec(np.std(Rwork_fin))
    Rfree_stdev = fourdec(np.std(Rfree_fin))
    Rwork_diff_sem = str(round(np.std(Rwork_diff)/sqrt(len(flag_sets)), 5))
    Rfree_diff_sem = str(round(np.std(Rfree_diff)/sqrt(len(flag_sets)), 5))
    # Rwork_diff_avg2 = fourdec(float(Rwork_fin_avg) - float(Rwork_init_avg))
    # Rfree_diff_avg2 = fourdec(float(Rfree_fin_avg) - float(Rfree_init_avg))

    csvfilename = project + "_R-values.csv"
    if not os.path.isfile(csvfilename):
        with open(csvfilename, "w") as csvfile:
            csvfile.write("# Shell      Rwork(init) Rwork(fin) Rwork(diff)"
                          "   Rfree(init) Rfree(fin) Rfree(diff)   "
                          "Rwork(StDev)   Rfree(StDev)  "
                          "Rwork(diff,SEM) Rfree(diff,SEM)\n")

    # Formating issues
    space_Rwork = ""
    space_Rfree = ""
    if float(Rwork_diff_avg) >= 0:
        space_Rwork = " "
    if float(Rfree_diff_avg) >= 0:
        space_Rfree = " "
    with open(csvfilename, "a") as csvfile:
        csvfile.write(twodec(shells[-2]) + "A->" + twodec(shells[-1]) + "A"
                      "      " + Rwork_init_avg + "     " + Rwork_fin_avg + ""
                      "     " + space_Rwork + Rwork_diff_avg + "        "
                      "" + Rfree_init_avg + "     " + Rfree_fin_avg + "     "
                      "" + space_Rfree + Rfree_diff_avg + "         "
                      "" + Rwork_stdev + "       " + Rfree_stdev + "     "
                      "" + Rwork_diff_sem + "       " + Rfree_diff_sem + "\n")

    # === Rgap ===
    Rgap_avg = fourdec(float(Rfree_fin_avg) - float(Rwork_fin_avg))
    csvfilename_gap = project + "_Rgap.csv"
    if not os.path.isfile(csvfilename_gap):
        with open(csvfilename_gap, "w") as csvfile:
            csvfile.write("# Resolution   Rwork   Rfree   Rfree-Rwork\n")
    with open(csvfilename_gap, "a") as csvfile:
        csvfile.write(twodec(shells[-1]) + "          " + Rwork_fin_avg + "   "
                      "" + Rfree_fin_avg + "   " + Rgap_avg + "\n")
    return csvfilename, csvfilename_gap
