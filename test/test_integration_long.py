import pytest
import os
import shutil
import filecmp
import warnings
from helper import run, config, tmp_environ
from pairef.preparation import which
from pairef import __version__


# Description of the integration tests
# ------------------------------------
# mdm2 - normal run using refmac, quite fast calculation
# CDO - normal run using refmac
# POLI - run using refmac with --TLSIN and --weight
# BO - run using refmac with --libin and --comfile
# TL - run with free flag No.2 with input mmcif - refmac
#                                               - phenix.refine
# hse11 - normal run using phenix.refine, quite fast calculation, low resolution ~2.9 A


if which("refmac5"):
    @pytest.mark.parametrize(  # res_i requires two decimals, others requires refin. programe to be specified
    "project, files, res_i, res_r, others",
    [("mdm2_0-10A", {"xyzin": "mdm2_1-60A.pdb", "hklin_unmerged": "mdm2_unmerged.mtz", "hklin": "mdm2_merged.mtz"}, "1.60", "1.5,1.4,1.3", "--refmac"),
     ("CDO_0-10A", {"xyzin": "CDO_2B5H_edit_refmac1.pdb", "hklin": "CDO_R.mtz", "hklin_unmerged": "CDO_XDS_ASCII.HKL"}, "2.00", "1.9,1.8,1.7,1.6,1.5,1.42", "--refmac"),
     ("POLI_TLS", {"xyzin": "poli67_edit12_refmac1.pdb", "hklin": "poli67_R.mtz", "hklin_unmerged": "poli67_XDS_ASCII.HKL", "tlsin": "poli67_edit12_refmac1_TLS+Biso.tlsin"}, "2.30", "2.2,2.1,2.0,1.9", "--refmac --weight 0.06 --tlsin poli67_edit12_refmac1_TLS+Biso.tlsin --TLS-ncyc 5"),
     ("BO_LIB", {"xyzin": "BO_edit94_refmac1.pdb", "hklin": "BO_R.mtz", "hklin_unmerged": "BO_XDS_ASCII.HKL", "libin": "BO_TRP-HIS_FC6.cif", "comin": "BO_setting.com"}, "2.59", "2.50", "--refmac --libin BO_TRP-HIS_FC6.cif --comfile BO_setting.com"),
     ("TL_cif_refmac", {"xyzin": "TL_3n21_edit05_refmac1_shaken.mmcif", "hklin": "TL_AUTOMATIC_DEFAULT_free_R.mtz"}, "1.80", "1.70,1.60,1.50", "--refmac --flag 2"),
     ("TL_cif_phenix", {"xyzin": "TL_3n21_edit05_refmac1_shaken.cif", "hklin": "TL_AUTOMATIC_DEFAULT_free_R.mtz", "defin": "TL_setting.def"}, "1.80", "1.70,1.60,1.50", "--phenix --flag 2 --def TL_setting.def"),
     ("hse11", {"xyzin": "hse11_2-90A.pdb", "hklin": "hse11_2-45A.mtz", "hklin_unmerged": "hse11_XDS_ASCII.HKL"}, "2.90", "2.70,2.45", "--phenix")],
     ids=["mdm2_0-10A", "CDO_0-10A", "POLI_TLS", "BO_LIB", "TL_cif_refmac", "TL_cif_phenix", "hse11"])
    def test_integration(project, files, res_i, res_r, others, tmp_environ):
        if "--phenix" in others and not which("phenix.refine"):
            warnings.warn(UserWarning("phenix.refine was not found in executables.\n"))
            pytest.skip("Skipping tests for phenix.refine.")
        print("\nPerforming integration tests, please wait...\n")
        dir_project = "pairef_" + project
        if os.path.isdir(dir_project):
            # Preparation - clean rests from previously examined test that
            # has been interrupted
            shutil.rmtree(dir_project)

        try:
            import urllib.request  # Python 3
            urllib_wget = urllib.request.urlretrieve
        except ImportError:  # Python 2
            import urllib
            urllib_wget = urllib.urlretrieve

        url_prefix = "https://raw.githubusercontent.com/MartinMalyMM/pairef_test_data/main/"
        for f in files.values():
            filename, log = urllib_wget(url_prefix + f, f)
            # try:
            #    filename, log = urllib_wget(url_prefix + f, f)
            #except urllib.error.HTTPError:        # Python 3
            #    filename, log = urllib_wget(url_prefix_alt + f, f)
            ## print(url_prefix + f)
            ## print(bool("Content-Length: 14" in str(log)))
            ## print(str(log))
            #if "Content-Length: 14" in str(log):  # Python 2
            #    print(url_prefix + f + " not found")
            #    print("trying alt")
            #    filename, log = urllib_wget(url_prefix_alt + f, f)
            assert os.path.isfile(f)
        if "TL_cif" in project:
            filename, log = urllib_wget(
                "https://pairef.fjfi.cvut.cz/docs/publication_examples/3-2_TL/pairef_TL_step0-10A/"
                "AUTOMATIC_DEFAULT_scaled_unmerged.mtz",
                "TL_AUTOMATIC_DEFAULT_scaled_unmerged.mtz")
            files["hklin_unmerged"] = "TL_AUTOMATIC_DEFAULT_scaled_unmerged.mtz"
            assert os.path.isfile("TL_AUTOMATIC_DEFAULT_scaled_unmerged.mtz")
        xyzin = files["xyzin"]
        hklin = files["hklin"]
        hklin_unmerged = files["hklin_unmerged"]

        arguments = "--HKLIN " + hklin + \
            " --XYZIN " + xyzin + \
            " -u " + hklin_unmerged + \
            " -i " + res_i + \
            " -r " + res_r + \
            " -p " + project + \
            " " + others
        cp = run(arguments)
        assert cp.returncode == 0
        # assert not cp.stderr
        if cp.stderr:
            print("STDERR: " + cp.stderr)

        expected_stdout = """
 __   _  ___ __  ___ ___
 )_) /_)  )  )_) )_  )_
/   / / _(_ / \ (__ (

automatic PAIRed REFinement protocol
version: """ + str(__version__) + """
run date and time: XX
user@host: XX

Please cite: "Paired refinement under the control of PAIREF"
M. Maly, K. Diederichs, J. Dohnalek, P. Kolenko (2020) IUCrJ 7

Command line arguments: """ + arguments + """

Program has been executed with following input parameters:
"""
        if "--phenix" in others:
            expected_stdout += " * Refinement software: phenix.refine\n"
        elif "--refmac" in others:
            expected_stdout += " * Refinement software: REFMAC5\n"
        expected_stdout += " * XYZIN: " + xyzin + "\n"
        expected_stdout += " * HKLIN: " + hklin + "\n"
        expected_stdout += " * HKLIN unmerged: " + hklin_unmerged + "\n"
        if "--libin" in others:
            expected_stdout += " * LIBIN: " + files["libin"] + "\n"
        if "--tlsin" in others:
            expected_stdout += " * TLSIN: " + files["tlsin"] + "\n"
        expected_stdout += " * Project name: " + project + "\n"
        expected_stdout += " * Resolution shells: " + res_r + "\n"
        if "--comfile" in others:
            expected_stdout += " * Com file for REFMAC5: " + files["comin"] + "\n"
        if "--def" in others:
            expected_stdout += " * Keyword file for phenix.refine: " + files["defin"] + "\n"
        if "--weight" in others:
            expected_stdout += " * Weight matrix: 0.06\n"  # POLI
        if "--TLS-ncyc" in others:                         # POLI
            expected_stdout += " * Number of number of cycles of TLS refinement: 5\n"

        stdout_all = cp.stdout.splitlines(True)
        assert "run date and time:" in stdout_all[7]
        stdout_all[7] = "run date and time: XX\n"
        assert "user@host:" in stdout_all[8]
        stdout_all[8] = "user@host: XX\n"
        flag = "00"
        j = 22  # usually
        for i, line in enumerate(stdout_all):
            if "Resolution of the merged diffraction data" in line:
                j = i - 1
            if "--flag" in others:
                if " * Data with FreeRflag set 2 will be excluded during refinement." in line:
                    flag = "02"  # TL
        if "--flag" in others:
            assert flag == "02"
        assert "".join(stdout_all[:j]) == expected_stdout

        assert os.path.isdir(dir_project)
        os.chdir(dir_project)

        with open("PAIREF_out.log", "r") as stdout_logfile:
            stdout_log = stdout_logfile.readlines()
        assert "run date and time:" in stdout_log[7]
        stdout_log[7] = "run date and time: XX\n"
        assert "user@host:" in stdout_log[8]
        stdout_log[8] = "user@host: XX\n"
        assert "".join(stdout_log[:j]) == expected_stdout

        calculation_ended = None
        for line in stdout_all:
            if "Calculation ended" in line:
                calculation_ended = True
        assert calculation_ended

        res_i_name = res_i.replace(".", "-")
        files_out = \
            [xyzin,
             hklin,
             hklin_unmerged,
             "PAIREF_cutoff.txt",
             "PAIREF_out.log",
             "PAIREF_" + project + ".html",
             "styles.css",
             project + "_CCfree.png",
             project + "_CC.png",
             project + "_CCwork.png",
             project + "_Comp_Mult.png",
             project + "_Intensities.png",
             project + "_merging_stats.csv",
             project + "_No_reflections.png",
             project + "_No_work_free_reflections.png",
             project + "_Optical_resolution.csv",
             project + "_Optical_resolution.png",
             project + "_R" + flag + "_Rgap.csv",
             project + "_R" + flag + "_R-values.csv",
             project + "_Rfree.png",
             project + "_Rgap.csv",
             project + "_Rgap.png",
             project + "_Rmerge_Rmeas_Rpim.png",
             project + "_R-values.csv",
             project + "_R-values.png",
             project + "_Rwork.png"]
        if "--refmac" in others:
            files_out += \
                [project + "_R" + flag + "_" + res_i_name + "A_comparison_at_" + res_i_name + "A.log",
                 project + "_R" + flag + "_" + res_i_name + "A_comparison_at_" + res_i_name + "A.mtz",
                 project + "_R" + flag + "_" + res_i_name + "A.csv",
                 project + "_R" + flag + "_" + res_i_name + "A.log",
                 project + "_R" + flag + "_" + res_i_name + "A.mmcif",
                 project + "_R" + flag + "_" + res_i_name + "A.mtz",
                 project + "_R" + flag + "_" + res_i_name + "A.pdb"]
        elif "--phenix" in others:
            files_out += \
                [project + "_R" + flag + "_" + res_i_name + "A.csv",
                 project + "_R" + flag + "_" + res_i_name + "A_001.log",
                 project + "_R" + flag + "_" + res_i_name + "A_001.cif",
                 project + "_R" + flag + "_" + res_i_name + "A_001.mtz",
                 project + "_R" + flag + "_" + res_i_name + "A_001.pdb"]
        for f in files_out:
            assert os.path.isfile(f)

        assert filecmp.cmp(project + "_R-values.csv",
                           project + "_R" + flag + "_R-values.csv")
        assert filecmp.cmp(project + "_Rgap.csv",
                           project + "_R" + flag + "_Rgap.csv")
        if "refmac" in others:
            assert filecmp.cmp(project + "_R" + flag + "_" + res_i_name + "A.log",
                               project + "_R" + flag + "_" + res_i_name + "A_comparison_at_" + res_i_name + "A.log")
            assert filecmp.cmp(project + "_R" + flag + "_" + res_i_name + "A.mtz",
                               project + "_R" + flag + "_" + res_i_name + "A_comparison_at_" + res_i_name + "A.mtz")
        with open("PAIREF_" + project + ".html", "r") as htmlfile:
            htmlfile_content = htmlfile.read()
        assert "Calculations are still in progress" not in htmlfile_content
        assert "</html>" in htmlfile_content

        # Clean up a bit
        os.chdir("..")
        for f in files.values():
            if os.path.isfile(f):
                os.remove(f)
        # shutil.rmtree(dir_project)

else:
    warnings.warn(UserWarning("Integration test could not be performed.\n"
                              "refmac5: Command not found."))
