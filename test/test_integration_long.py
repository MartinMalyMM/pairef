import pytest
import os
import shutil
import filecmp
import warnings
from helper import run, config, tmp_environ
from pairef.preparation import which
from pairef import __version__


if which("refmac5"):
    @pytest.mark.parametrize(  # res_i needs two decimals
    "project, xyzin, hklin, hklin_unmerged, res_i, res_r",
    [("mdm2_0-10A", "mdm2_1-60A.pdb", "mdm2_merged.mtz", "mdm2_unmerged.mtz", "1.60", "1.55,1.5,1.45,1.4,1.35,1.3"),
     ("CDO_0-10A", "CDO_2B5H_edit_refmac1.pdb", "CDO_R.mtz", "CDO_XDS_ASCII.HKL", "2.00", "1.9,1.8,1.7,1.6,1.5,1.42"),
     ("POLI_0-10A", "poli67_edit12_refmac1.pdb", "poli67_R.mtz", "poli67_XDS_ASCII.HKL", "2.30", "2.2,2.1,2.0,1.9")])
    def test_integration_0_10A_refmac(project, xyzin, hklin,
                                      hklin_unmerged, res_i, res_r, tmp_environ):
        print("\nPerforming integration tests, please wait, it can take "
              "sevelar hours...\n")
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
        urllib_wget(url_prefix + xyzin, xyzin)
        urllib_wget(url_prefix + hklin, hklin)
        urllib_wget(url_prefix + hklin_unmerged, hklin_unmerged)

        arguments = "--HKLIN " + hklin + \
            " --XYZIN " + xyzin + \
            " -u " + hklin_unmerged + \
            " -i " + res_i + \
            " -r " + res_r + \
            " -p " + project + \
            " --refmac"
        cp = run(arguments)
        assert cp.returncode == 0
        # assert not cp.stderr

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
 * Refinement software: REFMAC5
"""
        expected_stdout += " * XYZIN: " + xyzin + "\n"
        expected_stdout += " * HKLIN: " + hklin + "\n"
        expected_stdout += " * HKLIN unmerged: " + hklin_unmerged + "\n"
        expected_stdout += " * Project name: " + project + "\n"
        expected_stdout += " * Resolution shells: " + res_r + "\n"

        stdout_all = cp.stdout.splitlines(True)
        assert "run date and time:" in stdout_all[7]
        stdout_all[7] = "run date and time: XX\n"
        assert "user@host:" in stdout_all[8]
        stdout_all[8] = "user@host: XX\n"
        assert "".join(stdout_all[:22]) == expected_stdout

        assert os.path.isdir(dir_project)
        os.chdir(dir_project)

        with open("PAIREF_out.log", "r") as stdout_logfile:
            stdout_log = stdout_logfile.readlines()
        assert "run date and time:" in stdout_log[7]
        stdout_log[7] = "run date and time: XX\n"
        assert "user@host:" in stdout_log[8]
        stdout_log[8] = "user@host: XX\n"
        assert "".join(stdout_log[:22]) == expected_stdout

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
             project + "_R00_" + res_i_name + "A_comparison_at_" + res_i_name + "A.log",
             project + "_R00_" + res_i_name + "A_comparison_at_" + res_i_name + "A.mtz",
             project + "_R00_" + res_i_name + "A.csv",
             project + "_R00_" + res_i_name + "A.log",
             project + "_R00_" + res_i_name + "A.mmcif",
             project + "_R00_" + res_i_name + "A.mtz",
             project + "_R00_" + res_i_name + "A.pdb",
             project + "_R00_Rgap.csv",
             project + "_R00_R-values.csv",
             project + "_Rfree.png",
             project + "_Rgap.csv",
             project + "_Rgap.png",
             project + "_Rmerge_Rmeas_Rpim.png",
             project + "_R-values.csv",
             project + "_R-values.png",
             project + "_Rwork.png"]
        for f in files_out:
            assert os.path.isfile(f)

        assert filecmp.cmp(project + "_R-values.csv",
                           project + "_R00_R-values.csv")
        assert filecmp.cmp(project + "_Rgap.csv",
                           project + "_R00_Rgap.csv")
        assert filecmp.cmp(project + "_R00_" + res_i_name + "A.log",
                           project + "_R00_" + res_i_name + "A_comparison_at_" + res_i_name + "A.log")
        assert filecmp.cmp(project + "_R00_" + res_i_name + "A.mtz",
                           project + "_R00_" + res_i_name + "A_comparison_at_" + res_i_name + "A.mtz")
        with open("PAIREF_" + project + ".html", "r") as htmlfile:
            htmlfile_content = htmlfile.read()
        assert "Calculations are still in progress" not in htmlfile_content
        assert "</html>" in htmlfile_content

        # Clean up ?
        os.chdir("..")
        for input_file in ["mdm2_1-60A.pdb", "mdm2_merged.mtz", "mdm2_unmerged.mtz"]:
            if os.path.isfile(input_file):
                os.remove(input_file)
        # shutil.rmtree(dir_project)

else:
    warnings.warn(UserWarning("Integration test could not be performed.\n"
                              "refmac5: Command not found."))
