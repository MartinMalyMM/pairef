import pytest
import os
import sys
import shutil
import platform
import filecmp
import warnings
from helper import run, config, tmp_environ
from pairef.preparation import which
from pairef import __version__


if "CCP4" in os.environ:
    if platform.python_version_tuple()[0] == "3":  # CCP4 8
        import sysconfig
        site_packages = (sysconfig.get_path('purelib'))
        demo_data = os.path.join(site_packages, 'ccp4i2', 'demo_data')
        hklin_unmerged = os.path.join(demo_data, "mdm2", "mdm2_unmerged.mtz")
        # python_path = "python" + str(platform.python_version_tuple()[0]) + \
        #     "." + str(platform.python_version_tuple()[1])
        # hklin_unmerged = os.path.join(os.environ.get("CCP4"), "lib",
        #     python_path, "site-packages/ccp4i2/demo_data/mdm2/mdm2_unmerged.mtz")
    elif platform.python_version_tuple()[0] == "2":  # CCP4 7
        hklin_unmerged = os.path.join(os.environ.get("CCP4"),
            "share/ccp4i2/demo_data/mdm2/mdm2_unmerged.mtz")

dir_project = "pairef_test_quick_mdm2"
if which("refmac5") and os.path.isfile(hklin_unmerged):
    def test_integration_quick_mdm2(tmp_environ):
        print("\nPerforming a quick integration test using mdm2 demo data"
              ", please wait...\n")
        if os.path.isdir(dir_project):
            # Preparation - clean rests from previously examined test that
            # has been interrupted
            shutil.rmtree(dir_project)
        arguments = "--HKLIN " + str(config("mdm2_merged.mtz")) + \
            " --XYZIN " + str(config("mdm2_1-60A.pdb")) + \
            " -u " + hklin_unmerged + " -i 1.6 -r 1.55,1.5 -p test_quick_mdm2"
        cp = run(arguments)
        assert cp.returncode == 0
        # assert not cp.stderr

        expected_stdout_1 = """
 __   _  ___ __  ___ ___
 )_) /_)  )  )_) )_  )_
/   / / _(_ / \ (__ (

automatic PAIRed REFinement protocol
version: """ + str(__version__) + "\n"
        # run date and time: ...

        expected_stdout_2 = """
Please cite: "Paired refinement under the control of PAIREF"
M. Maly, K. Diederichs, J. Dohnalek, P. Kolenko (2020) IUCrJ 7

Command line arguments: """ + arguments + """

Program has been executed with following input parameters:
"""
        expected_stdout_2 += " * XYZIN: " + str(config("mdm2_1-60A.pdb")) + "\n"
        expected_stdout_2 += " * HKLIN: " + str(config("mdm2_merged.mtz")) + "\n"
        expected_stdout_2 += " * HKLIN unmerged: " + hklin_unmerged + "\n"
        expected_stdout_2 += """ * Project name: test_quick_mdm2
 * Resolution shells: 1.55,1.5

"""
        expected_stdout_2 += "Resolution of the merged diffraction data " + \
            str(config("mdm2_merged.mtz")) + ": 61.93-1.25 A\n"
        expected_stdout_2 += "Resolution of the unmerged diffraction data " + \
            hklin_unmerged + ": 61.90-1.24 A"
        expected_stdout_2 += """
Manual setting of initial high resolution limit will be used: 1.60 A.
High resolution diffraction limits: 1.55 A, 1.50 A
"""
        expected_stdout_2 += \
        " * Data with FreeRflag set 0 will be excluded during refinement.\n"

        expected_stdout_3 = """Refinement using REFMAC5:

   * Calculating initial statistics at 1.60 A resolution...
       Collecting statistics from logfiles...

   * Refining using data up to 1.55 A resolution...
       Running command:
       refmac5 HKLIN mdm2_merged.mtz XYZIN test_quick_mdm2_R00_1-60A.pdb HKLOUT test_quick_mdm2_R00_1-55A.mtz XYZOUT test_quick_mdm2_R00_1-55A.pdb LIBOUT test_quick_mdm2_R00_1-55A.cif
       Calculating statistics of the refined structure model... . . .
       Collecting statistics from logfiles...
       Updating graphs...
       Preliminary suggested cutoff: 1.XX A

   * Refining using data up to 1.50 A resolution...
       Running command:
       refmac5 HKLIN mdm2_merged.mtz XYZIN test_quick_mdm2_R00_1-55A.pdb HKLOUT test_quick_mdm2_R00_1-50A.mtz XYZOUT test_quick_mdm2_R00_1-50A.pdb LIBOUT test_quick_mdm2_R00_1-50A.cif
       Calculating statistics of the refined structure model... . . . .
       Collecting statistics from logfiles...
       Updating graphs...
       Preliminary suggested cutoff: 1.XX A

     * Calculating merging statistics... . . . . . . . .
       Using labels=I,SIGI

Suggested cutoff: 
"""
        stdout_all = cp.stdout.splitlines(True)
        assert "Preliminary suggested cutoff:" in stdout_all[43]
        stdout_all[43] = "       Preliminary suggested cutoff: 1.XX A\n"
        assert "Preliminary suggested cutoff:" in stdout_all[51]
        stdout_all[51] = "       Preliminary suggested cutoff: 1.XX A\n"
        stdout_1 = "".join(stdout_all[:7])
        stdout_2 = "".join(stdout_all[9:27])
        stdout_3 = "".join(stdout_all[32:57])
        assert stdout_1 == expected_stdout_1
        assert stdout_2 == expected_stdout_2
        assert stdout_3 == expected_stdout_3
        assert "1.60 A" in stdout_all[57]

        assert os.path.isdir(dir_project)
        os.chdir(dir_project)

        with open("PAIREF_out.log", "r") as stdout_logfile:
            stdout_log = stdout_logfile.readlines()
        assert "Preliminary suggested cutoff:" in stdout_log[43]
        stdout_log[43] = "       Preliminary suggested cutoff: 1.XX A\n"
        assert "Preliminary suggested cutoff:" in stdout_log[51]
        stdout_log[51] = "       Preliminary suggested cutoff: 1.XX A\n"
        print("STDOUT:\n" + cp.stdout)
        stdout_log_1 = "".join(stdout_log[:7])
        stdout_log_2 = "".join(stdout_log[9:27])
        stdout_log_3 = "".join(stdout_log[32:57])
        assert stdout_log_1 == expected_stdout_1
        assert stdout_log_2 == expected_stdout_2
        assert stdout_log_3 == expected_stdout_3
        assert "1.60 A" in stdout_log[57]

        files_out = \
            ["mdm2_1-60A.pdb",
             "mdm2_merged.mtz",
             "mdm2_unmerged.mtz",
             "PAIREF_cutoff.txt",
             "PAIREF_out.log",
             "PAIREF_test_quick_mdm2.html",
             "styles.css",
             "test_quick_mdm2_1-60A.pdb",
             "test_quick_mdm2_CCfree.png",
             "test_quick_mdm2_CC.png",
             "test_quick_mdm2_CCwork.png",
             "test_quick_mdm2_Comp_Mult.png",
             "test_quick_mdm2_Intensities.png",
             "test_quick_mdm2_merging_stats.csv",
             "test_quick_mdm2_No_reflections.png",
             "test_quick_mdm2_No_work_free_reflections.png",
             "test_quick_mdm2_Optical_resolution.csv",
             "test_quick_mdm2_Optical_resolution.png",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-50A.log",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-50A.mtz",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-55A.log",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-55A.mtz",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-55A_prev_pair.log",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-55A_prev_pair.mtz",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-60A.log",
             "test_quick_mdm2_R00_1-50A_comparison_at_1-60A.mtz",
             "test_quick_mdm2_R00_1-50A.csv",
             "test_quick_mdm2_R00_1-50A.log",
             "test_quick_mdm2_R00_1-50A.mmcif",
             "test_quick_mdm2_R00_1-50A.mtz",
             "test_quick_mdm2_R00_1-50A.pdb",
             # "test_quick_mdm2_R00_1-50A_sfcheck.log",
             # "test_quick_mdm2_R00_1-50A_sfcheck.out",
             # "test_quick_mdm2_R00_1-50A_sfcheck.ps",
             # "test_quick_mdm2_R00_1-50A_sfcheck.xml",
             "test_quick_mdm2_R00_1-50A_stats_vs_cycle.png",
             "test_quick_mdm2_R00_1-55A_comparison_at_1-55A.log",
             "test_quick_mdm2_R00_1-55A_comparison_at_1-55A.mtz",
             "test_quick_mdm2_R00_1-55A_comparison_at_1-60A.log",
             "test_quick_mdm2_R00_1-55A_comparison_at_1-60A.mtz",
             "test_quick_mdm2_R00_1-55A_comparison_at_1-60A_prev_pair.log",
             "test_quick_mdm2_R00_1-55A_comparison_at_1-60A_prev_pair.mtz",
             "test_quick_mdm2_R00_1-55A.csv",
             "test_quick_mdm2_R00_1-55A.log",
             "test_quick_mdm2_R00_1-55A.mmcif",
             "test_quick_mdm2_R00_1-55A.mtz",
             "test_quick_mdm2_R00_1-55A.pdb",
             # "test_quick_mdm2_R00_1-55A_sfcheck.log",
             # "test_quick_mdm2_R00_1-55A_sfcheck.out",
             # "test_quick_mdm2_R00_1-55A_sfcheck.ps",
             # "test_quick_mdm2_R00_1-55A_sfcheck.xml",
             "test_quick_mdm2_R00_1-55A_stats_vs_cycle.png",
             "test_quick_mdm2_R00_1-60A_comparison_at_1-60A.log",
             "test_quick_mdm2_R00_1-60A_comparison_at_1-60A.mtz",
             "test_quick_mdm2_R00_1-60A.csv",
             "test_quick_mdm2_R00_1-60A.log",
             "test_quick_mdm2_R00_1-60A.mmcif",
             "test_quick_mdm2_R00_1-60A.mtz",
             "test_quick_mdm2_R00_1-60A.pdb",
             # "test_quick_mdm2_R00_1-60A_sfcheck.log",
             # "test_quick_mdm2_R00_1-60A_sfcheck.out",
             # "test_quick_mdm2_R00_1-60A_sfcheck.ps",
             # "test_quick_mdm2_R00_1-60A_sfcheck.xml",
             "test_quick_mdm2_R00_Rgap.csv",
             "test_quick_mdm2_R00_R-values.csv",
             "test_quick_mdm2_Rfree.png",
             "test_quick_mdm2_Rgap.csv",
             "test_quick_mdm2_Rgap.png",
             "test_quick_mdm2_Rmerge_Rmeas_Rpim.png",
             "test_quick_mdm2_R-values.csv",
             "test_quick_mdm2_R-values.png",
             "test_quick_mdm2_Rwork.png"]
        for f in files_out:
            assert os.path.isfile(f)

        # Rfree decrease?
        # with open("test_quick_mdm2_R-values.csv", "r") as csvfile:
        #     csvfile_lines = csvfile.readlines()
        #     lines = [1, 2]
        #     for line in lines:
        #         assert float(csvfile_lines[line].split()[6]) < 0
        # But CC* is lower than CCwork so suggested cutoff should be 1.60 A
        with open("PAIREF_cutoff.txt", "r") as result:
            assert result.read() == "1.60"

        assert filecmp.cmp("test_quick_mdm2_R-values.csv",
                           "test_quick_mdm2_R00_R-values.csv")
        assert filecmp.cmp("test_quick_mdm2_Rgap.csv",
                           "test_quick_mdm2_R00_Rgap.csv")
        assert filecmp.cmp("test_quick_mdm2_R00_1-60A.log",
                           "test_quick_mdm2_R00_1-60A_comparison_at_1-60A.log")
        assert filecmp.cmp("test_quick_mdm2_R00_1-60A.mtz",
                           "test_quick_mdm2_R00_1-60A_comparison_at_1-60A.mtz")
        with open("PAIREF_test_quick_mdm2.html", "r") as htmlfile:
            htmlfile_content = htmlfile.read()
        assert "Calculations are still in progress" not in htmlfile_content
        assert "</html>" in htmlfile_content

        # Clean up
        os.chdir("..")
        shutil.rmtree(dir_project)
        
elif not which("refmac5"):
    warnings.warn(UserWarning("Integration test could not be performed.\n"
                              "refmac5: Command not found."))
else:
    warnings.warn(UserWarning(
        "Integration test could not be performed.\n"
        "Unmerged demo data ccp4i2/demo_data/mdm2_unmerged.mtz "
        "could not be found.\n"
        "For this test, the CCP4 Suite must be installed and the relating "
        "paths must be set correctly."))
