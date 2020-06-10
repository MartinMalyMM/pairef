import pytest
import platform
import sys
import os
import tempfile
import shutil
from helper import run, config, tmp_environ
from pairef.preparation import create_workdir, which


def test_python_version():
    assert int(platform.python_version_tuple()[0]) == 2
    assert int(platform.python_version_tuple()[1]) >= 7


RES_SHELLS_GOOD = "1.9,1.8,1.7,1.6,1.5,1.4"

prog = 'cctbx.python -m pairef'
help = \
"""usage: """ + prog + """ --XYZIN XYZIN --HKLIN HKLIN [-u HKLIN_UNMERGED]
                              [--LIBIN LIBIN] [--TLSIN TLSIN] [-c COMIN]
                              [-p PROJECT] [-r RES_SHELLS] [-n N_SHELLS]
                              [-s STEP] [-i RES_INIT] [-f FLAG] [-w WEIGHT]
                              [--ncyc NCYC] [--constant-grid] [--complete]
                              [--TLS-ncyc TLS_NCYC] [--TLSIN-keep] [-h]
                              [--prerefinement-ncyc PREREFINEMENT_NCYC]
                              [--prerefinement-reset-bfactor]
                              [--prerefinement-add-to-bfactor ADD_TO_BFACTOR]
                              [--prerefinement-set-bfactor SET_BFACTOR]
                              [--prerefinement-shake-sites [SHAKE_SITES]]
                              [--prerefinement-no-modification]"""


def test_no_args(tmp_environ):
    cp = run("")
    # assert not cp.stdout
    # assert cp.stderr == help + "\n" + prog + \
    #     ": error: argument --XYZIN/--xyzin is required\n"
    assert cp.stderr == "error: argument --XYZIN/--xyzin is required\n"


def test_file_not_exists(tmp_environ):
    cp = run("--HKLIN " + str(config("foo.mtz")) + ""
             " --XYZIN " + str(config("lysozyme_arp_model_2A.pdb")) + ""
             " -i 2" + ""
             " -r " + RES_SHELLS_GOOD + ""
             " -p file_not_exists")
    assert cp.returncode != 0
    # assert not cp.stdout
    # assert cp.stderr == help + "\ncctbx.python -m pairef: error: " \
    #     "The file " + str(config("foo.mtz")) + " does not exist!\n"
    assert cp.stderr == "error: " \
        "The file " + str(config("foo.mtz")) + " does not exist!\n"


error_res_shells_not_floats = "ERROR: Explicit definition of high " \
    "resolution shells (option -r) is not correct. " \
    "Values must be divided using commas without " \
    "any spaces (e.g. 2.1,2.0,1.9).\nAborting.\n"
# error_res_init_not_float = help + "\n" + prog + \
#     ": error: argument -i: invalid float value: 'abc'\n"
error_res_init_not_float = "error: argument -i: invalid float value: 'abc'\n"
error_only_one_value = "ERROR: Explicit definition of high " \
    "resolution shells (option -r) is not correct. " \
    "Values must be divided using commas without " \
    "any spaces (e.g. 2.1,2.0,1.9). " \
    "More values (than only one) must be provided as the first " \
    "value is the conservative diffraction limit.\nAborting.\n"
error_invalid_order = "ERROR: Explicit definition of high " \
    "resolution shells (option -r) is not correct. " \
    "Values must be set in the decreasing order " \
    "(e.g. 2.1,2.0,1.9).\nAborting.\n"
error_invalid_res_init1 = "ERROR: Explicit definition of high " \
    "resolution shells (option -r) is not correct. " \
    "Resolution of the first shell (2.20 A), " \
    "which was explicitely defined, is lower than or the same as " \
    "the resolution of data which were used for refinement of the " \
    "input structure model (2.10 A).\nAborting.\n"
error_invalid_res_init2 = "ERROR: Explicit definition of high " \
    "resolution shells (option -r) is not correct. " \
    "Resolution of the first shell (2.10 A), " \
    "which was explicitely defined, is lower than or the same as " \
    "the resolution of data which were used for refinement of the " \
    "input structure model (2.10 A).\nAborting.\n"


# This is testing the fuction def_res_shells
@pytest.mark.parametrize(
    ["res_init", "res_shells", "error_message"],
    # # [  # #("2", "abc", error_res_shells_not_floats),
    [("abc", "2.2,2,1.8", error_res_init_not_float)],
    # ("2.0", error_message_only_one_value),
    # #("2.1", "1.9,3.5,1.8", error_invalid_order),
    # #("2.1", "2.2,2,1.9", error_invalid_res_init1),
    # #("2.1", "2.1,2.0,1.9", error_invalid_res_init2)],
    # # ],
    # # ids=["res_shells_not_floats",
    ids=["res_init_not_float"])
    # #"invalid_order",
    # #"invalid_res_init1",
    # #"invalid_res_init2"
    # # ])
def test_invalid_res_shells_setting(tmp_environ,
                                    res_init, res_shells, error_message):
    if os.path.isdir("pairef_bad_res_shells_setting"):
        # Preparation - clean rests from previously examined test that
        # has been interrupted
        shutil.rmtree("pairef_bad_res_shells_setting")
    cp = run("--HKLIN " + str(config("lysozyme_l3600s_1-4A.mtz")) + ""
             " --XYZIN " + str(config("lysozyme_arp_model_2A.pdb")) + ""
             " -i " + res_init + ""
             " -r " + res_shells + ""
             " -p bad_res_shells_setting" + ""
             " -t")
    assert cp.returncode != 0
    # assert not cp.stdout
    assert cp.stderr == error_message


def test_valid_args(tmp_environ):
    if os.path.isdir("pairef_valid_args"):
        # Preparation - clean rests from previously examined test that
        # has been interrupted
        shutil.rmtree("pairef_valid_args")
    arguments = "--HKLIN " + str(config("lysozyme_l3600s_1-4A.mtz")) + \
        " --XYZIN " + str(config("lysozyme_arp_model_2A.pdb")) + \
        " -i 2" + \
        " -r " + RES_SHELLS_GOOD + \
        " -p valid_args" + \
        " -t"
    cp = run(arguments)
    expected_stdout_1 = """
 __   _  ___ __  ___ ___
 )_) /_)  )  )_) )_  )_
/   / / _(_ / \ (__ (

automatic PAIRed REFinement protocol
"""
# run date and time: ...
    expected_stdout_2 = """
Please reference: "Paired refinement under the control of PAIREF"
M. Maly, K. Diederichs, J. Dohnalek, P. Kolenko (2020) IUCrJ 7

Command line arguments: """ + arguments + """

Program has been executed with following input parameters:
"""
    expected_stdout_2 += " * XYZIN: " + str(os.path.dirname(__file__)) + "/" \
        "fixtures/lysozyme_arp_model_2A.pdb\n"
    expected_stdout_2 += " * HKLIN: " + str(os.path.dirname(__file__)) + "/" \
        "fixtures/lysozyme_l3600s_1-4A.mtz\n"
    expected_stdout_2 += """ * Project name: valid_args
 * Resolution shells: 1.9,1.8,1.7,1.6,1.5,1.4
 * Light-testing mode (REFMAC5 will not be executed).

"""
    expected_stdout_2 += "Resolution of the merged diffraction data " + \
        str(os.path.dirname(__file__)) + "/" \
        "fixtures/lysozyme_l3600s_1-4A.mtz: 39.35-1.40 A"
    expected_stdout_2 += """
Manual setting of initial high resolution limit will be used: 2.00 A.
High resolution diffraction limits: 1.90 A, 1.80 A, 1.70 A, 1.60 A, 1.50 A, 1.40 A
"""
    expected_stdout_2 += \
        " * Data with FreeRflag set 0 will be excluded during refinement.\n"
    expected_stdout_2 += "\n"
    expected_stdout_2 += "Current working directory: " \
        "" + tempfile.gettempdir() + "/pairef_valid_args\n"
    expected_stdout_2 += \
        "------> RESULTS AND THE CURRENT STATUS OF CALCULATIONS ARE LISTED" \
        " IN A HTML LOG FILE " + tempfile.gettempdir() + \
        "/pairef_valid_args/PAIREF_valid_args.html\n\n"
    stdout_all = cp.stdout.splitlines(True)
    stdout_1 = "".join(stdout_all[:6])
    stdout_2 = "".join(stdout_all[9:])
    assert stdout_1 == expected_stdout_1
    assert stdout_2 == expected_stdout_2
    assert os.path.isdir("pairef_valid_args")
    shutil.rmtree("pairef_valid_args")
    assert cp.returncode == 0
    assert not cp.stderr


def test_create_workdir(tmp_environ):
    # Preparation - clean rests from previously examined test that has been
    # interrupted
    if os.path.isdir("pairef_create_workdir"):
        os.rmdir("pairef_create_workdir")
    if os.path.isdir("pairef_create_workdir_new"):
        os.rmdir("pairef_create_workdir_new")
    # Test
    create_workdir("create_workdir")
    assert os.path.isdir("pairef_create_workdir")
    create_workdir("create_workdir")
    assert os.path.isdir("pairef_create_workdir_new")
    os.rmdir("pairef_create_workdir")
    os.rmdir("pairef_create_workdir_new")


def test_which():
    result = which("python")
    assert result
    result = which("ThisCommandShouldNotExist")
    assert not result
