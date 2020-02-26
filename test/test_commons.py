import pytest
import platform
import sys
import os
import tempfile
import shutil
from pairef.commons import twodec, twodecname, fourdec, extract_from_file
from helper import run, config


@pytest.mark.parametrize(["value_in", "value_twodec",
                          "value_twodecname", "value_fourdec"],
                         [(2, "2.00", "2-00", "2.0000"),
                          (0.617, "0.62", "0-62", "0.6170")])
def test_twodec_twodecname_fourdec(value_in, value_twodec, value_twodecname,
                                   value_fourdec):
    assert twodec(value_in) == value_twodec
    assert twodecname(value_in) == value_twodecname
    assert fourdec(value_in) == value_fourdec


filename = config("lysozyme_1-40A.log")
filename_fake = config("foo")
searched_stats = "Resolution limits"
searched_Rfree = "Free R factor"
table_stats_i = """Resolution limits                    =     39.355     1.402
Number of used reflections           =      21761
Percentage observed                  =    97.1541
Percentage of free reflections       =     5.0028
Overall R factor                     =     0.2103
Free R factor                        =     0.2236
Average Fourier shell correlation    =     0.9155
AverageFree Fourier shell correlation=     0.9072
Overall weighted R factor            =     0.2013
Free weighted R factor               =     0.2116
Overall weighted R2 factor           =     0.2598
Free weighted R2 factor              =     0.2707
Average correlation coefficient      =     0.8746
Overall correlation coefficient      =     0.9484
Free correlation coefficient         =     0.9443
Cruickshanks DPI for coordinate error=     0.0750
DPI based on free R factor           =     0.0714
Overall figure of merit              =     0.8166
ML based su of positional parameters =     0.0512
ML based su of thermal parameters    =     1.3238
"""
table_stats_f1 = """Resolution limits                    =     39.355     1.402
Number of used reflections           =      21761
Percentage observed                  =    97.1541
Percentage of free reflections       =     5.0028
"""
table_stats_f2 = """Overall R factor                     =     0.2094
Free R factor                        =     0.2220
Average Fourier shell correlation    =     0.9245
AverageFree Fourier shell correlation=     0.9165
Overall weighted R factor            =     0.2007
Free weighted R factor               =     0.2115
Overall weighted R2 factor           =     0.2580
Free weighted R2 factor              =     0.2746
Average correlation coefficient      =     0.8794
Overall correlation coefficient      =     0.9487
Free correlation coefficient         =     0.9438
Cruickshanks DPI for coordinate error=     0.0746
DPI based on free R factor           =     0.0709
Overall figure of merit              =     0.8337
ML based su of positional parameters =     0.0484
ML based su of thermal parameters    =     1.2421
"""
table_stats_f = table_stats_f1 + table_stats_f2
Rfree_i = "0.2236"
Rfree_f = "0.2220"
file_not_found = """ERROR: File """ + str(filename_fake) + """ was not found.
Aborting.
"""


@pytest.mark.parametrize(["filename", "searched", "skip_lines", "n_lines",
                          "nth_word", "not_found", "get_first", "text_test",
                          "returncode"],
                         [(filename, searched_stats, 0, 20,
                           False, "stop", False, table_stats_f, 0),
                          (filename, searched_stats, 0, 20,
                           False, "stop", True, table_stats_i, 0),
                          (filename, searched_stats, 4, 16,
                           False, "stop", False, table_stats_f2, 0),
                          (filename, searched_stats, 5, 1,
                           4, "stop", True, Rfree_i, 0),
                          (filename, searched_stats, 5, 1,
                           4, "stop", False, Rfree_f, 0),
                          (filename, searched_stats, 5, 1,
                           999, "N/A", False, "N/A", 0),
                          (filename, searched_Rfree, 0, 1,
                           4, "stop", True, Rfree_i, 0),
                          (filename, searched_Rfree, 0, 1,
                           4, "stop", False, Rfree_f, 0),
                          (filename, searched_Rfree, 0, 1,
                           999, "N/A", False, "N/A", 0),
                          # (filename_fake, searched_stats, 5, 1,
                          #  4, "stop", False, file_not_found, 1),
                          (filename, "foooooo", 0, 2,
                           False, "N/A", False, "N/A", 0),  # ??? ["N/A"]
                          (filename, "foooooo", 0, 2,
                           5, "N/A", False, "N/A", 0)],
                         ids=["table_stats_final_all",
                              "table_stats_init_all",
                              "table_stats_final_skip",
                              "table_stats_Rfree_init",
                              "table_stats_Rfree_final",
                              "table_stats_Rfree_final_N/A_index_error",
                              "Rfree_init",
                              "Rfree_final",
                              "Rfree_final_N/A_index_error",
                              # "file_not_found_abort",
                              "searched_not_found_N/A_lines",
                              "searched_not_found_N/A_word"])
def test_extract_from_file(filename, searched, skip_lines, n_lines,
                           nth_word, not_found, get_first,
                           text_test, returncode):
    text = extract_from_file(filename=filename, searched=searched,
                             skip_lines=skip_lines, n_lines=n_lines,
                             nth_word=nth_word, not_found=not_found,
                             get_first=get_first)
    if returncode == 0:
        assert "".join(text) == text_test
