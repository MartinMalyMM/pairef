import pytest
import os
import tempfile
import shutil
from pairef.graphs import matplotlib_bar, matplotlib_line, write_log_html
from helper import run, config, tmp_environ, AttrDict


def test_matplotlib_bar(tmp_environ):
    shutil.copy2(config("A_R-values.csv"), tempfile.gettempdir())
    args = AttrDict()
    args.project = "A"
    args.complete_cross_validation = False
    matplotlib_bar(args, values="R-values")
    assert os.path.isfile("A_R-values.png")
    os.remove("A_R-values.csv")
    os.remove("A_R-values.png")


@pytest.mark.parametrize(["statistics", "title", "filename_suffix", "pngfile"],
                         # [(["Rwork"], r"$\it{R}_{\mathrm{work}}$",
                         [(["Rwork"], "$\it{R}_{\mathrm{work}}$",
                          "Rwork", "NK_Rwork.png"),
                          # (["Rfree"], r"$\it{R}_{\mathrm{free}}$",
                          (["Rfree"], "$\it{R}_{\mathrm{free}}$",
                           "Rfree", "NK_Rfree.png")],
                         ids=["Rwork", "Rfree"])
def test_matplotlib_line_Rvalues(tmp_environ, statistics, title,
                                 filename_suffix, pngfile):
    files_in = ["NK_R00_1-80A.csv", "NK_R00_1-70A.csv", "NK_R00_1-60A.csv",
                "NK_R00_1-50A.csv"]
    for f in files_in:
        shutil.copy2(config(f), tempfile.gettempdir())
    matplotlib_line(shells=(1.8, 1.7, 1.6, 1.5), project="NK",
                    statistics=statistics, n_bins_low=10,
                    title=title, multiscale=False,
                    filename_suffix=filename_suffix)
    assert os.path.isfile(pngfile)
    files = files_in + [pngfile]
    for f in files:
        os.remove(f)


def test_matplotlib_line_Rgap(tmp_environ):
    shutil.copy2(config("NK_Rgap.csv"), tempfile.gettempdir())
    matplotlib_line(shells=(1.8), project="NK",
                    statistics="Rgap", n_bins_low=10,
                    title="$\it{R}_{\mathrm{free}}-"
                    "\it{R}_{\mathrm{work}}$", multiscale=False,
                    # title=r"$\it{R}_{\mathrm{free}}-"
                    # r"\it{R}_{\mathrm{work}}$", multiscale=False,
                    filename_suffix="Rgap")
    assert os.path.isfile("NK_Rgap.png")
    os.remove("NK_Rgap.csv")
    os.remove("NK_Rgap.png")


def test_write_log_html(tmp_environ):
    versions_dict = {"refmac_version": "N/A",
                     "pairef_version": "N/A"}
    args = AttrDict()
    args.xyzin = "fake.pdb"
    args.hklin = "fake.mtz"
    args.hklin_unmerged = "fake.HKL"
    args.libin = "fake.mtz"
    args.comin = "fake.com"
    args.tlsin = "fake.tls"
    args.defin = None
    args.phenix = None
    args.project = "A"
    args.ncyc = 11
    args.weight = None
    args.tlsin_keep = None
    args.tls_ncyc = None
    args.flag = 0
    args.add_to_bfactor = None
    args.complete_cross_validation = False
    args.no_modification = None
    args.prerefinement_ncyc = None
    args.reset_bfactor = None
    args.add_to_bfactor = None
    args.set_bfactor = None
    args.shake_sites = None
    args.constant_grid = None
    shutil.copy2(config("A_R-values.csv"), tempfile.gettempdir())
    shells_ready_with_res_init = (1.8, 1.7, 1.6)
    shells = (1.8, 1.7, 1.6, 1.5)
    flag_sets = range(20)
    htmlfilename_in_progress = write_log_html(
        shells, shells_ready_with_res_init, args, versions_dict, flag_sets)
    assert os.path.isfile(htmlfilename_in_progress)
    with open(htmlfilename_in_progress, "r") as htmlfile_in_progress:
        htmlfile_in_progress_content = htmlfile_in_progress.read()
    assert "Calculations are still in progress" in htmlfile_in_progress_content
    assert "</html>" in htmlfile_in_progress_content
    os.remove(htmlfilename_in_progress)

    htmlfilename_done = write_log_html(
        shells, shells, args, versions_dict, flag_sets, done=True)
    assert os.path.isfile(htmlfilename_done)
    with open(htmlfilename_done, "r") as htmlfile_done:
        htmlfile_done_content = htmlfile_done.read()
    assert "Calculations are still in progress" not in htmlfile_done_content
    assert "</html>" in htmlfile_done_content
    os.remove(htmlfilename_done)
    os.remove("A_R-values.csv")
    # os.remove("styles.css")
