# coding: utf-8
import os
import sys
import matplotlib
matplotlib.use('agg')  # TKinter makes problems, agg should work
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import numpy as np
from collections import namedtuple
import platform
import shutil
import cgi
import warnings
from .commons import twodec, twodecname, fourdec
from .preparation import which
from .settings import warning_dict, date_time


def xticklabels_compress(list, n_max=13, depth=1):
    """ If there are more than `n_max` bins, do not show all the labels
    in list.

    Args:
        list (list): containing labels for a graph
        n_max (int): maximal allowed number of values in list that are not ""
        depth (int)

    Returns:
        list: containing labels for a graph - compressed
    """
    if len(list) > 88:
        return list  # long lists cause problems
    if len([label for label in list if label is not ""]) >= n_max:
        if len(list) % 2 == 0:
            i = int(pow(2, depth - 1) + 1)  # i = 2 in case of depth == 1
            list[1] = ""
        else:
            i = int(pow(2, depth - 1))  # i = 1 in case of depth == 1
        erase = True
        while i < len(list) - 1:
            if erase:
                list[i] = ""
                erase = False
            else:
                erase = True
            i = i + 1 * depth
        list[- depth - 1] = ""
    # If it was not enough, do it again, recursively
    if len([label for label in list if label is not ""]) >= n_max:
        list = xticklabels_compress(list=list, n_max=n_max, depth=2*depth)
    return list


def matplotlib_bar(args, values="R-values", flag_sets=[], ready_shells=[]):
    """Plots and saves a bar chart using `matplotlib`.

    If `flag_sets` is an empty
    list, it is assumed that the values are saved in the file
    `args.project_values.csv`.

    In the other case, a chart showing results
    of complete cross-validation is ploted; the needed values
    are picked from files `project_RXX_values.csv` where `XX` is a number
    of a flag.

    Args:
        args (parser): Input arguments (including e. g. name of the project) \
                       parsed by `argparse` via function process_arguments()
        values (str): expected value: `"R-values"` (not ready yet: \
                      `"CC-values"`)
        flag_sets (list): List of free reflection flag sets (int)
        ready_shells (list)

    Returns:
        str: Name of the PNG file containing the chart.
    """
    # Graph setting
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["legend.loc"] = 'best'
    plt.rcParams["legend.fontsize"] = 'large'
    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')
    # fig.tight_layout()
    bar_width = 0.35
    opacity = 1
    ax.set_ylabel(r'$\Delta \it{R}$', rotation=0, fontsize=14)

    xticklabels_list = []
    values_work_list = []
    values_free_list = []
    errors_work_list = []
    errors_free_list = []

    def pick_work_free_from_csv_line(line, values_work_list, values_free_list,
                                     errors_work_list, errors_free_list,
                                     errors=False):
        if line.lstrip()[0] == "#":  # If it is a comment, do not load data
            continue_sign = True
            return values_work_list, values_free_list, errors_work_list, errors_free_list, continue_sign
        continue_sign = False
        try:
            values_work_list.append(float(line.split()[3]))
        except ValueError:
            values_work_list.append(None)
        try:
            values_free_list.append(float(line.split()[6]))
        except ValueError:
            values_free_list.append(None)
        if errors:  # return also standard error of mean
            try:
                errors_work_list.append(float(line.split()[9]))
            except ValueError:
                errors_work_list.append(float('nan'))            
            try:
                errors_free_list.append(float(line.split()[10]))
            except ValueError:
                errors_free_list.append(float('nan'))
        else:
            errors_work_list.append(float('nan'))
            errors_free_list.append(float('nan'))
        return values_work_list, values_free_list, errors_work_list, errors_free_list, continue_sign

    if flag_sets:
        # Chart for the complete cross-validation
        # (differences of statistics depending on to various FreeRflag sets)
        pngfilename = args.project + "_" + values + "_complete_" \
            "" + twodecname(ready_shells[-1]) + "A.png"
        # Load data - statistics relating to individual free refl. sets
        for flag in flag_sets:
            xticklabel = str(flag)
            xticklabels_list.append(xticklabel)
            csvfilename = args.project + "_R" + str(flag).zfill(2) + "_" \
                "" + values + ".csv"
            with open(csvfilename, "r") as csvfile:
                last_line = csvfile.readlines()[-1]
            values_work_list, values_free_list, errors_work_list, \
                errors_free_list, continue_sign = \
                pick_work_free_from_csv_line(
                    last_line, values_work_list, values_free_list, 
                    errors_work_list, errors_free_list)
        xticklabels_list = xticklabels_compress(xticklabels_list, n_max=21)
        # Count numbers of increases and decreases
        values_work_positive = sum(1 for i in values_work_list if float(i) > 0)
        values_work_negative = sum(1 for i in values_work_list if float(i) < 0)
        values_work_zero = sum(1 for i in values_work_list if float(i) == 0)
        values_free_positive = sum(1 for i in values_free_list if float(i) > 0)
        values_free_negative = sum(1 for i in values_free_list if float(i) < 0)
        values_free_zero = sum(1 for i in values_free_list if float(i) == 0)
        # Load data - average statistics
        csvfilename = args.project + "_" + values + ".csv"
        with open(csvfilename, "r") as csvfile:
            last_line = csvfile.readlines()[-1]
        values_work_list, values_free_list, errors_work_list, \
            errors_free_list, continue_sign = pick_work_free_from_csv_line(
                last_line, values_work_list, values_free_list,
                errors_work_list, errors_free_list, errors=True)
        xticklabels_list.append("avrg")
        # Prepare lists of colors
        color1 = ["#0065BD"] * len(flag_sets) + ["#156570"]
        color2 = ["#6AADE4"] * len(flag_sets) + ["#00B2A9"]
        # Define labels and graph title
        ax.set_xlabel(r'$\mathrm{Free\ reflection\ set}$', fontsize=14)
        if values == "R-values":
            values_work_label = r'$\it{R}_{\mathrm{work}} ( ' \
                '' + str(values_work_positive) + r'\uparrow , ' \
                '' + str(values_work_negative) + r'\downarrow , ' \
                '' + str(values_work_zero) + r'\rightarrow' \
                ')$'
            values_free_label = r'$\it{R}_{\mathrm{free}} ( ' \
                '' + str(values_free_positive) + r'\uparrow , ' \
                '' + str(values_free_negative) + r'\downarrow , ' \
                '' + str(values_free_zero) + r'\rightarrow' \
                ')$'
            values_abb = "R"
        elif values == "CC-values":
            values_work_label = r'CC$_\mathrm{work}$'
            values_free_label = r'CC$_\mathrm{free}$'
            values_abb = "CC"
        graph_title = r"$\mathrm{" + twodec(ready_shells[-2]) + r"\AA}" \
            r"\rightarrow \mathrm{" + twodec(ready_shells[-1]) + r"\AA}$"

    else:
        # Define labels
        if values == "R-values":
            values_work_label = r'$\it{R}_{\mathrm{work}}$'
            values_free_label = r'$\it{R}_{\mathrm{free}}$'
            values_abb = "R"
        elif values == "CC-values":
            values_work_label = r'CC$_\mathrm{work}$'
            values_free_label = r'CC$_\mathrm{free}$'
            values_abb = "CC"
        ax.set_xlabel(r'$\mathrm{Resolution\ step\ (\AA})$', fontsize=14)
        # Chart showing differences of statistics depending on resolution
        csvfilename = args.project + "_" + values + ".csv"
        pngfilename = args.project + "_" + values + ".png"
        # Define graph title and bar colors
        if args.complete_cross_validation:
            graph_title = r'$\mathrm{Differences\ of\ overall\ } ' \
                '' + values_abb + r'$-$\mathrm{values\ averaged\ over\ free\ sets}$'
            color1 = "#156570"
            color2 = "#00B2A9"
            errors = True
        else:
            graph_title = r'$\mathrm{Differences\ of\ overall\ } ' \
                '' + values_abb + r'$-$\mathrm{values}$'
            color1 = "#0065BD"
            color2 = "#6AADE4"
            errors = False
        # Load data
        with open(csvfilename, "r") as csvfile:
            for line in csvfile.readlines():
                values_work_list, values_free_list, errors_work_list, \
                    errors_free_list, continue_sign = \
                    pick_work_free_from_csv_line(
                        line, values_work_list, values_free_list,
                        errors_work_list, errors_free_list, errors)
                if continue_sign:
                    continue
                xticklabel = line.split()[0]
                # xticklabel = xticklabel.replace("A", r"\AA")
                # angstroem units are mentioned in the x-axis label
                xticklabel = xticklabel.replace("A", "")
                xticklabel = xticklabel.replace("->", r"\rightarrow")
                xticklabel = r"$\mathrm{" + xticklabel + "}$"
                xticklabels_list.append(xticklabel)
        xticklabels_list = xticklabels_compress(xticklabels_list, n_max=21)

    n_groups = len(xticklabels_list)
    index = np.arange(n_groups)  # NumPy
    # Plot the chart
    bar1 = ax.bar(index, values_work_list, bar_width,
                  alpha=opacity, color=color1, linewidth=0,
                  label=values_work_label,
                  yerr=errors_work_list, capsize=5, ecolor='orange')
    bar2 = ax.bar(index + bar_width, values_free_list, bar_width,
                  alpha=opacity, color=color2, linewidth=0,
                  label=values_free_label,
                  yerr=errors_free_list, capsize=5, ecolor='orange')  # , hatch="/")
    ax.set_title(graph_title, fontsize=16, y=1.04)
    ax.set_xticks(index + bar_width/2)
    ax.set_xticklabels(xticklabels_list)
    ax.legend()
    if len(xticklabels_list) < 4:
        xticklabels_rotation = 0
    elif len(xticklabels_list) < 9:
        xticklabels_rotation = 30
    else:
        xticklabels_rotation = 90
    plt.setp(ax.get_xticklabels(), rotation=xticklabels_rotation,
             horizontalalignment='center', fontsize=14)
    with warnings.catch_warnings():
        # There is a bug in matplotlib 1.x.x
        # https://github.com/matplotlib/matplotlib/issues/5209
        warnings.simplefilter(action='ignore', category=FutureWarning)
        plt.savefig(pngfilename, bbox_inches="tight", dpi=96)
    plt.clf()
    plt.close('all')
    return pngfilename


def matplotlib_line(shells, project, statistics, n_bins_low, title, flag=0,
                    multiscale=False, filename_suffix="", refinement="refmac"):
    """Plots statistics values (choice by `statistics`)
    `project+"_"+twodecname(shells[*])+"A.csv"`.
    Generate and save plot `project+"_"+statistic+".png` using `matplotlib`.

    Args:
        shells (list): containing `float`
        project (str): Name of the project
        statistics (list): List of names of statistics to be plotted (`str`)
        n_bins_low (int)
        title (str)
        flag (int)
        multiscale (bool): Use 2 different y-axis for data lines
        filename_suffix (str)
        refinement (str): "refmac" or "phenix"

    Returns:
        str:
            Name of the generated file with plot
            (`project+"_"+statistic+".png`)
    """
    # Graph setting
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "cm"
    # Legend setting (Rgap and multiscale graph have another setting)
    legend_loc = "center left"
    plt.rcParams["legend.fontsize"] = 'large'
    legend_bbox_to_anchor = (1.04, 0.5)
    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')
    dpi = 96
    ax.set_xlabel(r'Resolution ($\mathrm{\AA}$)')
    ax.set_ylabel('')
    ax.grid(color='#DCDCDC')
    xticklabels_rotation = 0

    # Load xtickslabels from a model refined up to the highest resolution
    xshell_list = []
    xticklabels_list = []

    # Set filename of output PNG and graph title
    if filename_suffix:
        pngfilename = project + "_" + filename_suffix + ".png"
    else:
        pngfilename = project + "_" + title + ".png"
    graph_title = title

    # Set markers and colors
    markers = Line2D.filled_markers
    markers_cycle = 0
    # Colors are inspired by https://davidmathlogic.com/colorblind/
    # Gray color was added and blue colors were modified using
    # CTU graphical manual
    # colors = ["#0065BD", "#9B9B9B", "#0B672A", "#44AA99", "#6AADE4",
    #           "#DDCC77", "#CC6677", "#AA4499", "#882255", "#000000"]
    colors = ["#808080", "#6AADE4", "#DC267F", "#FE6100", "#0065BD", "#FFB000"]
    colors_cycle = 0
    lns = []  # List that will contain lists of datalines, needed for
              # multiscale graphs

    for statistic in statistics:
        # if "Rwork" in statistics or "Rfree" in statistics
        #    or "CCwork" in statistics or "CCfree" in statistics:
        if statistic == "Rwork" or statistic == "Rfree" or \
                statistic == "CCwork" or statistic == "CCfree":
            values_list_list = []
            if statistic == "Rwork":
                column = 6
                ax.axhline(y=0.42, color='r', linestyle='-')
            if statistic == "Rfree":
                column = 7
                ax.axhline(y=0.42, color='r', linestyle='-')
            if statistic == "CCwork":
                column = 8
            if statistic == "CCfree":
                column = 9
            # graph_title = '$\it{' + statistic.replace(statistic[-4:], "") +
            # '}_{\mathrm{' + statistic[-4:] + '}}$'

            csvfilename = project + "_R" + str(flag).zfill(2) + "_" \
                "" + twodecname(shells[-1]) + "A.csv"
            with open(csvfilename, "r") as csvfile:
                for line in csvfile.readlines():
                    if line.lstrip()[0] == "#":  # If it is a comment,
                        continue                 # do not load data
                    xshell_list.append(int(line.split()[0]))
                    xticklabels_list.append(twodec(float(line.split()[3])))
            xticklabels_list = xticklabels_compress(xticklabels_list)

            # Load statistic relating to the structure models
            for i in range(len(shells)):
                csvfilename = project + "_R" + str(flag).zfill(2) + "_" \
                    "" + twodecname(shells[i]) + "A.csv"
                values_list_list.append([])
                with open(csvfilename, "r") as csvfile:
                    for line in csvfile.readlines():
                        if line.lstrip()[0] == "#":  # If it is a comment,
                            continue                 # do not load data
                        try:
                            values_list_list[i].append(
                                float(line.split()[column]))
                        except ValueError:
                            values_list_list[i].append(None)
                # Make the `markers` and `colors` lists cycled inf.
                if i + 1 - markers_cycle > len(markers):
                    markers_cycle = markers_cycle + len(markers)
                if i + 1 - colors_cycle > len(colors):
                            colors_cycle = colors_cycle + len(colors)
                # values_label = r'~$\mathrm{' + twodec(shells[i]) + r'\ \AA}$'
                values_label = twodec(shells[i]) + r' $\mathrm{\AA}$'

                # Insert `None` into missing values to obey
                # "ValueError: x and y must have same first dimension"
                missing_values = \
                    len(xticklabels_list) - len(values_list_list[i])
                for j in range(missing_values):
                    values_list_list[i].append(None)

                # Vertical line - the conservative high resolution limit
                plt.axvline(x=n_bins_low, color='#9B9B9B', linestyle='--')

                lns += ax.plot(xshell_list, values_list_list[i],
                               label=values_label,
                               marker=markers[i - markers_cycle], markersize=4,
                               markeredgewidth=0,
                               color=colors[i - colors_cycle])
            statistics.remove(statistic)

        elif statistic == "Rgap" or statistic == "res_opt":
            values_list = []
            if statistic == "Rgap":
                csvfilename = project + "_R" + str(flag).zfill(2) + "_Rgap.csv"
                values_column = 3
            else:  # statistic == "res_opt"
                csvfilename = project + "_Optical_resolution.csv"
                values_column = 1
            with open(csvfilename, "r") as csvfile:
                for line in csvfile.readlines():
                    if line.lstrip()[0] == "#":  # If it is a comment,
                        continue                 # do not load data
                    xticklabels_list.append(twodec(float(line.split()[0])))
                    values_list.append(float(line.split()[values_column]))
            xticklabels_list = xticklabels_compress(xticklabels_list)
            xshell_list = range(len(values_list))
            values_label = title
            legend_loc = "best"
            legend_bbox_to_anchor = None
            lns += ax.plot(xshell_list, values_list, label=values_label,
                           marker="o", markersize=4,
                           markeredgewidth=0, color=colors[1])
            statistics.remove(statistic)

        elif statistic == "Rwork_cyc" or statistic == "Rfree_cyc":
            values_list = []
            if statistic == "Rfree_cyc":
                values_column = 2  # for refmac
                values_range = (20, 26)  # for phenix
                values_label = r'$\it{R}_\mathrm{free}$'
                color = "#6AADE4"
            else:  # Rwork
                values_column = 1  # for refmac
                values_range = (27, 33)  # for phenix
                values_label = r'$\it{R}_\mathrm{work}$'
                color = "#0065BD"
            prefix = project + "_R" + str(flag).zfill(2) + "_" + \
                    twodecname(shells[-1]) + "A"
            if refinement == "refmac":
                logfilename = prefix + ".log"
                with open(logfilename, "r") as logfile:
                    lines = logfile.readlines()
                for i in range(len(lines)):
                    if "    Ncyc    Rfact    Rfree     FOM      -LL     " \
                            "-LLfree  rmsBOND  zBOND rmsANGL  zANGL rmsCHIRAL $$" \
                            in lines[i]:
                        j = i + 2
                while lines[j].split()[0] != '$$':
                    values_list.append(float(lines[j].split()[values_column]))
                    if not len(xticklabels_list) == len(xshell_list):
                        xticklabels_list.append(lines[j].split()[0])
                    j = j + 1
            elif refinement == "phenix":
                xticklabels_rotation = 90
                pdbfilename = prefix + "_001.pdb"
                with open(pdbfilename, "r") as pdbfile:
                    lines = pdbfile.readlines()
                for i in range(len(lines)):
                    if "REMARK  stage r-work r-free bonds angles " \
                            "b_min b_max b_ave n_water shift" \
                            in lines[i]:
                        j = i + 1
                while lines[j].split()[-1][-1] != '-':  # until hline ---------
                    values_list.append(float(
                        lines[j][values_range[0]:values_range[1]]))
                    if not len(xticklabels_list) == len(xshell_list):
                        xticklabel = lines[j][7:18]
                        xticklabel.strip()
                        xticklabels_list.append(xticklabel)
                    j = j + 1
            xticklabels_list = xticklabels_compress(xticklabels_list, n_max=21)
            xshell_list = range(len(values_list))
            legend_loc = "best"
            legend_bbox_to_anchor = None
            lns += ax.plot(xshell_list, values_list, label=values_label,
                           marker="o", markersize=5,
                           markeredgewidth=0, color=color)
            ax.set_xlabel(r'Cycle')
            # dpi=64

    if (os.path.isfile(project + "_merging_stats.csv") and statistics) \
            or "n_work" in statistics or "n_free" in statistics:
        if "graph_title" not in locals():
            graph_title = title
        if "pngfilename" not in locals():
            pngfilename = project + "_" + title + ".png"
        if "n_work" in statistics or "n_work" in statistics:
            csvfilename = project + "_R" + str(flag).zfill(2) + "_" + \
                twodecname(shells[-1]) + "A.csv"
            xticklabels_column = 3
        elif os.path.isfile(project + "_merging_stats.csv"):
            csvfilename = project + "_merging_stats.csv"
            xticklabels_column = 2
        if not xshell_list and not xticklabels_list:
            with open(csvfilename, "r") as csvfile:
                for line in csvfile.readlines():
                    if line.lstrip()[0] == "#":  # If it is a comment,
                        continue                 # do not load data
                    xshell_list.append(int(line.split()[0]))
                    xticklabels_list.append(twodec(float(
                        line.split()[xticklabels_column])))
            xticklabels_list = xticklabels_compress(xticklabels_list)
        for i, statistic in enumerate(statistics):
            if statistic == "CC*" and os.path.isfile(csvfilename):
                values_label = r'CC$^*$'
                values_column = 14
                marker = "*"
                color = "g"
                linestyle = '--'
            else:
                marker = markers[i - markers_cycle]
                color = colors[i - colors_cycle]
                linestyle = '-'
                if statistic == "n_obs" and os.path.isfile(csvfilename):
                    values_label = 'No. observed refl.'
                    values_column = 3
                elif statistic == "n_unique" and os.path.isfile(csvfilename):
                    values_label = 'No. unique refl.'
                    values_column = 4
                elif statistic == "n_work" and os.path.isfile(csvfilename):
                    values_label = 'No. work refl.'
                    values_column = 4
                elif statistic == "n_free" and os.path.isfile(csvfilename):
                    values_label = 'No. free refl.'
                    values_column = 5
                elif statistic == "Multiplicity" \
                        and os.path.isfile(csvfilename):
                    values_label = statistic
                    values_column = 5
                elif statistic == "Completeness" \
                        and os.path.isfile(csvfilename):
                    values_label = statistic
                    values_column = 6
                elif statistic == "<I>" and os.path.isfile(csvfilename):
                    values_label = r'<$\it{I}$>'
                    values_column = 7
                elif statistic == "<I/sI>" and os.path.isfile(csvfilename):
                    values_label = r'<$I/\sigma(I)$>'
                    values_column = 8
                elif statistic == "Rmerge" and os.path.isfile(csvfilename):
                    values_label = r'$\it{R}_\mathrm{merge}$'
                    values_column = 9
                elif statistic == "Rmeas" and os.path.isfile(csvfilename):
                    values_label = r'$\it{R}_\mathrm{meas}$'
                    values_column = 10
                elif statistic == "Rpim" and os.path.isfile(csvfilename):
                    values_label = r'$\it{R}_\mathrm{pim}$'
                    values_column = 11
                elif statistic == "CChalf" and os.path.isfile(csvfilename):
                    values_label = r'CC$_\mathrm{1/2}$'
                    values_column = 12
                else:
                    break
            values_list = []
            with open(csvfilename, "r") as csvfile:
                for line in csvfile.readlines():
                    if line.lstrip()[0] == "#":  # If it is a comment
                        continue                 # do not load data
                    try:
                        values_list.append(float(line.split()[values_column]))
                    except ValueError:
                        values_list.append(None)
            values_list = values_list[:len(xshell_list)]
            if len(values_list) < len(xticklabels_list):
                values_list = values_list + \
                    [None]*(len(xticklabels_list) - len(values_list))
            if multiscale and i == 1:
                legend_bbox_to_anchor = (1.12, 0.5)
                ax2 = ax.twinx()
                ax2.tick_params('y', colors=color)
                lns += ax2.plot(xshell_list, values_list, label=values_label,
                                marker=marker, markersize=6,
                                linestyle=linestyle, markeredgewidth=0,
                                color=color)
            else:
                lns += ax.plot(xshell_list, values_list, label=values_label,
                               marker=marker, markersize=6,
                               linestyle=linestyle, markeredgewidth=0,
                               color=color)
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=legend_loc,
              bbox_to_anchor=legend_bbox_to_anchor)
    ax.set_title(graph_title, fontsize=16, y=1.03)
    plt.xticks(xshell_list)
    ax.set_xticklabels(xticklabels_list)
    plt.setp(ax.get_xticklabels(), horizontalalignment='center',
             rotation=xticklabels_rotation)
    plt.savefig(pngfilename, bbox_inches="tight", dpi=dpi)
    plt.clf()
    plt.close('all')
    return pngfilename


def write_log_html(shells, ready_shells, args, versions_dict, flag_sets,
                   res_cur=0, ready_merging_statistics=False, done=False):
    """Created html output log.

    Args:
        shells (list)
        ready_shells (list)
        args
        versions_dict (dict): Dictionary containing keys "refmac_version" and
                              "pairef_version"
        flag_sets (list)
        res_cur (float)
        ready_merging_statistics (bool)
        done (bool)

    Returns:
        str: Name of the created HTML file
    """

    def warning_orangebox(warning_keys, page):
        """If any of warning in `warning_keys` appears, show it on the page
        in yellow box.

        Args:
            warning_keys (list)
            page (str)

        Retruns:
            str
        """
        for warning_key in warning_keys:
            if warning_key in warning_dict:
                page += '\t\t<div class="orangebox">\n\t\t\t'
                page += warning_dict[warning_key]
                page += '\n\t\t</div>\n'
        return page

    import getpass
    import socket


    page = """<!DOCTYPE html>
<head>
    <title>PAIREF - results """ + args.project + """</title>
    <link rel="stylesheet" type="text/css" href="styles.css">\n"""
    if not done:
        page += """\t<meta http-equiv="refresh" content="15">\n"""
    page += """</head>
<body>
    <div id="header">
        <div id="title">PAIREF</div>
        <div id="subtitle">Automatic PAIRed REFinement protocol</div>
    </div>

    <h1>Results - """ + args.project + """</h1>\n"""

    if done:
        page += """\t<div class="greenbox">Calculations ended""" \
            """ successfully.</div>\n"""
    else:
        page += """\t<div class="yellowbox">\n"""
        page += """\t\t<span class="progress">""" \
            """Calculations are still in progress...</span> """
        if len(ready_shells) >= 2:
            page += "Displaying results for the " \
                "following high resolution diffraction limits: "
            for shell in ready_shells:
                page += twodec(shell) + " &#8491;, "
            page = page[:-2]  # erase the last comma and space
        page += """\n\t\t<span class="reload">""" \
            """<a href="javascript:window.location.reload(true)">""" \
            """REFRESH</a></span>\n"""
        page += "\t</div>\n"

    page += "\t<h2>Input parameters</h2>\n"
    page += "\t\t<table>\n"
    page += "\t\t<tr><td>Project:</td><td>" + args.project + "</td></tr>\n"
    page += "\t\t<tr><td>High resolution diffraction limits:</td><td>"
    for shell in shells:
        page += twodec(shell) + " &#8491;, "
    page = page[:-2]  # erase the last comma and space
    page += "</td></tr>\n"
    page += "\t\t<tr><td>Initial structure model:</td><td>" \
        "" + args.xyzin + "</td></tr>\n"
    page += "\t\t<tr><td>Merged diffraction data:</td><td>" \
        "" + args.hklin + "</td></tr>\n"
    if args.hklin_unmerged:
        page += "\t\t<tr><td>Unmerged diffraction data:</td>" \
            "<td>" + args.hklin_unmerged + "</td></tr>\n"
    if args.libin:
        page += "\t\t<tr><td>Crystallographic restraints:</td>" \
            "<td>" + args.libin + "</td></tr>\n"
    if args.tlsin:
        page += "\t\t<tr><td>Input TLS file:</td><td>" + args.tlsin
        if args.tlsin_keep:
            page += "(keep using the same TLS input file " \
                "in all the refinement runs)"
        page += "</td></tr>\n"
    if args.comin:
        page += "\t\t<tr><td>Keywords for REFMAC5:</td>" \
            "<td>" + args.comin + "</td></tr>\n"
    if args.defin:
        page += "\t\t<tr><td>Keywords for phenix.refine:</td>" \
            "<td>" + args.defin + "</td></tr>\n"
    if args.weight:
        page += "\t\t<tr><td>Weight matrix for REFMAC5:</td>" \
            "<td>" + fourdec(args.weight) + "</td></tr>\n"
    if args.tls_ncyc:
        page += "\t\t<tr><td>Number of cycles of TLS refinement:</td>" \
            "<td>" + str(args.tls_ncyc) + "</td></tr>\n"
    if args.ncyc:
        page += "\t\t<tr><td>Number of refinement cycles that " \
            "will be performed in every resolution step:</td>" \
            "<td>" + str(args.ncyc) + "</td></tr>\n"
    if args.prerefinement_ncyc:
        page += "\t\t<tr><td>Number of pre-refinement cycles that will be " \
            "performed before the paired refinement protocol:</td>" \
            "<td>" + str(args.prerefinement_ncyc) + "</td></tr>\n"
    if args.complete_cross_validation:
        page += "\t\t<tr><td><strong>" \
            "Complete cross-validation is performed.</strong></td>" \
            "<td>(across " + str(len(flag_sets)) + " " \
            "free reflection sets)</td></tr>\n"
    else:
        page += "\t\t<tr><td>Data with FreeRflag set " + str(args.flag) + " " \
            "will be excluded during refinement.</td>" \
            "<td></td></tr>\n"

    if (args.complete_cross_validation or args.no_modification or
            args.reset_bfactor or args.add_to_bfactor or args.set_bfactor or
            args.shake_sites):
        page += "\t\t<tr><td>Modification of the input structure model:\n"
        page += "\t\t\t<ul>\n"
    if args.reset_bfactor:
        page += "\t\t\t<li>Reset B-factors to the mean value</li>\n"
    if args.add_to_bfactor:
        page += "\t\t\t<li>Add value to B-factors: " + \
            twodec(args.add_to_bfactor) + "</li>\n"
    if args.set_bfactor:
        page += "\t\t\t<li>Set B-factors to the value: " + \
            twodec(args.set_bfactor) + "</li>\n"
    if args.shake_sites:
        page += "\t\t\t<li>Randomize coordinates with the " + \
            "given mean error value: " + twodec(args.shake_sites) + "</li>\n"
    if args.no_modification:
        page += "\t\t\t<li>No modification</li>\n"
    if (args.complete_cross_validation or args.no_modification or
            args.reset_bfactor or args.add_to_bfactor or args.set_bfactor or
            args.shake_sites):
        page += "\t\t\t</ul>\n"

    if args.constant_grid:
        page += "\t\t<tr><td>The same FFT grid will be kept through the " \
            "whole paired refinement.</td><td></td></tr>\n"
    page += "\t\t</table>\n"
    
    page += "\t\t<h2>Run details and program versions</h2>\n"
    page += "\t\t<table>\n"
    page += "\t\t<tr><td>Working directory:</td>" \
        "<td>" + os.getcwd() + "</td></tr>\n"
    page += "\t\t<tr><td>Run date and time:</td>" \
        "<td>" + date_time + "</td></tr>\n"
    page += "\t\t<tr><td>user@host:</td><td>" \
        "" + getpass.getuser() + "@" + socket.gethostname() + "</td></tr>\n"
    page += "\t\t<tr><td><i>PAIREF</i> version:</td>" \
        "<td>" + versions_dict["pairef_version"] + "</td></tr>\n"
    if args.phenix:
        page += "\t\t<tr><td><i>phenix.refine</i> version:</td>" \
            "<td>" + versions_dict["phenix_version"] + "</td></tr>\n"
    else:
        page += "\t\t<tr><td><i>REFMAC5</i> version:</td>" \
            "<td>" + versions_dict["refmac_version"] + "</td></tr>\n"
    page += "\t\t<tr><td>Python version and executable:</td>" \
        "<td>" + platform.python_version() + " " + \
        sys.executable + "</td></tr>\n"
    page += "\t\t<tr><td>matplotlib version:</td>" \
        "<td>" + matplotlib.__version__ + "</td></tr>\n"
    page += "\t\t</table>\n"

    # Show general warnings if there are some
    warning_keys = ["workdir", "low_res", "refinement_version_mismatch",
                    "refinement_not_before", "no_modification",
                    "mtzdump", "binning", "flags", "amb_labels"]
    page = warning_orangebox(warning_keys, page)
    if ready_shells:
        if len(ready_shells) >= 2:
            page += '\t<h2>Overall values</h2>\n'
            page = warning_orangebox(["overall_R"], page)
            page += '\t\t<div class="wrap">\n'
            page += '\t\t\t<div class="column">\n'
            page += '\t\t\t\t<a href="' + args.project + '_R-values.png">'
            page += '<img src="' + args.project + '_R-values.png' \
                '?shell=' + twodecname(ready_shells[-1]) + '" ' \
                'alt="R-values"></a><br />\n'
            page += '\t\t\t\tRaw data: <a href="' + args.project + '' \
                '_R-values.csv">' + args.project + '_R-values.csv' \
                '</a><br />\n'
            page += '\t\t\t\t<pre>\n'
            with open(args.project + "_R-values.csv", "r") as csvfile:
                page += csvfile.read()
            page += '</pre>'
            page += '\n\t\t\t\t<p class="note">Note: For each incremental ' \
                'step of resolution from X->Y, the <i>R</i>-values were ' \
                'calculated at resolution X.'
            if args.complete_cross_validation:
                page += ' Standard error of mean is shown in orange.'
            page += '</p>\n'
            page += '\t\t\t</div>\n'
            page += '\t\t\t<div class="column">\n'
            page += '\t\t\t\t<a href="' + args.project + '_Rgap.png">' \
                '<img src="' + args.project + '_Rgap.png' \
                '?shell=' + twodecname(ready_shells[-1]) + '" ' \
                'alt="Rgap"></a><br />\n'
            page += '\t\t\t\tRaw data: <a href="' \
                '' + args.project + '_Rgap.csv">' + args.project + '' \
                '_Rgap.csv</a>\n'
            page += '\t\t\t\t<p class="note">Note: The R-values were ' \
                'calculated at ' + twodec(shells[0]) + ' &#8491; resolution.' \
                '</p>\n'
            page += '\t\t\t</div>\n'
            page += '\t\t</div>\n'

            if args.complete_cross_validation:
                for shell in ready_shells[1:]:  # exclude the initial difr. l.
                    page += '\t\t<a href="' + args.project + '_R-values_com' \
                        'plete_' + twodecname(shell) + 'A.png">' \
                        '<img src="' + args.project + '_R-values_complete' \
                        '_' + twodecname(shell) + 'A.png' \
                        '?shell=' + twodecname(shell) + '' \
                        '" alt="R-values (complete cross-validation, ' \
                        '' + twodec(shell) + ' A)"></a>\n'

        # "Rfree", "CCfree", "Rwork", "CCwork", "No_work_free_refl." graphs
        if not args.complete_cross_validation:
            page += '\t<h2>Dependence on resolution</h2>\n'
            # Show relating warnings if there are some
            warning_keys = ["Nfree", "high_R", "low_R", "low_CC"]
            page = warning_orangebox(warning_keys, page)
            # Show "Rfree", "CCfree", "Rwork", "CCwork" graphs
            graphs = ["Rfree", "CCfree", "Rwork", "CCwork"]
            for i, graph in enumerate(graphs):
                pngfilename = args.project + "_" + graph + ".png"
                if os.path.isfile(pngfilename):
                    page += '\t\t<a href="' + args.project + '_' + graph + '' \
                        '.png">' \
                        '<img src="' + args.project + '_' + graph + '.png' \
                        '?shell=' + twodecname(ready_shells[-1]) + '' \
                        '" alt="' + graph + '"></a>\n'
                if i % 2 == 1:  # Display max. 2 graphs in row
                    page += "\t\t<br />\n"
            # Link to raw data
            page += '\t\tRaw data:\n'
            page += '\t\t<ul>\n'
            for shell in ready_shells:
                page += '\t\t\t<li><a href="' + args.project + '_' \
                    'R' + str(args.flag).zfill(2) + '_' \
                    '' + twodecname(shell) + 'A.csv">' + args.project + '_' \
                    'R' + str(args.flag).zfill(2) + '_' \
                    '' + twodecname(shell) + 'A.csv</a></li>\n'
            page += '\t\t</ul>\n'
            # Show No. work free reflections graph
            graph = "No_work_free_reflections"
            pngfilename = args.project + "_" + graph + ".png"
            if os.path.isfile(pngfilename):
                page += '\t\t<a href="' + args.project + '_' + graph + '' \
                    '.png">' \
                    '<img src="' + args.project + '_' + graph + '.png' \
                    '?shell=' + twodecname(ready_shells[-1]) + '' \
                    '" alt="' + graph + '"></a>\n'

    # Optical resolution
    pngfilename = args.project + "_Optical_resolution.png"
    csvfilename = args.project + "_Optical_resolution.csv"
    if os.path.isfile(pngfilename):
        page += "\t<h2>Optical resolution</h2>\n"
        page += '\t\t<a href="' + pngfilename + '">' \
            '<img src="' + pngfilename + \
            '?shell=' + twodecname(ready_shells[-1]) + '" ' \
            'alt="Optical resolution"></a><br />\n'
        page += '\t\tRaw data: <a href="' + csvfilename + '">' + \
            csvfilename + '</a>\n'

    csvfilename = args.project + "_merging_stats.csv"
    if args.hklin_unmerged and ready_merging_statistics \
            and os.path.isfile(csvfilename):
        graphs = ["CC", "Intensities", "Comp_Mult", "Rmerge_Rmeas_Rpim",
                  "No_reflections"]
        page += "\t<h2>Merging statistics</h2>\n"
        warning_keys = ["merging_stats", "CC*"]
        page = warning_orangebox(warning_keys, page)
        for graph in graphs:
            page += '\t\t<a href="' + args.project + '_' + graph + '' \
                '.png">' \
                '<img src="' + args.project + '_' + graph + '.png" alt="' \
                '' + graph + '"></a>\n'
        page += '\t\t<p>Raw data: <a href="' + args.project + '' \
            '_merging_stats.csv">' \
            '' + args.project + '_merging_stats.csv</a></p>\n'
        page += '\t\t<pre>\n'
        with open(args.project + "_merging_stats.csv", "r") as csvfile:
            page += cgi.escape(csvfile.read())
        page += '\t\t</pre>\n'

    # Statistics vs. cycle
    if ((args.complete_cross_validation or args.prerefinement_ncyc) and
            (ready_shells or res_cur)) or \
            (not args.complete_cross_validation and len(ready_shells) >= 2):
        if res_cur:
            ready_shells = ready_shells + [res_cur]
        page += "\t<h2>Statistics vs. cycle</h2>\n"
        for flag in flag_sets:
            page += '\t\t<div class="scroll">\n'
            page += '\t\t<div class="wrap">\n'
            for i, shell in enumerate(ready_shells):
                prefix = args.project + '_R' + str(flag).zfill(2) + \
                    '_' + twodecname(shell) + 'A'
                pngfilename = prefix + '_stats_vs_cycle.png'
                if args.phenix:
                    logfilename = prefix + '_001.log'
                    pdbfilename = prefix + '_001.pdb'
                    ciffilename = prefix + '_001.cif'
                else:  # refmac
                    logfilename = prefix + '.log'
                    pdbfilename = prefix + '.pdb'
                    ciffilename = prefix + '.mmcif'
                if os.path.isfile(pngfilename):
                    page += '\t\t\t<div class="column">\n'
                    page += '\t\t\t\t<a href="' + pngfilename + '">' \
                        '<img src="' + pngfilename + '" alt="' \
                        '' + twodec(shell) + 'A statistics vs. cycle, ' \
                        'flag ' + str(flag).zfill(2) + '"></a><br />\n'
                    page += '\t\t\t\t<a href="' + logfilename + '">' \
                        'Log file</a> from refinement at ' \
                        '' + twodec(shell) + ' &#8491;<br />\n'
                    page += '\t\t\t\t' \
                        'Structure model refined at ' + twodec(shell) + ' ' \
                        '&#8491; <a href="' + pdbfilename + '">PDB</a> ' \
                        '<a href="' + ciffilename + '">mmCIF</a>\n'
                    page += '\t\t\t</div>\n'
            page += '\t\t</div>\n'
            page += '\t\t</div>\n'
    page += """
\t<h2>References</h2>
\t\tPlease reference the used software:
\t\t<ul>
\t\t<li>Paired refinement under the control of <i>PAIREF</i>. M. Maly, K. Diederichs, J. Dohnalek, P. Kolenko (2020) <i>IUCrJ</i> <b>7</b></li>
"""
    if not args.phenix or "ccp4" in sys.executable.lower():
        page += "\t\t<li>Overview of the <i>CCP</i>4 suite and current developments. Collaborative Computational Project, Number 4 (2011) <i>Acta Cryst.</i> D<b>67</b>:235–242</li>\n"
    if args.phenix:
        page += "\t\t<li>Macromolecular structure determination using X-rays, neutrons and electrons: recent developments in <i>Phenix</i>. D. Liebschner, P.V. Afonine, M.L. Baker, G. Bunkóczi, V.B. Chen, T.I. Croll, B. Hintze, L.W. Hung, S. Jain, A.J. McCoy, N.W. Moriarty, R.D. Oeffner, B.K. Poon, M.G. Prisant, R.J. Read, J.S. Richardson, D.C. Richardson, M.D. Sammito, O.V. Sobolev, D.H. Stockwell, T.C. Terwilliger, A.G. Urzhumtsev, L.L. Videau, C.J. Williams, P.D. Adams (2019) <i>Acta Cryst.</i> D<b>75</b>:861-877</li>\n"
        page += "\t\t<li>Towards automated crystallographic structure refinement with <i>phenix.refine</i>. P.V. Afonine, R.W. Grosse-Kunstleve, N. Echols, J.J. Headd, N.W. Moriarty, M. Mustyakimov, T.C. Terwilliger, A. Urzhumtsev, P.H. Zwart, P.D. Adams (2012) <i>Acta Cryst.</i> D<b>68</b>:352-67</li>\n"
    else:
        page += "\t\t<li><i>REFMAC</i>5 for the refinement of macromolecular ""crystal structures. G.N. Murshudov, P. Skubak, A.A. Lebedev, N.S. Pannu, R.A. Steiner, R.A. Nicholls, M.D. Winn, F. Long, A.A. Vagin (2011) <i>Acta Cryst.</i> D<b>67</b>:355–367</li>\n"
    if which("sfcheck"):
        page += "\t\t<li><i>SFCHECK</i>: a unified set of procedures for evaluating the quality of macromolecular structure-factor data and their agreement with the atomic model. A.A. Vaguine, J. Richelle, S.J. Wodak (1999) <i>Acta Cryst.</i> D<b>55</b>:191–205</li>\n"
    page += """
\t\t<li>The <i>Computational Crystallography Toolbox</i>: crystallographic algorithms in a reusable software framework. R.W. Grosse-Kunstleve, N.K. Sauter, N.W. Moriarty, P.D. Adams (2002) <i>J. Appl. Crystallogr.</i> <b>35</b>:126–136</li>
\t\t</ul>
\t\t&nbsp;
\t\t<ul>
\t\t<li>Linking crystallographic model and data quality.  P.A. Karplus & K. Diederichs (2012) <i>Science</i> <b>336</b>:1030–3</li>
\t\t<li>Assessing and maximizing data quality in macromolecular crystallography. P.A. Karplus & K. Diederichs (2015) <i>Cur. Op. in Str. Biology</i> <b>34</b>:60–68</li>
\t\t<li>Better models by discarding data? P.A. Karplus & K. Diederichs (2013) <i>Acta Cryst.</i> D<b>59</b>:1215–1222</li>
\t\t</ul>


</body>
</html>"""

    htmlfilename = "PAIREF_" + args.project + ".html"
    with open(htmlfilename, "w") as htmlfile:
        htmlfile.write(page)

    # Styles
    cssfilepath = str(os.path.dirname(os.path.abspath(__file__))) + \
        "/static/styles.css"
    if os.path.isfile(cssfilepath):
        shutil.copy2(cssfilepath, ".")

    return htmlfilename
