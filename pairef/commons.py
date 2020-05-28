# coding: utf-8
import sys
from .settings import warning_dict


def twodec(var):
    """Returns number with 2 decimals as a string.
    If a float is not given, it returns a string.

    Args:
        var (float)

    Returns:
        str
    """
    try:
        i = ('{:.{prec}f}'.format(var, prec=2))
    except ValueError:
        i = str(var)
    return(i)


def twodecname(var):
    """Returns number with 2 decimals but intead of decimal point is used "-"

    Args:
        var (float)

    Returns:
        str
    """
    res = twodec(var)
    out = res.replace(".", "-")
    return(out)


def fourdec(var):
    """Returns number with 4 decimals as a string.
    If a float is not given, it returns a string.

    Args:
        var (float)

    Returns:
        str
    """
    try:
        i = ('{:.{prec}f}'.format(var, prec=4))
    except ValueError:
        i = str(var)
    return(i)


def extract_from_file(filename, searched, skip_lines, n_lines,
                      nth_word=False, not_found="stop", get_first=False):
    """Returns line(s) or word relating to the search based on
    `searched` string in the file `filename`.

    If the `filename` is not found, abort (always).
    Default behavior: The last case of `searched` string match is used
    and if the `searched` string is not found in the file `filename`,
    the program is stopped with an error message.

    Args:
        filename (str): Name of the file
        searched (str): Searched string
        skip_lines (int): Number of lines to be skipped
        n_lines (int): Number of lines that should be returned
        nth_word (bool or int): `False` if lines should be returned
            or a order in a line of the word that should be picked
        not_found (str): If the `searched` string is not found, write an error
            message and exit (if `not_found="stop"`, default option)
            or (if `not_found="N/A"`) return `"N/A"` or `["N/A"]`.
        get_first (bool): Use the first case of `searched` string match.

    Returns:
        list or str:
            List containing number of lines (controlled by `n_lines`)
            with the `skip_line` offset (if `nth_word=False`) or
            picked word (`nth_word`-th word) in the `skip_lines`-th
            following line (if `nth_word=True`)
    """

    # Python 2 and 3 compatibility
    try:
        FileNotFoundError
    except NameError:
        FileNotFoundError = IOError

    try:
        # str() is needed as "TypeError: coercing to Unicode:
        #                     need string or buffer, PosixPath found"
        with open(str(filename), "r") as f:
            file_lines = f.readlines()
    except FileNotFoundError:
        sys.stderr.write("ERROR: File " + str(filename) + " was not found.\n"
                         "Aborting.\n")
        sys.exit(1)
    for i in range(len(file_lines)):
        if searched in file_lines[i]:
            j = i
            if get_first:
                break
    if "j" not in locals():
        if not_found == "stop":
            sys.stderr.write("ERROR: File " + str(filename) + " is not in a "
                             "proper format. Statistics could not be found.\n"
                             "Aborting.\n")
            sys.exit(1)
        elif not_found == "N/A":
            if not nth_word:
                lines_array = ["N/A"]
                return lines_array
            else:
                word = "N/A"
                return word
    if not nth_word:
        lines_array = file_lines[j + skip_lines:j + skip_lines + n_lines]
        return lines_array
    else:
        line = file_lines[j + skip_lines]
        try:
            word = line.split()[nth_word]
        except IndexError:
            if not_found == "stop":
                sys.stderr.write("ERROR: File " + str(filename) + " is not in "
                                 "a proper format. "
                                 "Statistics could not be found.\n"
                                 "Aborting.\n")
                sys.exit(1)
            elif not_found == "N/A":
                word = "N/A"
        return word


def warning_my(key, message):
    message = "WARNING: " + message
    if key in warning_dict:
        if not message in warning_dict[key]:
            sys.stderr.write(message + "\n")
            warning_dict[key] = warning_dict[key] + "<br />\n" + message
        # else (if message in warning_dict[key]) do nothing
    else:
        sys.stderr.write(message + "\n")
        warning_dict[key] = message
    return True


def try_symlink(src, dst):
    """Make new symlink to `src` if the `dst` file does not exist yet. If it is
    not possible to make symlinks (difficulties on Windows), just make a copy
    of the file even if exists already.

    Args:
        src (str): File name of the source
        dst (str): File name of the destination

    Returns:
        bool: True"""

    import os
    if hasattr(os, "symlink"):
        # New symlink only if it has not been made previously
        if not os.path.isfile(dst):
            os.symlink(src, dst)
    else:  # Windows
        import shutil
        shutil.copy2(src, dst)
    return True
