#!/usr/bin/env python3
#
# pylint: disable=line-too-long

"""
This module provides helpers to extract the time step value and base name of a ECOGEN solution file and
to list all files with same basename related to the detected time step.

For example, if 'resultB64_CPU2_TIME10.vtu' file is given as input and test case has been partitionned on
4 partitions, files 'resultB64_CPU0_TIME10.vtu', 'resultB64_CPU1_TIME10.vtu', 'resultB64_CPU2_TIME10.vtu'
and 'resultB64_CPU3_TIME10.vtu' are listed (with directory path).
"""

import os
import re
from pathlib import Path

def extract_name_and_timestep(filename):
    """
    Extract base name and time step from filename.

    Returns:
        tuple: (basename, timestep)
               basename is 'result_CPU', 'resultB64_CPU', or 'result'
               time_step is an integer
    """
    # Detect time in files named either "result_CPU<n>_TIME<m>.<ext>", resultB64_CPU<n>_TIME<m>.<ext> or "result_TIME<m>.<ext>"
    match = re.search(r"(result(?:B64)?_CPU|result)(?:(?:\d|\*)+)?_TIME(\d+)\.(.*)", filename)
    if match:
        return match.group(1),int(match.group(2)),match.group(3)

    raise ValueError(f"Filename {filename} does not match expected pattern.")


def get_files_for_time(directory, basename, timestep, extension):
    """
    Find all files matching a given base name and time step.

    Args:
        directory: Path to search in
        basename: 'result_CPU', 'resultB64_CPU', or 'result'
        timestep: Time step value

    Returns:
        list: Sorted list of matching filenames
    """
    files = []

    if not os.path.isdir(directory) :
        raise ValueError(f"Directory {directory} does not exist.")

    for file in Path(directory).glob("*"):
        if not file.is_file():
            continue

        base, time, ext = extract_name_and_timestep(file.name)

        if base == basename and time == timestep and ext == extension:
            files.append(directory+"/"+file.name)

    if not files:
        raise ValueError(f"No file with name '{basename}.*_TIME{timestep}.*' in directory '{directory}'.")

    return sorted(files)  # Sort for consistent order
