#!/usr/bin/env python3
#
# pylint: disable=line-too-long
#
"""
This script allows to compare solutions fields of two ECOGEN output files at GNU format (.out extension) and
with same partitioning:
  - file path, base names, time steps and extensions are extracted from filenames provided as input
  - if solution is saved under a distributed form (file per process output), all partitions are recovered and compared (even if a partition index is explicitely provided).

Usage and help can be obtained with the `-h` command line option:
    `./scripts/tools/compare_dataset/compare_gnu.py -h`

To compare results of directories 'dir1' and 'dir2' at timestep 3:
    `./scripts/tools/compare_dataset/compare_gnu.py --epsilon 1e-9 <dir1>/datasets/result_CPU*_TIME3.out <dir2>/datasets/result_CPU*_TIME3.out`

or, with same result:
    `./scripts/tools/compare_dataset/compare_gnu.py --epsilon 1e-9 <dir1>/datasets/result_CPU1_TIME3.out <dir2>/datasets/result_CPU2_TIME3.out`
"""

import argparse
import sys
import os
import numpy as np
import list_parallel_files as lpf


#########################################################################################################
# tabs are lists of arrays of 2 dimensions: We have one 2D array for each input file associated to a CPU
# and each array stores the solution fields in columns
# (so array[i,j] is the value of the jth solution field at the ith cell)
#########################################################################################################
def compare_dataarrays(tab1, tab2, epsilon):
    """Compare solution fields partition by partition at a tolerance normalized by the field magnitude computed over the entire domain"""
    all_equal = True

    if len(tab1) != len(tab2):
        return False

    # compute magnitude of each solution field
    ## Max magnitude of each solution field
    all_rows1 = np.vstack(tab1)
    all_rows2 = np.vstack(tab2)
    if all_rows1.shape != all_rows2.shape:
        return False

    mag1 = np.max(np.fabs(all_rows1),axis=0)
    mag2 = np.max(np.fabs(all_rows2),axis=0)

    ## Use the max magnitude find using the 2 output directories
    mag = np.fmax(mag1,mag2)

    for f,ar1 in enumerate(tab1):
        ar2 = tab2[f]
        for i in range(ar1.shape[1]):
            allclose = np.allclose(ar1[:,i], ar2[:,i], atol=epsilon*mag[i])
            if not allclose:
                maxerr = np.max(np.abs(ar1[:,i]-ar2[:,i]))
                print(f"Files CPU'{f}' - Columns '{i+1}' differ (max magnitude: '{mag[i]}', max diff: '{maxerr}').\n")
                all_equal = False

    return all_equal

##########################################################################################################
# Remark: not tested on large cases and not optimized at all (may be very slow and memory expensive)
#########################################################################################################
def main():
    """Compare coordinates and solution values between two sets of solutions files at GNU output file format
    and for a given time step. Files have to use the same partitionning and all partitions are compared."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Compare coordinates and solution values between two GNU files.")

    parser.add_argument("file1", help="First GNU file")
    parser.add_argument("file2", help="Second GNU file")
    parser.add_argument('-e',"--epsilon", type=float, default=1e-12, help="Maximal tolerance for comparison (default: 1e-12)")
    args = parser.parse_args()

    # Extract directory and filename
    dir1, fname1 = os.path.split(args.file1)
    dir2, fname2 = os.path.split(args.file2)

    # Extract basename and timestep value from the filename
    basename1,time_value1,ext1 = lpf.extract_name_and_timestep(fname1)
    basename2,time_value2,ext2 = lpf.extract_name_and_timestep(fname2)

    # Get all files with this TIME value in each directory
    files_in_dir1 = lpf.get_files_for_time(dir1,basename1, time_value1, ext1)
    files_in_dir2 = lpf.get_files_for_time(dir2,basename2, time_value2, ext2)

    # Load all parallel files of a given time step and add outputs to a list
    tab1 = []
    tab2 = []

    for f1 in files_in_dir1:
        tab1.append(np.loadtxt(f1))

    for f2 in files_in_dir2:
        tab2.append(np.loadtxt(f2))

    print("\nFile parsing OK\n")

    # Compare solution fields partition by partition but normalize the tolerance by the maximal
    # magnitude of the field over the whole domain
    cell_equal = compare_dataarrays(tab1, tab2, args.epsilon)

    if cell_equal:
        print("\nAll data are equal within epsilon.\n")
        sys.exit(0)
    else:
        print("\nSome data differ.\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
