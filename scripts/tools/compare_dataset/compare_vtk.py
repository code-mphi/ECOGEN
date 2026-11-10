#!/usr/bin/env python3
#
# pylint: disable=line-too-long,invalid-name
#

"""
This script allows to compare solutions fields of two ECOGEN output files at VTK or VTKHDF format (.vtu, .vtr or .vtkhdf extension) and
with same partitioning:
  - file path, base names, time steps and extensions are extracted from filenames provided as input
  - if solution is saved under a distributed form (file per process output), all partitions are recovered and compared (even if a partition index is explicitely provided).

Comparison is possible between:
  - unstructured files with same partitioning (.vtu and .vtu file, unstructured VTKHDF file and .vtu file, unstructured VTKHDF file and unstructured VTKHDF file).
  - Two .vtr files containing Cartesian grids with same partitioning.
  - Cartesian VTKHDF files and/or serial .vtr files containing a Cartesian grid (because Cartesian VTKHDF file is saved under its serial form).

The `CPUrank` field of VTKHDF file is ignored when the file is compared to VTK file format and for Cartesian grids.

Usage and help can be obtained with the `-h` command line option:
    `./scripts/tools/compare_dataset/compare_vtk.py -h`

To compare results of directories 'dir1' and 'dir2' at timestep 3:
    `./scripts/tools/compare_dataset/compare_vtk.py --epsilon 1e-9 <dir1>/datasets/result_TIME3.vtkhdf <dir2>/datasets/result_CPU*_TIME3.vtu`

or, with same result:
    `./scripts/tools/compare_dataset/compare_vtk.py --epsilon 1e-9 <dir1>/datasets/result_TIME3.vtkhdf <dir2>/datasets/result_CPU2_TIME3.vtu`
"""

import argparse
import sys
import os
import vtk
import h5py
import numpy as np
from colorama import Fore
import list_parallel_files as lpf

from vtk.util import numpy_support

####################################################################################################
# Store Cell Data objects of one VTK partition under the for of a dictionnary whose key is the
# solution field name
####################################################################################################
def get_partitiondata(data):
    """ Store cell data objects of a VTK partition as a dictionnary"""
    array_dict = {} # empty dictionnary to store arrays of current partition

    for i in range(data.GetNumberOfArrays()):
        name = data.GetArrayName(i) # get field name
        vtk_array = data.GetArray(i) # get VTK array (type VTKarray)
        np_array = vtk_to_numpy(vtk_array)  # Convert VTK array to NumPy array

        # Keep track of partition to which each array belong in order to print which CPU file has too large errors
        array_dict[name] = np_array

    # Remark: dictionnaries are passed by references so a copy of dictionnary linked to the partition
    # is needed. Here it is automatically done because the object return by the function is copied.
    return array_dict


####################################################################################################
# travel a list of vtk objects and append data fields in a list of dictionnaries using the data name
####################################################################################################
# vtk_objects contains the object readed in the file (so it is either a list/array of unstructured
# meshes of the CPU files or 1 unique partitionned mesh composed of multiple partitions)
def get_celldataarrays(vtk_objects):
    """Travel a list of vtk objects and append data fields in a list of dictionnaries using the data name"""
    arrays = [] # list of dictionnaries.
    # Each list item will contain a dictionnary where
    # the name of the solution field is associated to the array of solution values
    for data_object in vtk_objects:

        class_name = data_object.GetClassName()
        if class_name == "vtkPartitionedDataSet" : # Unstructured grid saved at VTKHDF file format
            for i in range(data_object.GetNumberOfPartitions()):
                data = data_object.GetPartition(i).GetCellData()
                arrays.append(get_partitiondata(data))

        elif class_name == "vtkUnstructuredGrid" : # Unstructured grid saved at VTU file format.
            data = data_object.GetCellData()
            arrays.append(get_partitiondata(data))

        elif class_name == "vtkRectilinearGrid" : # Rectilinear grid saved at custom VTKHDF or at VTR file format
            data = data_object.GetCellData()
            arrays.append(get_partitiondata(data))

        else :
            print( Fore.RED + "\n  -- Error: unexpected grid type: " + class_name + "\n" + Fore.RESET )
            sys.exit(1)

    # Remark: a new item (each item being a dictionary) is added to the list for each partition
    # (so either for each CPU file or for each "partition" of vtkhdf file)
    return arrays

########################################################################
# Convert vtkArray into a numpy one: each vtk array is an array of tuple
########################################################################
def vtk_to_numpy(vtk_array):
    """Convert a vtkArray in a numpy array"""

    return np.array([vtk_array.GetTuple(i) for i in range(vtk_array.GetNumberOfTuples())])

##########################################
# Compute magnitude of each solution field
##########################################
def magnitude(arrays1,arrays2):
    """Compute maximal magnitude of each solution field over the whole domain"""

    mag = {} # dictionnary to store the global maximum magnitude of each field (computed using data of all parititons)
    for f,ar1 in enumerate(arrays1):
        ar2 = arrays2[f]

        common = sorted(set(ar1.keys()) & set(ar2.keys())) # Sort set for predictable results

        for name in common:
            mag1 = np.max(np.fabs(ar1[name]))
            mag2 = np.max(np.fabs(ar2[name]))
            if name in mag:
                mag[name] = max(mag[name],mag1,mag2) # maximum magnitude over the 2 compared files for and all partitions already seen
            else:
                mag[name] = max(mag1,mag2) # First partition seen: Initialization with maximum magnitude over the 2 compared file

    return mag

######################
# Compare array values
######################
def compare_arrays(arrays1, arrays2, epsilon, mag, ignore_CPUrank):
    """Compare solution fields partition by partition at a tolerance normalized by the field magnitude computed over the entire domain"""
    all_equal = True

    for f,ar1 in enumerate(arrays1):
        ar2 = arrays2[f]

        common = sorted(set(ar1.keys()) & set(ar2.keys())) # Sort set for predictable results

        for name in common:

            # In some cases we want to ignore the CPUrank solution field
            if name == "CPUrank" and ignore_CPUrank:
                continue

            if ar1[name].shape != ar2[name].shape:

                print(Fore.RED + f"\n Error: Partition {f} - Different shapes of arrays for array '{name}' of the 2 compared"
                      f" files ({ar1[name].shape} and {ar2[name].shape})\n" + Fore.RESET)
                print(Fore.RED + "        Non-comparable arrays.\n" + Fore.RESET)
                return False

            allclose = np.allclose(ar1[name], ar2[name], atol=epsilon*mag[name])

            if not allclose:
                maxerr = np.max(np.abs(ar1[name]-ar2[name]))
                print(Fore.RED + f"\nPartition {f} - CellDataArray '{name}' differ (max magnitude: '{mag[name]}', max diff: '{maxerr}').\n" + Fore.RESET)
                all_equal = False

        for name in set(ar1.keys()) - set(ar2.keys()):
            # In some cases we want to ignore the CPUrank solution field
            if name == "CPUrank" and ignore_CPUrank:
                continue

            print(Fore.RED + f"\nPartition {f} - CellDataArray '{name}' only in first file.\n" + Fore.RESET)
            all_equal = False

        for name in set(ar2.keys()) - set(ar1.keys()):
            # In some cases we want to ignore the CPUrank solution field
            if name == "CPUrank" and ignore_CPUrank:
                continue

            print(Fore.RED + f"\nPartition {f} - CellDataArray '{name}' only in second file.\n" + Fore.RESET)
            all_equal = False

    return all_equal

####################################################################################################
# Check if data fields are the same at a given tolerance normalized by the maximal field magnitude.
# tab1 and 2 are lists of VTK objects (each list item stores the solution of one partition)
####################################################################################################
def compare_celldataarrays(tab1, tab2, epsilon, ext1, ext2):
    """Compute maximal magnitude of each solution field over the entire domain and compare fields partition by partition"""

    arrays1 = get_celldataarrays(tab1)
    arrays2 = get_celldataarrays(tab2)

    # Detect VTKHDF file
    file1_isVTKHDF = ext1=="vtkhdf"
    file2_isVTKHDF = ext2=="vtkhdf"

    # Detect cartesian VTKHDF file
    file1_isCartesianVTKHDF = (file1_isVTKHDF and tab1[0].GetClassName() == "vtkRectilinearGrid")
    file2_isCartesianVTKHDF = (file2_isVTKHDF and tab2[0].GetClassName() == "vtkRectilinearGrid")

    if len(arrays1) != len(arrays2): # Files dosen't have the same number of partitions
        print(Fore.RED + "\n  -- Error: differences for the number of partitions in each directory.")

        if file1_isCartesianVTKHDF or file2_isCartesianVTKHDF :
            print( "            Comparison of partitionned data to cartesian VTKHDF output files is not yet implemented.")
            print( "            Please, compare cartesian VTKHDF file to another cartesian VTKHDF file or to a serial file.\n" + Fore.RESET)
        else :
            print( "            Only files with identical partitionning can be compared.\n"  + Fore.RESET )

        sys.exit(1)

    # Compute solution magnitude for each field
    mag = magnitude(arrays1,arrays2)

    # In some cases, we want to ignore the 'CPUrank' field and to compare only other fields:
    #   - in case of the comparison of 2 VTKHDF cartesian outputs: file outputs have to be the same whatever the number of MPI process used, at the exception of the CPUrank array
    #   - when comparing a VTKHDF output (unstructured grid or rectilinear grid) and a non-HDF VTK output (.vtu or .vtr extension): the CPUrank array is not saved in .vtu and .vtr files
    #   (because the information is retrieved automatically from the file partitioning). Remark, when parsing unstructured VTKHDF files, we get an object of type 'vtkPartitionedDataset'
    ignore_CPUrank = False
    if (file1_isCartesianVTKHDF and file2_isCartesianVTKHDF) or (file1_isVTKHDF and not file2_isVTKHDF) or ( file2_isVTKHDF and not file1_isVTKHDF):
        ignore_CPUrank = True


    # Array comparison
    all_equal = compare_arrays(arrays1,arrays2,epsilon,mag, ignore_CPUrank)
    return all_equal

def flatten_first_three_dims(arr):
    """Flatten first 3 dimensions, preserve 4th if it exists."""
    if arr.ndim == 3:
        return arr.reshape(-1)

    if arr.ndim == 4:
        return arr.reshape(-1, arr.shape[3])

    raise ValueError("Array must be 3D or 4D")

def read_vtk_file(filename):
    """Read a VTK file using suitable reader based on file extension."""
    # Try different VTK readers based on file extension
    if filename.endswith(".vtr"):
        reader = vtk.vtkXMLRectilinearGridReader()
    elif filename.endswith(".vtu"):
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif filename.endswith(".vtk"):
        reader = vtk.vtkDataSetReader()
    elif filename.endswith(".vtkhdf"):
        # For now, Rectilinear cartesian grid cannot be readed by VTK module:
        # check mesh type and use VTK reader if possible, custom reader otherwise
        f = h5py.File(filename,'r')
        if f['VTKHDF'].attrs['Type'].decode() == "RectilinearGrid" :
            # Parse Cell Data field and store it under the form of a VTK object
            grid = vtk.vtkRectilinearGrid()
            if 'CellData' in f['VTKHDF']:
                for array_name in f['VTKHDF']['CellData'].keys():
                    array_data = f['VTKHDF']['CellData'][array_name][:]
                    # Flattening array is required to pass it to VTK (because we have too many dimensions)
                    # We want to flatten the geometric direction: for scalar data, we obtain an array of range 1 and dim ncells.
                    # For vectorial data, we obtain an array of range 2 and dims (ncells,3)
                    flat_array = flatten_first_three_dims(array_data)
                    vtk_array = numpy_support.numpy_to_vtk(flat_array)
                    vtk_array.SetName(array_name)
                    grid.GetCellData().AddArray(vtk_array)
            return grid

        # 'RectilinearGrid' keyword is not found: Standard VTKHDF file
        reader = vtk.vtkHDFReader()
    else:
        print(Fore.RED + f"\n  -- Error: Unexpected file extension for file {filename}.\n" + Fore.RESET)
        sys.exit(1)

    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

##########################################################################################################
# Remark: not tested on large cases and not optimized at all (may be very slow and memory expensive)
#########################################################################################################
def main():
    """Compare solution values between two sets of solutions files at VTK or VTKHDF file format
    (.vtk, .vtr or .vtkhdf extension) and for a given time step.
    Files have to use the same partitionning and all partitions are compared. Mesh is ignored"""

    parser = argparse.ArgumentParser(description="Compare all Solution values between two VTK or VTKHDF files. Mesh coordinates are not compared.")

    parser.add_argument("file1", help="First VTK or VTKHDF file.")
    parser.add_argument("file2", help="Second VTK or VTKHDF file.")
    parser.add_argument('-e', '--epsilon', type=float, default=1e-12, help="Maximale tolerance for comparison (default: 1e-12)")
    args = parser.parse_args()

    # Extract directory and filename
    dir1, fname1 = os.path.split(args.file1)
    dir2, fname2 = os.path.split(args.file2)

    # Extract basename and timestep value from the filename
    basename1,time_value1, ext1 = lpf.extract_name_and_timestep(fname1)
    basename2,time_value2, ext2 = lpf.extract_name_and_timestep(fname2)

    # Get all files with this TIME value in each directory
    files_in_dir1 = lpf.get_files_for_time(dir1,basename1, time_value1, ext1)
    files_in_dir2 = lpf.get_files_for_time(dir2,basename2, time_value2, ext2)

    # Load all parallel files of a given time step
    tab1 = []
    tab2 = []

    for f1 in files_in_dir1:
        tab1.append(read_vtk_file(f1))

    for f2 in files_in_dir2:
        tab2.append(read_vtk_file(f2))

    # Extension of reference file of each directory will be used in 'compare_celldataarrays' to detect VTKHDF files
    cell_equal = compare_celldataarrays(tab1, tab2, args.epsilon, ext1, ext2)

    if cell_equal:
        print("\nAll CellDataArrays are equal within epsilon.\n")
        sys.exit(0)
    else:
        print("\nSome CellDataArrays differ.\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
