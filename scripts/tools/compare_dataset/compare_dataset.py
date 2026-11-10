#!/usr/bin/env python3
#
# pylint: disable=line-too-long
#

"""
process automatically 'result_validation' and 'result_reference' directories created by ECOGEN non regression script or CI.
 It compares all pair of files of both directories.
 It outputs the list of test cases for which differences are greater than the wanted tolerance.

Get help:
    ```
     ./scripts/tools/compare_dataset/compare_dataset.py -h
    ```

Usage example:
    ```
     ./scripts/tools/compare_dataset/compare_dataset.py --epsilon 1e-9 -r <results_dir> -v <validation_dir>
    ```
"""

import os
import sys
import argparse
import subprocess
import list_parallel_files as lpf

def usage():
    """Print help for command line options"""
    script_name = os.path.basename(sys.argv[0])
    print(f"""
Usage: {script_name} [-h] [-e,--epsilon epsilon] [-r,--ref reference_dir] [-v,--val validation_dir]

       Process 2 directories containing a list of test cases results (typically the 'results_validation' and 'results_reference'
       directories created by the nonreg script) and compare results data at each time step ('results_CPU*_TIME<n>' files).
       All CPUs data of a given time step are processed at the same time to avoid normalization issues due to the partial vue of the data of each partition.

       Remark: Mesh coordinates are not compared for VTK output files.

        -e, --epsilon <epsilon> maximale tolerance for comparison. Defaults to 1e-12.
        -r, --ref <reference_dir> first directory. Defaults to 'results_reference'.
        -v, --val <validation_dir> second directory. Defaults to 'results_validation'.
""")


def run_diff_command(ref_dir, val_dir):
    """Run diff command and return the output lines"""
    try:
        cmd = [
            'diff', '-rq',
            '-x', '*.xml',
            '-x', 'infoCalcul.out',
            '-x', 'mainOutput.out',
            ref_dir, val_dir
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.stdout.splitlines()
    except subprocess.SubprocessError as e:
        print(f"Error running diff command: {e}", file=sys.stderr)
        return []


def parse_diff_output(diff_lines):
    """Parse diff output to extract time steps having file pairs that differ on at least one partition"""
    timesteps = {}

    for line in diff_lines:
        if line.startswith('Only in'):
            print(f"    Warning: {line}")

        elif line.startswith('Files') and 'differ' in line and 'datasets/result' in line:
            # Extract file paths from diff output
            # Format: "Files ref_file and val_file differ"
            parts = line.split()
            if len(parts) >= 4:
                ref_file = parts[1]
                val_file = parts[3]

                # The key has to use the test case output dir name and the time step
                ref_parent_dir = os.path.dirname(ref_file).split('/')[-2] # name of directory parent to dataset one
                ref_path =  os.path.join(*os.path.dirname(ref_file).split('/')[0:-2]) # path of this directory

                val_parent_dir = os.path.dirname(val_file).split('/')[-2] # name of directory parent to dataset one
                val_path =  os.path.join(*os.path.dirname(val_file).split('/')[0:-2]) # path of this directory
                assert ref_parent_dir == val_parent_dir

                # Extract time step
                time_step_ref = lpf.extract_name_and_timestep(ref_file)[1]
                time_step_val = lpf.extract_name_and_timestep(val_file)[1]
                assert time_step_ref == time_step_val

                # Extract
                if time_step_ref:
                    timesteps[f"'.*TIME{time_step_ref}' files of '{ref_parent_dir}'" ] = (ref_file, val_file,ref_path,val_path)

    return timesteps


def run_comparison(ref_file, val_file, epsilon, exec_path):
    """Based on file extension, run the gnu or vtk comparison script on each time step for which some partitions differ"""
    filename = os.path.basename(ref_file)
    ext = filename.split('.')[-1]

    if ext == 'out':
        compare_script = "compare_gnu.py"
    elif ext in ['vtk', 'vtr', 'vtu', 'vtkhdf']:
        compare_script = "compare_vtk.py"
    else:
        print(f"Unexpected file extension: {ext} !!!")
        return False, f"Unsupported file extension: {ext}"

    try:
        cmd = [
            'python3',
            os.path.join(exec_path, compare_script),
            ref_file, val_file,
            f'--epsilon={epsilon}'
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return result.returncode == 0, result.stdout + result.stderr
    except subprocess.SubprocessError as e:
        return False, f"Error running comparison script: {e}"


def main():
    """process automatically 'result_validation' and 'result_reference' directories created by ECOGEN non regression
    script or CI. It compares all pair of files of both directories.
    It outputs the list of test cases for which differences are greater than the wanted tolerance."""

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Compare result directories', add_help=False)
    parser.add_argument('-h', '--help', action='store_true', help='Show help message')
    parser.add_argument('-e', '--epsilon', type=float, default=1e-12, help='Maximum tolerance for comparison')
    parser.add_argument('-r', '--ref', type=str, default='results_reference', help='Reference directory')
    parser.add_argument('-v', '--val', type=str, default='results_validation', help='Validation directory')

    args = parser.parse_args()

    if args.help:
        usage()
        sys.exit(0)

    epsilon = args.epsilon
    ref_dir = args.ref
    val_dir = args.val

    # Get script directory
    exec_path = os.path.dirname(os.path.abspath(__file__))

    # Check if directories exist
    if not os.path.isdir(ref_dir):
        print(f"Error: Reference directory '{ref_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(val_dir):
        print(f"Error: Validation directory '{val_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Run diff command and parse output
    print("Running diff to find differing files...")
    diff_lines = run_diff_command(ref_dir, val_dir)
    timesteps = parse_diff_output(diff_lines)

    if not timesteps:
        print("No differences found in result files.\n")
        print("\n==================== Final summary =======================")
        print("Results match at wanted tolerance !!")
        print("==========================================================\n")
        return

    print(f"    Common subdirs: Found {len(timesteps)} timestep(s) with differences.\n")

    # Process each timestep
    list_diff = []

    for time_step, (ref_file, val_file, ref_dir, val_dir) in timesteps.items():

        print(f"Comparing: {time_step} directory in '{ref_dir}' and '{val_dir}' directories")

        success, output = run_comparison(ref_file, val_file, epsilon, exec_path)

        if not success:
            print("\nComparison failed for:")
            print(f"  {ref_file}  -- {val_file}")
            print(f"  {output}")
            print("-----------------------\n")

            list_diff.extend([os.path.dirname(ref_file), os.path.dirname(val_file)])

    # Print final summary
    print("\n==================== Final summary =======================")
    if not list_diff:
        print("Results match at wanted tolerance !!")
    else:
        for i in range(0, len(list_diff), 2):
            print(f"Results {list_diff[i]} and {list_diff[i+1]} differ")
    print("==========================================================\n")


if __name__ == "__main__":
    main()
