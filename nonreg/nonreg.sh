#!/bin/bash
#set -e # Exit with nonzero exit code if anything fails

# -------------------------------------------------------------------------
#
# Configuration for non regression tests
#
# Arguments: see error handling
#
# Assigns global variables:
#   TEST_OUTPUTS - output filename for `run.sh` script (terminal test outputs)
#   TEST_STATUSES - report filename for `run.sh` script (test exit statuses)
#   TMP_TEST_LIST - filename for temporary file that will store the complete test list
#   REPORT_DIR - name of folder for report storage
#   MAIN_WRKTREE_DIR - path of the main worktree
#
# Creates:
#   REPORT_DIR - path of folder to store reports (relative to current ecogen directory)
#
# Moves:
#   results -> results_save folder
#
# Deletes:
#   nonreg/results_reference - reference results from a previous non regression run
#   nonreg/results_validation - validation results from a previous non regression run
#
# Returns:  None
#
# Examples:
#   configure_nonreg libTests/tests.list origin/devel origin/feature/my-new-feature
#
# -------------------------------------------------------------------------
configure_nonreg() {

  local n=2

  if [ $# != $n ]
  then
    echo ""
    echo "configure_nonreg: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "         <reference_branch> - branch to use as reference"
    echo "         <validation_branch> - branch to validate"
    echo ""
    return 1
  fi

  local reference_branch validation_branch
  reference_branch=$1
  validation_branch=$2

  #Settings
  #--------
  TEST_OUTPUTS="mainOutput.out"
  TEST_STATUSES="report.out"
  TMP_TEST_LIST="ECOGEN_nonReg.tmp.list"

  echo "Branch to validate: $validation_branch"
  echo "Validation against branch: $reference_branch"

  local validationBranchSimplified
  validationBranchSimplified=$(echo "$validation_branch" | sed -e "s/\//_/g")

  # Stores the absolute path of the main worktree directory
  MAIN_WRKTREE_DIR=$(pwd)

  #Creating folder for report storage
  #----------------------------------
  # Wrning: relative path for REPORT_DIR: we will have to ensure to save all outputs in the main
  # worktree
  REPORT_DIR="./nonreg/reports/$(date +%Y%m%d_%H%M)_$validationBranchSimplified"
  mkdir -p "$REPORT_DIR"
  #Cleaning results
  #----------------
  run_cmd rm -Rf "$REPORT_DIR"/results_reference "$REPORT_DIR"/results_validation
  run_cmd echo "nonreg/results_reference and nonreg/results_validation deleted"

  # Moves results directory to preserve any pre-existing results (otherwise move
  # the new empty 'results' directory)
  if [ -d results ]
  then
    if ! mv -v results results_save
    then
      return 1
    fi
  fi

  #Check if Gmsh is installed
  #--------------------------

  #On Ubuntu Gmsh redirects to stdout while on MacOS
  #it redirects to stderr (even if success)
  if ! gmsh --version &> /dev/null
  then
    echo "WARNING: Gmsh is not installed or not avalaible in your PATH."
    echo "All unstructured nonreg test cases will not tested."
  fi
}

# -------------------------------------------------------------------------
#
# Concatenate files containing list of tests with possible modification of
# the expected return value for all the tests of a given file.
#
# Arguments: see error handling
#
# Uses global variables:
#   REPORT_DIR - name of directory to save regression results
#
# Returns:  0 for success, 1 for failure
#
# Remark: passing arrays as function args is very sensitive to shell type and version so this
# function may be unstable.
#
# Example:
#   Create a file containing all test cases of 'ECOGEN_nonReg_full.list' without
#   modifications of expected return values and all test cases of  'ECOGEN_nonReg_vtkhdf.list'
#   with 2 as expected return value.
#
#   > files=(./nonreg/ECOGEN_nonReg_full.list ./nonreg/ECOGEN_nonReg_vtkhdf.list)
#   > error_code=('' 2)
#   > concat_test_lists files error_code tmp.list
#
# -------------------------------------------------------------------------
concat_test_lists() {
  local n=3
  if [ $# != $n ]
  then
    echo "concat_test_lists: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <test_lists> - array containing names of files with lists of non regression tests  "
    echo "          <err_codes> - array of integers or empty values: if 'err_codes[i]' is not empty,"\
      "it is added as expected return value to all the test cases listed in file 'test_lists[i]'."
    echo "          <output_file> - name of output file that will contain the complete list of test cases."
    echo ""
    return 1
  fi

  declare input_files # avoid warning on possibly referenced but unset variable
  eval "input_files=(\"\${$1[@]}\")"
  eval "err=(\"\${$2[@]}\")"
  local output_file=$3

  n="${#input_files[@]}" # '${!input_files[@]}' syntax is not compatible with zsh so we don't use it (in ordre to have maximal portability)

  for ((i=0; i<n; i++));do
    # fill 'error' by the wanted error code or by '$3' if we want to keep the third column of the file
    local tmp="${err[*]:$i:1}"
    local error="${tmp:-\$3}"

    # Now remove lines starting by '#' from test list and print either the 3 first columns of non empty lines of the file
    # or the 2 first columns and the wanted error code as third column.
    # Again, syntax '${array[@]:$i:1}' is used to ensure compatibility with
    # zsh whose arrays are 1-index based while bash ones are 0-index based. The syntax means that we look at then entire array ('@'),
    # with a offset of '$i' and we take '1' element
    touch "$output_file" # Avoid errors with `>>` sometime not able to create the file
    grep -v '\#'  "${input_files[@]:$i:1}" 2> /dev/null | awk 'NF {print $1,$2,'"$error"'}' >> "$output_file"
  done
  return 0
}

# -------------------------------------------------------------------------
#
# Build current branch.
#
# Arguments: see error handling
#
# Assigns:
#   error - set error to 1 if build fails
#
# Returns:  0 for success, 1 for failure
#
# Examples:
#   > build_branch test_statuses_file.out my_branch_name
#
# Warning:
#   Can be called from a different worktree than the ecogen repository to test:
#   in this case, depending on the path in the '$output_file' variable, results
#   can be saved in the main worktree (MAIN_WRKTREE_DIR) or in the temporary one
#   (../ecogen_devel)
# -------------------------------------------------------------------------
build_branch() {

  local n=2

  if [ $# != $n ]
  then
    echo ""
    echo "build_branch: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "         <output_file> - name and path of output compilation report file"
    echo "         <branch> - name of current branch (used for output messages only)"
    echo ""
    return 1
  fi

  local output_file branch
  output_file="$1"
  branch="$2"

  echo "Cleaning and compiling using Makefile..."
  make clean > "$output_file"

  if make relwithassert -j 4 >> "$output_file" 2>&1
  then
    echo "Compilation OK"
    return 0
  else
    local sha
    sha="$(git rev-parse --short HEAD)"
    echo "Compilation failed. Running on $branch (commit $sha) abort."
    error=1
    return $error
  fi
}

# -------------------------------------------------------------------------
#
# Test current branch.
#
# Arguments: see error handling
#
# Uses global variables:
#   TEST_OUTPUTS - output filename for `run.sh` script (terminal test outputs)
#   TEST_STATUSES - report filename for `run.sh` script (test exit statuses)
#
# Moves:
#   results -> <REPORT_DIR>/<results_dir>
#
# Assigns:
#   error - set error to 1 if run.sh fails (i.e. contains at least 1 failing test)
#
# Returns:  0 for success, 1 for failure
#
# Examples:
#   test_branch my_test_list.list my_result_dir branch_name
#
# Warning:
#   Can be called from a different worktree than the ecogen repository to test:
#   in this case, depending on the path in the '$results_dir' variable, results
#   can be saved in the main worktree (MAIN_WRKTREE_DIR) or in the temporary one
#   (../ecogen_devel)
# -------------------------------------------------------------------------
test_branch() {

  local n=4
  if [ $# != $n ]
  then
    echo ""
    echo "test_branch: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <test_list> - file with the list of non regression tests"
    echo "          <results_dir> - name and path of results directory"
    echo "          <branch> - name of branch"
    echo "          <ecogen_exec> - path and name of ECOGEN binary"
    echo ""
    return 1
  fi

  local sha
  sha="$(git rev-parse --short HEAD)"

  local test_list results_dir branch

  test_list=$1
  results_dir=$2
  branch=$3
  ecogen=$4

  echo "Running non-regression tests for $branch branch (commit $sha)... More details in main output file ./results/$TEST_OUTPUTS and report file ./results/$TEST_STATUSES"
  chmod u+x ./scripts/run.sh

  # Mark the script as "failing"
  local err=0
  if ! ./scripts/run.sh -l "$test_list" -o ./results/"$TEST_OUTPUTS" -r ./results/"$TEST_STATUSES" -e "$ecogen"
  then
    echo "Error: some of the non regression tests are failing..."
    printf '\n%s\n' "-- Step status: Failed"
    error=1
    err=1
  else
   printf '\n%s\n' "-- Step status: Succeed"
  fi

  mv -v results "$results_dir"

  return "$err"
}

# -------------------------------------------------------------------------
#
# Eval differences in git tree based on git diff arguments stored in
# <diff_args> variable and prune the output to keep only added or
# deleted lines depending on the sign ('+' or '-') passed as argument.
#
# Examples:
#   prune_diff '+'
#   prune_diff '-'
#
# -------------------------------------------------------------------------
prune_diff() {
  local n=1
  if [ $# != $n ]
  then
    echo ""
    echo "prune_diff: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <sign> - '+' or '-' to search for added or removed diffs"
    echo ""
    return 1
  fi

  local sign="$1"

  eval git diff "$diff_args" | grep "^$sign"  | grep -E '[[:alnum:]]' | grep -vE "^\\$sign\\$sign\\$sign"  | grep -vE "^\\$sign\#"
}


# -------------------------------------------------------------------------
#
# Count diff lines based on git diff arguments stored in <diff_args>
# variable and the given sign ('+' or '-') passed as argument.
#
# Examples:
#   count_line_diff '+'
#   count_line_diff '-'
#
# -------------------------------------------------------------------------
count_line_diffs() {
  local n=1
  if [ $# != $n ]
  then
    echo ""
    echo "count_line_diff: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <sign> - '+' or '-' to search for added or removed diffs"
    echo ""
    return 1
  fi

  prune_diff "$1" | wc -l
}

# -------------------------------------------------------------------------
#
# Store diff lines based on git diff arguments stored in <diff_args>
# variable and the given sign ('+' or '-') passed as argument.
#
# Examples:
#   store_diff_list '+'
#   store_diff_list '-'
#
# -------------------------------------------------------------------------
store_diff_list() {
  local n=1
  if [ $# != $n ]
  then
    echo ""
    echo "store_diff_list: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <sign> - '+' or '-' to search for added or removed diffs"
    echo ""
    return 1
  fi

  local sign="$1"
  prune_diff "$sign" | sed s/'^'"$sign"// | awk NF=1 | sort | uniq
}

# -------------------------------------------------------------------------
#
# Write in DIFF_FILTER and FIND_FILTER the command line arguments that have
# to be passed to, respectively, the diff and find commands, in order to ignore:
#   - the output directories of:
#      - tests whose inputs have been modified
#      - tests added in the non regression test list
#      - tests deleted from the non regression test list
#   - the lines containing tests paths added/deleted from the non-reg test
#     list (only if the second argument is 1)
#
# The list of test cases paths has to be passed as first argument. If the
# second argument is 0, ignore only the ouptut test directories, if the
# second argument is 1, ignore both the output test directories and lines
# containing test paths.
#
# Modifies global variable:
#   DIFF_FILTER: concatenation of the diff exclusion commands
#   FIND_FILTER: concatenation of the find exclusion commands
#
#
# Example:
#   exclude_from_diff "./dir1/dir2/ dir3/dir4" 0
#   exclude_from_diff "./dir5" 1
#
# -------------------------------------------------------------------------
exclude_from_diff() {
  local n=2
  if [ $# != $n ]
  then
    echo ""
    echo "exclude_from_diff: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <test_cases_paths> - list of test cases paths"
    echo "          <ignore_test_paths> - 1 to ignore lines containing"
    echo "                                the test paths"
    echo ""
    return 1
  fi

  while IFS= read -r filepath
  do
    # Find the main.xml file in the given path and extract the run names
    # Remark: In case of removed test cases (test cases existing in reference branch but not validation one),
    # the filepath doesn't exist so test case is not added to list of ignored tests
    if [ -d "$filepath" ]
    then
      while IFS= read -r run_name
      do
        DIFF_FILTER+=" -x $run_name" # Output directories to exclude from final diff call
        FIND_FILTER+=" -path */$run_name/* -prune -o" # Same output directories to exclude from 'find' command in final comparison of vtkhdf files.
                                                     # Remark: path should not be quoted for the command to work
        if [ "$2" -eq "1" ]
        then
          # Test cases that will be ignored when comparing final 'report.out' files (lines 'Test case: '$filepath' ...'  will be ignored)
          DIFF_FILTER+=" --ignore-matching-lines=$filepath"
        fi
      done < <(find "$filepath" -type f -name "main.xml" -exec sed -n 's:.*<run>\(.*\)</run>.*:\1:p' {} +)
    fi
  done <<< "$1"
}

# -------------------------------------------------------------------------
#
# Check non regression test cases modifications.
# Ignore whitespace changes and blank lines.
#
# Arguments: see error handling
#
# Assigns global variables:
#   INPUT_DIFFS: list of non regression test cases with modified inputs
#   ADDED_LIST: list of non regression tests added in <test_list>
#   REMOVED_LIST: list of non regression tests deleted from <test_list>
#   DIFF_FILTER: computed diff exclusion arguments
#
# Returns:
#    0 for success (no diff, news tests added to the test list or modification
#      of the inputs of existing tests),
#    1 if tests have been deleted from the test list
#    2 for unexpected failure
#
# Examples:
#   eval_diff my_test_list.list origin/devel ""
#   eval_diff my_test_list.list origin/devel my_branch
#
#
# -------------------------------------------------------------------------
eval_diff() {

  local n=3
  if [ $# != $n ]
  then
    echo ""
    echo "count_diff_lines: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <test_list> - array of files with list of non-regression tests"
    echo "          <reference_branch> - name of reference branch"
    echo "          <validation_branch> - name of validation branch"
    echo "                                An empty string can be passed to check uncommited changes in validation branch."
    echo ""
    return 2
  fi

  local test_list reference_branch validation_branch

  if declare -p "$1" >/dev/null 2>&1; then
    # argument passed by variable name
    eval "test_list=(\${$1[@]})"
  else
    # argument is a literal value
    test_list="$1"
  fi

  reference_branch="$2"
  validation_branch="$3"

  # First: check for modifications in existing test cases inputs (directory add or removal are ignored)
  local err=0

  local files
  files=$(git diff --ignore-blank-lines --ignore-space-change --diff-filter=M --stat "$reference_branch"..${validation_branch:+"$validation_branch"} -- nonreg/nonregTests | sed '$d' | awk '{print $1}' )

  if [ -n "$files" ]
  then
    echo "Warning: Modified test cases ↵"

    INPUT_DIFFS=$(echo "$files" | xargs -n 1 dirname | sort | uniq)

    while IFS= read -r line; do
      echo "         $line"
    done <<< "$INPUT_DIFFS"
  fi

  # For now we still want to check the diffs in the output results of
  # modified tests. Uncomment the next line to change this behaviour.
  #exclude_from_diff "$INPUT_DIFFS" 0

  # Second: check for added/removed test cases in the non regression test list
  local diff_args="--ignore-blank-lines --ignore-space-change $reference_branch..$validation_branch -- ${test_list[*]}"

  # The 'git diff --numstat' command allow to find easily stat infos
  # on number of lines changed between commits but it counts blank
  # lines added just before a new line (it is view as a unique
  # modification).  Thus we cannot directly use its output to give
  # informations on added and deleted lines.
  #
  # Filtering the git diff output works as soon as 'git diff' is not
  # empty.  In consequence, we check for empty diff with 'git diff
  # --numstat', then we treat non empty diffs with combinations of
  # diff and grep commands.
  local has_diff
  has_diff=$(eval git diff --numstat "$diff_args" | wc -l)

  if [[ "$has_diff" -eq "0" ]]
  then
    return "$err"

  else
    echo "Warning: modifications in one or several test lists: ${test_list[*]}"

    diff_args="-U0 $diff_args"

    local added
    added="$(count_line_diffs '+')"
    if [[ "$added" -gt "0"  ]]
    then
      echo "         Added   $added line(s) of test case(s)."
      ADDED_LIST="$(store_diff_list '+')"
      exclude_from_diff "$ADDED_LIST" 1

    fi

    local removed
    removed="$(count_line_diffs '-')"
    if [[ "$removed" -gt "0"  ]]
    then
      echo "         Removed $removed line(s) of test case(s)."
      echo "         => pipeline marked as failed!"
      # remark: this command is useless in fact: as directory doesn't exist, exclude_from_diff command will not be able to
      # extract output directory name so it cannot be added to the list of exclude directories.
      # A visual check by user is needed.
      REMOVED_LIST="$(store_diff_list '-')"
      exclude_from_diff "$REMOVED_LIST" 1
      err=1
    fi

    return "$err"
  fi

}

# -------------------------------------------------------------------------
#
# Find vtkhdf files in directory while ignoring directories stored in
# 'FIND_FILTER'
#
# Arguments: see error handling
#
# Returns:
#    0 if command is run with suitable number of args
#    1 if wrong number of args
#
# Examples:
# FIND_FILTER="-path dir1 -prune -o -path dir2 -prune -o" find_vtkhdf_with_filter result_validation
#
# -------------------------------------------------------------------------

find_vtkhdf_with_filter() {

  local n=1
  if [ $# != $n ]
  then
    echo ""
    echo "find_vtkhdf_with_filter: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <dir> - directory in which listing files"
    echo ""
    return 1
  fi

  # echo command is needed on zsh to avoid error ' unknown primary or operator' when FIND_FILTER is filled and
  # 'No such file or directory' when empty. $(printf ...) cannot be quoted without errors.
  find "$1" $(printf '%s ' "$FIND_FILTER") -type f -name '*.vtkhdf'  -print | awk -F'/' '{ print $(NF-2)"/"$(NF-1)"/"$NF }' | sort

}

# -------------------------------------------------------------------------
#
# Check if all vtkhdf files present in 'dir1' exists in 'dir2'
#
# Arguments: see error handling
#
# Returns:
#    0 if command is run with suitable number of args
#    1 if wrong number of args
#
# Examples:
#   compare_vtkhdf_file_list_with_ref results_reference results_validation
#
# -------------------------------------------------------------------------
compare_vtkhdf_file_list_with_ref() {

  local n=2
  if [ $# != $n ]
  then
    echo ""
    echo "compare_vtkhdf_file_list_with_ref: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <dir1> - reference directory"
    echo "          <dir2> - directory in which we check presence of vtkhdf files of 'dir1'"
    echo ""
    return 1
  fi

  cmd="comm -23 <(find_vtkhdf_with_filter $1) <(find_vtkhdf_with_filter $2) | sed -En \"s|(.*)|File '\1' only in '$1' |p\""

  if [ "$DRY_RUN" = "true" ]; then
    echo "dry-run: $cmd"
  else
    eval "$cmd"
  fi

  return 0
}

# -------------------------------------------------------------------------
#
# Check if 2 directories contains the same vtkhdf file list
#
# Returns:
#    0 if number of args is ok
#    1 if wrong number of args
#
# Examples:
#   compare_vtkhdf_file_list results_reference results_validation
#
# -------------------------------------------------------------------------
compare_vtkhdf_file_list() {

  local n=2
  if [ $# != $n ]
  then
    echo ""
    echo "compare_vtkhdf_file_list: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <dir1> - first directory"
    echo "          <dir2> - second directory"
    echo ""
    return 1
  fi

  compare_vtkhdf_file_list_with_ref "$1" "$2" && compare_vtkhdf_file_list_with_ref "$2" "$1"

  return 0
}

# -------------------------------------------------------------------------
#
# Check if 2 directories contains the same vtkhdf files
#
# Uses global variables:
#   FIND_FILTER - find exclusion arguments
#
# Returns:
#    0 if directories have no differences
#    1 if wrong number of args
#    2 if directories contain different list of vtkhdf files
#
# Examples:
#   compare_vtkhdf results_reference results_validation
#
# -------------------------------------------------------------------------
compare_vtkhdf() {

  local n=2
  if [ $# != $n ]
  then
    echo ""
    echo "compare_vtkhdf: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <dir1> - first directory"
    echo "          <dir2> - second directory"
    echo ""
    return 1
  fi

  error=0

  # Set directories from arguments
  dir1="$1"
  dir2="$2"

  # 'FIND_FILTER' is used here
  out=$(compare_vtkhdf_file_list "$dir1" "$dir2")

  if [ -n "$out" ]; then
    echo "$out"
    if [ "$DRY_RUN" != "true" ]; then
      error=2
    fi
  fi

  # Remark: next command is asymmetrical: if a file in 'dir1' has no counterpart in 'dir2'
  # a message is printed while it is not true for files in 'dir2' but not in 'dir1'.
  # It doesn't matter as file lists have previously been compared.
  if [ "$DRY_RUN" != "true" ]; then # find will fail in dry_run mode
    # Loop through all .vtkhdf files in the first directory recursively
    while IFS= read -r file; do
      # Get the relative path of the file
      rel_path="${file#"$dir1"/}"

      # Construct the corresponding path in the second directory
      counterpart="$dir2/$rel_path"

      # Check if the counterpart file exists
      if [ -f "$counterpart" ]; then
        # Compare the files using h5diff
        if ! run_cmd h5diff -q --enable-error-stack "$file" "$counterpart"; then
          echo "Differences found in: '$rel_path'"
          error=2
        fi
      else
        echo "Counterpart not found for: '$rel_path'"
        error=2
      fi
      # In next 'find' call echo command is needed on zsh to avoid error ' unknown primary or operator'
      # when FIND_FILTER is filled and 'No such file or directory' when empty. $(printf...) command cannot
      # be quoted without error
    done < <(find "$dir1" $(printf '%s ' "$FIND_FILTER") -type f -name '*.vtkhdf' -print)
  fi

  return "$error"
}

# -------------------------------------------------------------------------
#
# Check if 2 directories contains the same vtkhdf files
#
# Uses global variables:
#   TEST_OUTPUTS - output filename for `run.sh` script (terminal test outputs)
#   DIFF_FILTER - diff exclusion arguments
#   FIND_FILTER - find exclusion arguments
#
#
# Returns:
#    0 if called with suitable number of args
#    1 if wrong number of args
#
# Examples:
#   compare_results results_reference results_validation
#
# -------------------------------------------------------------------------
compare_results() {

  local n=3
  if [ $# != $n ]
  then
    echo ""
    echo "compare_results: Wrong number of arguments."
    echo "    Expects $n args:"
    echo "          <dir1> - first directory"
    echo "          <dir2> - second directory"
    echo "          <output_file> - path and name of file in which saving the output"
    echo ""
    return 1
  fi

  # Set directories from arguments
  local dir1="$1"
  local dir2="$2"
  local output_file="$3"

  echo "-----------------------------------------"
  echo "Computing differences between reference and validation branches..."

  echo "
  Exclude list for diff: $DIFF_FILTER
  "
  echo "
  Exclude list for find: $FIND_FILTER
  "

  local cmd1
  cmd1="diff --ignore-blank-lines --ignore-space-change $DIFF_FILTER -x '*.vtkhdf' -x $TEST_OUTPUTS -qr -x 'infoCalcul.out' -x '.DS_Store' $dir1 $dir2 | tee -a $output_file 2>&1 || true"
  local cmd2
  # 'FIND_FILTER' global variable is used inside 'compare_vtkhdf'
  cmd2="compare_vtkhdf $dir1 $dir2 | tee -a $output_file 2>&1 || true"

  # Here the run_cmd script fails for normal mode due to the combination of output redirection and
  # protected arguments passed to 'eval' ('\$*')
  if [ "$DRY_RUN" = "true" ]; then
    echo "dry-run: $cmd1"
    echo "dry-run: $cmd2"
  else
    eval "$cmd1"
    eval "$cmd2"
  fi

}

# -------------------------------------------------------------------------
#
# Formatted print of the branch name passed as argument
#
# -------------------------------------------------------------------------
print_message() {
  mess="$1"

  echo ""
  printf '%*s\n' ${#mess} | tr ' ' '*'
  echo "$mess"
  printf '%*s\n' ${#mess} | tr ' ' '*'
}

# -------------------------------------------------------------------------
#
# Formatted print of the branch name passed as argument of the non
# regression message
#
# -------------------------------------------------------------------------
print_nonreg_message() {

  print_message "------- RUNNING $1 BRANCH NONREG TESTS -------"

}


# -------------------------------------------------------------------------
#
# Clean worktree at exit
#
# -------------------------------------------------------------------------
function cleanup {
  echo ""
  echo "Exit. Clean temporary worktree."
  git worktree remove -f ecogen_reference 2>/dev/null
  git branch -D reference 2>/dev/null
}

# -------------------------------------------------------------------------
#
# Check differences between reference and validation branch
#
# Moves:
#   results_save -> results folder (if exist)
#
# Uses global variables:
#   TEST_OUTPUTS - output filename for `run.sh` script (terminal test outputs)
#   REPORT_DIR - name of folder for report storage
#   DIFF_FILTER - diff exclusion arguments
#
# Deletes:
#   nonreg/results_reference - reference results from a previous non regression run
#   nonreg/results_validation - validation results from a previous non regression run
#
# Returns: 0 for success, 1 for failure
#
# Example:
#   check_and_clean
#
# -------------------------------------------------------------------------
check_and_clean() {

  echo ""
  echo "************************************************************************"
  echo "--------   FINISHING NONREG TESTS (CLEANING/COMPUTING DIFF)  -----------"
  echo "************************************************************************"
  echo "Cleaning..."

  cd "$MAIN_WRKTREE_DIR" ||
    (echo "Unable to restore results directory and compute diffs" && return 1 )

  if [ -d results_save ]
  then
    mv -v results_save results
  fi

  out=$(compare_results "$REPORT_DIR/results_reference" "$REPORT_DIR/results_validation" "$REPORT_DIR/diff_nonreg.out")

  if [ -s "$REPORT_DIR"/diff_nonreg.out ] # we can use a relative path for REPORT_DIR because we are in MAIN_WRKTREE_DIR
  then
    echo "Differences exist: details available in $REPORT_DIR/diff_nonreg.out file."
    printf '\n%s\n' "-- Step status: Failed"
    error=1
  else
    echo "Non-regression tests complete: tests common to both branches present no difference!"
    printf '\n%s\n' "-- Step status: Succeed"
  fi
}

# -------------------------------------------------------------------------
#
# Dry-run mode: echo commands in dry-run mode,
# run them otherwise
#
# -------------------------------------------------------------------------
run_cmd() {
  if [ "$DRY_RUN" = "true" ]; then
    echo "dry-run: $*"
  else
    # $* is escaped to avoid special char evaluation before command execution which may lead to missing argument
    # when the name of evaluation branch starts by # for example (# is evaluated as a comment). See PR #119.
    eval "\$*"
  fi
}

# -------------------------------------------------------------------------
#
# Print help for main program.
#
# -------------------------------------------------------------------------
usage() {
  echo ""
  echo "Usage: $(basename "$0") [-d] [-r reference_branch] [-l test_list [err_code]] [-l test_list2 [err_code2]...] "
  echo ""
  echo "    Build ecogen using the repository 'Makefile', run non-regression tests and validate results"\
       " against a reference branch."
  echo ""
  echo "        -d dry-run mode."
  echo "        -r <reference_branch> branch to use as reference branch."\
       " Defaults to 'origin/devel'."
  echo "        -l <test_lists> [n] path and name of file with test list "\
       "(can be used multiple times to concatenate tests)."\
       "If an integer 'err_code' is provided, 'err_code' is added as expected error value for all the tests of the file."\
       " Defaults to './nonreg/ECOGEN_nonReg_full.list'."
  echo ""
}

# -------------------------------------------------------------------------
#
# Main non regression script when called from command line.
#
# Arguments:
#   [-r <reference_branch>] - name of branch used as reference
#   [-l <test_list>] - file containing the list of tests to run
#
# -------------------------------------------------------------------------
main() {
  trap cleanup EXIT # Execute cleanup function at exit command

  error=0
  DRY_RUN=false

  test_list_files=()
  file_error_codes=()

  # Parse command line options
  while getopts ':hdr:l:' OPTION; do
    case "$OPTION" in
      l)
        test_list_files+=("$OPTARG")
        # Check if the next argument exists and is a positive integer: OPTIND is the index of next arguments and
        # ${!OPTIND} allow to get the value of the argument at index OPTIND
        local next_arg
        next_arg="${!OPTIND}"
        if [[ -n "$next_arg" && "$next_arg" =~ ^[1-9][0-9]*$ ]]
        then
          file_error_codes+=("$next_arg")
          # Advance OPTIND
          ((OPTIND++))
        else
          file_error_codes+=("")
        fi

        ;;
      r)
        reference_branch_name=$OPTARG
        ;;
      d)
        DRY_RUN=true
        ;;
      h)
        usage
        exit 0
        ;;
      \?)
        echo ""
        echo "Error: Illegal option -$OPTARG"
        usage
        exit 1
        ;;
      :)
        echo ""
        echo "Error: Option -$OPTARG requires an argument"
        usage
        ;;
    esac
  done

  # Assign default values if parameters are unset
  test_list_files=("${test_list_files[@]:-./nonreg/ECOGEN_nonReg_full.list}")
  reference_branch_name=${reference_branch_name:=origin/devel}

  validation_branch_name="$(git rev-parse --abbrev-ref HEAD)"

  # Parse command line options and configure non regression tests.
  # Defines global variables for output dirs and files and creates
  # report directories
  if ! configure_nonreg "$reference_branch_name" "$validation_branch_name"
  then
    # Either the function is called with wrong arguments, or old
    # 'results' directory cannot be moved in 'results_save'. Quit the
    # process to avoid loting results
    echo "Error when configuring non-regression checks."
    exit 1
  fi

  print_message "------- EVALUATING TEST CASES DIFFERENCES -------"
  # Check if non regression tests has been modified. Here
  # we don't want to pass the name of the validation branch because it
  # will ignore possible uncommited modifications.
  if ! eval_diff test_list_files "$reference_branch_name" ""
  then
    if [ "$?" -eq "2" ]
    then
      echo "Error when checking modification in non regression test list."
      exit 1
    else
      # $? == 1 => Error: test deletion in the test list!
      error=1
      printf '\n%s\n' "-- Step status: Failed"
    fi
  else
    printf '\n%s\n' "-- Step status: Succeed"
  fi

  # Concatenate the input list files with wanted error codes
  concat_test_lists test_list_files file_error_codes "$REPORT_DIR/$TMP_TEST_LIST"

  # Build and test the branch to validate
  print_nonreg_message "$validation_branch_name"
  if run_cmd build_branch "$REPORT_DIR/output_compile.out" "$validation_branch_name"
  then
    run_cmd test_branch "$REPORT_DIR/$TMP_TEST_LIST" "$REPORT_DIR/results_validation" "$validation_branch_name" "./ECOGEN"
  fi

  # Clone, build and test the reference branch
  print_nonreg_message "$reference_branch_name"

  echo "Fetching and checking out $reference_branch_name branch..."

  if ! git fetch > "$REPORT_DIR"/output_fetchAndCheckout.out
  then
    echo "Warning: unable to fetch the remote repository."
    echo "Attempt to run using last commit locally available on $reference_branch_name."
  fi

  count=$(git worktree list | \grep  -c ecogen_reference)
  if [ "$count" ]
  then
    # unclean worktree: delete the old tree
    git worktree remove -f ecogen_reference >> "$REPORT_DIR"/output_fetchAndCheckout.out 2>&1
    git branch -D reference >> "$REPORT_DIR"/output_fetchAndCheckout.out 2>&1
  fi

  # Checkout the reference branch in a new tree
  if ! git worktree add ../ecogen_reference -B reference "$reference_branch_name" >> "$REPORT_DIR"/output_fetchAndCheckout.out
  then
    echo "Error during checkout on $reference_branch_name. Running on $reference_branch_name abort."
    error=1
  else
    cd ../ecogen_reference ||
      (echo "Unable to move in 'ecogen_reference' directory. Exit." && exit 1)

    # Concatenate the input list files with wanted error codes (file will be moved in final result dir later)
    concat_test_lists test_list_files file_error_codes "$TMP_TEST_LIST"

    # Here, we are not anymore in the ecogen repo we want to test (MAIN_WRKTREE_DIR): provide
    # suitable path to the output files and directories (allows to directly save 'reference' results in MAIN_WRKTREE_DIR
    # without having to move them)
    if run_cmd build_branch "$MAIN_WRKTREE_DIR/$REPORT_DIR/output_compile_ref.out" "$reference_branch_name"
    then
      run_cmd test_branch "$TMP_TEST_LIST" "$MAIN_WRKTREE_DIR/$REPORT_DIR/results_reference" "$reference_branch_name"  "./ECOGEN"
    fi
  fi

  # keep track of temporary test list of ref branch
  mv "$TMP_TEST_LIST" "$MAIN_WRKTREE_DIR/$REPORT_DIR/$TMP_TEST_LIST-ref"

  cd "$MAIN_WRKTREE_DIR" ||
    (echo "Unable to move in '$MAIN_WRKTREE_DIR' directory. Next step may fail.")

  # Check test differences and final cleaning
  if ! check_and_clean
  then
    error=1
  fi

  # Sending error code if something went wrong
  exit "$error"
}

# -------------------------------------------------------------------------
#
# Definition of two usages for current script file:
#
#     - the file can be sourced to load the defined functions and allow
#       individual function call
#
#     - the file can be executed as a script to run the non regression
#       pipeline
#
# Examples:
#   source "./nonreg/nonreg.sh"
#   source "./nonreg/nonreg.sh" && usage
#   ./nonreg/nonreg.sh
#   ./nonreg/nonreg.sh -h
#   ./nonreg/nonreg.sh -r branch_toto -l ./libTests/tests.list
#
# -------------------------------------------------------------------------

# Do not call the main function when we source the .sh file
if [ "$0" = "${BASH_SOURCE[0]}" ]
then
  # Script call from command line: execute the main func
  main "$@"
fi
