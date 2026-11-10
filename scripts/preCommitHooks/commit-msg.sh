#!/bin/bash
#
# A hook script to check that the commit title starts by the branch name.
#
# Called by "git commit" with one argument, the name of the file
# that has the commit message.  The hook should exit with non-zero
# status after issuing an appropriate message if it wants to stop the
# commit.  The hook is allowed to edit the commit message file.
#

# Get the first line (title) of the commit message
commit_title=$(head -n 1 "$1")

# Get the current branch name
if [[ $CI_COMMIT_REF_NAME ]]; then
  # Inside CI, branch is always checked out as detached HEAD but in
  # this case we want to check the commit messages.
  branch_name=$CI_COMMIT_REF_NAME
  commit_sha=$CI_COMMIT_SHORT_SHA
else
  # Avoid issues during interactive rebase or merge conflicts
  # resolution. Do not check commit messages in this case.
  branch_name=$(git rev-parse --abbrev-ref HEAD)
  commit_sha=$(git rev-parse --short=6 --ref HEAD)
fi

# Check if the branch name is not HEAD (which happens during rebases)
# and don't check the commit title in that case
if [ "$branch_name" != "HEAD" ]; then
  # Check if the commit title starts with the branch name followed by a colon
  if [[ ! $commit_title =~ ^$branch_name: && ! "$commit_title" =~ ^"Merge" ]]; then
    echo "Error: Commit '$commit_sha: $commit_title' message title must start with the branch name '$branch_name:'"
    exit 1
  else
    echo "Info: Commit '$commit_sha: $commit_title' has a valid title format."
  fi
else
  echo "Warning: Commit '$commit_sha: $commit_title' skip check of the commit title format."
fi
