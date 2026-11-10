#!/bin/sh

# Get the current branch name
branch_name=$(git rev-parse --abbrev-ref HEAD)

# Get the commit message file path
commit_msg_file=$1
# Get the commit source (see https://git-scm.com/docs/githooks#_hooks)
commit_source=$2

# Read the current commit message
commit_msg=$(cat "$commit_msg_file")

# Help detect merge commits (sometimes, commit_source is empty, for
# example when conflicts have been solved)
if git rev-parse -q --verify MERGE_HEAD >/dev/null 2>&1
then
  commit_source="merge"
fi

# Only run for certain commit types (e.g., not during rebase with 'reword' or 'edit')
if [ "$commit_source" != "merge" ] && [ "$commit_source" != "squash" ] && [ "$commit_source" != "commit" ]; then

  # Commit source is empty when -m and -t options are not provided or in rebase cases.
  # Discriminate rebase cases by checking that the branch is not HEAD.
  if [ "$branch_name" != "HEAD" ]; then
    # Check if the commit message starts with the branch name followed by a colon
    if ! echo "$commit_msg" | grep -q "^$branch_name:"; then
      # If not, prepend the branch name followed by a colon to the commit message
      echo "$branch_name: $commit_msg" > "$commit_msg_file"

      # Print a message to notify the user that the commit title was modified
      echo "commit-msg: Commit title was automatically updated to prefix it with the branch name."
    fi
  fi
fi
