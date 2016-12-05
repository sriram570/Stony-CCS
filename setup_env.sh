#!/bin/bash
repo_dir=`dirname $0`
repo_path=$(python -c 'import sys,os; \
           print os.path.realpath(os.path.expanduser(sys.argv[1]))' "$repo_dir")
export PATH=$PATH:$repo_path/external/poaV2
unset repo_dir
unset repo_path
