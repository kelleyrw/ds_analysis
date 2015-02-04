#!/bin/bash

# get the project dir
tools_dir=$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -P $tools_dir/.. && pwd)

# add to env 
export DSA=${project_dir}
export PATH=$PATH:$DSA/bin/release
export PYTHONPATH=/usr/local/lib/python2.7/site-packages/gurobipy:$MR/gurobi_study/python:$DSA/bin/release:$RS/utilities:/Users/rkelley/Development/rovi/Research/projects/python:$DSA/bin/temp/external_libs/boost_numpy/darwin-4.2.1/release/address-model-64/architecture-ia64/:$DSA/bin/temp/programs/reach_study/darwin-4.2.1/release/address-model-64/architecture-ia64/
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$DSA/bin/release
