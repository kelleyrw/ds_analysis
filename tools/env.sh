#!/bin/bash

# get the project dir
tools_dir=$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -P $tools_dir/.. && pwd)

# add to env 
export DSA=${project_dir}
export PATH=$PATH:$DSA/bin/release
export PYTHONPATH=$PYTHONPATH:$DSA/bin/release
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$DSA/bin/release
