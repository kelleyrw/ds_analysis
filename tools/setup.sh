#!/bin/bash

# get the project dir
tools_dir=$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -P $tools_dir/.. && pwd)

# setup external libs
pushd ${project_dir}/external_libs

echo -e "\nInstalling PyUblas..."
./setup_pyublas.sh

echo -e "\n\nInstalling Boost.Numpy"
./setup_boost_numpy.sh

popd

