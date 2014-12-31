#!/bin/bash

# setup external libs
pushd external_libs

echo -e "\nInstalling PyUblas..."
./setup_pyublas.sh

echo -e "\n\nInstalling Boost.Numpy"
./setup_boost_numpy.sh

popd

