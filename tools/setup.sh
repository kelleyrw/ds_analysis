#!/bin/bash

# setup external libs
pushd external_libs
./setup_pyublas.sh
./setup_boost_numpy.sh
popd

