#!/bin/bash

version=v2013.1
git clone http://git.tiker.net/trees/pyublas.git pyublas/${version} 
pushd pyublas/${version}
git checkout ${version}
popd
