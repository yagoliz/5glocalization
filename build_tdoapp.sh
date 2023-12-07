#!/bin/bash

rm -rf third_party/tdoapp
git clone https://github.com/yagoliz/tdoapp.git third_party/tdoapp

cd third_party/tdoapp
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKS=On -G Ninja
cmake --build . -j 8