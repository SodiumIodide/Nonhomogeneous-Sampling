#!/usr/bin/env bash

rm -rf ./build/*
cd ./build/
cmake ../src
make
./NonHomog
