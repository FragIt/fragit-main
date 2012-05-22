#!/usr/bin/env bash

PATH_TO_SRC=../src
OLD_PYTHONPATH=$PYTHONPATH

export PYTHONPATH=$PYTHONPATH:$PATH_TO_SRC

python test_all.py

export PYTHONPATH=$OLD_PYTHONPATH
rm *.pyc
