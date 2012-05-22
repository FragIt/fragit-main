#!/usr/bin/env bash

PATH_TO_SRC=../src
OLD_PYTHONPATH=$PYTHONPATH
COVERAGE="python /usr/local/lib/python2.6/dist-packages/coverage"

export PYTHONPATH=$PYTHONPATH:$PATH_TO_SRC
#python test_util.py

$COVERAGE run test_Config.py
$COVERAGE report #--omit="test_Config.py"
$COVERAGE html -d html #--omit="test_Config.py"

export PYTHONPATH=$OLD_PYTHONPATH
