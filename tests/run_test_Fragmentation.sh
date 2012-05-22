#!/usr/bin/env bash

PATH_TO_SRC=../src
OLD_PYTHONPATH=$PYTHONPATH
COVERAGE="python /usr/local/lib/python2.6/dist-packages/coverage"

export PYTHONPATH=$PYTHONPATH:$PATH_TO_SRC

$COVERAGE run test_Fragmentation.py
$COVERAGE report --omit="test_Fragmentation.py,/usr*"
$COVERAGE html -d html --omit="test_Fragmentation.py,/usr*"

export PYTHONPATH=$OLD_PYTHONPATH
