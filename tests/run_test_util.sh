#!/usr/bin/env bash

PATH_TO_SRC=../src
OLD_PYTHONPATH=$PYTHONPATH
COVERAGE="python /usr/local/lib/python2.6/dist-packages/coverage"

export PYTHONPATH=$PYTHONPATH:$PATH_TO_SRC
#python test_util.py

$COVERAGE run test_util.py
$COVERAGE report --omit="test_util.py,/usr*"
$COVERAGE html -d html --omit="test_util.py,/usr*"

export PYTHONPATH=$OLD_PYTHONPATH
