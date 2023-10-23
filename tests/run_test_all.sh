#!/usr/bin/env bash

PARAM=
if [ "$#" == "1" ]
then
  PARAM=$1
fi

if [ -e htmlcov ]
then
  rm -rf htmlcov
fi

PATH_TO_SRC=../fragit
OLD_PYTHONPATH=$PYTHONPATH

export PYTHONPATH=$PYTHONPATH:$PATH_TO_SRC

coverage run --include=../fragit/*.py test_all.py $PARAM

coverage html

export PYTHONPATH=$OLD_PYTHONPATH
rm *.pyc
