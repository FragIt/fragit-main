#!/usr/bin/env bash

all:
	python setup.py build
	python setup.py install --prefix=/home/cstein/Programs

test:
	cd tests; $(MAKE) $(MFLAGS) test;

test-full:
	cd tests; $(MAKE) $(MFLAGS) test-full;
