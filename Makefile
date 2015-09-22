#!/usr/bin/env bash

all:
	python setup.py build
	python setup.py install --prefix=/home/cstein/Programs/fragit-1.4

test:
	cd tests; $(MAKE) $(MFLAGS) test;

test-full:
	cd tests; $(MAKE) $(MFLAGS) test-full;

tag:
	@echo "git tag fragit-1.x.x -m \"FragIt 1.x.x Release\""
