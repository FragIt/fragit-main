#!/usr/bin/env bash

all:
	python setup.py build
	python setup.py install --prefix=/home/cstein/Programs

test:
	cd tests; $(MAKE) $(MFLAGS) test;

test-full:
	cd tests; $(MAKE) $(MFLAGS) test-full;

tag:
	@echo "git tag fragit-1.3.6 -m \"FragIt 1.3.6 Release\""
