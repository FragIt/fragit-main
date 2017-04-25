#!/usr/bin/env python
# This file is part of FragIt (http://github.com/FragIt/fragit-main),
# a library to fragment molecules for use in fragment based methods
# in quantum chemistry.
#
# Copyright (C) 2012-2016, Casper Steinmann
#
#
import sys
from distutils.core import setup

# use the source code to get version information
from src.strings import version_str

__doc__="""FragIt: a tool to fragment molecules for fragment based methods.

FragIt is a python based tool that allows you to quickly fragment
a molecule and use the output as an input file in quantum chemistry
programs that supports such fragment based methods.

Currently, only the Fragment Molecular Orbital (FMO) method in GAMESS
is supported, but FragIt has been developed to easily allow for other
output writers to be added quickly."""


# Chosen from http://www.python.org/pypi?:action=list_classifiers
classifiers = """\
Development Status :: 5 - Stable
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GPL2 or later
Natural Language :: English
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering :: Chemistry
Topic :: Software Development :: Libraries :: Python Modules
"""

def setup_fragit():
    doclines = __doc__.split("\n")

    setup(name="fragit",
        version=version_str,
        url = "https://github.com/FragIt/fragit-main",
        author = "Casper Steinmann",
        author_email = "casper.steinmann@gmail.com",
        maintainer = "Casper Steinmann",
        maintainer_email = "casper.steinmann@gmail.com",
        license = "GPL2 or later",
        description = doclines[0],
        long_description = "\n".join(doclines[2:]),      
        classifiers = filter(None, classifiers.split("\n")),
        platforms = ["Any."],
        package_dir={'fragit': 'src'},
        packages=['fragit'],
        scripts=['scripts/fragit', 'scripts/fragit-conf'],
        package_data = {'': ['pymol_template','jmol_template']}, # relative to 'packages' specified above
        data_files=[
            ('', ['README.md','LICENSE', 'CHANGES.md']),
            ('share', ['share/README.md']),
            ('share/templates', [
                'share/templates/pymol.py', 'share/templates/jmol.py'
            ]),
            ('share/hmo', [
                'share/hmo/STO-3G', 'share/hmo/3-21G', 'share/hmo/6-31G',
                'share/hmo/6-31G(d)', 'share/hmo/6-31G*', 'share/hmo/cc-pVDZ',
                'share/hmo/pcseg-0', 'share/hmo/pcseg-1'
            ])
        ]
  )

if __name__ == '__main__':
  setup_fragit()
