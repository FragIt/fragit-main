#!/usr/bin/env python
import sys
from distutils.core import setup

# use the source code to get version information
from src.strings import version_str
from src.strings import doc_str

__author__ = "Casper Steinmann"
__copyright__ = "Copyright 2017"
__license__ = 'GPL2 or later'
__version__ = version_str
__email__ = "casper.steinmann@gmail.com"
__url__ = "https://github.com/FragIt/fragit-main"
__doc__= doc_str


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
        version=__version__,
        url = __url__,
        author = __author__,
        author_email = __email__,
        maintainer = __author__,
        maintainer_email = __email__,
        license = __license__,
        description = doclines[0],
        long_description = "\n".join(doclines[2:]),      
        classifiers = filter(None, classifiers.split("\n")),
        platforms = "Any",
        package_dir={'fragit': 'src'},
        packages=['fragit'],
        scripts=['scripts/fragit', 'scripts/fragit-conf'],
        data_files=[
            ('', ['README.md', 'LICENSE', 'CHANGES.md']),
            ('share', ['share/README.md']),
            ('share/templates',
              [
                'share/templates/pymol.py', 'share/templates/jmol.py'
              ]
            ),
            ('share/hmo',
              [
                'share/hmo/STO-3G', 'share/hmo/3-21G', 'share/hmo/6-31G',
                'share/hmo/6-31G(d)', 'share/hmo/6-31G*', 'share/hmo/cc-pVDZ',
                'share/hmo/pcseg-0', 'share/hmo/pcseg-1'
              ]
            )
        ]
  )

if __name__ == '__main__':
    setup_fragit()
