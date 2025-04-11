"""
Copyright (C) 2011-2025 Casper Steinmann
"""

version = ("1", "9", "2")
version_str = ".".join(version)

doc_str = """FragIt is a tool to fragment large molecules for use in fragment based quantum chemistry methods.
The output of FragIt is typically an input file which can be run in a quantum chemistry program.
The preferable way to use FragIt is to use configuration files (see option -c below) based on
output from the fragit-conf program also included with FragIt.
For further usage information, see the wiki on the github repository.
"""

__doc__ = doc_str
