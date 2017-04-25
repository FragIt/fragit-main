"""
Copyright (C) 2011-2016 Casper Steinmann
"""

import sys
import string
import os

from .util import substitute_file

# filenames for templates
filenames = {'pymol':'pymol.py',
             'jmol' :'jmol.py'}

extension = {'pymol':'py',
             'jmol' :'jmol'}

class Template(object):
    def __init__(self, directories, infile, outfile):
        self.template_directory = os.path.join(directories['share'], 'templates')
        self.infile = infile
        self.outfile = outfile
        self.data = list() # holds the template
        self.fragments_data = list()
        self.buffer_data = list()
        self.active_data = list()
        self.backbone_data = list()
        self.template_type = None # must be set later
        self.template_filename = "" # must be set later
        self.load_structure_string = None
        self.replacements = {}
        self.fragment_charges = list()

    def _setTemplateType(self,value):
        if not isinstance(value, str): raise ValueError("Template type is a string value.")
        if value not in filenames: raise ValueError("Template type '%s' is not valid." % value)
        self.template_type = value
        self.template_filename = os.path.join(self.template_directory, filenames[value])

    def _setLoadStructureString(self,value):
        """ Command used in template to load a structure.

            Arguments:
            ----------
            value -- the command the program must use to load a structure
        """
        if not isinstance(value, str): raise ValueError("Load structure string is a string value.")
        self.load_structure_string = value + "\n"

    def setFragmentsData(self, value):
        self.fragments_data = value

    def setBufferData(self, value):
        self.buffer_data = value

    def setActiveData(self, value):
        self.active_data = value

    def setBackboneData(self, value):
        self.backbone_data = value

    def setPairData(self,value):
        self.fragpairs = value

    def setFragmentCharges(self, value):
        self.fragment_charges = value

    def formatSingleFragment(self, fragment_data):
        raise NotImplementedError

    def formatFragments(self, fragments):
        raise NotImplementedError

    def formatBackbone(self):
        raise NotImplementedError

    def formatBuffer(self):
        raise NotImplementedError

    def formatActive(self):
        raise NotImplementedError

    def formatBreakPoints(self):
        raise NotImplementedError

    def formatFragmentCharges(self):
        raise NotImplementedError

    def override(self):
        if self.load_structure_string is None: raise ValueError("Load structure string not defined")

        self.replacements['FRAGMENTS'] = self.formatFragments()
        self.replacements['BUFFER'] = self.formatBuffer()
        self.replacements['ACTIVE'] = self.formatActive()
        self.replacements['BACKBONE'] = self.formatBackbone()
        self.replacements['BREAKPOINTS'] = self.formatBreakPoints()
        self.replacements['LOADCOMMAND'] = self.load_structure_string % (self.infile)
        self.replacements['FRAGMENTQ'] = self.formatFragmentCharges()

    def write(self):
        outfile = "{0:s}.{1:s}".format(self.outfile, extension[self.template_type])
        substitute_file(self.template_filename, outfile, self.replacements)
