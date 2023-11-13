"""
Copyright (C) 2011-2023 Casper Steinmann
"""

import os

from fragit.util import substitute_file

# filenames for templates
filenames = {"pymol": "pymol.py",
             "jmol": "jmol.spt"}

extension = {"pymol": "py",
             "jmol": "spt"}


class Template(object):
    def __init__(self, directories, infile, outfile):
        self.template_directory = os.path.join(directories["share"], "templates")
        self.infile = infile
        self.outfile = outfile
        self.data = list()  # holds the template
        self.fragments_data = list()
        self.buffer_data = list()
        self.active_data = list()
        self.backbone_data = list()
        self.template_type = None  # must be set later
        self.template_filename = ""  # must be set later
        self.load_structure_string = None
        self.replacements = {}
        self.fragment_charges = list()

    def _set_template_type(self, value):
        if not isinstance(value, str):
            raise ValueError("Template type is a string value.")
        if value not in filenames:
            raise ValueError("Template type '%s' is not valid." % value)
        self.template_type = value
        self.template_filename = os.path.join(self.template_directory, filenames[value])

    def _set_load_structure_string(self, value):
        """ Command used in template to load a structure.

            Arguments:
            ----------
            value -- the command the program must use to load a structure
        """
        if not isinstance(value, str):
            raise ValueError("Load structure string is a string value.")
        self.load_structure_string = value + "\n"

    def set_fragments_data(self, value):
        self.fragments_data = value

    def set_buffer_data(self, value):
        self.buffer_data = value

    def set_active_data(self, value):
        self.active_data = value

    def set_backbone_data(self, value):
        self.backbone_data = value

    def set_pair_data(self, value):
        self.fragpairs = value

    def set_fragment_charges(self, value):
        self.fragment_charges = value

    def format_single_fragment(self, fragment_data):
        raise NotImplementedError

    def format_fragments(self):
        raise NotImplementedError

    def format_backbone(self):
        raise NotImplementedError

    def format_buffer(self):
        raise NotImplementedError

    def format_active(self):
        raise NotImplementedError

    def format_break_points(self):
        raise NotImplementedError

    def format_fragment_charges(self):
        raise NotImplementedError

    def override(self):
        if self.load_structure_string is None:
            raise ValueError("Load structure string not defined")

        self.replacements['FRAGMENTS'] = self.format_fragments()
        self.replacements['BUFFER'] = self.format_buffer()
        self.replacements['ACTIVE'] = self.format_active()
        self.replacements['BACKBONE'] = self.format_backbone()
        self.replacements['BREAKPOINTS'] = self.format_break_points()
        self.replacements['LOADCOMMAND'] = self.load_structure_string % (self.infile)
        self.replacements['FRAGMENTQ'] = self.format_fragment_charges()

    def write(self):
        outfile = "{0:s}.{1:s}".format(self.outfile, extension[self.template_type])
        substitute_file(self.template_filename, outfile, self.replacements)
