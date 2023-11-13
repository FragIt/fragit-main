"""
Copyright (C) 2013-2017 Casper Steinmann
"""
from typing import List

import numpy as np
from openbabel import openbabel

from fragit.writer import Standard
from fragit.util import get_filename_and_extension, Z2LABEL
from fragit.template import Template
from fragit.jmol import JmolTemplate
from fragit.pymol import PymolTemplate


class XYZ(Standard):
    def __init__(self, fragmentation, directories):
        Standard.__init__(self, fragmentation, directories)

    def setup(self):
        self._setup_active_fragments_information()
        if self._do_pymol:
            self._dump_pymol()
        if self._do_jmol:
            self._dump_jmol()

    def _setup_active_fragments_information(self):
        self._active_atoms = []

    def _dump_pymol(self):
        pt = PymolTemplate(self._directories, self._input_filename, self._output_filename)
        self._set_template_data(pt)
        self._write_template_file(pt)

    def _dump_jmol(self):
        pt = JmolTemplate(self._directories, self._input_filename, self._output_filename)
        self._set_template_data(pt)
        self._write_template_file(pt)

    def _set_template_data(self, template: Template):
        template.set_fragments_data(self._fragmentation.get_fragments())
        template.set_buffer_data(self._fragment_layers)
        template.set_active_data(self._active_atoms)
        template.set_backbone_data(self._fragmentation.get_backbone_atoms())
        template.set_pair_data(self._fragmentation.get_explicitly_break_atom_pairs())

    @staticmethod
    def _write_template_file(template):
        template.override()
        template.write()

    def build_single_fragment(self, fragment: List[int]):
        """
            fragment -- atom idx of the current fragment
        """
        output_atoms = [self._fragmentation.get_ob_atom(index) for index in fragment]
        output_types = [atom.GetAtomicNum() for atom in output_atoms]

        return output_atoms, output_types

    @staticmethod
    def write_fragment_xyz(atms: List[openbabel.OBAtom], nuclear_charges: List[int]) -> str:
        """ Generates the xyz file format for a single fragment

            Arguments:
            atms -- list of openbabel atoms
            types -- list of nuclear charges
        """
        xyz_line = "{0:s} {1:20.12f} {2:20.12f} {3:20.12f}\n"
        n = len(atms)
        s = "{0:d}\n{1:s}\n".format(n, "")
        for nucz, _obatom in zip(nuclear_charges, atms):
            (x, y, z) = (_obatom.GetX(), _obatom.GetY(), _obatom.GetZ())
            s += xyz_line.format(Z2LABEL[nucz], x, y, z)

        return s

    def write_file(self, filename):
        """ Dumps all fragments to individual .xyz files

            Arguments:
            filename -- base for the filename with extension
        """
        ff, ext = get_filename_and_extension(filename)
        filename_template = "{0:s}_{1:s}_{2:03d}{3:s}"

        # dump all fragments
        for ifg, fragment in enumerate(self._fragmentation.get_fragments(), start=1):
            (atoms, nuclear_charges) = self.build_single_fragment(fragment)
            ss = self.write_fragment_xyz(atoms, nuclear_charges)
            with open(filename_template.format(ff, "fragment", ifg, ext), "w") as f:
                f.write(ss)
