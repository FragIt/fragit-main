"""
Copyright (C) 2013-2017 Casper Steinmann
"""
import numpy

from .writer import Standard
from .util import getFilenameAndExtension, Z2LABEL

class XYZ(Standard):
    def __init__(self, fragmentation, directories):
        Standard.__init__(self,fragmentation, directories)

    def setup(self):
        self._setupLayeredInformation()
        self._setupActiveFragmentsInformation()
        if self._do_pymol: self._dump_pymol()
        if self._do_jmol: self._dump_jmol()

    def _setupLayeredInformation(self):
        self._fragment_layers = self._getFragmentLayersFromFragment()

    def _getFragmentLayersFromFragment(self):
        fragments = self._fragmentation.getFragments()
        return numpy.array([1 for i in fragments])

    def _setupActiveFragmentsInformation(self):
        self._active_atoms = []

    def _dump_pymol(self):
        from pymol import PymolTemplate
        pt = PymolTemplate(self._directories, self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _dump_jmol(self):
        from jmol import JmolTemplate
        pt = JmolTemplate(self._directories, self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _setTemplateData(self, template):
        template.setFragmentsData(self._fragmentation.getFragments())
        template.setBufferData(self._fragment_layers)
        template.setActiveData(self._active_atoms)
        template.setBackboneData(self._fragmentation.getBackboneAtoms())
        template.setPairData(self._fragmentation.getExplicitlyBreakAtomPairs())

    def _writeTemplateFile(self, template):
        template.override()
        template.write()

    def _build_single_fragment(self, fragment):
        """
            fragment -- atom idx of the current fragment
        """
        output_atoms = [self._fragmentation.mol.GetAtom(id) for id in fragment]
        output_types = [atom.GetAtomicNum() for atom in output_atoms]

        return output_atoms, output_types

    def fragment_xyz(self, atms, types):
        """ Generates the xyz file format for a single fragment

            Arguments:
            atms -- list of openbabel atoms
            types -- list of nuclear charges
        """
        xyz_line = "{0:s} {1:20.12f} {2:20.12f} {3:20.12f}\n"
        n = len(atms)
        s = "{0:d}\n{1:s}\n".format(n,"")
        for nucz, _obatom in zip(types, atms):
            (x,y,z) = (_obatom.GetX(), _obatom.GetY(), _obatom.GetZ())
            s += xyz_line.format(Z2LABEL[nucz], x, y, z)

        return s

    def writeFile(self, filename):
        """ Dumps all fragments to individual .xyz files

            Arguments:
            filename -- base for the filename with extension
        """
        ff, ext = getFilenameAndExtension(filename)
        filename_template = "{0:s}_{1:s}_{2:03d}{3:s}"

        # dump all fragments
        for ifg, fragment in enumerate(self._fragmentation.getFragments(), start=1):
            (atms, types) = self._build_single_fragment(fragment)
            ss = self.fragment_xyz(atms, types)
            with open(filename_template.format(ff, "fragment", ifg, ext), "w") as f:
                f.write(ss)
