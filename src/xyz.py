"""
Copyright (C) 2013-2016 Casper Steinmann
"""
from numpy import sqrt, dot, where, array
from writer import Standard
from util import WriteStringToFile

from util import file_extension,is_list,listTo2D,join2D,is_int
from util import listToRanges,listOfRangesToString,Uniqify,ravel2D
from util import deepLength,listDiff,intlistToString
from util import getFilenameAndExtension

class XYZ(Standard):
    def __init__(self, fragmentation):
        Standard.__init__(self,fragmentation)

    def setup(self):
        self._setupLayeredInformation()
        self._setupActiveFragmentsInformation()
        #self._validateMultiLayerInformation()
        if self._do_pymol: self._dump_pymol()
        if self._do_jmol: self._dump_jmol()

    def _setupLayeredInformation(self):
        self._fragment_layers = self._getFragmentLayersFromFragment()

    def _getFragmentLayersFromFragment(self):
        fragments = self._fragmentation.getFragments()
        return array([1 for i in fragments])

    def _setupActiveFragmentsInformation(self):
        self._active_atoms = []

    def _dump_pymol(self):
        from pymol import PymolTemplate
        pt = PymolTemplate(self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _dump_jmol(self):
        from jmol import JmolTemplate
        pt = JmolTemplate(self._input_filename, self._output_filename)
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

    def _build_single_fragment(self,fragment):
        """
                fragment: atom idx of the current fragment
                pairs     : pairs of breaking points
        """
        output_atoms = [self._fragmentation.mol.GetAtom(id) for id in fragment]
        output_types = [atom.GetAtomicNum() for atom in output_atoms]

        return output_atoms, output_types

    def fragment_xyz(self, atms, types):
        """Generates the xyz file format based on the atoms, types,
           ids and neighbours of each fragment
        """
        n = len(atms)
        s = "%i\n%s\n" % (n,"")
        for id, (type, atom) in enumerate(zip(types,atms)):
            (x,y,z) = (atom.GetX(), atom.GetY(), atom.GetZ())
            s += "%s %20.12f %20.12f %20.12f\n" % (self._elements.GetSymbol(type),
                                                   x, y, z)
        return s

    def writeFile(self, filename):
        """Dumps all fragments to individual
             .xyz files.
        """
        ff,ext = getFilenameAndExtension(filename)
        filename_template = "%s_%s_%03i%s"

        # first we dump all capped fragments
        for ifg,fragment in enumerate(self._fragmentation.getFragments()):
            (atms, types) = self._build_single_fragment(fragment)
            ss = self.fragment_xyz(atms, types)
            with open(filename_template % (ff,"fragment",ifg+1,ext), "w") as f:
                f.write(ss)
