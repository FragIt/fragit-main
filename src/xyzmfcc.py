"""
Copyright (C) 2013-2017 Casper Steinmann
"""
import numpy

from .mfcc import MFCC, Cap
from .pymol import PymolTemplate
from .jmol import JmolTemplate
from .writer import Standard
from .util import getFilenameAndExtension
from .util import shares_elements, calculate_hydrogen_position


class XYZMFCC(Standard):
    def __init__(self, fragmentation, directories):
        Standard.__init__(self,fragmentation, directories)
        self._mfcc = MFCC(fragmentation)

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
        return numpy.array([1 for i in fragments])

    def _setupActiveFragmentsInformation(self):
        self._active_atoms = []

    def _dump_pymol(self):
        pt = PymolTemplate(self._directories, self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _dump_jmol(self):
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

    def _build_single_fragment(self, fragment, caps):
        atomnames = ["" for i in fragment]
        if -1 in fragment:
            atoms = [None for i in fragment]
            nucz = [0 for a in atoms]
            neighbours = [-1 for a in atoms]
            ids = [-1 for a in atoms]
        else:
            atoms = [self._fragmentation.getOBAtom(i) for i in fragment]
            if self._fragmentation.hasAtomNames():
                names = self._fragmentation.getAtomNames()
                atomnames = []
                for i in fragment:
                    try:
                        atomnames.append(names[i-1])
                    except IndexError:
                        print("Warning: FragIt could not correctly name atom {0:d}.".format(i))
                        print("         The problem could be with your PDB file.")
                        atomnames.append("X")
            nucz  = [a.GetAtomicNum() for a in atoms]
            neighbours = [-1 for a in atoms]
            ids = [i for i in fragment]

        if caps is not None:
            for icap,cap in enumerate(caps):
                if shares_elements( fragment, cap.getAtomIDs() ):
                    for id,atom,atomname,z,nbr in zip(cap.getAtomIDs(), cap.getAtoms(), cap.getAtomNames(), cap.getNuclearCharges(), cap.getNeighbourList() ):
                        if id not in fragment:
                            atoms.append( atom )
                            atomnames.append( atomname )
                            nucz.append( z )
                            neighbours.append( nbr )
                            ids.append( id )

        return Cap(atoms, atomnames, ids, nucz, neighbours)

    def getCaps(self):
        return self._mfcc.getCaps()

    def BuildCappedFragment(self, fragment):
        return self._build_single_fragment(fragment, self.getCaps())

    def BuildFragment(self, fragment):
        return self._build_single_fragment(fragment, None)

    def _fragment_xyz(self, fragment ):
        """Generates the xyz file format based on the atoms, types,
           ids and neighbours of each fragment
        """
        # NB! the word fragment here is actually of type Cap. Just to be sure
        # nobody is actually doing something utterly wrong, check that here.
        if not isinstance(fragment, Cap):
            raise TypeError("_fragment_xyz expected an object of type Cap.")
        atoms = fragment.getAtoms()
        nuczs = fragment.getNuclearCharges()
        nbrls = fragment.getNeighbourList()

        n = len(atoms)
        s = "%i\n%s\n" % (n,"")
        for id, (atom, nucz, neighbour) in enumerate(zip(atoms,nuczs,nbrls)):
            (x,y,z) = (atom.GetX(), atom.GetY(), atom.GetZ())
            if atom.GetAtomicNum() != nucz:
                # atom is the light atom and it is connected to the nbrs[id] atom
                heavy_atom = self._fragmentation.getOBAtom( neighbour )
                (x,y,z) = calculate_hydrogen_position( heavy_atom, atom )
            s += "%s %20.12f %20.12f %20.12f\n" % (self._elements.GetSymbol(nucz),
                                                   x, y, z)
        return s

    def writeFile(self, filename):
        """Dumps all caps and capped fragments to individual files
        """
        ff,ext = getFilenameAndExtension(filename)
        filename_template = "{0}_{1}_{2:03d}{3}"

        # these are the capped fragments
        for ifg, fragment in enumerate(self._fragmentation.getFragments(), start=1):
            capped_fragment = self.BuildCappedFragment( fragment )
            ss = self._fragment_xyz( capped_fragment )
            with open( filename_template.format(ff, "fragment", ifg, ext), 'w' ) as f:
                f.write(ss)

        # these are the caps
        for icap, cap in enumerate(self.getCaps(), start=1):
            ss = self._fragment_xyz(cap)
            with open( filename_template.format(ff, "cap", icap, ext), 'w' ) as f:
                f.write(ss)
