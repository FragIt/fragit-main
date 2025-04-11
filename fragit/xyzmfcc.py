"""
Copyright (C) 2013-2023 Casper Steinmann
"""
import numpy as np
from typing import List, Optional

from openbabel import openbabel

from fragit.mfcc import MFCC, Cap
from fragit.pymol import PymolTemplate
from fragit.jmol import JmolTemplate
from fragit.writer import Standard
from fragit.util import get_filename_and_extension, Z2LABEL
from fragit.util import shares_elements, calculate_hydrogen_position


class XYZMFCC(Standard):
    def __init__(self, fragmentation, directories):
        Standard.__init__(self, fragmentation, directories)
        self._mfcc = MFCC(fragmentation)

    def setup(self):
        self._setup_layered_information()
        self._setup_active_fragments_information()
        if self._do_pymol: self._dump_pymol()
        if self._do_jmol: self._dump_jmol()

    def _setup_layered_information(self):
        self._fragment_layers = self._get_fragment_layers_from_fragment()

    def _get_fragment_layers_from_fragment(self):
        fragments = self._fragmentation.get_fragments()
        return np.array([1 for i in fragments])

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

    def _set_template_data(self, template):
        template.set_fragments_data(self._fragmentation.get_fragments())
        template.set_buffer_data(self._fragment_layers)
        template.set_active_data(self._active_atoms)
        template.set_backbone_data(self._fragmentation.get_backbone_atoms())
        template.set_pair_data(self._fragmentation.get_explicitly_break_atom_pairs())

    def _write_template_file(self, template):
        template.override()
        template.write()

    def _build_single_fragment(self, fragment, caps: Optional[List[Cap]]):
        atomnames = ["" for _ in fragment]
        atoms: List[Optional[openbabel.OBAtom]]
        if -1 in fragment:
            atoms = [None for _ in fragment]
            nucz = [0 for _ in atoms]
            neighbours = [-1 for _ in atoms]
            ids = [-1 for _ in atoms]
        else:
            atoms = [self._fragmentation.get_ob_atom(i) for i in fragment]
            if self._fragmentation.has_atom_names():
                names = self._fragmentation.get_atom_names()
                atomnames = []
                for i in fragment:
                    try:
                        atomnames.append(names[i-1])
                    except IndexError:
                        print("Warning: FragIt could not correctly name atom {0:d}.".format(i))
                        print("         The problem could be with your PDB file.")
                        atomnames.append("X")
            nucz = [a.GetAtomicNum() for a in atoms if a is not None]
            neighbours = [-1 for _ in atoms]
            ids = [i for i in fragment]

        if caps is not None:
            for icap, cap in enumerate(caps):
                if shares_elements(fragment, cap.get_atom_ids()):
                    for index, atom, atomname, z, nbr in zip(cap.get_atom_ids(), cap.get_atoms(), cap.get_atom_names(), cap.get_nuclear_charges(), cap.get_neighbour_list()):
                        if index not in fragment:
                            atoms.append(atom)
                            atomnames.append(atomname)
                            nucz.append(z)
                            neighbours.append(nbr)
                            ids.append(index)

        return Cap(atoms, atomnames, ids, nucz, neighbours)

    def get_caps(self):
        return self._mfcc.get_caps()

    def build_capped_fragment(self, fragment):
        return self._build_single_fragment(fragment, self.get_caps())

    def build_fragment(self, fragment):
        return self._build_single_fragment(fragment, None)

    def _fragment_xyz(self, fragment: Cap):
        """Generates the xyz file format based on the atoms, types,
           ids and neighbours of each fragment
        """
        # NB! the word fragment here is actually of type Cap. Just to be sure
        # nobody is actually doing something utterly wrong, check that here.
        if not isinstance(fragment, Cap):
            raise TypeError("_fragment_xyz expected an object of type Cap.")
        atoms = fragment.get_atoms()
        nuczs = fragment.get_nuclear_charges()
        nbrls = fragment.get_neighbour_list()

        n = len(atoms)
        s = "%i\n%s\n" % (n,"")
        for id, (atom, nucz, neighbour) in enumerate(zip(atoms, nuczs, nbrls)):
            (x, y, z) = (atom.GetX(), atom.GetY(), atom.GetZ())
            if atom.GetAtomicNum() != nucz:
                # atom is the light atom and it is connected to the nbrs[id] atom
                heavy_atom = self._fragmentation.get_ob_atom(neighbour)
                (x, y, z) = calculate_hydrogen_position(heavy_atom, atom)
            s += "%s %20.12f %20.12f %20.12f\n" % (Z2LABEL[nucz],
                                                   x, y, z)
        return s

    def write_file(self, filename):
        """Dumps all caps and capped fragments to individual files
        """
        ff,ext = get_filename_and_extension(filename)
        filename_template = "{0}_{1}_{2:03d}{3}"

        # these are the capped fragments
        for ifg, fragment in enumerate(self._fragmentation.get_fragments(), start=1):
            capped_fragment = self.build_capped_fragment(fragment)
            ss = self._fragment_xyz(capped_fragment)
            with open(filename_template.format(ff, "fragment", ifg, ext), "w") as f:
                f.write(ss)

        # these are the caps
        for icap, cap in enumerate(self.get_caps(), start=1):
            ss = self._fragment_xyz(cap)
            with open(filename_template.format(ff, "cap", icap, ext), "w" ) as f:
                f.write(ss)
