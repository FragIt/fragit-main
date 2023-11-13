"""
Copyright (C) 2013-2023 Casper Steinmann
"""

from fragit.fragit_exceptions import OBNotFoundException

try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")


class Cap(object):
    def __init__(self, atoms, atomnames, ids, nucz, nbrs):
        self._atoms = atoms[:]
        self._atom_names = atomnames[:]
        self._ids = ids[:]
        self._nucz = nucz[:]
        self._nbrs = nbrs[:]
        self._charge = 0
        self._recalculate = False
        self._ignore = False

    def get_atoms(self):
        return self._atoms

    def get_charge(self):
        return self._charge

    def get_nuclear_charges(self):
        return self._nucz

    def get_neighbour_list(self):
        return self._nbrs

    def get_atom_ids(self):
        return self._ids

    def get_atom_names(self):
        return self._atom_names

    def set_charge(self, value):
        self._charge = value

    def do_recalculation(self):
        self._recalculate = True

    def get_recalculation_state(self):
        return self._recalculate

    def do_ignore(self, value):
        self._ignore = value

    def get_ignore(self):
        return self._ignore


class MFCC(object):

    def __init__(self, fragmentation):
        self._fragmentation = fragmentation
        self._caps = []
        self._identify_caps()

    def _identify_caps(self):
        self._mfcc_order = 0
        if self._fragmentation.get_output_format() == 'XYZ-MFCC':
            self._mfcc_order = self._fragmentation.get_mfcc_order()
            if self._mfcc_order <= 0:
                raise ValueError("You must specify the order of capping.")
            self._build_caps()

    def get_caps(self):
        return self._caps

    def has_caps(self):
        return len(self._caps) > 0

    def _build_caps(self):
        for pair in self._fragmentation.get_explicitly_break_atom_pairs():
            self._caps.append(self._build_cap(pair))

    def _build_cap(self, pair):
        """ Builds a cap around a pair of fragmentation points.

            a cap contains the following information:
              a list of OBAtoms that makes out the cap
              IDs of the neighboring OBAtoms
              Type (or nuclear charge) of the atom
              Neighbour ID of the atom - this is used to identify hydrogens that needs to be repositioned later.
        """
        cap_atm = [self._fragmentation.get_ob_atom(index) for index in pair]
        cap_atmnam = ["" for _ in pair]
        if self._fragmentation.has_atom_names():
            cap_atmnam = [self._fragmentation._atom_names[index-1] for index in pair]

        cap_ids = [a.GetIdx() for a in cap_atm]
        cap_typ = [a.GetAtomicNum() for a in cap_atm]
        cap_nbs = [-1 for _ in cap_atm]
        order = 0
        while order < self._mfcc_order:
            order += 1
            cap_atm, cap_ids, cap_typ, cap_nbs, cap_atmnam = self._extend_cap(cap_atm, cap_ids, cap_typ, cap_nbs, cap_atmnam, order == self._mfcc_order)

        return Cap(cap_atm, cap_atmnam, cap_ids, cap_typ, cap_nbs)

    def _extend_cap(self, atms, ids, typs, nbs, nams, is_final_cap):
        """ Extends the current cap with neighboring atom IDs.

            If called with is_final_cap == True then atoms replaced with
            hydrogens and will OPTIONALLY be translated later.

            Arguments:
            atms -- the OpenBabel atoms currently in the cap
            ids -- the IDs of the atoms in the atms list
            typs -- the atom types of the atoms in the atms list
            nbs -- the indices of the neighbours of the cap
            nams -- the atom names of the atoms in the atms list
            is_final_cap -- a boolean that signals that hydrogens should be added if True

            Returns:
            updated atms, ids, typs, nbs, nams
        """
        atms_out = atms[:]
        atm_namout = nams[:]
        ids_out = ids[:]
        typs_out = typs[:]
        nbs_out = nbs[:]
        for atom in atms:
            for atomext in openbabel.OBAtomAtomIter(atom):
                if atomext in atms:
                    continue
                atms_out.append(atomext)
                ids_out.append(atomext.GetIdx())
                nbs_out.append(atom.GetIdx())
                if self._fragmentation.has_atom_names():
                    if not is_final_cap:
                        atm_namout.append(self._fragmentation._atom_names[atomext.GetIdx() - 1])
                    else:
                        # we are sure that when it is a final cap, it is
                        # a hydrogen.
                        atm_namout.append(" H  ")
                else:
                    atm_namout.append("   ")

                if is_final_cap:
                    typs_out.append(1)
                else:
                    typs_out.append(atomext.GetAtomicNum())
        return atms_out[:], ids_out[:], typs_out[:], nbs_out[:], atm_namout[:]
