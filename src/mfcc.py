"""
Copyright (C) 2013-2016 Casper Steinmann
"""

try:
    import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")
import numpy

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

    def getAtoms(self):
        return self._atoms

    def getCharge(self):
        return self._charge

    def getNuclearCharges(self):
        return self._nucz

    def getNeighbourList(self):
        return self._nbrs

    def getAtomIDs(self):
        return self._ids

    def getAtomNames(self):
        return self._atom_names

    def setCharge(self, value):
        self._charge = value

    def doRecalculation(self):
        self._recalculate = True

    def getRecalculationState(self):
        return self._recalculate

    def doIgnore(self, value):
        self._ignore = value

    def getIgnore(self):
        return self._ignore

    def getAtomNames(self):
        return self._atom_names

class MFCC(object):

    def __init__(self, fragmentation):
        self._fragmentation = fragmentation
        self._caps = []
        self._identifyCaps()

    def _identifyCaps(self):
        self._mfcc_order = 0
        if self._fragmentation.getOutputFormat() == 'XYZ-MFCC':
            self._mfcc_order = self._fragmentation.getMFCCOrder()
            if self._mfcc_order <= 0:
                raise ValueError("You must specify the order of capping.")
            self._build_caps()

    def getCaps(self):
        return self._caps

    def hasCaps(self):
        return len(self._caps) > 0

    def _build_caps(self):
        for pair in self._fragmentation.getExplicitlyBreakAtomPairs():
            self._caps.append( self._build_cap(pair) )

    def _build_cap(self, pair):
        """ Builds a cap around a pair of fragmentation points.

            a cap contains the following information:
              a list of OBAtoms that makes out the cap
              IDs of the neighboring OBAtoms
              Type (or nuclear charge) of the atom
              Neighbour ID of the atom - this is used to identify hydrogens that needs to be repositioned later.
        """
        cap_atm = [self._fragmentation.getOBAtom(id) for id in pair]
        cap_atmnam = ["" for id in pair]
        if self._fragmentation.hasAtomNames():
            cap_atmnam = [self._fragmentation._atom_names[id-1] for id in pair]

        cap_ids = [a.GetIdx() for a in cap_atm]
        cap_typ = [a.GetAtomicNum() for a in cap_atm]
        cap_nbs = [-1 for a in cap_atm]
        order = 0
        while order < self._mfcc_order:
            order += 1
            cap_atm, cap_ids, cap_typ, cap_nbs, cap_atmnam = self._extend_cap(cap_atm, cap_ids, cap_typ, cap_nbs, cap_atmnam, order == self._mfcc_order)

        return Cap( cap_atm, cap_atmnam, cap_ids, cap_typ, cap_nbs)


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
        ids_out  = ids[:]
        typs_out = typs[:]
        nbs_out = nbs[:]
        for atom in atms:
            for atomext in openbabel.OBAtomAtomIter(atom):
                if atomext in atms: continue
                atms_out.append(atomext)
                ids_out.append(atomext.GetIdx())
                nbs_out.append(atom.GetIdx())
                if self._fragmentation.hasAtomNames():
                    if not is_final_cap:
                        atm_namout.append( self._fragmentation._atom_names[atomext.GetIdx() -1] )
                    else:
                        # we are sure that when it is a final cap, it is
                        # a hydrogen.
                        atm_namout.append( " H  " )
                else:
                    atm_namout.append( "   " )

                if is_final_cap:
                    typs_out.append(1)
                else:
                    typs_out.append(atomext.GetAtomicNum())
        return atms_out[:], ids_out[:], typs_out[:], nbs_out[:], atm_namout[:]
