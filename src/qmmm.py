"""
**********************************************************************
qmmm.py

Copyright (C) 2013 Casper Steinmann

This file is part of the FragIt project.

FragIt is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

FragIt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.
***********************************************************************/
"""
#import os
#import sys
#import logging

import openbabel
import numpy

from util import calculate_hydrogen_position
from util import ravel2D, listDiff, Uniqify
#from config import FragItConfig

#from fragmentation import Fragmentation

class QMMM(object):
    """ Actions related to performing additional refinement specific to QM/MM codes
    """

    def __init__(self, fragmentation, qmlist):
        self._fragmentation = fragmentation
        self._qmfrags = map(int, qmlist)
        self._fd = FragmentDistances(fragmentation)
        retract = lambda L: [l-1 for l in L]
        self._qmfrags = retract(self._qmfrags)
        self._qmfrags.sort()
        self._qmfrags.reverse()
        self._qm_charges = fragmentation.getFragmentCharges()

    def pop_qm_fragment(self):
        """ Remove the qm fragments from the fragmentation. Adds hydrogens to both the
            qm-fragment that is returned and to the neighbouring fragments if bonds were cut.
        """

        # at this point, fragments have been generated
        fragments = self._fragmentation.getFragments()

        # investigate which fragments should be included in the QM-region
        qmfrags = self._qmfrags[:]
        distance = self._fragmentation.getActiveAtomsDistance()
        qmfrags = Uniqify(self._addQM(qmfrags))
        qmfrags.sort()
        qmfrags.reverse()

        qm_region_charge = 0

        fragments_for_qm_no_hydrogens = []

        # instead of removing the qm-fragment, we rather adopt a quite novel approach
        # in which we signal to the user of the API that the fragment is not to be used
        # further by setting its original atom numbers to -1
        for idx in qmfrags:
            qm_region_charge += self._qm_charges[idx]
            old_fragment = fragments.pop(idx) 
            fragments.insert(idx, [-1 for i in old_fragment])
            fragments_for_qm_no_hydrogens.insert(0,old_fragment[:])

        # for simplicity, let us just squash the qm fragments into one big fragment
        fragments_for_qm_no_hydrogens = ravel2D(fragments_for_qm_no_hydrogens)

        breaks = self._fragmentation.getExplicitlyBreakAtomPairs()
        if len(breaks) == 0: return fragments_for_qm_no_hydrogens

        # below here: add hydrogens to qm-fragment and to the rest of the capped structure
        fragment_for_qm = fragments_for_qm_no_hydrogens[:]

        # first, fix the QM-fragment, removing any bond-breaks that reside
        # inside (or bordering) the qm-fragment. if breaks are bordering,
        # add appropriate hydrogen atoms.
        lenqmfrags = len(fragments_for_qm_no_hydrogens)
        remove_breaks = []
        for ibreak, bbreak in enumerate(breaks):
            lendiff = len(listDiff(fragments_for_qm_no_hydrogens, list(bbreak)))
            difflen = lenqmfrags - lendiff

            # the break is not present in the qm-region, leave it
            if difflen == 0: continue
            # mark the ibreak'th item for removal. no hydrogens to be added
            if difflen == 2: remove_breaks.append(ibreak)
            # the break is bordering the qm-region and the mm-region
            if difflen == 1:

                for iibreak in bbreak:
                    # fix the qm-fragment first
                    if iibreak in fragments_for_qm_no_hydrogens:
                        new_atoms = self.satisfyValency(fragments_for_qm_no_hydrogens, iibreak, bbreak)
                        if len(new_atoms) > 0: fragment_for_qm.extend(new_atoms)

                    # then fix the fragments themselves
                    # INFO/WARNING: this is a lists of lists thing. BE CAREFULL
                    for ifragment, fragment in enumerate(fragments):
                        if iibreak in fragment:
                            new_atoms = self.satisfyValency(fragment, iibreak, bbreak)
                            fragments[ifragment].extend(new_atoms)

                # also mark the ibreak'th item for removal
                remove_breaks.append(ibreak)

        for ibreak in remove_breaks:
            self._fragmentation.popExplicitlyBreakAtomPairs(breaks[ibreak])

        #print("FRAGIT: [qm] {0}".format(fragment_for_qm))
        return (fragment_for_qm, qm_region_charge)

    def satisfyValency(self, fragment, iheavy, bbreak):
        """ Satisfies the valency of atom number iheavy in the supplied fragment. Returns a new fragment with all atoms in the correct place.
        """

        # heavy is the heavy atom that wants a hydrogen
        heavy = self._fragmentation.getOBAtom(iheavy)
        ilight = 0
        if bbreak.index(iheavy) == 0: ilight = 1
        light = self._fragmentation.getOBAtom(bbreak[ilight])
        ival = heavy.GetImplicitValence()
        rval = heavy.GetValence()
        new_atoms = []
        if ival != rval:
            if self._fragmentation.mol.AddHydrogens(heavy):
                for nbatom in openbabel.OBAtomAtomIter(heavy):
                    if nbatom.GetIdx() not in fragment:
                        x,y,z = calculate_hydrogen_position(heavy, light)
                        nbatom.SetVector(x,y,z)
                        new_atoms.append(nbatom.GetIdx())
        return new_atoms

    def _addQM(self, qmfragments):
        qmfrags = qmfragments[:]
        for qmfragment in qmfragments:
            qmfrags.extend( self._fd.getHydrogenBoundFragments(qmfragment) )
            qmfrags.extend( self._fd.getCovalentlyBoundFragments(qmfragment) )
        return qmfrags

class FragmentDistances(object):
    def __init__(self, fragmentation):
        """ Fragments   : array of arrays of atomic numbers [[a1,a2,a3],[a4,a5,a6], ... ]
        """
        self._fragmentation = fragmentation
        self._fragments = fragmentation.getFragments()
        self._fragment_donors = [self._donors_from_fragment(i) for i in range(len(self._fragments))]
        self._fragment_acceptors = [self._acceptors_from_fragment(i) for i in range(len(self._fragments))]

        #print len(self._fragment_acceptors), len(self._fragments)

    def getHydrogenBoundFragments(self, idx):
        """ Obtains all hydrogen bound fragments to the idx'th fragment
        """
        hb_fragments = []
        donors = []
        acceptors = []

        if self._fragmentation.doQMMMHydrogenBondDonors():
            donors = self._fragment_donors[idx]
        if self._fragmentation.doQMMMHydrogenBondAcceptors():
            acceptors = self._fragment_acceptors[idx]

        # first we will find any donor (current fragment) -> acceptor (whole system) pairs
        isDonor = False
        for ifg, _acceptors in enumerate(self._fragment_acceptors):
            if ifg == idx: continue
            for A in _acceptors:
                if isDonor: break
                for H,D in donors:
                    isDonor = self._isHydrogenBond(D, H, A)
                    if isDonor:
                        hb_fragments.append( ifg )
                        break
            isDonor = False

        # then we will find acceptors (current fragment) -> donor (whole system)
        isAcceptor = False
        for ifg, _donors in enumerate(self._fragment_donors):
            if ifg == idx: continue
            for H,D in _donors:
                if isAcceptor: break
                for A in acceptors:
                    isAcceptor = self._isHydrogenBond(D, H, A)
                    if isAcceptor:
                        hb_fragments.append( ifg )
                        break
            isAcceptor = False

        return hb_fragments

    def getCovalentlyBoundFragments(self, idx):
        atomids = self._fragments[idx]
        bonds = self._fragmentation.getExplicitlyBreakAtomPairs()

        other_fragment_atomids = []
        other_fragments = []

        # let us find the nearby fragments that are covalently connected
        # currently, this only works with nearest neighbours
        for (l,r) in bonds:
            if l in atomids:
                other_fragment_atomids.append(r)
            if r in atomids:
                other_fragment_atomids.append(l)

        # lets find the associated fragments
        for ifg,atoms in enumerate(self._fragments):
            if len(atoms) > len(listDiff(atoms, other_fragment_atomids)):
                other_fragments.append(ifg)

        return other_fragments

    def _isHydrogenBond(self, D, H, A):
        """ angle > 110 is according to:
            http://pac.iupac.org/publications/pac/pdf/2011/pdf/8308x1637.pdf

            D:   donor
            H:   hydrogen
            A:   acceptor
        """
        RAH = A.GetDistance(H)
        RAD = A.GetDistance(D)
        if RAH < self._fragmentation.getHBondDistanceMin() and RAD < self._fragmentation.getHBondDistanceMax():
            ANG = D.GetAngle(H,A) # in degrees
            if ANG > self._fragmentation.getHBondAngle():
                return True
        return False

    def _donors_from_fragment(self, idx):
        """ Searches a for a hydrogen-bond donor (aha a hydrogen). we consider only
            X --- H-N
            X --- H-O
            as potential donors. the atoms will be returned as a list of tuples (H, ?)
        """
        obatoms = [self._fragmentation.getOBAtom(i) for i in self._fragments[idx]]
        # lets locate all the hydrogen atoms connected to a nitrogen or oxygen
        donors = []
        for obatom in obatoms:
            if obatom.IsHydrogen() and obatom.IsHbondDonorH():
                for otheratom in obatoms:
                    if obatom == otheratom: continue
                    if obatom.IsConnected(otheratom) and otheratom.IsHbondDonor():
                        donors.append((obatom, otheratom))
                        break

        return donors

    def _acceptors_from_fragment(self, idx):
        """ Searches a for a hydrogen-bond acceptor (aha a carboxyl oxygen). we consider only
            X --- O=C
            as potential donors. the atoms will be returned
        """
        obatoms = [self._fragmentation.getOBAtom(i) for i in self._fragments[idx]]
        acceptors = []
        for obatom in obatoms:
            if obatom.IsHbondAcceptor():
                acceptors.append(obatom)
        return acceptors

    def _is_doner(self, atom):
        return atom.IsNitrogen() or atom.IsOxygen()
