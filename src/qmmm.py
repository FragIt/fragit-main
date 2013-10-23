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
from util import ravel2D, listDiff
#from config import FragItConfig

#from fragmentation import Fragmentation

class QMMM(object):
    """ Actions related to performing additional refinement specific to QM/MM codes
    """

    def __init__(self, fragmentation, qmlist):
        self._fragmentation = fragmentation
        self._qmfrags = map(int, qmlist)
        retract = lambda L: [l-1 for l in L]
        self._qmfrags = retract(self._qmfrags)
        self._qmfrags.sort()
        self._qmfrags.reverse()

    def pop_qm_fragment(self):
        """ Remove the qm fragments from the fragmentation. Adds hydrogens to both the
            qm-fragment that is returned and to the neighbouring fragments if bonds were cut.
        """

        # at this point, fragments have been generated
        fragments = self._fragmentation.getFragments()

        fragments_for_qm_no_hydrogens = []
        # instead of removing the qm-fragment, we rather adopt a quite novel approach
        # in which we signal to the user of the API that the fragment is not to be used
        # further by setting its original atom numbers to -1
        for idx in self._qmfrags:
            #print("FRAGIT: removing fragment {0}".format(idx))
            old_fragment = fragments.pop(idx) 
            fragments.insert(idx, [-1 for i in old_fragment])
            fragments_for_qm_no_hydrogens.insert(0,old_fragment[:])

        # for simplicity, let us just squash the qm fragments into one big fragment
        fragments_for_qm_no_hydrogens = ravel2D(fragments_for_qm_no_hydrogens)

        breaks = self._fragmentation.getExplicitlyBreakAtomPairs()
        if len(breaks) == 0: return fragments_for_qm_no_hydrogens

        # below here: add hydrogens to qm-fragment and to the rest of the capped structure
        fragment_for_qm = fragments_for_qm_no_hydrogens[:]

        #print("FRAGIT: [qm] {0}".format(fragment_for_qm))

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
        return fragment_for_qm

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
