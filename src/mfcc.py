"""
**********************************************************************
mfcc.py

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

#from util import *
#from config import FragItConfig

#from fragmentation import Fragmentation

class Cap(object):
    def __init__(self, atoms, ids, nucz, nbrs):
        self._atoms = atoms[:]
        self._ids = ids[:]
        self._nucz = nucz[:]
        self._nbrs = nbrs[:]
        self._charge = 0

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
            OBAtoms
            IDs of the OBAtoms
            Type (or nuclear charge) of the atom
            Neighbour ID of the atom - this is used to identify hydrogens that needs to be repositioned.
        """
        cap_atm = [self._fragmentation.getOBAtom(id) for id in pair]
        cap_ids = [a.GetIdx() for a in cap_atm]
        cap_typ = [a.GetAtomicNum() for a in cap_atm]
        cap_nbs = [-1 for a in cap_atm]
        order = 0
        while order < self._mfcc_order:
            order += 1
            cap_atm, cap_ids, cap_typ, cap_nbs = self._extend_cap(cap_atm, cap_ids, cap_typ, cap_nbs, order == self._mfcc_order)

        return Cap( cap_atm, cap_ids, cap_typ, cap_nbs)

        #return (cap_atm, cap_ids, cap_typ, 0, cap_nbs)

    def _extend_cap(self, atms, ids, typs, nbs, is_final_cap):
        """Extends the current cap with neighboring atoms.
           if this is_final_cap then atoms are hydrogens. they will
           OPTIONALLY be translated later.
        """
        atms_out = atms[:]
        ids_out  = ids[:]
        typs_out = typs[:]
        nbs_out = nbs[:]
        for atom in atms:
            for atomext in openbabel.OBAtomAtomIter(atom):
                if atomext in atms: continue
                atms_out.append(atomext)
                ids_out.append(atomext.GetIdx())
                nbs_out.append(atom.GetIdx())
                if is_final_cap:
                    typs_out.append(1)
                else:
                    typs_out.append(atomext.GetAtomicNum())
        return atms_out[:], ids_out[:], typs_out[:], nbs_out[:]
