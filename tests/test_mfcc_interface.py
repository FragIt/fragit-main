"""
Copyright (C) 2017 Casper Steinmann
"""
import copy
import os
import sys
import unittest

from src.config import FragItDataPE
from src.fragmentation import Fragmentation
from src.util import fileToMol
from src.xyzmfcc import XYZMFCC

class TestMFCCModule(unittest.TestCase):
    def setUp(self):
        pass

    def delete_file(self,filename):
        try:
            f = open(filename)
        except IOError:
            return
        finally:
            f.close()
            os.remove(filename)

    def test_mfcc_no_covalent_bonds(self):
        """ Tests MFCC without covalent bonds """
        filename = "temp.inp"
        molecule = fileToMol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()

        ref_fragments = fragmentation.getFragments()
        self.assertEqual(len(ref_fragments),4)

        mfcc = XYZMFCC(fragmentation, {})
        fragments = [mfcc.BuildFragment(fragment) for fragment in fragmentation.getFragments()]
        caps = mfcc.getCaps()
        capped_fragments = [mfcc.BuildCappedFragment(fragment) for fragment in fragmentation.getFragments()]

        self.assertEqual(len(fragments), len(capped_fragments))
        self.assertEqual(len(caps), 0) # no caps when fragmenting water
        for frag, cfrag in zip(fragments, capped_fragments):
            fset = set(frag.getAtomIDs())
            cset = set(cfrag.getAtomIDs())
            self.assertEqual(len(cset.difference(fset)), 0)

    def test_mfcc_with_covalent_bonds(self):
        """ Tests MFCC with covalent bonds """
        filename = "temp.inp"
        molecule = fileToMol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()

        ref_fragments = fragmentation.getFragments()
        self.assertEqual(len(ref_fragments),5)

        mfcc = XYZMFCC(fragmentation, {})
        fragments = [mfcc.BuildFragment(fragment) for fragment in fragmentation.getFragments()]
        caps = mfcc.getCaps()
        capped_fragments = [mfcc.BuildCappedFragment(fragment) for fragment in fragmentation.getFragments()]

        self.assertEqual(len(fragments), len(capped_fragments))
        self.assertEqual(len(caps), 4) # 4 internal caps in A-A-A-A because 4 bonds broken
        results = [6, 12, 12, 12, 6] # number of additional atoms per fragment when capped
        for i, (frag, cfrag) in enumerate(zip(fragments, capped_fragments)):
            fset = set(frag.getAtomIDs())
            cset = set(cfrag.getAtomIDs())
            self.assertEqual(len(cset.difference(fset)), results[i])

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestMFCCModule))
