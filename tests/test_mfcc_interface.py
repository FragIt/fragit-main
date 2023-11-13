"""
Copyright (C) 2017 Casper Steinmann
"""
import copy
import os
import sys
import unittest

from fragit.config import FragItDataPE
from fragit.fragmentation import Fragmentation
from fragit.util import file_to_mol
from fragit.xyzmfcc import XYZMFCC

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
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        ref_fragments = fragmentation.get_fragments()
        self.assertEqual(len(ref_fragments),4)

        mfcc = XYZMFCC(fragmentation, {})
        fragments = [mfcc.build_fragment(fragment) for fragment in fragmentation.get_fragments()]
        caps = mfcc.get_caps()
        capped_fragments = [mfcc.build_capped_fragment(fragment) for fragment in fragmentation.get_fragments()]

        self.assertEqual(len(fragments), len(capped_fragments))
        self.assertEqual(len(caps), 0) # no caps when fragmenting water
        for frag, cfrag in zip(fragments, capped_fragments):
            fset = set(frag.get_atom_ids())
            cset = set(cfrag.get_atom_ids())
            self.assertEqual(len(cset.difference(fset)), 0)

    def test_mfcc_with_covalent_bonds(self):
        """ Tests MFCC with covalent bonds """
        filename = "temp.inp"
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        ref_fragments = fragmentation.get_fragments()
        self.assertEqual(len(ref_fragments),5)

        mfcc = XYZMFCC(fragmentation, {})
        fragments = [mfcc.build_fragment(fragment) for fragment in fragmentation.get_fragments()]
        caps = mfcc.get_caps()
        capped_fragments = [mfcc.build_capped_fragment(fragment) for fragment in fragmentation.get_fragments()]

        self.assertEqual(len(fragments), len(capped_fragments))
        self.assertEqual(len(caps), 4) # 4 internal caps in A-A-A-A because 4 bonds broken
        results = [6, 12, 12, 12, 6] # number of additional atoms per fragment when capped
        for i, (frag, cfrag) in enumerate(zip(fragments, capped_fragments)):
            fset = set(frag.get_atom_ids())
            cset = set(cfrag.get_atom_ids())
            self.assertEqual(len(cset.difference(fset)), results[i])

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestMFCCModule))
