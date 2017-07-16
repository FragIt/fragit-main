"""
Copyright (C) 2017 Casper Steinmann
"""
import os
import sys
import unittest

from src.config import FragItDataPE
from src.fragmentation import Fragmentation
from src.util import fileToMol
from src.qmmm import QMMM

class TestQMMMModule(unittest.TestCase):

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


    def test_qmmm_simple_0(self):
        """ qm/mm interface faulty parameters """
        self.assertRaises(TypeError, QMMM, None, 1)
        self.assertRaises(TypeError, QMMM, None, 1.0)
        self.assertRaises(TypeError, QMMM, None, "test")

        # clean break with zero elements in QM region
        self.assertRaises(ValueError, QMMM, None, [])


    def test_qmmm_simple_1(self):
        """ 1 qm fragment in qm/mm interface without bond breaking """
        filename = "temp.inp"
        molecule = fileToMol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()

        ref_fragments = fragmentation.getFragments()

        qmmm = QMMM(fragmentation, [1]) # extracts first fragment (for user this is #1)
        self.assertEqual(len(ref_fragments),4)
        qmfrag, qmcharge = qmmm.pop_qm_fragment()

        new_fragments = fragmentation.getFragments()

        self.assertEqual(qmfrag, [1,2,3]) # extracts first fragment
        self.assertEqual(len(new_fragments),4) # must leave all fragments intact
        self.assertEqual(new_fragments[0], [-1, -1, -1]) # -1 signals "DO NOT USE"


    def test_qmmm_simple_2(self):
        """ 2 qm fragments in qm/mm interface without bond breaking """
        filename = "temp.inp"
        molecule = fileToMol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()

        ref_fragments = fragmentation.getFragments()
        self.assertEqual(len(ref_fragments),4)

        qmmm = QMMM(fragmentation, [1,3]) # extracts fragments 1 and 3
        qmfrag, qmcharge = qmmm.pop_qm_fragment()

        new_fragments = fragmentation.getFragments()

        self.assertEqual(qmfrag, [1,2,3,7,8,9])
        self.assertEqual(len(new_fragments),4)
        self.assertEqual(new_fragments[0], [-1, -1, -1])
        self.assertEqual(new_fragments[2], [-1, -1, -1])

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestQMMModule))
