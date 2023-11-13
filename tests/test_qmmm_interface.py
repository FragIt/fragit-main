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
from fragit.qmmm import QMMM

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
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        ref_fragments = fragmentation.get_fragments()
        self.assertEqual(len(ref_fragments),4)

        qmmm = QMMM(fragmentation, [1]) # extracts first fragment (for user this is #1)
        qmfrag, qmcharge = qmmm.pop_qm_fragment()

        new_fragments = fragmentation.get_fragments()

        self.assertEqual(qmfrag, [1,2,3]) # extracts first fragment
        self.assertEqual(len(new_fragments),4) # must leave all fragments intact
        self.assertEqual(new_fragments[0], [-1, -1, -1]) # -1 signals "DO NOT USE"


    def test_qmmm_simple_2(self):
        """ 2 qm fragments in qm/mm interface without bond breaking """
        filename = "temp.inp"
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        ref_fragments = fragmentation.get_fragments()
        self.assertEqual(len(ref_fragments),4)

        qmmm = QMMM(fragmentation, [1,3]) # extracts fragments 1 and 3
        qmfrag, qmcharge = qmmm.pop_qm_fragment()
        new_fragments = fragmentation.get_fragments()

        self.assertEqual(qmfrag, [1,2,3,7,8,9])
        self.assertEqual(len(new_fragments),4)
        self.assertEqual(new_fragments[0], [-1, -1, -1])
        self.assertEqual(new_fragments[2], [-1, -1, -1])

    def test_qmmm_advanced_1(self):
        """ 1 qm fragment in qm/mm interface with bond breaking """
        filename = "temp.inp"
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        ref_fragments = copy.deepcopy(fragmentation.get_fragments()[:])
        self.assertEqual(len(ref_fragments),5)

        qmmm = QMMM(fragmentation, [1]) # extracts first fragment (for user this is #1)
        qmfrag, qmcharge = qmmm.pop_qm_fragment()
        new_fragments = fragmentation.get_fragments()[:]

        self.assertEqual(len(new_fragments), 5)
        self.assertEqual(len(qmfrag), len(ref_fragments[0])+1) # one atom added to QM region
        self.assertEqual(len(new_fragments[1]), len(ref_fragments[1])+1) # one atom added to MM region

    def test_qmmm_advanced_2(self):
        """ 2 qm fragments in qm/mm interface with bond breaking """
        filename = "temp.inp"
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        ref_fragments = copy.deepcopy(fragmentation.get_fragments()[:])
        self.assertEqual(len(ref_fragments),5)

        qmmm = QMMM(fragmentation, [2,3])
        qmfrag, qmcharge = qmmm.pop_qm_fragment()
        new_fragments = copy.deepcopy(fragmentation.get_fragments()[:])

        self.assertEqual(len(new_fragments), 5)
        self.assertEqual(len(qmfrag), len(ref_fragments[1]) + len(ref_fragments[2])+2) # two atoms added to QM region from two fragments
        self.assertEqual(len(new_fragments[0]), len(ref_fragments[0])+1) # one atom added to MM region for one fragment
        self.assertEqual(len(new_fragments[3]), len(ref_fragments[3])+1) # one atom added to MM region for one fragment

    def test_qmmm_advanced_3(self):
        """ 1 qm fragment + cov bound in qm/mm interface with bond breaking """
        filename = "temp.inp"
        m1 = file_to_mol("tests/5ala.xyz")
        frag1 = Fragmentation(m1, defaults=FragItDataPE)
        frag1.begin_fragmentation()
        frag1.do_fragmentation()
        frag1.finish_fragmentation()
        qm1 = QMMM(frag1, [2,3,4])
        qmf1, qmq1 = qm1.pop_qm_fragment()

        m2 = file_to_mol("tests/5ala.xyz")
        frag2 = Fragmentation(m2, defaults=FragItDataPE)
        frag2.values['qmmm']['includecovalent'] = True
        frag2.begin_fragmentation()
        frag2.do_fragmentation()
        frag2.finish_fragmentation()
        qm2 = QMMM(frag2, [3])
        qmf2, qmq2 = qm2.pop_qm_fragment()

        # must produce identical results
        self.assertEqual(len(qmf1), len(qmf2))
        self.assertEqual(frag1.get_fragments(), frag2.get_fragments())

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestQMMModule))
