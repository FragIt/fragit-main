"""
Copyright (C) 2011-2017 Casper Steinmann
"""
import os
import unittest
import openbabel

from fragit.fragmentation import Fragmentation
from fragit.util import fileToMol
from fragit.config import FragItDataFMO

class TestFragmentationModule(unittest.TestCase):

    def setUp(self):
        self.filename_pdb = "tests/1UAO.pdb"

        # for testing, use OpenBabel functionality directly
        self.molecule = fileToMol(self.filename_pdb)
        self.fragmentation = Fragmentation(self.molecule, defaults=FragItDataFMO)

    def tearDown(self):
        pass

    def delete_file(self,filename):
        try:
            f = open(filename)
        except IOError:
            return
        finally:
            f.close()
            os.remove(filename)

    def test_FragmentationDefaultParameters(self):
        frg = self.fragmentation
        self.assertEqual(frg.mol != None, True)

    def test_FragmentationSetActiveFragments(self):
        self.fragmentation.setActiveFragments([1, 2, 3])
        self.assertEqual(self.fragmentation.active_fragments, [1, 2, 3])

    def test_FragmentationApplySmartProtectPatterns(self):
        self.fragmentation.applySmartProtectPatterns()
        self.assertEqual(self.fragmentation.getExplicitlyProtectedAtoms(), [1, 2, 3, 4, 10])

    def test_FragmentationDetermineFormalCharges(self):
        #self.fragmentation.determineFormalCharges()
        self.assertAlmostEqual(sum(self.fragmentation.formalCharges), -2)

    def test_FragmentationGetProtectedAtoms(self):
        self.assertEqual(self.fragmentation.getExplicitlyProtectedAtoms(), [])

    def test_FragmentationAddProtectedAtomsAfterProtect(self):
        self.fragmentation.setProtectedAtoms()
        self.fragmentation.addExplicitlyProtectedAtoms([44, 55, 67])
        self.assertEqual(self.fragmentation.getExplicitlyProtectedAtoms(), [1, 2, 3, 4, 10, 44, 55, 67])

    def test_FragmentationAddBrokenBond(self):
        pass

    def test_FragmentationIsBondProtected(self):
        bond_pair = (2, 3)
        self.assertEqual(self.fragmentation.isBondProtected(bond_pair), False)
        self.fragmentation.addExplicitlyProtectedAtoms([2])
        self.assertEqual(self.fragmentation.isBondProtected(bond_pair), True)

    def test_FragmentationRealBondBreakerNoProtect(self):
        bond_atoms = (2, 3)
        self.fragmentation.addFragmentationAtomPair(bond_atoms)
        self.assertEqual(self.fragmentation.getExplicitlyBreakAtomPairs(), [(2, 3)])

    def test_FragmentationIsValidExplicitBond(self):
        self.assertRaises(ValueError, self.fragmentation.isValidExplicitBond, (1, 1))
        self.assertRaises(ValueError, self.fragmentation.isValidExplicitBond, (2, 4))

    def test_FragmentationBreakBondsWithNoProtect(self):
        self.fragmentation.breakBonds()
        self.assertEqual(len(self.fragmentation.getExplicitlyBreakAtomPairs()), 9)

    def test_FragmentationBreakBondsExplcitWithNoProtect(self):
        self.fragmentation.addExplicitlyBreakAtomPairs([(111, 112)])
        self.fragmentation.breakBonds()
        self.assertEqual(len(self.fragmentation.getExplicitlyBreakAtomPairs()), 10)

    def test_FragmentationBreakBondsWithProtect(self):
        self.fragmentation.setProtectedAtoms()
        self.fragmentation.breakBonds()
        self.assertEqual(len(self.fragmentation.getExplicitlyBreakAtomPairs()), 8)

    def test_FragmentationBreakBondsExplcitWithProtect(self):
        self.fragmentation.setProtectedAtoms()
        self.fragmentation.addExplicitlyBreakAtomPairs([(111, 112)])
        self.fragmentation.breakBonds()
        self.assertEqual(len(self.fragmentation.getExplicitlyBreakAtomPairs()), 9)

    def test_FragmentationDetermineFragmentsNoBreaking(self):
        self.assertRaises(ValueError, self.fragmentation.determineFragments)

    def test_FragmentationDetermineFragmentsWithBreaking(self):
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        for fragment, expected_size in zip(self.fragmentation.getFragments(), [7, 21, 12, 14, 15, 14, 7, 14, 24, 10]):
            self.assertEqual(len(fragment), expected_size)

    def test_FragmentationDetermineFragmentsWithBreakingAndGrouping(self):
        self.fragmentation.setFragmentGroupCount(2)
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.doFragmentGrouping()
        self.assertEqual(len(self.fragmentation.getFragments()), 5)
        for fragment, expected_size in zip(self.fragmentation.getFragments(), [28, 26, 29, 21, 34]):
            self.assertEqual(len(fragment), expected_size)

    def test_FragmentationDetermineFragmentsWithBreakingAndGroupingTriple(self):
        self.fragmentation.setFragmentGroupCount(3)
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.doFragmentGrouping()
        self.assertEqual(len(self.fragmentation.getFragments()), 4)
        for fragment, expected_size in zip(self.fragmentation.getFragments(), [40, 43, 45, 10]):
            self.assertEqual(len(fragment), expected_size)

    def test_FragmentationFindFragmentsCastErrors(self):
        self.assertRaises(ValueError, self.fragmentation.getAtomsInSameFragment, 1, 1)
        self.assertRaises(ValueError, self.fragmentation.getAtomsInSameFragment, "")

    def test_FragmentationFindFragmentsNoGrouping(self):
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        res1 = [12, 13, 31, 32]
        res1.extend(range(35, 43))
        self.assertEqual(self.fragmentation.getAtomsInSameFragment(11, 0), [3, 4, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])
        self.assertEqual(self.fragmentation.getAtomsInSameFragment(12, 0), res1)

    def test_FragmentationFindFragmentsNoGroupingWithProtect(self):
        self.fragmentation.setProtectedAtoms()
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        res1 = list(range(1,12))
        res1.extend(list(range(14,31)))
        res2 = [12, 13, 31, 32]
        res2.extend(list(range(35, 43)))
        self.assertEqual(self.fragmentation.getAtomsInSameFragment(11, 0), res1)
        self.assertEqual(self.fragmentation.getAtomsInSameFragment(12, 0), res2)

    def test_FragmentationFragmentChargeAfterFragment(self):
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.determineFragmentCharges()
        self.assertEqual(self.fragmentation.getFragmentCharges(), [ 1, 0, -1, 0, -1, 0, 0, 0, 0, -1])

    def test_FragmentationFragmentChargeAfterProtectAndFragment(self):
        self.fragmentation.setProtectedAtoms()
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.determineFragmentCharges()
        self.assertEqual(self.fragmentation.getFragmentCharges(), [1, -1, 0, -1, 0, 0, 0, 0, -1])

    def test_FragmentationFragmentChargeAfterFragmentAndGroup(self):
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.setFragmentGroupCount(2)
        self.fragmentation.doFragmentGrouping()
        self.fragmentation.determineFragmentCharges()
        self.assertEqual(self.fragmentation.getFragmentCharges(), [1, -1, -1, 0, -1])

    def test_FragmentationFragmentChargeAfterProtectFragmentAndGroup(self):
        self.fragmentation.setProtectedAtoms()
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.setFragmentGroupCount(2)
        self.fragmentation.doFragmentGrouping()
        self.fragmentation.determineFragmentCharges()
        self.assertEqual(self.fragmentation.getFragmentCharges(), [0, -1, 0, 0, -1])

    def test_FragmentationGetOBAtom(self):
        test_atom = self.fragmentation.getOBAtom(1)
        self.assertEqual(type(test_atom), type(openbabel.OBAtom()))

    def test_FragmentationNameFragments(self):
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.nameFragments()
        self.assertEqual(self.fragmentation.getFragmentNames(), ['GLY', 'GLY', 'TYR', 'ASP', 'PRO', 'GLU', 'THR', 'GLY', 'THR', 'TRP'])

    def test_FragmentationNameFragmentsProtect(self):
        self.fragmentation.setProtectedAtoms()
        self.fragmentation.breakBonds()
        self.fragmentation.determineFragments()
        self.fragmentation.nameFragments()
        self.assertEqual(self.fragmentation.getFragmentNames(), ['GLY', 'TYR', 'ASP', 'PRO', 'GLU', 'THR', 'GLY', 'THR', 'TRP'])

    def test_FragmentationNameFragmentsGroupByTwo(self):
        self.fragmentation.breakBonds()
        self.fragmentation.setFragmentGroupCount(2)
        self.fragmentation.determineFragments()
        self.fragmentation.doFragmentGrouping()
        self.fragmentation.nameFragments()
        self.assertEqual(self.fragmentation.getFragmentNames(), ['GLY', 'TYR', 'PRO', 'THR', 'THR'])

    def test_writereadconfiguration_basic(self):
        filename = "temp.cfg"
        handle = open(filename, 'w')
        self.fragmentation.writeConfigurationToFile(handle)
        handle.close()
        otherfrag = Fragmentation(self.molecule)
        otherfrag.readConfigurationFromFile(filename)
        for key in otherfrag.values.keys():
          for key2 in otherfrag.values[key].keys():
            self.assertEqual(self.fragmentation.values[key][key2], otherfrag.values[key][key2])
        self.delete_file(filename)

def suite():
    s = unittest.TestSuite()
    s.addTest(unittest.makeSuite(TestFragmentationModule))
    return s

if __name__ == '__main__':
    unittest.main()
