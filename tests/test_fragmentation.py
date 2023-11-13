"""
Copyright (C) 2011-2017 Casper Steinmann
"""
import os
import unittest
from openbabel import openbabel

from fragit.fragmentation import Fragmentation
from fragit.util import file_to_mol
from fragit.config import FragItDataFMO

class TestFragmentationModule(unittest.TestCase):

    def setUp(self):
        self.filename_pdb = "tests/1UAO.pdb"

        # for testing, use OpenBabel functionality directly
        self.molecule = file_to_mol(self.filename_pdb)
        self.fragmentation = Fragmentation(self.molecule, defaults=FragItDataFMO)

    def tearDown(self):
        pass

    def delete_file(self, filename):
        try:
            f = open(filename)
        except IOError:
            return
        finally:
            f.close()
            os.remove(filename)

    def test_FragmentationDefaultParameters(self):
        frg = self.fragmentation
        self.assertEqual(frg.mol is not None, True)

    def test_FragmentationApplySmartProtectPatterns(self):
        self.fragmentation.apply_smart_protect_patterns()
        self.assertEqual(self.fragmentation.get_explicitly_protected_atoms(), [1, 2, 3, 4, 10])

    def test_FragmentationDetermineFormalCharges(self):
        #self.fragmentation.determineFormalCharges()
        self.assertAlmostEqual(sum(self.fragmentation.formalCharges), -2)

    def test_FragmentationGetProtectedAtoms(self):
        self.assertEqual(self.fragmentation.get_explicitly_protected_atoms(), [])

    def test_FragmentationAddProtectedAtomsAfterProtect(self):
        self.fragmentation.set_protected_atoms()
        self.fragmentation.add_explicitly_protected_atoms([44, 55, 67])
        self.assertEqual(self.fragmentation.get_explicitly_protected_atoms(), [1, 2, 3, 4, 10, 44, 55, 67])

    def test_FragmentationAddBrokenBond(self):
        pass

    def test_FragmentationIsBondProtected(self):
        bond_pair = (2, 3)
        self.assertEqual(self.fragmentation.is_bond_protected(bond_pair), False)
        self.fragmentation.add_explicitly_protected_atoms([2])
        self.assertEqual(self.fragmentation.is_bond_protected(bond_pair), True)

    def test_FragmentationRealBondBreakerNoProtect(self):
        bond_atoms = (2, 3)
        self.fragmentation.add_fragmentation_atom_pair(bond_atoms)
        self.assertEqual(self.fragmentation.get_explicitly_break_atom_pairs(), [(2, 3)])

    def test_FragmentationIsValidExplicitBond(self):
        self.assertRaises(ValueError, self.fragmentation.is_valid_explicit_bond, (1, 1))
        self.assertRaises(ValueError, self.fragmentation.is_valid_explicit_bond, (2, 4))

    def test_FragmentationBreakBondsWithNoProtect(self):
        self.fragmentation.break_bonds()
        self.assertEqual(len(self.fragmentation.get_explicitly_break_atom_pairs()), 9)

    def test_FragmentationBreakBondsExplcitWithNoProtect(self):
        self.fragmentation.add_explicitly_break_atom_pairs([(111, 112)])
        self.fragmentation.break_bonds()
        self.assertEqual(len(self.fragmentation.get_explicitly_break_atom_pairs()), 10)

    def test_FragmentationBreakBondsWithProtect(self):
        self.fragmentation.set_protected_atoms()
        self.fragmentation.break_bonds()
        self.assertEqual(len(self.fragmentation.get_explicitly_break_atom_pairs()), 8)

    def test_FragmentationBreakBondsExplcitWithProtect(self):
        self.fragmentation.set_protected_atoms()
        self.fragmentation.add_explicitly_break_atom_pairs([(111, 112)])
        self.fragmentation.break_bonds()
        self.assertEqual(len(self.fragmentation.get_explicitly_break_atom_pairs()), 9)

    def test_FragmentationDetermineFragmentsNoBreaking(self):
        self.assertRaises(ValueError, self.fragmentation.build_fragments)

    def test_FragmentationDetermineFragmentsWithBreaking(self):
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        for fragment, expected_size in zip(self.fragmentation.get_fragments(), [7, 21, 12, 14, 15, 14, 7, 14, 24, 10]):
            self.assertEqual(len(fragment), expected_size)

    def test_FragmentationDetermineFragmentsWithBreakingAndGrouping(self):
        self.fragmentation.set_fragment_group_count(2)
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.do_fragment_grouping()
        self.assertEqual(len(self.fragmentation.get_fragments()), 5)
        for fragment, expected_size in zip(self.fragmentation.get_fragments(), [28, 26, 29, 21, 34]):
            self.assertEqual(len(fragment), expected_size)

    def test_FragmentationDetermineFragmentsWithBreakingAndGroupingTriple(self):
        self.fragmentation.set_fragment_group_count(3)
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.do_fragment_grouping()
        self.assertEqual(len(self.fragmentation.get_fragments()), 4)
        for fragment, expected_size in zip(self.fragmentation.get_fragments(), [40, 43, 45, 10]):
            self.assertEqual(len(fragment), expected_size)

    def test_FragmentationFindFragmentsCastErrors(self):
        self.assertRaises(ValueError, self.fragmentation.get_atoms_in_same_fragment, 1, 1)
        self.assertRaises(ValueError, self.fragmentation.get_atoms_in_same_fragment, "")

    def test_FragmentationFindFragmentsNoGrouping(self):
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        res1 = [12, 13, 31, 32]
        res1.extend(range(35, 43))
        self.assertEqual(self.fragmentation.get_atoms_in_same_fragment(11, 0), [3, 4, 10, 11, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30])
        self.assertEqual(self.fragmentation.get_atoms_in_same_fragment(12, 0), res1)

    def test_FragmentationFindFragmentsNoGroupingWithProtect(self):
        self.fragmentation.set_protected_atoms()
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        res1 = list(range(1,12))
        res1.extend(list(range(14,31)))
        res2 = [12, 13, 31, 32]
        res2.extend(list(range(35, 43)))
        self.assertEqual(self.fragmentation.get_atoms_in_same_fragment(11, 0), res1)
        self.assertEqual(self.fragmentation.get_atoms_in_same_fragment(12, 0), res2)

    def test_FragmentationFragmentChargeAfterFragment(self):
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.determine_fragment_charges()
        self.assertEqual(self.fragmentation.get_fragment_charges(), [1, 0, -1, 0, -1, 0, 0, 0, 0, -1])

    def test_FragmentationFragmentChargeAfterProtectAndFragment(self):
        self.fragmentation.set_protected_atoms()
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.determine_fragment_charges()
        self.assertEqual(self.fragmentation.get_fragment_charges(), [1, -1, 0, -1, 0, 0, 0, 0, -1])

    def test_FragmentationFragmentChargeAfterFragmentAndGroup(self):
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.set_fragment_group_count(2)
        self.fragmentation.do_fragment_grouping()
        self.fragmentation.determine_fragment_charges()
        self.assertEqual(self.fragmentation.get_fragment_charges(), [1, -1, -1, 0, -1])

    def test_FragmentationFragmentChargeAfterProtectFragmentAndGroup(self):
        self.fragmentation.set_protected_atoms()
        #self.fragmentation.determineFormalCharges()
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.set_fragment_group_count(2)
        self.fragmentation.do_fragment_grouping()
        self.fragmentation.determine_fragment_charges()
        self.assertEqual(self.fragmentation.get_fragment_charges(), [0, -1, 0, 0, -1])

    def test_FragmentationGetOBAtom(self):
        test_atom = self.fragmentation.get_ob_atom(1)
        self.assertEqual(type(test_atom), type(openbabel.OBAtom()))

    def test_FragmentationNameFragments(self):
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.name_fragments()
        self.assertEqual(self.fragmentation.get_fragment_names(), ['GLY', 'GLY', 'TYR', 'ASP', 'PRO', 'GLU', 'THR', 'GLY', 'THR', 'TRP'])

    def test_FragmentationNameFragmentsProtect(self):
        self.fragmentation.set_protected_atoms()
        self.fragmentation.break_bonds()
        self.fragmentation.build_fragments()
        self.fragmentation.name_fragments()
        self.assertEqual(self.fragmentation.get_fragment_names(), ['GLY', 'TYR', 'ASP', 'PRO', 'GLU', 'THR', 'GLY', 'THR', 'TRP'])

    def test_FragmentationNameFragmentsGroupByTwo(self):
        self.fragmentation.break_bonds()
        self.fragmentation.set_fragment_group_count(2)
        self.fragmentation.build_fragments()
        self.fragmentation.do_fragment_grouping()
        self.fragmentation.name_fragments()
        self.assertEqual(self.fragmentation.get_fragment_names(), ['GLY', 'TYR', 'PRO', 'THR', 'THR'])

    def test_writereadconfiguration_basic(self):
        filename = "temp.cfg"
        handle = open(filename, 'w')
        self.fragmentation.write_configuration_to_file(handle)
        handle.close()
        otherfrag = Fragmentation(self.molecule)
        otherfrag.read_configuration_from_file(filename)
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
