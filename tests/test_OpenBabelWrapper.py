"""
**********************************************************************
tests/test_OpenBabelWrapper.py - test functionality of openbabel

Copyright (C) 2011 Casper Steinmann

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
import os
import unittest
import openbabel

import OpenBabelWrapper as OBW

def write_file(filename):
    f = open(filename,'w')
    f.write("9\n")
    f.write("\n")
    f.write("C         -3.97094        3.73817       -0.07959\n")
    f.write("C         -2.66836        3.55462        0.16705\n")
    f.write("C         -4.75529        2.94754       -1.07449\n")
    f.write("H         -5.16763        3.61369       -1.83856\n")
    f.write("H         -5.58647        2.43984       -0.57539\n")
    f.write("H         -4.14386        2.19010       -1.57513\n")
    f.write("H         -4.51083        4.50680        0.46887\n")
    f.write("H         -2.08836        2.80010       -0.35577\n")
    f.write("H         -2.15162        4.16209        0.90394")
    f.close()

def remove_file(filename):
    os.remove(filename)

class TestOpenBabelWrapperModule(unittest.TestCase):

    def setUp(self):
        self.filename_xyz = "test.xyz"
        write_file(self.filename_xyz)

        # for testing OpenBabel functionality directly
        self.molecule = openbabel.OBMol()
        self.conversion = openbabel.OBConversion()
        self.conversion.SetInFormat("xyz")
        self.conversion.ReadFile(self.molecule, self.filename_xyz)

        # for testing OpenBabelWrapper functionality
        self.obwmolecule = OBW.Molecule(self.filename_xyz)

    def tearDown(self):
        remove_file(self.filename_xyz)

    def test_OpenBabelAPIMoleculeFromFileSimple(self):
        self.assertEqual(self.molecule.NumAtoms(), 9)

    def test_OBWMoleculeAtomCount(self):
        self.assertEqual(self.obwmolecule.getAtomCount(), 9)

    def test_OBWMoleculeIsOK(self):
        self.assertTrue(self.obwmolecule.isOK())

    def test_OBWMoleculeTotalCharge(self):
        self.assertEqual(self.obwmolecule.getTotalCharge(), 0)

    def test_OBWMoleculeElementSymbol(self):
        elements = ["H","C","N","O","S","Cl"]
        element_id = [1,6,7,8,16,17]
        for i in element_id:
            idx = element_id.index(i)
            self.assertEqual(self.obwmolecule.getElementSymbol(i), elements[idx])

    def test_OBWMoleculeMatchPatternSimple(self):
        self.assertEqual(self.obwmolecule.MatchPattern("C"), [(1,),(2,),(3,)])
        self.assertEqual(self.obwmolecule.MatchPattern("CC=C"), [(3,1,2)])
        self.assertEqual(self.obwmolecule.MatchPattern("C=C"), [(1,2)])

    def test_OBWMoleculeMatchPatternAdvanced(self):
	self.assertEqual(self.obwmolecule.MatchPattern("[$(C=C)][$(CC)]"), [(1,3)])

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestOpenBabelWrapperModule))
  return s

if __name__ == '__main__':
    unittest.main()
