"""
**********************************************************************
tests/test_StandardWriter.py - test cases for Gamess Writer

Copyright (C) 2012 Casper Steinmann

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
from writer import Standard
from util import fileToMol
from fragmentation import Fragmentation

class TestStandardWriterModule(unittest.TestCase):

    def setUp(self):
      self.molecule = fileToMol("1UAO.pdb")
      self.fragmentation = Fragmentation(self.molecule)
      self.standardwriter = Standard(Fragmentation)

    def test_defaults(self):
      self.assertEqual(self.standardwriter._nlayers, 1)

    def test_notimplementederror(self):
      self.assertRaises(NotImplementedError, self.standardwriter.writeFile)
      self.assertRaises(NotImplementedError, self.standardwriter.setup)

    def test_setcentralfragmentid(self):
      self.assertRaises(ValueError, self.standardwriter.setCentralFragmentID, -1)

    def test_setbuffermaxdistance(self):
      self.assertRaises(ValueError, self.standardwriter.setBufferMaxDistance, -0.1)

    def test_setboundariesfromstring_two_layers(self):
      layerstring = "1.0"
      nlayers = 2
      self.standardwriter.setBoundariesFromString(layerstring)
      self.assertEqual(self.standardwriter._nlayers, nlayers)
      self.assertEqual(self.standardwriter._boundaries, [1.0])

    def test_setboundariesfromstring_three_layers(self):
      layerstring = "1.0,2.0"
      nlayers = 3
      self.standardwriter.setBoundariesFromString(layerstring)
      self.assertEqual(self.standardwriter._nlayers, nlayers)
      self.assertEqual(self.standardwriter._boundaries, [1.0,2.0])

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestStandardWriterModule))
  return s

if __name__ == '__main__':
    unittest.main()
