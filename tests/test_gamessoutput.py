"""
**********************************************************************
tests/test_GamessOutput.py - test cases for Gamess Writer

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
import sys
sys.path.append('../src')
import unittest
from gamess import Gamess

from util import fileToMol, ReadStringListFromFile
from fragmentation import Fragmentation

class TestGamessOutputModule(unittest.TestCase):

    def setUp(self):
      self.molecule = fileToMol("watercluster4.xyz")
      self.fragmentation = Fragmentation(self.molecule)
      self.fixtures = 'gamess-fmo-fixtures'

    def delete_file(self,filename):
        try:
                f = open(filename)
        except IOError:
                return
        finally:
                f.close()
                os.remove(filename)

    def test_water_1(self):
      filename = "temp.inp"
      otherfile = self.fixtures + "/water_1.fixture"
      gamess = Gamess(self.fragmentation)
      gamess.setup()
      gamess.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

    def test_water_2(self):
      filename = "temp.inp"
      otherfile = self.fixtures + "/water_2.fixture"
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      gamess = Gamess(self.fragmentation)
      gamess.setup()
      gamess.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

    def test_water_3(self):
      filename = "temp.inp"
      otherfile = self.fixtures + "/water_3.fixture"
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      gamess = Gamess(self.fragmentation)
      gamess.setCentralFragmentID(1)
      gamess.setup()
      gamess.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

    def test_water_4(self):
      filename = "temp.inp"
      otherfile = self.fixtures + "/water_4.fixture"
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      gamess = Gamess(self.fragmentation)
      gamess.setCentralFragmentID(1)
      gamess.setBoundariesFromString("1.0")
      gamess.setup()
      gamess.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

    def test_water_5(self):
      filename = "temp.inp"
      otherfile = self.fixtures + "/water_5.fixture"
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      gamess = Gamess(self.fragmentation)
      gamess.setCentralFragmentID(1)
      gamess.setBoundariesFromString("1.0")
      gamess.setActiveAtomsDistance(1.0)
      gamess.setBufferMaxDistance(1.0)
      gamess.setup()
      gamess.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestGamessOutputModule))
  return s

if __name__ == '__main__':
    unittest.main()
