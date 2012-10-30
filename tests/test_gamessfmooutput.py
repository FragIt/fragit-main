"""
**********************************************************************
tests/test_GamessFMOOutput.py - test cases for GamessFMO Writer

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
from gamessfmo import GamessFMO

from util import fileToMol, ReadStringListFromFile
from fragmentation import Fragmentation

class TestGamessFMOOutputModule(unittest.TestCase):

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
      gamessfmo = GamessFMO(self.fragmentation)
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
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
      gamessfmo = GamessFMO(self.fragmentation)
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
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
      gamessfmo = GamessFMO(self.fragmentation)
      gamessfmo.setCentralFragmentID(1)
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
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
      gamessfmo = GamessFMO(self.fragmentation)
      gamessfmo.setCentralFragmentID(1)
      gamessfmo.setBoundariesFromString("1.0")
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
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
      gamessfmo = GamessFMO(self.fragmentation)
      gamessfmo.setCentralFragmentID(1)
      gamessfmo.setBoundariesFromString("1.0")
      gamessfmo.setActiveAtomsDistance(1.0)
      gamessfmo.setBufferMaxDistance(1.0)
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestGamessFMOOutputModule))
  return s

if __name__ == '__main__':
    unittest.main()
