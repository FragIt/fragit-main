"""
Copyright (C) 2012-2016 Casper Steinmann
"""
import os
import sys
import unittest
from src.gamessfmo import GamessFMO

from src.util import fileToMol, ReadStringListFromFile
from src.fragmentation import Fragmentation

class TestGamessFMOOutputModule(unittest.TestCase):

    def setUp(self):
      self.molecule = fileToMol("tests/watercluster4.xyz")
      self.fragmentation = Fragmentation(self.molecule)
      self.fixtures = 'tests/gamess-fmo-fixtures'

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
      directories = {'share':''}
      gamessfmo = GamessFMO(self.fragmentation, directories)
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
      directories = {'share':''}
      gamessfmo = GamessFMO(self.fragmentation, directories)
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
      directories = {'share':''}
      gamessfmo = GamessFMO(self.fragmentation, directories)
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
      directories = {'share':''}
      gamessfmo = GamessFMO(self.fragmentation, directories)
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
      directories = {'share':''}
      gamessfmo = GamessFMO(self.fragmentation, directories)
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
