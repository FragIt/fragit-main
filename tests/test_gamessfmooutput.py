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
      molecule = fileToMol("tests/watercluster4.xyz")
      fragmentation = Fragmentation(molecule)
      directories = {'share':''}
      gamessfmo = GamessFMO(fragmentation, directories)
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
      molecule = fileToMol("tests/watercluster4.xyz")
      fragmentation = Fragmentation(molecule)
      directories = {'share':''}
      gamessfmo = GamessFMO(fragmentation, directories)
      fragmentation.beginFragmentation()
      fragmentation.doFragmentation()
      fragmentation.finishFragmentation()
      directories = {'share':''}
      gamessfmo = GamessFMO(fragmentation, directories)
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
      molecule = fileToMol("tests/watercluster4.xyz")
      fragmentation = Fragmentation(molecule)
      fragmentation.beginFragmentation()
      fragmentation.doFragmentation()
      fragmentation.finishFragmentation()
      directories = {'share':''}
      gamessfmo = GamessFMO(fragmentation, directories)
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
      molecule = fileToMol("tests/watercluster4.xyz")
      fragmentation = Fragmentation(molecule)
      fragmentation.beginFragmentation()
      fragmentation.doFragmentation()
      fragmentation.finishFragmentation()
      directories = {'share':''}
      gamessfmo = GamessFMO(fragmentation, directories)
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
      molecule = fileToMol("tests/watercluster4.xyz")
      fragmentation = Fragmentation(molecule)
      fragmentation.beginFragmentation()
      fragmentation.doFragmentation()
      fragmentation.finishFragmentation()
      directories = {'share':''}
      gamessfmo = GamessFMO(fragmentation, directories)
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

      # self.delete_file(filename)

    # -----------
    # 5 ala tests
    # -----------
    def test_5ala_1_afo(self):
      filename = "temp.inp"
      otherfile = self.fixtures + "/5ala_1_afo.fixture"
      molecule = fileToMol("tests/5ala.xyz")
      fragmentation = Fragmentation(molecule)
      fragmentation.setFMOAFOFragmentation()
      fragmentation.beginFragmentation()
      fragmentation.doFragmentation()
      fragmentation.finishFragmentation()
      directories = {'share':''}
      gamessfmo = GamessFMO(fragmentation, directories)
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

    def test_5ala_1_hop(self):
      filename = "temp_lol.inp"
      otherfile = self.fixtures + "/5ala_1_hop.fixture"
      molecule = fileToMol("tests/5ala.xyz")
      fragmentation = Fragmentation(molecule)
      fragmentation.beginFragmentation()
      fragmentation.doFragmentation()
      fragmentation.finishFragmentation()
      directories = {'share':'share'}
      gamessfmo = GamessFMO(fragmentation, directories)
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      self.delete_file(filename)

    def test_5ala_2_hop(self):
      filename = "temp_lol.inp"
      otherfile = self.fixtures + "/5ala_2_hop.fixture"
      molecule = fileToMol("tests/5ala.xyz")
      fragmentation = Fragmentation(molecule)
      fragmentation.setQMBasis('3-21G:6-31G*')
      fragmentation.beginFragmentation()
      fragmentation.doFragmentation()
      fragmentation.finishFragmentation()
      directories = {'share':'share'}
      gamessfmo = GamessFMO(fragmentation, directories)
      gamessfmo.setCentralFragmentID(1)
      gamessfmo.setBoundariesFromString("1.0")
      gamessfmo.setup()
      gamessfmo.writeFile(filename)
      generated = ReadStringListFromFile(filename)
      fixture = ReadStringListFromFile(otherfile)

      self.assertEqual(len(generated), len(fixture))
      for i in range(len(fixture)):
        self.assertEqual(generated[i], fixture[i])

      # self.delete_file(filename)

    #def test_5ala_3_hop(self):
    #  """ Regression test to make sure if only one basis
    #      set is specified we get the correct answer for HOP
    #  """
    #  filename = "temp_lol.inp"
    #  otherfile = self.fixtures + "/5ala_2_hop.fixture"
    #  molecule = fileToMol("tests/5ala.xyz")
    #  fragmentation = Fragmentation(molecule)
    #  fragmentation.setQMBasis('3-21G:6-31G(d)')
    #  fragmentation.beginFragmentation()
    #  fragmentation.doFragmentation()
    #  fragmentation.finishFragmentation()
    #  directories = {'share':'share'}
    #  gamessfmo = GamessFMO(fragmentation, directories)
    #  gamessfmo.setCentralFragmentID(1)
    #  gamessfmo.setBoundariesFromString("1.0")
    #  gamessfmo.setup()
    #  gamessfmo.writeFile(filename)
    #  generated = ReadStringListFromFile(filename)
    #  fixture = ReadStringListFromFile(otherfile)

    #  self.assertEqual(len(generated), len(fixture))
    #  for i in range(len(fixture)):
    #    self.assertEqual(generated[i], fixture[i])

    #  # self.delete_file(filename)

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestGamessFMOOutputModule))
  return s

if __name__ == '__main__':
    unittest.main()
