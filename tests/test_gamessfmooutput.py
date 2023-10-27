"""
Copyright (C) 2012-2017 Casper Steinmann
"""
import os
import sys
import unittest

from fragit.config import FragItDataFMO
from fragit.gamessfmo import GamessFMO
from fragit.fragmentation import Fragmentation
from fragit.util import fileToMol, ReadStringListFromFile

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
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
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
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
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
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
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
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
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
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
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

        self.delete_file(filename)

    # -----------
    # 5 ala tests
    # -----------
    def test_5ala_1_afo(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_1_afo.fixture"
        molecule = fileToMol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
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

    def test_5ala_2_afo(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_2_afo.fixture"
        molecule = fileToMol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.setQMBasis('3-21G:6-31G(d)')
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

        self.delete_file(filename)


    ## HOP TESTS ##
    def test_5ala_1_hop(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_1_hop.fixture"
        molecule = fileToMol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.setFMOHOPFragmentation()
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
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_2_hop.fixture"
        molecule = fileToMol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.setQMBasis('3-21G:6-31G*')
        fragmentation.setFMOHOPFragmentation()
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

        self.delete_file(filename)

    def test_5ala_3_hop(self):
        """ Correct input for HOP with one basis """
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_3_hop.fixture"
        molecule = fileToMol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.setFMOHOPFragmentation()
        fragmentation.setQMBasis('3-21G')
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

        self.delete_file(filename)

    ## Test FMO EFP writer ##
    # basic test, all waters replaced by H2ORHF waters.
    def test_2form8wat_1(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/2form8wat_1.fixture"
        molecule = fileToMol("tests/2form8wat.pdb")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()
        fragmentation.setFMOEFPWatersFromLayer(1)
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.setup()
        gamessfmo.writeFile(filename)
        generated = ReadStringListFromFile(filename)
        fixture = ReadStringListFromFile(otherfile)

        self.assertEqual(len(generated), len(fixture))

        ignoring = False
        for i in range(len(fixture)):
            if "EFRAG" in generated[i] or "EFRAG" in fixture[i]:
                ignoring = True

            if ignoring:
                if "END" in generated[i] or "END" in fixture[i]:
                    ignoring = False

            if not ignoring:
                self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    def test_2form8wat_2(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/2form8wat_2.fixture"
        molecule = fileToMol("tests/2form8wat.pdb")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()
        fragmentation.setFMOEFPWatersFromLayer(1)
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.setCentralFragmentID(1)
        gamessfmo.setBoundariesFromString("3.0")
        gamessfmo.setup()
        gamessfmo.writeFile(filename)
        generated = ReadStringListFromFile(filename)
        fixture = ReadStringListFromFile(otherfile)

        self.assertEqual(len(generated), len(fixture))

        ignoring = False
        for i in range(len(fixture)):
            if "EFRAG" in generated[i] or "EFRAG" in fixture[i]:
                ignoring = True

            if ignoring:
                if "END" in generated[i] or "END" in fixture[i]:
                    ignoring = False

            if not ignoring:
                self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestGamessFMOOutputModule))
  return s

if __name__ == '__main__':
    unittest.main()
