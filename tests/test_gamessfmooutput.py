"""
Copyright (C) 2012-2017 Casper Steinmann
"""
import os
import sys
import unittest

from fragit.config import FragItDataFMO
from fragit.gamessfmo import GamessFMO
from fragit.fragmentation import Fragmentation
from fragit.util import file_to_mol, read_string_list_from_file

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
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    def test_water_2(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/water_2.fixture"
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    def test_water_3(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/water_3.fixture"
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.set_central_fragment_id(1)
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    def test_water_4(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/water_4.fixture"
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.set_central_fragment_id(1)
        gamessfmo.set_boundaries_from_string("1.0")
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    def test_water_5(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/water_5.fixture"
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.set_central_fragment_id(1)
        gamessfmo.set_boundaries_from_string("1.0")
        gamessfmo.set_active_atoms_distance(1.0)
        gamessfmo.set_buffer_max_distance(1.0)
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

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
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.set_fmoafo_fragmentation()
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    def test_5ala_2_afo(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_2_afo.fixture"
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.set_qm_basis('3-21G:6-31G(d)')
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':'share'}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.set_central_fragment_id(1)
        gamessfmo.set_boundaries_from_string("1.0")
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)


    ## HOP TESTS ##
    def test_5ala_1_hop(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_1_hop.fixture"
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.set_fmohop_fragmentation()
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':'share'}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)


    def test_5ala_2_hop(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_2_hop.fixture"
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.set_qm_basis('3-21G:6-31G*')
        fragmentation.set_fmohop_fragmentation()
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':'share'}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.set_central_fragment_id(1)
        gamessfmo.set_boundaries_from_string("1.0")
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    def test_5ala_3_hop(self):
        """ Correct input for HOP with one basis """
        filename = "temp.inp"
        otherfile = self.fixtures + "/5ala_3_hop.fixture"
        molecule = file_to_mol("tests/5ala.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.set_fmohop_fragmentation()
        fragmentation.set_qm_basis('3-21G')
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        directories = {'share':'share'}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.set_central_fragment_id(1)
        gamessfmo.set_boundaries_from_string("1.0")
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

        self.assertEqual(len(generated), len(fixture))
        for i in range(len(fixture)):
            self.assertEqual(generated[i], fixture[i])

        self.delete_file(filename)

    ## Test FMO EFP writer ##
    # basic test, all waters replaced by H2ORHF waters.
    def test_2form8wat_1(self):
        filename = "temp.inp"
        otherfile = self.fixtures + "/2form8wat_1.fixture"
        molecule = file_to_mol("tests/2form8wat.pdb")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        fragmentation.set_fmoefp_waters_from_layer(1)
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

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
        molecule = file_to_mol("tests/2form8wat.pdb")
        fragmentation = Fragmentation(molecule, defaults=FragItDataFMO)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()
        fragmentation.set_fmoefp_waters_from_layer(1)
        directories = {'share':''}
        gamessfmo = GamessFMO(fragmentation, directories)
        gamessfmo.set_central_fragment_id(1)
        gamessfmo.set_boundaries_from_string("3.0")
        gamessfmo.setup()
        gamessfmo.write_file(filename)
        generated = read_string_list_from_file(filename)
        fixture = read_string_list_from_file(otherfile)

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
