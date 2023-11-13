"""
Copyright (C) 2012-2017 Casper Steinmann
"""
import os
import unittest

from fragit.config import FragItDataPE
from fragit.fragmentation import Fragmentation
from fragit.util import file_to_mol
from fragit.writer import Standard
from fragit.xyz import XYZ

class TestStandardWriterModule(unittest.TestCase):

    def setUp(self):
      directories = {'share':''}
      self.molecule = file_to_mol("tests/1UAO.pdb")
      self.fragmentation = Fragmentation(self.molecule)
      self.standardwriter = Standard(self.fragmentation, directories)

    def test_defaults(self):
      self.assertEqual(self.standardwriter._nlayers, 1)

    def test_notimplementederror(self):
      self.assertRaises(NotImplementedError, self.standardwriter.write_file, "")
      self.assertRaises(NotImplementedError, self.standardwriter.setup)

    def test_setcentralfragmentid(self):
      self.assertRaises(ValueError, self.standardwriter.set_central_fragment_id, -1)

    def test_setbuffermaxdistance(self):
      self.assertRaises(ValueError, self.standardwriter.set_buffer_max_distance, -0.1)

    def test_setboundariesfromstring_two_layers(self):
      layerstring = "1.0"
      nlayers = 2
      self.standardwriter.set_boundaries_from_string(layerstring)
      self.assertEqual(self.standardwriter._nlayers, nlayers)
      self.assertEqual(self.standardwriter._boundaries, [1.0])

    def test_setboundariesfromstring_three_layers(self):
      layerstring = "1.0,2.0"
      nlayers = 3
      self.standardwriter.set_boundaries_from_string(layerstring)
      self.assertEqual(self.standardwriter._nlayers, nlayers)
      self.assertEqual(self.standardwriter._boundaries, [1.0,2.0])

class TestXYZWriterModule(unittest.TestCase):
    def setUp(self):
        pass

    def delete_file(self,filename):
        try:
            f = open(filename)
        except IOError:
            return
        finally:
            f.close()
            os.remove(filename)

    def test_write_xyz(self):
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        xyzwriter = XYZ(fragmentation, {})
        os.chdir("tests")
        xyzwriter.write_file("temp_water.xyz")

        # we do not test the written coordinates but merely number of atoms written
        nat_written = 0
        for i in range(len(fragmentation.get_fragments())):
            filename = "temp_water_fragment_{0:03d}.xyz".format(i+1)
            m2 = file_to_mol(filename)
            nat_written += m2.NumAtoms()
            self.delete_file(filename)
        self.assertEqual(nat_written, molecule.NumAtoms())

        os.chdir("..")


class TestXYZMFCCWriterModule(unittest.TestCase):
    def setUp(self):
        pass

    def delete_file(self,filename):
        try:
            f = open(filename)
        except IOError:
            return
        finally:
            f.close()
            os.remove(filename)

    def test_write_xyz(self):
        molecule = file_to_mol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.begin_fragmentation()
        fragmentation.do_fragmentation()
        fragmentation.finish_fragmentation()

        xyzwriter = XYZ(fragmentation, {})
        os.chdir("tests")
        xyzwriter.write_file("temp_water.xyz")

        # we do not test the written coordinates but merely number of atoms written
        nat_written = 0
        for i in range(len(fragmentation.get_fragments())):
            filename = "temp_water_fragment_{0:03d}.xyz".format(i+1)
            m2 = file_to_mol(filename)
            nat_written += m2.NumAtoms()
            self.delete_file(filename)
        self.assertEqual(nat_written, molecule.NumAtoms())

        os.chdir("..")

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestStandardWriterModule))
  s.addTest(unittest.makeSuite(TestXYZWriterModule))
  return s
