"""
Copyright (C) 2012-2017 Casper Steinmann
"""
import os
import unittest

import openbabel

from fragit.config import FragItDataPE
from fragit.fragmentation import Fragmentation
from fragit.util import fileToMol
from fragit.writer import Standard
from fragit.xyz import XYZ

class TestStandardWriterModule(unittest.TestCase):

    def setUp(self):
      directories = {'share':''}
      self.molecule = fileToMol("tests/1UAO.pdb")
      self.fragmentation = Fragmentation(self.molecule)
      self.standardwriter = Standard(self.fragmentation, directories)

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
        molecule = fileToMol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()

        xyzwriter = XYZ(fragmentation, {})
        os.chdir("tests")
        xyzwriter.writeFile("temp_water.xyz")

        # we do not test the written coordinates but merely number of atoms written
        nat_written = 0
        for i in range(len(fragmentation.getFragments())):
            filename = "temp_water_fragment_{0:03d}.xyz".format(i+1)
            m2 = fileToMol(filename)
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
        molecule = fileToMol("tests/watercluster4.xyz")
        fragmentation = Fragmentation(molecule, defaults=FragItDataPE)
        fragmentation.beginFragmentation()
        fragmentation.doFragmentation()
        fragmentation.finishFragmentation()

        xyzwriter = XYZ(fragmentation, {})
        os.chdir("tests")
        xyzwriter.writeFile("temp_water.xyz")

        # we do not test the written coordinates but merely number of atoms written
        nat_written = 0
        for i in range(len(fragmentation.getFragments())):
            filename = "temp_water_fragment_{0:03d}.xyz".format(i+1)
            m2 = fileToMol(filename)
            nat_written += m2.NumAtoms()
            self.delete_file(filename)
        self.assertEqual(nat_written, molecule.NumAtoms())

        os.chdir("..")

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestStandardWriterModule))
  s.addTest(unittest.makeSuite(TestXYZWriterModule))
  return s
