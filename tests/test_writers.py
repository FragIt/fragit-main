"""
Copyright (C) 2012-2016 Casper Steinmann
"""
import os
import unittest
from src.writer import Standard
from src.util import fileToMol
from src.fragmentation import Fragmentation

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

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestStandardWriterModule))
  s.addTest(unittest.makeSuite(TestXYZWriterModule))
  return s
