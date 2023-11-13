"""
Copyright (C) 2011-2017 Casper Steinmann
"""
import os
import unittest
from fragit.gamessfmo import GamessFMO
from fragit.xyzmfcc import XYZMFCC
from fragit.xyz import XYZ
from fragit.outputformats import *

class TestOutputFormatsModule(unittest.TestCase):

    def setUp(self):
        pass

    def test_getwriterandextension(self):
        self.assertRaises(ValueError, get_writer_and_extension, "bogus")
        self.assertEqual(get_writer_and_extension("GAMESS-FMO"), (GamessFMO,".inp"))


    def test_supported_output_formats(self):
        self.assertEqual( supported_output_formats(), {'GAMESS-FMO':GamessFMO, 'XYZ-MFCC': XYZMFCC, 'XYZ': XYZ} )

    def test_supported_output_fileexts(self):
        self.assertEqual( supported_output_fileexts(), {'GAMESS-FMO':'.inp', 'XYZ-MFCC': '.xyz', 'XYZ': '.xyz'} )

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestOutputFormatsModule))
  return s

if __name__ == '__main__':
    unittest.main()
