"""
**********************************************************************
tests/test_OutputFormats.py - test cases for OutputFormats

Copyright (C) 2011 Casper Steinmann

This file is part of the FragIt project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
"""
import os
import unittest
from gamess import Gamess
import outputformats

class TestOutputFormatsModule(unittest.TestCase):

    def setUp(self):
        pass

    def test_getwriterandextension(self):
        self.assertRaises(ValueError, outputformats.get_writer_and_extension, "bogus")
        self.assertEqual(outputformats.get_writer_and_extension("GAMESS"), (Gamess,".inp"))


    def test_supported_output_formats(self):
        self.assertEqual( outputformats.supported_output_formats(), {'GAMESS':Gamess} )

    def test_supported_output_fileexts(self):
        self.assertEqual( outputformats.supported_output_fileexts(), {'GAMESS':'.inp'} )

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestOutputFormatsModule))
  return s

if __name__ == '__main__':
    unittest.main()
