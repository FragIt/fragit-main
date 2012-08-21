"""
**********************************************************************
tests/test_GamessOutput.py - test cases for Gamess Writer

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
import unittest
from gamess import Gamess

from util import fileToMol
from fragmentation import Fragmentation

class TestGamessOutputModule(unittest.TestCase):

    def setUp(self):
      self.molecule = fileToMol("1UAO.pdb")
      self.fragmentation = Fragmentation(self.molecule)
      self.gamess = Gamess(Fragmentation)

    #
    # STATIC string tests. 
    #

    def test_SYSTEMgroup(self):
      STRING=" $SYSTEM MWORDS=125 $END"
      self.assertEqual(self.gamess.SYSTEMgroup(), STRING)

    def test_GDDIgroup(self):
      STRING=" $GDDI NGROUP=1 $END"
      self.assertEqual(self.gamess.GDDIgroup(), STRING)

    def test_SCFgroup(self):
      STRING=" $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END"
      self.assertEqual(self.gamess.SCFgroup(), STRING)

    def test_BASISgroup(self):
      STRING=" $BASIS GBASIS=N21 NGAUSS=3 $END"
      self.assertEqual(self.gamess.BASISgroup(), STRING)

    def test_CONTRLgroup(self):
      STRING = " $CONTRL NPRINT=-5 ISPHER=1 LOCAL=BOYS\n         RUNTYP=ENERGY\n $END"
      self.assertEqual(self.gamess.CONTRLgroup(), STRING)

    #
    # integration tests
    #
    def test_DATAgroup(self):
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      self.gamess = Gamess(self.fragmentation)
      self.gamess.setup()
      STRING=" $DATA\n\nc1\nH-1 1\nC-1 6\nN-1 7\nO-1 8\n $END"
      self.assertEqual(self.gamess.DATAgroup(), STRING)

    #
    # heavy integration tests. Possible to skip these
    #

    @unittest.skipUnless('-full' in sys.argv,'Skipping time-consuming integration tests')
    def test_IntegrationFMOPRPgroupSingleLayer(self):
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      self.gamess = Gamess(self.fragmentation)
      self.gamess.setup()
      STRING=" $FMOPRP NPRINT=9 NGUESS=6 $END"
      self.assertEqual(self.gamess.FMOPRPgroup(), STRING)

    @unittest.skipUnless('-full' in sys.argv,'Skipping time-consuming integration tests')
    def test_IntegrationFMOPRPgroupMultipleLayer(self):
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      self.gamess = Gamess(self.fragmentation)
      self.gamess.setCentralFragmentID(1)
      self.gamess.setBoundariesFromString("1.0")
      self.gamess.setup()
      STRING=" $FMOPRP NPRINT=9 NGUESS=6 $END"
      self.assertEqual(self.gamess.FMOPRPgroup(), STRING)

    @unittest.skipUnless('-full' in sys.argv,'Skipping time-consuming integration tests')
    def test_IntegrationFMOPRPgroupMultipleLayerFDOptimization(self):
      self.fragmentation.beginFragmentation()
      self.fragmentation.doFragmentation()
      self.fragmentation.finishFragmentation()
      self.gamess = Gamess(self.fragmentation)
      self.gamess.setCentralFragmentID(1)
      self.gamess.setBoundariesFromString("1.0")
      self.gamess.setActiveAtomsDistance("2.0")
      self.gamess.setBufferMaxDistance("2.0")
      self.gamess.setup()
      STRING=" $FMOPRP NPRINT=9 NGUESS=134 $END"
      self.assertEqual(self.gamess.FMOPRPgroup(), STRING)

def suite():
  s = unittest.TestSuite()
  s.addTest(unittest.makeSuite(TestGamessOutputModule))
  return s

if __name__ == '__main__':
    unittest.main()
