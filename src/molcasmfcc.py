"""
**********************************************************************
molcasmfcc.py

Copyright (C) 2013 Casper Steinmann

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
from numpy import sqrt, dot, where, array
from writer import Standard
from util import WriteStringToFile

from util import file_extension,is_list,listTo2D,join2D,is_int
from util import listToRanges,listOfRangesToString,Uniqify,ravel2D
from util import deepLength,listDiff

class MolcasMFCC(Standard):
  def __init__(self, fragmentation):
    Standard.__init__(self,fragmentation)

  def setup(self):
    #self._setupLayeredInformation()
    #self._setupActiveFragmentsInformation()
    #self._validateMultiLayerInformation()
    if self._do_pymol: self._dump_pymol()
    if self._do_jmol: self._dump_jmol()

  def _dump_pymol(self):
    from pymol import PymolTemplate
    pt = PymolTemplate(self._input_filename, self._output_filename)
    self._setTemplateData(pt)
    self._writeTemplateFile(pt)

  def _dump_jmol(self):
    from jmol import JmolTemplate
    pt = JmolTemplate(self._input_filename, self._output_filename)
    self._setTemplateData(pt)
    self._writeTemplateFile(pt)

  def _setTemplateData(self, template):
    template.setFragmentsData(self._fragmentation.getFragments())
    #template.setBufferData(self._fragment_layers)
    #template.setActiveData(self._active_atoms)
    #template.setBackboneData(self._fragmentation.getBackboneAtoms())
    #template.setPairData(self._fragmentation.getExplicitlyBreakAtomPairs())

  def _writeTemplateFile(self, template):
    template.override()
    template.write()

  def _build_single_fragment(self,fragment,pairs):
    """
        fragment: atom idx of the current fragment
        pairs   : pairs of breaking points
    """

    output_fragment = fragment[:]

    # locate pairs which we would also like to include
    for pair in pairs:
      if pair[0] in output_fragment or pair[1] in output_fragment:
        output_fragment.extend(pair)
    output_fragment = Uniqify(output_fragment)
    caps = listDiff(output_fragment, fragment)
    print fragment
    print caps

  def writeFile(self, filename):
    """Dumps all fragments to individual
       input files.

       initially, we dump xyz-files. this will change later
    """
    filename_template = "%s%i.xyz"
    # first we dump all fragments
    self._build_single_fragment(self._fragmentation.getFragments()[2], self._fragmentation.getExplicitlyBreakAtomPairs())
