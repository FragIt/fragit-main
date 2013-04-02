"""
**********************************************************************
xyzmfcc.py

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
from util import deepLength,listDiff,intlistToString
from util import getFilenameAndExtension

class XYZMFCC(Standard):
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

  def _build_single_fragment(self,fragments,ifg,caps):
    """
        fragment: atom idx of the current fragment
        pairs   : pairs of breaking points
    """
    output_fragment = fragments[ifg][:]
    output_atoms = [self._fragmentation.mol.GetAtom(id) for id in output_fragment]
    output_types = [atom.GetAtomicNum() for atom in output_atoms]

    lc = None
    rc = None
    if ifg > 0:
      lc = caps[ifg-1]
    if ifg < len(fragments)-1:
      rc = caps[ifg]

    #if lc is not None: print "lc:", lc[1]
    #if rc is not None: print "rc:", rc[1]

    if lc is not None:
      for id,atom_id in enumerate(lc[1]):
        if atom_id in output_fragment: continue
        output_fragment.append(atom_id)
        output_atoms.append(lc[0][id])
        output_types.append(lc[2][id])

    if rc is not None:
      for id,atom_id in enumerate(rc[1]):
        if atom_id in output_fragment: continue
        output_fragment.append(atom_id)
        output_atoms.append(rc[0][id])
        output_types.append(rc[2][id])

    return output_atoms, output_fragment, output_types

  def fragment_xyz(self, atms, ids, types):

    n = len(atms)
    s = "%i\n%s\n" % (n,"")#intlistToString(ids))
    for id, atom in enumerate(atms):
      s += "%s %20.12f %20.12f %20.12f\n" % (self._elements.GetSymbol(types[id]),
                                             atom.GetX(),
                                             atom.GetY(),
                                             atom.GetZ())

    return s

  def writeFile(self, filename):
    """Dumps all fragments to individual
       .xyz files.
    """
    ff,ext = getFilenameAndExtension(filename)
    filename_template = "%s_%s_%03i.%s"
    # first we dump all fragments
    for ifg in range(len(self._fragmentation.getFragments())):
      (atms, ids, types) = self._build_single_fragment(self._fragmentation.getFragments(), ifg, self._fragmentation._caps)
      ss = self.fragment_xyz(atms, ids, types)
      with open(filename_template % (ff,"FRAGMENT",ifg+1,ext), "w") as f:
        f.write(ss)

    for ifg,(atms, ids, types) in enumerate(self._fragmentation._caps):
      ss = self.fragment_xyz(atms, ids, types)
      with open(filename_template % (ff,"CAP",ifg+1,ext), "w") as f:
        f.write(ss)
