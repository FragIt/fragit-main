"""
**********************************************************************
writer.py

Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2012 Casper Steinmann

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

import openbabel

from util import floatlistFromString, intlistFromString

class Standard(object):
	def __init__(self, fragmentation):
		self._fragmentation = fragmentation
		self._elements = openbabel.OBElementTable()
		self._nlayers = 1
		self._boundaries = []
		self._title = ""
		self._active_fragments = []
		self._central_fragment = 0
		self._active_atoms_distance = 0.0
		self._freeze_backbone = False
		self._input_filename = None
		self._do_jmol = False
		self._do_pymol = False

	def writeFile(self):
		raise NotImplementedError

	def setup(self):
		raise NotImplementedError

	def setBoundariesFromString(self, value):
		self._boundaries = floatlistFromString(value)
		self._nlayers = len(self._boundaries)+1

	def setCentralFragmentID(self, value):
		if value < 0: raise ValueError
		self._central_fragment = value

	def setActiveFragments(self, value):
		self._active_fragments = intlistFromString(value)

	def setActiveAtomsDistance(self, value):
		self._active_atoms_distance = value

	def setBufferMaxDistance(self, value):
		if value < 0.0: raise ValueError
		self._buffer_maximum_distance = value

	def setFreezeBackbone(self):
		self._freeze_backbone = True

	def setJmolOutput(self, infile, outfile):
		self._input_filename = infile
		self._output_filename = outfile
		self._do_jmol = True

	def setPymolOutput(self, infile, outfile):
		self._input_filename = infile
		self._output_filename = outfile
		self._do_pymol = True
