"""
**********************************************************************
openbabelwrapper.py

Copyright (C) 2011-2016 Casper Steinmann

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

try:
    import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")
from util import file_extension

class Molecule(object):

	def __init__(self, filename):
		self._setup(filename)

	def _setup(self, filename):
		self._molecule = openbabel.OBMol()
		self._loadMolecule(filename)
		self._computePartialAtomicCharges()
		self._setupSmartsPattern()

	def _loadMolecule(self, filename):
		file_format = self._get_format_from_filename(filename)
		conversion = openbabel.OBConversion()
		conversion.SetInFormat(file_format)
		conversion.ReadFile(self._molecule, filename)

	def _get_format_from_filename(self,filename):
		return file_extension(filename)[1:]

	def _computePartialAtomicCharges(self):
		self._charge_model = openbabel.OBChargeModel.FindType("mmff94")
		self._charge_model.ComputeCharges(self._molecule)

	def _setupSmartsPattern(self):
		self._pattern = openbabel.OBSmartsPattern()

	def isOK(self):
		value = self.getAtomCount() > 1
		return value

	def MatchPattern(self, pattern_to_match):
		self._pattern.Init(pattern_to_match)
		self._pattern.Match(self._molecule)
		match = [m for m in self._pattern.GetUMapList()]
		return match

	def getPartialAtomCharges(self):
		return self._charge_model.GetPartialCharges()

	def getTotalCharge(self):
		return int(sum(self.getPartialAtomCharges()))

	def getElementSymbol(self,atom_index):
		table = openbabel.OBElementTable()
		return table.GetSymbol(atom_index)

	def getAtomCount(self):
		return self._molecule.NumAtoms()
