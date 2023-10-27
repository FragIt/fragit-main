"""
Copyright (C) 2011-2016 Casper Steinmann
"""

from .fragit_exceptions import OBNotFoundException

try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")
from .util import file_extension, Z2LABEL

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
        return Z2LABEL[atom_index]

    def getAtomCount(self):
        return self._molecule.NumAtoms()
