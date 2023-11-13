"""
Copyright (C) 2011-2016 Casper Steinmann
"""

from fragit.fragit_exceptions import OBNotFoundException

try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")
from fragit.util import file_extension, Z2LABEL


class Molecule(object):

    def __init__(self, filename):
        self._setup(filename)

    def _setup(self, filename):
        self._molecule = openbabel.OBMol()
        self._load_molecule(filename)
        self._compute_partial_atomic_charges()
        self._setup_smarts_pattern()

    def _load_molecule(self, filename):
        file_format = self._get_format_from_filename(filename)
        conversion = openbabel.OBConversion()
        conversion.SetInFormat(file_format)
        conversion.ReadFile(self._molecule, filename)

    @staticmethod
    def _get_format_from_filename(filename):
        return file_extension(filename)[1:]

    def _compute_partial_atomic_charges(self):
        self._charge_model = openbabel.OBChargeModel.FindType("mmff94")
        self._charge_model.ComputeCharges(self._molecule)

    def _setup_smarts_pattern(self):
        self._pattern = openbabel.OBSmartsPattern()

    def is_ok(self):
        value = self.get_atom_count() > 1
        return value

    def match_pattern(self, pattern_to_match):
        self._pattern.Init(pattern_to_match)
        self._pattern.Match(self._molecule)
        match = [m for m in self._pattern.GetUMapList()]
        return match

    def get_partial_atom_charges(self):
        return self._charge_model.GetPartialCharges()

    def get_total_charge(self) -> int:
        return int(sum(self.get_partial_atom_charges()))

    @staticmethod
    def get_element_symbol(atom_index: int) -> str:
        return Z2LABEL[atom_index]

    def get_atom_count(self):
        return self._molecule.NumAtoms()
