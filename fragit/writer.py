"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2017 Casper Steinmann
"""
from typing import List

from .fragit_exceptions import OBNotFoundException
try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")

from fragit.util import float_list_from_string, int_list_from_string


class Standard(object):
    def __init__(self, fragmentation, directories):
        self._fragmentation = fragmentation
        self._fragment_layers: List[int] = [1 for _ in fragmentation.get_fragments()]
        self._nlayers: int = 1
        self._boundaries = []
        self._title = ""
        self._active_fragments = []
        self._central_fragment = 0
        self._active_atoms_distance = 0.0
        self._freeze_backbone = False
        self._input_filename = None
        self._output_filename = None
        self._do_jmol = False
        self._do_pymol = False
        self._verbose = self._fragmentation.get_verbose()
        self._directories = directories
        self._buffer_maximum_distance = 0.0

    def write_file(self, filename):
        raise NotImplementedError

    def setup(self):
        raise NotImplementedError

    def set_boundaries_from_string(self, value: str):
        self._boundaries = float_list_from_string(value)
        self._nlayers = len(self._boundaries)+1

    def set_central_fragment_id(self, value: int):
        if value < 0:
            raise ValueError
        self._central_fragment = value

    def set_active_fragments(self, value: str):
        self._active_fragments = int_list_from_string(value)

    def set_active_atoms_distance(self, value: float):
        self._active_atoms_distance = value

    def set_buffer_max_distance(self, value: float):
        if value < 0.0:
            raise ValueError
        self._buffer_maximum_distance = value

    def set_freeze_backbone(self):
        self._freeze_backbone = True

    def set_jmol_output(self, infile: str, outfile: str):
        self._input_filename = infile
        self._output_filename = outfile
        self._do_jmol = True

    def set_pymol_output(self, infile, outfile):
        self._input_filename = infile
        self._output_filename = outfile
        self._do_pymol = True
