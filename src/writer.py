"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2016 Casper Steinmann
"""

try:
    import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")

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
        self._verbose = self._fragmentation.getVerbose()

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
