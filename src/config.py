"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2016 Casper Steinmann
"""
import sys
try:
    # python 2.X
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser

from .util import *

class FragItDataBase(dict):
    """default data for FragIt"""
    def __init__(self, *args):
        dict.__init__(self, args)
        self.data_types=dict()
        self.data_types['maxfragsize'] = int
        self.data_types['writer'] = str
        self.data_types['groupcount']=int
        self.data_types['boundaries']=str
        self.data_types['buffer']=float
        self.data_types['active']=float
        self.data_types['freezebackbone']=bool
        self.data_types['writepymol']=bool
        self.data_types['writejmol']=bool
        self.data_types['centralfragment']=int
        self.data_types['pairs']=str
        self.data_types['atomids']=str
        self.data_types['chargemodel']=str
        self.data_types['order']=int
        self.data_types['useatomnames']=bool
        self.data_types['includehbonddonors']=bool
        self.data_types['includehbondacceptors']=bool
        self.data_types['hbondangle']=float
        self.data_types['hbonddistancemin']=float
        self.data_types['hbonddistancemax']=float
        self.data_types['includecovalent']=bool
        self.data_types['includeallwithin']=float
        self.data_types['verbose']=bool

        # items here are complex values that need
        # specific parsing later on
        self.data_types['peptide']=str
        self.data_types['a-d-pyranose']=str
        self.data_types['nterminal']=str
        self.data_types['pairs']=str
        self.data_types['atomids']=str
        self.data_types['combinefragments'] = str

        self['fragmentation'] = dict()
        self['fragmentation']['maxfragsize']=100
        self['fragmentation']['writer']="XYZ"
        self['fragmentation']['groupcount']=1
        self['fragmentation']['chargemodel']="MMFF94"
        self['fragmentation']['combinefragments'] = "" # list of integers

        self['output'] = dict()
        self['output']['verbose']=True
        self['output']['boundaries']=""
        self['output']['buffer']=0.0
        self['output']['active']=0.0
        self['output']['freezebackbone']=False
        self['output']['writepymol']=False
        self['output']['writejmol']=False
        self['output']['centralfragment']=0
        self['output']['useatomnames'] = False

        # Fragmentation patterns are set in the individual settings below
        self['fragmentpatterns'] = dict()

        # Protection patterns are set in the individual settings below
        self['protectpatterns'] = dict()

        self['mergepatterns'] = dict()
        self['mergepatterns']['glycine']="" # do not merge by default

        self['explicitfragmentpairs'] = dict()
        self['explicitfragmentpairs']['pairs']="" # semi-colon separated list, i.e. 11,12;32,33;44,45

        self['explicitprotectatoms'] = dict()
        self['explicitprotectatoms']['atomids']="" # list of integers

        self['mfcc'] = dict()
        self['mfcc']['order'] = 0

        # options to control QM/MM refinement
        self['qmmm'] = dict()
        self['qmmm']['includehbonddonors'] = False
        self['qmmm']['includehbondacceptors'] = False
        # angle > 110 is according to:
        # http://pac.iupac.org/publications/pac/pdf/2011/pdf/8308x1637.pdf
        self['qmmm']['hbondangle'] = 110.0
        self['qmmm']['hbonddistancemin'] = 2.5
        self['qmmm']['hbonddistancemax'] = 3.9
        self['qmmm']['includecovalent'] = False
        self['qmmm']['includeallwithin'] = 0.0

    def getType(self, option, section):
        if "pattern" in section:
            return str

        if option not in self.data_types:
            raise ValueError("Option '%s' is not recognized." % option)

        return self.data_types[option]

class FragItDataFMO(FragItDataBase):
    """ Initializes FragIt with options which are applicable to the
        fragment molecular orbital (FMO) and related methods.
        Some options set in FragItDataBase will be overwritten.
    """
    def __init__(self):
        FragItDataBase.__init__(self)

        self['fragmentation']['writer']="GAMESS-FMO"

        # fragmentation patterns for FMO as discussed in the original PLoS ONE publication
        # DOI: 10.1371/journal.pone.0044480
        self['fragmentpatterns']['peptide']="[$(CN)][$(C(=O)NCC(=O))]"
        self['fragmentpatterns']['a-d-pyranose']="[$(C1C(CO)OC(O)C(O)C1(O))][$(OC1C(O)C(O)CC(CO)O1)]"
        self['fragmentpatterns']['dnabackbone'] = "[$(CCOP)][$(CC1OCCC1)]"

        # fragmentation patterns for FMO as discussed in the original PLoS ONE publication
        # DOI: 10.1371/journal.pone.0044480

        # protection patterns are needed to remove small fragments
        self['protectpatterns']['nterminal']="[$([NH2]),$([NH3])]CC(=O)[$(NCC=O)]"

class FragItDataPE(FragItDataBase):
    """ Initializes FragIt with options which are applicable to the
        polarizable embedding (PE) approach. This is mostly tuned
        for potential generation through the polarizable embedding
        assistant script (PEAS).
        Some options set in FragItDataBase will be overwritten.
    """
    def __init__(self):
        FragItDataBase.__init__(self)

        self['fragmentation']['writer']="XYZ-MFCC"

        # DOI: xx
        self['fragmentpatterns']['peptide']="[$([CX3](=[OX1])[NX3][CX4])][$([NX3][CX3][CX4])]"

        # DOI: xx
        self['fragmentpatterns']['dnabackbone'] = "[$(POCC)][$(OC1COCC1)]"

        # utilize the MFCC principle. Standard is cap-order 2 (for peptides)
        self['mfcc']['order'] = 2

        # use atom names when using PE.
        self['output']['useatomnames'] = True


# export all config settings so they can be
# loaded at a later time.
ConfigSettings = {'BARE': FragItDataBase, 'FMO': FragItDataFMO, 'PE': FragItDataPE}


class FragItConfig(object):
    def __init__(self, defaults=FragItDataFMO, **kwargs):
        filename = kwargs.get('filename', None)
        verbose = kwargs.get('verbose', False)
        self.cfg = RawConfigParser()
        self.values = defaults()
        self._addSections()

        if filename is not None:
            if verbose:
                print("reading from '{}'".format(filename))
            self.readConfigurationFromFile(filename)

    def _addSections(self):
        """Updates the RawRawConfigParser with values from the data array
        """
        for section in self.values.keys():
            if not self.cfg.has_section(section): self.cfg.add_section(section)
            for key in self.values[section].keys():
                value = self.values[section][key]
                if "atomids" == key or "pairs" == key:
                    value = ""
                self.cfg.set(section,key,value)

    def readConfigurationFromFile(self, filename):
        try:
                with open(filename,'r') as f: pass
        except IOError:
                print("The configuration file '{}' does not exist. Aborting.".format(filename))
                sys.exit()
        self.cfg.read(filename)

        # code to parse data from sections into values
        # do section check and sanity checks here too
        for section in self.cfg.sections():
            if section not in self.values:
                raise KeyError("Section '%s' is not recognized." % section)

            for key in self.cfg.options(section):
                if key not in self.values[section] and "pattern" not in section: # dubious hack to make custom patterns writable.
                    raise KeyError("Option '%s' in '%s' is not recognized." % (key,section))

                format = self.values.getType(key,section)
                value = format(self.cfg.get(section,key))
                if format == type(True) and type(self.cfg.get(section,key)) == type(""):
                    value = (self.cfg.get(section,key)).lower() == "true"
                self.values[section][key] = value

    def writeConfigurationToFile(self,file):
        if isinstance(file, str):
            raise ValueError("Deprecated: File parameter currently only accepts a file handle, not filename.")

        self._addSections()
        self.cfg.write(file)

    def setMaximumFragmentSize(self, value):
        if not isinstance(value, int):
            raise TypeError("Expected an integer to define maximum fragment size.")
        if value <= 0:
            raise ValueError("Maximum fragment sizes must be positive.")
        self.values['fragmentation']['maxfragsize'] = value

    def getMaximumFragmentSize(self):
        return self.values['fragmentation']['maxfragsize']

    def setMinimumFragmentSize(self, value):
        if not isinstance(value, int):
            raise TypeError("Expected an integer to define minimum fragment size.")
        if value <= 0:
            value = -1
        self.values['fragmentation']['minfragsize'] = value

    def getMinimumFragmentSize(self):
        return self.values['fragmentation']['minfragsize']

    def getChargeModel(self):
        return self.values['fragmentation']['chargemodel']

    def setChargeModel(self, value):
        #if(value.upper()) not in ["MMFF94", "GASTEIGER"]: raise ValueError("Only 'MMFF94' and 'GASTEIGER' charge-models supported.")
        self.values['fragmentation']['chargemodel'] = value.upper()

    def setFragmentGroupCount(self, value):
        if not isinstance(value, int):
            raise TypeError("Expected integer input in fragment group count.")
        if value <= 0: value = 1
        self.values['fragmentation']['groupcount'] = value

    def getFragmentGroupCount(self):
        return self.values['fragmentation']['groupcount']

    def setWriter(self,value):
        if not isinstance(value, str):
            raise TypeError
        self.values['fragmentation']['writer'] = value

    def getWriter(self):
        return self.values['fragmentation']['writer']

    def getBreakPatterns(self):
        return self.values['fragmentpatterns']

    def setBreakPatterns(self,value):
        if not isinstance(value, dict):
            raise TypeError
        self.values['fragmentpatterns'] = value

    def getProtectPatterns(self):
        return self.values['protectpatterns']

    def setProtectPatterns(self,value):
        if not isinstance(value, dict):
            raise TypeError
        self.values['protectpatterns'] = value

    def clearProtectPatterns(self):
        self.values['protectpatterns'] = dict()

    def getCombineFragments(self):
        values = self.values['fragmentation']['combinefragments']
        if len(values) > 0:
            list_of_ids = values.split(",")
            return list(map(int, list_of_ids))
        return []

    def setCombineFragments(self, value):
        if not isinstance(value, str):
            raise TypeError
        if len(value) == 0:
            return
        self.values['fragmentation']['combinefragments'] = value

    def getExplicitlyProtectedAtoms(self):
        values = self.values['explicitprotectatoms']['atomids']
        if len(values) > 0:
            list_of_ids = values.split(",")
            return list(map(int, list_of_ids))
        return []

    def addExplicitlyProtectedAtoms(self,value):
        if not isinstance(value, list): raise TypeError
        list_of_ids = self.getExplicitlyProtectedAtoms()
        list_of_ids.extend(value)
        list_of_ids = Uniqify(list_of_ids)
        list_of_ids.sort()
        list_of_ids = map(str, list_of_ids)
        self.values['explicitprotectatoms']['atomids'] = ",".join(list_of_ids)

    def getExplicitlyBreakAtomPairs(self):
        values = self.values['explicitfragmentpairs']['pairs']
        if len(values) > 0:
            if values[-1] == ";": values = values[:-1]
            list_of_ids = values.split(";")
            return list(map(self._pair_to_tuple, list_of_ids))
        return []

    def _pair_to_tuple(self,value):
        values = value.split(",")
        return tuple(map(int, values))

    def _pair_from_tuple(self,value):
        if not isinstance(value, tuple) and len(value) != 2:
            raise ValueError
        return "%i,%i" % (value[0],value[1])

    def addExplicitlyBreakAtomPairs(self,value):
        values = self.getExplicitlyBreakAtomPairs()
        if not isinstance(values, list):
            raise ValueError("Error: Expected list in addExplicitlyBreakAtomPairs. Got: '{}'".format(type(values)))
        if value not in values: values.extend(value)
        values = Uniqify(values)
        values.sort()
        values_str = map(self._pair_from_tuple, values)
        self.values['explicitfragmentpairs']['pairs'] = ";".join(values_str)

    def popExplicitlyBreakAtomPairs(self,value):
        values = self.getExplicitlyBreakAtomPairs()
        if value in values:
            values.remove(value)
        values_str = map(self._pair_from_tuple, values)
        self.values['explicitfragmentpairs']['pairs'] = ";".join(values_str)

    def getOutputFormat(self):
        return self.values['fragmentation']['writer']

    def setOutputFormat(self, value):
        if not isinstance(value, str): raise TypeError
        self.values['fragmentation']['writer'] = value

    def enableMergeGlycinePattern(self):
        self.values['mergepatterns']['glycine'] = "O=CN[CX4H2]" # use to match glycine to get fragment indices

    def getMergePatterns(self):
        return self.values['mergepatterns']

    # output options
    def getBoundaries(self):
        return self.values['output']['boundaries']

    def setBoundaries(self, value):
        self.values['output']['boundaries'] = value

    def getCentralFragmentID(self):
        return self.values['output']['centralfragment']

    def setCentralFragmentID(self, value):
        if not isinstance(value, int):
            raise TypeError
        self.values['output']['centralfragment'] = value

    def getWriteJmolScript(self):
        return self.values['output']['writejmol']

    def setWriteJmolScript(self, value):
        if not isinstance(value, bool):
            raise TypeError
        self.values['output']['writejmol'] = value

    def getWritePymolScript(self):
        return self.values['output']['writepymol']

    def setWritePymolScript(self, value):
        if not isinstance(value, bool):
            raise TypeError
        self.values['output']['writepymol'] = value

    def getFreezeBackbone(self):
        return self.values['output']['freezebackbone']

    def getBufferDistance(self):
        return self.values['output']['buffer']

    def getActiveAtomsDistance(self):
        return self.values['output']['active']

    def useAtomNames(self):
        return self.values['output']['useatomnames']

    def getVerbose(self):
        return self.values['output']['verbose']

    def getMFCCOrder(self):
        return self.values['mfcc']['order']

    # options for QM/MM
    def getHBondAngle(self):
        return self.values['qmmm']['hbondangle']

    def getHBondDistanceMin(self):
        return self.values['qmmm']['hbonddistancemin']

    def getHBondDistanceMax(self):
        return self.values['qmmm']['hbonddistancemax']

    def doQMMMHydrogenBondDonors(self):
        return self.values['qmmm']['includehbonddonors']

    def doQMMMHydrogenBondAcceptors(self):
        return self.values['qmmm']['includehbondacceptors']

    def doQMMMIncludeCovalent(self):
        return self.values['qmmm']['includecovalent']

    def doQMMMIncludeAllWithin(self):
        return abs(self.values['qmmm']['includeallwithin']) > 0.0

    def getQMMMIncludeAllWithinDistance(self):
        return self.values['qmmm']['includeallwithin']

