"""
**********************************************************************
config.py

Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2013 Casper Steinmann

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
import sys

from util import *
from ConfigParser import RawConfigParser

class FragItData(dict):
  """default data for FragIt"""
  def __init__(self, *args):
    dict.__init__(self, args)
    self.data_types=dict()
    self.data_types['maxfragsize'] = int
    #self.data_types['minfragsize'] = int
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

    # items here are complex values that need
    # specific parsing later on
    self.data_types['peptide']=str
    self.data_types['a-d-pyranose']=str
    self.data_types['nterminal']=str
    self.data_types['pairs']=str
    self.data_types['atomids']=str
    self.data_types['combinefragments'] = str
    self.data_types['lcap']=str
    self.data_types['rcap']=str

    self['fragmentation'] = dict()
    self['fragmentation']['maxfragsize']=50
    self['fragmentation']['writer']="GAMESS-FMO"
    self['fragmentation']['groupcount']=1
    self['fragmentation']['chargemodel']="MMFF94"
    self['fragmentation']['combinefragments'] = "" # list of integers

    self['output'] = dict()
    self['output']['boundaries']=""
    self['output']['buffer']=0.0
    self['output']['active']=0.0
    self['output']['freezebackbone']=False
    self['output']['writepymol']=False
    self['output']['writejmol']=False
    self['output']['centralfragment']=0
    self['output']['useatomnames'] = False

    self['fragmentpatterns'] = dict()
    self['fragmentpatterns']['peptide']="[$(CN)][$(C(=O)NCC(=O))]"
    self['fragmentpatterns']['a-d-pyranose']="[$(C1C(CO)OC(O)C(O)C1(O))][$(OC1C(O)C(O)CC(CO)O1)]"
    self['fragmentpatterns']['dnabackbone'] = "[$(CCOP)][$(CC1OCCC1)]"

    self['protectpatterns'] = dict()
    self['protectpatterns']['nterminal']="[$([NH2]),$([NH3])]CC(=O)[$(NCC=O)]"

    self['mergepatterns'] = dict()
    self['mergepatterns']['glycine']="" # do not merge by default

    self['explicitfragmentpairs'] = dict()
    self['explicitfragmentpairs']['pairs']="" # semi-colon separated list, i.e. 11,12;32,33;44,45

    self['explicitprotectatoms'] = dict()
    self['explicitprotectatoms']['atomids']="" # list of integers

    self['mfcc'] = dict()
    #self['mfcc']['lcap'] = "" # string, no default capping to the left
    #self['mfcc']['rcap'] = "" # string, no default capping to the right
    self['mfcc']['order'] = 0

  def getType(self, option, section):
    if "pattern" in section: return str
    if not self.data_types.has_key(option):
      raise ValueError("Option '%s' is not recognized." % option)
    return self.data_types[option]

class FragItConfig(object):
  def __init__(self, defaults=FragItData):
    self.cfg = RawConfigParser()
    self.values = defaults()
    self._addSections()

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
        print "The configuration file '%s' does not exist. Aborting." % filename
        sys.exit()
    self.cfg.read(filename)

    # code to parse data from sections into values
    # do section check and sanity checks here too
    for section in self.cfg.sections():
      if not self.values.has_key(section):
        raise KeyError("Section '%s' is not recognized." % section)

      for key in self.cfg.options(section):
        if not self.values[section].has_key(key) and "pattern" not in section: # dubious hack to make custom patterns writable.
          raise KeyError("Option '%s' in '%s' is not recognized." % (key,section))

        format = self.values.getType(key,section)
        value = format(self.cfg.get(section,key))
        if format == type(True):
          value = (self.cfg.get(section,key)).lower() == "true"
        self.values[section][key] = value

  def writeConfigurationToFile(self,file):
    f = open(file,"w")
    self._addSections()
    self.cfg.write(f)
    f.close()

  def setMaximumFragmentSize(self, value):
    if not is_int(value): raise TypeError
    if value <= 0: raise ValueError("Fragment sizes cannot be zero or negative")
    self.values['fragmentation']['maxfragsize'] = value

  def getMaximumFragmentSize(self):
    return self.values['fragmentation']['maxfragsize']

  def setMinimumFragmentSize(self, value):
    if not is_int(value): raise TypeError
    if value <= 0: value = -1
    self.values['fragmentation']['minfragsize'] = value

  def getMinimumFragmentSize(self):
    return self.values['fragmentation']['minfragsize']

  def getChargeModel(self):
    return self.values['fragmentation']['chargemodel']

  def setChargeModel(self, value):
    #if(value.upper()) not in ["MMFF94", "GASTEIGER"]: raise ValueError("Only 'MMFF94' and 'GASTEIGER' charge-models supported.")
    self.values['fragmentation']['chargemodel'] = value.upper()

  def setFragmentGroupCount(self, value):
    if not is_int(value): raise TypeError
    if value <= 0: value = 1
    self.values['fragmentation']['groupcount'] = value

  def getFragmentGroupCount(self):
    return self.values['fragmentation']['groupcount']

  def setWriter(self,value):
    if not is_string(value): raise TypeError
    self.values['fragmentation']['writer'] = value

  def getWriter(self):
    return self.values['fragmentation']['writer']

  def getBreakPatterns(self):
    return self.values['fragmentpatterns']

  def setBreakPatterns(self,value):
    if not is_dict(value): raise TypeError
    self.values['fragmentpatterns'] = value

  def getProtectPatterns(self):
    return self.values['protectpatterns']

  def setProtectPatterns(self,value):
    if not is_dict(value): raise TypeError
    self.values['protectpatterns'] = value

  def clearProtectPatterns(self):
    self.values['protectpatterns'] = dict()

  def getCombineFragments(self):
    values = self.values['fragmentation']['combinefragments']
    if len(values) > 0:
      list_of_ids = values.split(",")
      return map(int, list_of_ids)
    return []

  def setCombineFragments(self, value):
    if not is_string(value): raise TypeError
    if len(value) == 0: return
    self.values['fragmentation']['combinefragments'] = value

  def getExplicitlyProtectedAtoms(self):
    values = self.values['explicitprotectatoms']['atomids']
    if len(values) > 0:
      list_of_ids = values.split(",")
      return map(int, list_of_ids)
    return []

  def addExplicitlyProtectedAtoms(self,value):
    if not is_list(value): raise TypeError
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
      return map(self._pair_to_tuple, list_of_ids)
    return []

  def _pair_to_tuple(self,value):
    values = value.split(",")
    return tuple(map(int, values))

  def _pair_from_tuple(self,value):
    if not is_tuple(value) and len(value) != 2: raise ValueError
    return "%i,%i" % (value[0],value[1])

  def addExplicitlyBreakAtomPairs(self,value):
    values = self.getExplicitlyBreakAtomPairs()
    if value not in values: values.extend(value)
    values = Uniqify(values)
    values.sort()
    values_str = map(self._pair_from_tuple, values)
    self.values['explicitfragmentpairs']['pairs'] = ";".join(values_str)

  def popExplicitlyBreakAtomPairs(self,value):
    values = self.getExplicitlyBreakAtomPairs()
    if value in values: values.remove(value)
    values_str = map(self._pair_from_tuple, values)
    self.values['explicitfragmentpairs']['pairs'] = ";".join(values_str)

  def getOutputFormat(self):
    return self.values['fragmentation']['writer']

  def setOutputFormat(self, value):
    if not is_string(value): raise TypeError
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
    if not is_int(value): raise TypeError
    self.values['output']['centralfragment'] = value

  def getWriteJmolScript(self):
    return self.values['output']['writejmol']

  def setWriteJmolScript(self, value):
    if not is_bool(value): raise TypeError
    self.values['output']['writejmol'] = value

  def getWritePymolScript(self):
    return self.values['output']['writepymol']

  def setWritePymolScript(self, value):
    if not is_bool(value): raise TypeError
    self.values['output']['writepymol'] = value

  def getFreezeBackbone(self):
    return self.values['output']['freezebackbone']

  def getBufferDistance(self):
    return self.values['output']['buffer']

  def getActiveAtomsDistance(self):
    return self.values['output']['active']

  def useAtomNames(self):
    return self.values['output']['useatomnames']


#  def getMFCCLeftCap(self):
#    return self.values['mfcc']['lcap']
#
#  def getMFCCRightCap(self):
#    return self.values['mfcc']['rcap']

  def getMFCCOrder(self):
    return self.values['mfcc']['order']

if __name__ == '__main__':
  cfg = FragItConfig()
  cfg.readConfigurationFromFile("my1.conf")
  cfg.writeConfigurationToFile("my2.conf")
