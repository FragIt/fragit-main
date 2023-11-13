"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2017 Casper Steinmann
"""
from configparser import RawConfigParser
import sys
from typing import List, Tuple, Dict

from fragit.util import remove_duplicates


class FragItDataBase(dict):
    """default data for FragIt"""

    def __init__(self, *args):
        dict.__init__(self, args)
        self.data_types = dict()
        self.data_types["maxfragsize"] = int
        self.data_types["writer"] = str
        self.data_types["groupcount"] = int
        self.data_types["boundaries"] = str
        self.data_types["buffer"] = float
        self.data_types["active"] = float
        self.data_types["freezebackbone"] = bool
        self.data_types["writepymol"] = bool
        self.data_types["writejmol"] = bool
        self.data_types["centralfragment"] = int
        self.data_types["pairs"] = str
        self.data_types["atomids"] = str
        self.data_types["chargemodel"] = str
        self.data_types["order"] = int
        self.data_types["useatomnames"] = bool
        self.data_types["includehbonddonors"] = bool
        self.data_types["includehbondacceptors"] = bool
        self.data_types["hbondangle"] = float
        self.data_types["hbonddistancemin"] = float
        self.data_types["hbonddistancemax"] = float
        self.data_types["includecovalent"] = bool
        self.data_types["includeallwithin"] = float
        self.data_types["verbose"] = bool
        self.data_types["dohop"] = bool
        self.data_types["efpwaters"] = int

        # items here are complex values that need
        # specific parsing later on
        self.data_types["peptide"] = str
        self.data_types["a-d-pyranose"] = str
        self.data_types["nterminal"] = str
        self.data_types["pairs"] = str
        self.data_types["atomids"] = str
        self.data_types["combinefragments"] = str
        self.data_types["basis"] = str

        self["fragmentation"] = dict()
        self["fragmentation"]["maxfragsize"] = 100
        self["fragmentation"]["writer"] = "XYZ"
        self["fragmentation"]["groupcount"] = 1
        self["fragmentation"]["chargemodel"] = "MMFF94"
        self["fragmentation"]["combinefragments"] = ""  # list of integers

        self["output"] = dict()
        self["output"]["verbose"] = True
        self["output"]["boundaries"] = ""
        self["output"]["buffer"] = 0.0
        self["output"]["active"] = 0.0
        self["output"]["freezebackbone"] = False
        self["output"]["writepymol"] = False
        self["output"]["writejmol"] = False
        self["output"]["centralfragment"] = 0
        self["output"]["useatomnames"] = False

        # Fragmentation patterns are set in the individual settings below
        self["fragmentpatterns"] = dict()

        # Protection patterns are set in the individual settings below
        self["protectpatterns"] = dict()

        self["mergepatterns"] = dict()
        self["mergepatterns"]["glycine"] = ""  # do not merge by default

        self["explicitfragmentpairs"] = dict()
        self["explicitfragmentpairs"]["pairs"] = ""  # semi-colon separated list, i.e. 11,12;32,33;44,45

        self["explicitprotectatoms"] = dict()
        self["explicitprotectatoms"]["atomids"] = ""  # list of integers

        self["mfcc"] = dict()
        self["mfcc"]["order"] = 0

        # options to control QM/MM refinement
        self["qmmm"] = dict()
        self["qmmm"]["includehbonddonors"] = False
        self["qmmm"]["includehbondacceptors"] = False
        # angle > 110 is according to:
        # http://pac.iupac.org/publications/pac/pdf/2011/pdf/8308x1637.pdf
        self["qmmm"]["hbondangle"] = 110.0
        self["qmmm"]["hbonddistancemin"] = 2.5
        self["qmmm"]["hbonddistancemax"] = 3.9
        self["qmmm"]["includecovalent"] = False
        self["qmmm"]["includeallwithin"] = 0.0

        # options for QM level of theory
        self["qm"] = dict()
        self["qm"]["basis"] = ""  # colon separated list for multi-layer runs

        # fmo specific options
        self["fmo"] = dict()
        self["fmo"]["dohop"] = False
        self["fmo"]["efpwaters"] = 0  # disable EFP waters

    def get_type(self, option, section):
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

        self["fragmentation"]["writer"] = "GAMESS-FMO"

        # fragmentation patterns for FMO as discussed in the original PLoS ONE publication
        # DOI: 10.1371/journal.pone.0044480
        self["fragmentpatterns"]["peptide"] = "[$(CN)][$(C(=O)NCC(=O))]"
        self["fragmentpatterns"]["a-d-pyranose"] = "[$(C1C(CO)OC(O)C(O)C1(O))][$(OC1C(O)C(O)CC(CO)O1)]"
        self["fragmentpatterns"]["dnabackbone"] = "[$(CCOP)][$(CC1OCCC1)]"

        # fragmentation patterns for FMO as discussed in the original PLoS ONE publication
        # DOI: 10.1371/journal.pone.0044480

        # protection patterns are needed to remove small fragments
        self["protectpatterns"]["nterminal"] = "[$([NH2]),$([NH3])]CC(=O)[$(NCC=O)]"

        # don"t use atom names when using FMO. This can cause
        # annoying errors in GAMESS
        self["output"]["useatomnames"] = False

        # basis set for QM calculations
        self["qm"]["basis"] = "3-21G"


class FragItDataPE(FragItDataBase):
    """ Initializes FragIt with options which are applicable to the
        polarizable embedding (PE) approach. This is mostly tuned
        for potential generation through the polarizable embedding
        assistant script (PEAS).
        Some options set in FragItDataBase will be overwritten.
    """

    def __init__(self):
        FragItDataBase.__init__(self)

        self["fragmentation"]["writer"] = "XYZ-MFCC"

        # DOI: xx
        self["fragmentpatterns"]["peptide"] = "[$([CX3](=[OX1])[NX3][CX4])][$([NX3][CX3][CX4])]"

        # DOI: xx
        self["fragmentpatterns"]["dnabackbone"] = "[$(POCC)][$(OC1COCC1)]"

        # utilize the MFCC principle. Standard is cap-order 2 (for peptides)
        self["mfcc"]["order"] = 2

        # use atom names when using PE.
        self["output"]["useatomnames"] = True


# export all config settings so they can be
# loaded at a later time.
ConfigSettings = {'BARE': FragItDataBase, 'FMO': FragItDataFMO, 'PE': FragItDataPE}


class FragItConfig(object):
    def __init__(self, defaults=FragItDataBase, **kwargs):
        filename = kwargs.get("filename", None)
        verbose = kwargs.get("verbose", False)
        self.cfg = RawConfigParser()
        self.values = defaults()
        self._add_sections()

        if filename is not None:
            # this code needs to be fixed but the verbosity flag should probably
            # be on by default (lots of stuff printed) and then we should probably
            # have a quench flag if one wants to remove all the output (for some reason)
            print("Info: FragIt [CONFIG] reading config from '{}':".format(filename))
            print("")
            print("   Contents of the config file")
            print(" -------------------------------")
            with open(filename, 'r') as f:
                print(f.read())
            print(" -------------------------------")

            self.read_configuration_from_file(filename)

    def _add_sections(self):
        """Updates the RawRawConfigParser with values from the data array
        """
        for section in self.values.keys():
            if not self.cfg.has_section(section):
                self.cfg.add_section(section)
            for key in self.values[section].keys():
                value = self.values[section][key]
                if "atomids" == key or "pairs" == key:
                    value = ""
                self.cfg.set(section, key, value)

    def read_configuration_from_file(self, filename):
        try:
            with open(filename, 'r') as f:
                pass
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
                if key not in self.values[section] and "pattern" not in section:  # dubious hack to make custom patterns writable.
                    raise KeyError("Option '%s' in '%s' is not recognized." % (key, section))

                fmt = self.values.get_type(key, section)
                value = fmt(self.cfg.get(section, key))
                if fmt == type(True) and type(self.cfg.get(section, key)) == type(""):
                    value = (self.cfg.get(section, key)).lower() == "true"
                self.values[section][key] = value

    def write_configuration_to_file(self, file):
        if isinstance(file, str):
            raise ValueError("Deprecated: File parameter currently only accepts a file handle, not filename.")

        self._add_sections()
        self.cfg.write(file)

    def set_maximum_fragment_size(self, value: int):
        if not isinstance(value, int):
            raise TypeError("Expected an integer to define maximum fragment size.")
        if value <= 0:
            raise ValueError("Maximum fragment sizes must be positive.")
        self.values["fragmentation"]["maxfragsize"] = value

    def get_maximum_fragment_size(self) -> int:
        return int(self.values["fragmentation"]["maxfragsize"])

    def set_minimum_fragment_size(self, value: int):
        if not isinstance(value, int):
            raise TypeError("Expected an integer to define minimum fragment size.")
        if value <= 0:
            value = -1
        self.values["fragmentation"]["minfragsize"] = value

    def get_minimum_fragment_size(self) -> int:
        return int(self.values["fragmentation"]["minfragsize"])

    def get_charge_model(self) -> str:
        return str(self.values["fragmentation"]["chargemodel"])

    def set_charge_model(self, value: str):
        self.values["fragmentation"]["chargemodel"] = value.upper()

    def set_fragment_group_count(self, value: int):
        if not isinstance(value, int):
            raise TypeError("Expected integer input in fragment group count.")
        if value <= 0:
            value = 1
        self.values["fragmentation"]["groupcount"] = value

    def get_fragment_group_count(self) -> int:
        return int(self.values["fragmentation"]["groupcount"])

    def set_writer(self, value: str):
        if not isinstance(value, str):
            raise TypeError
        self.values["fragmentation"]["writer"] = value

    def get_writer(self) -> str:
        return str(self.values["fragmentation"]["writer"])

    def get_break_patterns(self) -> Dict[str, str]:
        output = dict()
        for key in self.values["fragmentpatterns"]:
            output[key] = self.values["fragmentpatterns"][key]
        return output

    def set_break_patterns(self, value: Dict[str, str]):
        if not isinstance(value, dict):
            raise TypeError
        self.values["fragmentpatterns"] = value

    def get_protect_patterns(self) -> Dict[str, str]:
        output = dict()
        for key in self.values["protectpatterns"]:
            output[key] = self.values["protectpatterns"][key]
        return output

    def set_protect_patterns(self, value: Dict[str, str]):
        if not isinstance(value, dict):
            raise TypeError
        self.values["protectpatterns"] = value

    def clear_protect_patterns(self):
        self.values["protectpatterns"] = dict()

    def get_combine_fragments(self) -> List[int]:
        values = self.values["fragmentation"]["combinefragments"]
        if len(values) > 0:
            list_of_ids = values.split(",")
            return list(map(int, list_of_ids))
        return []

    def set_combine_fragments(self, value: str):
        if not isinstance(value, str):
            raise TypeError
        if len(value) == 0:
            return
        self.values["fragmentation"]["combinefragments"] = value

    def get_explicitly_protected_atoms(self) -> List[int]:
        values = self.values["explicitprotectatoms"]["atomids"]
        if len(values) > 0:
            list_of_ids = values.split(",")
            return list(map(int, list_of_ids))
        return []

    def add_explicitly_protected_atoms(self, value: List[int]):
        if not isinstance(value, list):
            raise TypeError
        list_of_ids = self.get_explicitly_protected_atoms()
        list_of_ids.extend(value)
        list_of_ids = remove_duplicates(list_of_ids)
        list_of_ids.sort()
        output = map(str, list_of_ids)
        self.values["explicitprotectatoms"]["atomids"] = ",".join(output)

    def get_explicitly_break_atom_pairs(self) -> List[Tuple[int, int]]:
        values = self.values["explicitfragmentpairs"]["pairs"]
        if len(values) > 0:
            if values[-1] == ";":
                values = values[:-1]
            list_of_ids = values.split(";")
            return [self._pair_to_tuple(s) for s in list_of_ids]
        return []

    @staticmethod
    def _pair_to_tuple(value: str) -> Tuple[int, int]:
        values: List[str] = value.split(",")
        return int(values[0]), int(values[1])

    @staticmethod
    def _pair_from_tuple(value: Tuple[int, int]) -> str:
        if not isinstance(value, tuple) and len(value) != 2:
            raise ValueError
        return f"{value[0]:d},{value[1]:d}"

    def add_explicitly_break_atom_pairs(self, value):
        values = self.get_explicitly_break_atom_pairs()
        if not isinstance(values, list):
            raise ValueError("Error: Expected list in addExplicitlyBreakAtomPairs. Got: '{}'".format(type(values)))
        if value not in values:
            values.extend(value)
        values = remove_duplicates(values)
        values.sort()
        values_str = map(self._pair_from_tuple, values)
        self.values["explicitfragmentpairs"]["pairs"] = ";".join(values_str)

    def pop_explicitly_break_atom_pairs(self, value):
        values = self.get_explicitly_break_atom_pairs()
        if value in values:
            values.remove(value)
        values_str = map(self._pair_from_tuple, values)
        self.values["explicitfragmentpairs"]["pairs"] = ";".join(values_str)

    def get_output_format(self):
        return self.values["fragmentation"]["writer"]

    def set_output_format(self, value):
        if not isinstance(value, str): raise TypeError
        self.values["fragmentation"]["writer"] = value

    def enable_merge_glycine_pattern(self):
        self.values["mergepatterns"]["glycine"] = "O=CN[CX4H2]"  # use to match glycine to get fragment indices

    def get_merge_patterns(self):
        return self.values["mergepatterns"]

    # output options
    def get_boundaries(self) -> str:
        return str(self.values["output"]["boundaries"])

    def set_boundaries(self, value: str):
        self.values["output"]["boundaries"] = value

    def get_central_fragment_id(self):
        return self.values["output"]["centralfragment"]

    def set_central_fragment_id(self, value):
        if not isinstance(value, int):
            raise TypeError
        self.values["output"]["centralfragment"] = value

    def get_write_jmol_script(self):
        return self.values["output"]["writejmol"]

    def set_write_jmol_script(self, value):
        if not isinstance(value, bool):
            raise TypeError
        self.values["output"]["writejmol"] = value

    def get_write_pymol_script(self):
        return self.values["output"]["writepymol"]

    def set_write_pymol_script(self, value):
        if not isinstance(value, bool):
            raise TypeError
        self.values["output"]["writepymol"] = value

    def get_freeze_backbone(self):
        return self.values["output"]["freezebackbone"]

    def get_buffer_distance(self):
        return self.values["output"]["buffer"]

    def get_active_atoms_distance(self):
        return self.values["output"]["active"]

    def use_atom_names(self):
        return self.values["output"]["useatomnames"]

    def get_verbose(self):
        return self.values["output"]["verbose"]

    def set_verbose(self, value):
        if not isinstance(value, bool):
            raise TypeError
        self.values["output"]["verbose"] = value

    def get_mfcc_order(self):
        return self.values["mfcc"]["order"]

    # options for QM/MM
    def get_h_bond_angle(self) -> float:
        return float(self.values["qmmm"]["hbondangle"])

    def get_h_bond_distance_min(self) -> float:
        return float(self.values["qmmm"]["hbonddistancemin"])

    def get_h_bond_distance_max(self) -> float:
        return float(self.values["qmmm"]["hbonddistancemax"])

    def do_qmmm_hydrogen_bond_donors(self) -> bool:
        return bool(self.values["qmmm"]["includehbonddonors"])

    def do_qmmm_hydrogen_bond_acceptors(self) -> bool:
        return bool(self.values["qmmm"]["includehbondacceptors"])

    def do_qmmm_include_covalent(self) -> bool:
        return bool(self.values["qmmm"]["includecovalent"])

    def do_qmmm_include_all_within(self) -> bool:
        return bool(abs(self.values["qmmm"]["includeallwithin"]) > 0.0)

    def get_qmmm_include_all_within_distance(self) -> float:
        return float(self.values["qmmm"]["includeallwithin"])

    # options for QM
    def get_qm_basis(self) -> List[str]:
        return list(self.values["qm"]["basis"].split(":"))

    def set_qm_basis(self, value: str):
        self.values["qm"]["basis"] = value

    # options for FMO
    def do_fmohop_fragmentation(self) -> bool:
        return bool(self.values["fmo"]["dohop"])

    def set_fmoafo_fragmentation(self):
        self.values["fmo"]["dohop"] = False

    def set_fmohop_fragmentation(self):
        self.values["fmo"]["dohop"] = True

    def do_fmoefp_waters(self) -> bool:
        return bool(self.get_fmoefp_waters_from_layer() > 0)

    def get_fmoefp_waters_from_layer(self) -> int:
        return int(self.values["fmo"]["efpwaters"])

    def set_fmoefp_waters_from_layer(self, value: int):
        if not isinstance(value, int):
            raise TypeError
        self.values["fmo"]["efpwaters"] = value
