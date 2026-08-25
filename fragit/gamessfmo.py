"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2023 Casper Steinmann
"""
import os
from typing import Dict, List, Tuple

# from numpy import sqrt, dot, where, array
import numpy as np

from openbabel import openbabel

from fragit.writer import Standard
from fragit.util import write_string_to_file
from fragit.util import list_to_2d, list_2d_to_str
from fragit.util import list_to_ranges, list_of_ranges_to_string, remove_duplicates, flatten
from fragit.util import Z2LABEL, LABEL2Z

# $BASIS set data depending on basis set
GAMESS_BASIS_GROUP = dict()
GAMESS_BASIS_GROUP['STO-3G'] = "GBASIS=STO NGAUSS=3"
GAMESS_BASIS_GROUP['3-21G'] = "GBASIS=N21 NGAUSS=3"
GAMESS_BASIS_GROUP['6-31G'] = "GBASIS=N31 NGAUSS=6"
GAMESS_BASIS_GROUP['6-31G*'] = "GBASIS=N31 NGAUSS=6 NDFUNC=1"
GAMESS_BASIS_GROUP['6-31G(d)'] = "GBASIS=N31 NGAUSS=6 NDFUNC=1"
GAMESS_BASIS_GROUP['6-31+G*'] = "GBASIS=N31 NGAUSS=6 NDFUNC=1 DIFFSP=.T."
GAMESS_BASIS_GROUP['6-31+G(d)'] = "GBASIS=N31 NGAUSS=6 NDFUNC=1 DIFFSP=.T."
GAMESS_BASIS_GROUP['cc-pVDZ'] = "GBASIS=CCD"
GAMESS_BASIS_GROUP['cc-pVTZ'] = "GBASIS=CCT"
GAMESS_BASIS_GROUP['aug-cc-pVDZ'] = "GBASIS=ACCD"
GAMESS_BASIS_GROUP['aug-cc-pVTZ'] = "GBASIS=ACCT"
GAMESS_BASIS_GROUP['DFTB-C'] = "GBASIS=DFTB"

# basis set data for atoms in $DATA group
GAMESS_DATA_BASIS: Dict[str, Dict[str, str]] = dict()
GAMESS_DATA_BASIS['STO-3G'] = dict()
GAMESS_DATA_BASIS['STO-3G']['H'] = 'STO 3'
GAMESS_DATA_BASIS['STO-3G']['C'] = 'STO 3'
GAMESS_DATA_BASIS['STO-3G']['N'] = 'STO 3'
GAMESS_DATA_BASIS['STO-3G']['O'] = 'STO 3'
GAMESS_DATA_BASIS['STO-3G']['S'] = 'STO 3'

GAMESS_DATA_BASIS['3-21G'] = dict()
GAMESS_DATA_BASIS['3-21G']['H'] = 'N21 3'
GAMESS_DATA_BASIS['3-21G']['C'] = 'N21 3'
GAMESS_DATA_BASIS['3-21G']['N'] = 'N21 3'
GAMESS_DATA_BASIS['3-21G']['O'] = 'N21 3'
GAMESS_DATA_BASIS['3-21G']['S'] = 'N21 3'

GAMESS_DATA_BASIS['6-31G'] = dict()
GAMESS_DATA_BASIS['6-31G']['H'] = 'N31 6'
GAMESS_DATA_BASIS['6-31G']['C'] = 'N31 6'
GAMESS_DATA_BASIS['6-31G']['N'] = 'N31 6'
GAMESS_DATA_BASIS['6-31G']['O'] = 'N31 6'
GAMESS_DATA_BASIS['6-31G']['S'] = 'N31 6'

# see emsl database for specifics
GAMESS_DATA_BASIS['6-31G(d)'] = dict()
GAMESS_DATA_BASIS['6-31G(d)']['H'] = 'N31 6'
GAMESS_DATA_BASIS['6-31G(d)']['C'] = 'N31 6\n  D   1\n    1      0.8000000              1.0000000'
GAMESS_DATA_BASIS['6-31G(d)']['N'] = 'N31 6\n  D   1\n    1      0.8000000              1.0000000'
GAMESS_DATA_BASIS['6-31G(d)']['O'] = 'N31 6\n  D   1\n    1      0.8000000              1.0000000'
GAMESS_DATA_BASIS['6-31G(d)']['S'] = 'N31 6\n  D   1\n    1      0.6500000              1.0000000'

GAMESS_DATA_BASIS['6-31G*'] = dict()
GAMESS_DATA_BASIS['6-31G*']['H'] = 'N31 6'
GAMESS_DATA_BASIS['6-31G*']['C'] = 'N31 6\n  D   1\n    1      0.8000000              1.0000000'
GAMESS_DATA_BASIS['6-31G*']['N'] = 'N31 6\n  D   1\n    1      0.8000000              1.0000000'
GAMESS_DATA_BASIS['6-31G*']['O'] = 'N31 6\n  D   1\n    1      0.8000000              1.0000000'
GAMESS_DATA_BASIS['6-31G*']['S'] = 'N31 6\n  D   1\n    1      0.6500000              1.0000000'

GAMESS_DATA_BASIS['cc-pVDZ'] = dict()
GAMESS_DATA_BASIS['cc-pVDZ']['H'] = 'CCD'
GAMESS_DATA_BASIS['cc-pVDZ']['C'] = 'CCD'
GAMESS_DATA_BASIS['cc-pVDZ']['N'] = 'CCD'
GAMESS_DATA_BASIS['cc-pVDZ']['O'] = 'CCD'
GAMESS_DATA_BASIS['cc-pVDZ']['S'] = 'CCD'

GAMESS_DATA_BASIS['pcseg-0'] = dict()
GAMESS_DATA_BASIS['pcseg-0']['H'] = 'PCseg-0'
GAMESS_DATA_BASIS['pcseg-0']['C'] = 'PCseg-0'
GAMESS_DATA_BASIS['pcseg-0']['N'] = 'PCseg-0'
GAMESS_DATA_BASIS['pcseg-0']['O'] = 'PCseg-0'
GAMESS_DATA_BASIS['pcseg-0']['S'] = 'PCseg-0'

GAMESS_DATA_BASIS['pcseg-1'] = dict()
GAMESS_DATA_BASIS['pcseg-1']['H'] = 'PCseg-1'
GAMESS_DATA_BASIS['pcseg-1']['C'] = 'PCseg-1'
GAMESS_DATA_BASIS['pcseg-1']['N'] = 'PCseg-1'
GAMESS_DATA_BASIS['pcseg-1']['O'] = 'PCseg-1'
GAMESS_DATA_BASIS['pcseg-1']['S'] = 'PCseg-1'


class GamessFMO(Standard):
    def __init__(self, fragmentation, directories):
        Standard.__init__(self, fragmentation, directories)

        self._water_fragments = []
        self._active_atoms: List[int] = []

    def setup(self):
        self.setup_layer_information()
        self._setup_active_fragments_information()
        self.validate_multi_layer_information()
        if self._do_pymol:
            self._dump_pymol()
        if self._do_jmol:
            self._dump_jmol()

    def setup_layer_information(self):
        self._fragment_layers = self.compute_fragment_layers()

    def compute_fragment_layers(self) -> List[int]:
        fragments = self._fragmentation.get_fragments()
        layers = [1 for _ in fragments]

        if self._fragmentation.do_fmoefp_waters():
            #  - figure out which fragments these oxygens belong to.
            #    If this is not a multilayer calculation, ALL
            #    waters are promoted (demoted?) to EFP waters
            water_oxygen = self._fragmentation.get_water_molecules()
            self._water_fragments = list(sorted([self._get_fragment_from_atom(i) for i in water_oxygen]))

        if self._central_fragment == 0:
            return layers
        other_fragment = fragments[self._central_fragment-1]
        distances = self.get_fragment_distances_vector(other_fragment)
        layers = self.get_layers_from_distances(distances)

        if self._fragmentation.do_fmoefp_waters():
            #  - if it is a multilayer calculation, all waters in
            #    the layer specified by getFMOEFPWatersFromLayer
            #    are promoted (demoted?) to EFP waters and multilayer
            #    is switched off.
            #    self._fragment_layers[active_fragment_id] = 2
            water_fragments = []
            for idx, water_fragment in enumerate(sorted(self._water_fragments)):
                if layers[water_fragment] == self._fragmentation.get_fmoefp_waters_from_layer():
                    water_fragments.append(water_fragment)
            self._water_fragments = water_fragments[:]
            layers = [1 for _ in fragments]
            self._nlayers = 1
        return layers

    def _dump_pymol(self):
        from fragit.pymol import PymolTemplate
        pt = PymolTemplate(self._directories, self._input_filename, self._output_filename)
        self._set_template_data(pt)
        self._write_template_file(pt)

    def _dump_jmol(self):
        from fragit.jmol import JmolTemplate
        pt = JmolTemplate(self._directories, self._input_filename, self._output_filename)
        self._set_template_data(pt)
        self._write_template_file(pt)

    def _set_template_data(self, template):
        template.set_fragments_data(self._fragmentation.get_fragments())
        template.set_buffer_data(self._fragment_layers)
        template.set_active_data(self._active_atoms)
        template.set_backbone_data(self._fragmentation.get_backbone_atoms())
        template.set_pair_data(self._fragmentation.get_explicitly_break_atom_pairs())
        template.set_fragment_charges(self._fragmentation.get_fragment_charges())

    @staticmethod
    def _write_template_file(template):
        template.override()
        template.write()

    def _setup_active_fragments_information(self):
        active_atoms = self._get_active_atoms_from_fragments()
        if self._active_atoms_distance > 0.0:
            active_atoms = self._get_active_atoms_from_distance()
            if self._verbose:
                print("Info: FragIt [GAMESS-FMO] found {0:d} atoms which should be active".format(len(active_atoms)))

        self._active_atoms = active_atoms[:]
        atoms = self._active_atoms[:]

        # if no fragment is specified as central, just
        # keep current multilayer stuff (if any) and exit
        if self._central_fragment == 0:
            self._active_frags = []
            return

        # 1) If there are active atoms, we must find the associated fragments
        # 2) The associated fragments are labelled active
        # 3) The active fragments have their atoms made flexible
        # we now have region A
        if len(self._active_atoms) > 0:
            active_frags = list()
            active_frags.append(self._central_fragment - 1)  # central must also be active
            for atom in atoms:
                ifrg = self._get_fragment_from_atom(atom)
                active_frags.extend([ifrg])
            active_frags = remove_duplicates(active_frags)
            active_frags = sorted(active_frags)
            self._active_fragments = active_frags[:]

            # promote active fragments to layer 2
            for active_fragment_id in active_frags:
                self._fragment_layers[active_fragment_id] = 2
    
            # add active fragment atoms to list of active atoms
            fragments = self._fragmentation.get_fragments()
            for frag in active_frags:
                atoms.extend(fragments[frag])
            atoms = remove_duplicates(atoms)
            atoms = sorted(atoms)
            if self._verbose and len(atoms) != len(self._active_atoms):
                print("Info: FragIt [GAMESS-FMO] active region is now {0:d} atoms ({1:d} fragments) ".format(len(active_atoms), len(active_frags)))

        # Optionally freeze backbone atoms in the active region
        if self._freeze_backbone:
            for item in self._fragmentation.get_backbone_atoms():
                if item in atoms:
                    atoms.remove(item)
                    continue

        atoms = remove_duplicates(atoms)
        atoms = sorted(atoms)
        self._active_atoms = atoms[:]

    def _get_fragment_from_atom(self, atom: int) -> int:
        for i, fragment in enumerate(self._fragmentation.get_fragments()):
            if atom in fragment:
                return i
        raise ValueError(f"Could not find atom {atom:d} in any fragment.")

    def validate_multi_layer_information(self):
        """ Validates multilayer (and FD) information and attempts
            to fix any issues should they be present.

            This method makes sure that the buffer region, b, around the active region A
            does not have close contacts between A and F.
        """
        active_atoms = self._active_atoms[:]
        active_fragments = self._active_fragments[:]
        fragments = self._fragmentation.get_fragments()
        for atom in active_atoms:  # this is for specifying something lame on the command line
            if not self.is_atom_in_active_layer(atom):
                self._active_atoms.remove(atom)

        # Here we make sure that fragments in A are not physically close
        # to fragments in F.
        # We do this by extending B (which includes A) with a buffer region, b,
        # with the buffer distance that a user wants.
        # Technically we promote the fragments in b to layer 2
        if len(active_atoms) > 0 and len(active_fragments) > 0:
            for active_fragment_index in active_fragments:
                active_fragment = fragments[active_fragment_index]
                distances = self.get_fragment_distances_vector(active_fragment)
                selection = np.where(np.asarray(distances) < self._buffer_maximum_distance)
                for index in selection[0]:
                    self._fragment_layers[index] = 2

            if self._verbose:
                frags = list()
                for i, k in enumerate(self._fragment_layers):
                    if k == 2:
                        frags.append(i)

                atms = list()
                for frag in frags:
                    atms.extend(fragments[frag])

                print("Info: FragIt [GAMESS-FMO] region B is {0:d} atoms ({1:d} fragments)".format(len(atms), len(frags)))

        self._active_atoms = active_atoms[:]

    def is_atom_in_active_layer(self, atom: int) -> bool:
        frags = self._fragmentation.get_fragments()
        fragment_layers = self._fragment_layers[:]
        for i, fragment in enumerate(frags):
            if (atom in fragment) and (fragment_layers[i] == 2):
                return True
        return False

    def write_file(self, filename):
        out_string_template = "%s%s%s%s%s%s%s%s%s%s%s%s%s"
        out_string = out_string_template % (self.write_gamess_system_group(),
                                            self.write_gamess_gddi_group(),
                                            self.write_gamess_scf_group(),
                                            self.write_gamess_contrl_group(),
                                            self.write_gamess_basis_group(),
                                            self.write_fmo_fmoprp_group(),
                                            self.write_fmo_fmo_group(),
                                            self.write_fmo_fmobnd_group(),
                                            self.write_gamess_data_group(),
                                            self.write_fmo_hyb_group(),
                                            self.write_fmo_fmoxyz_group(),
                                            self.write_fmo_fmoefp_group(),
                                            self.write_gamess_efrag_group())
        write_string_to_file(filename, out_string)

    def write_fmo_hyb_group(self) -> str:
        """ Generates the FMOHYB input group """
        s = ""
        nbonds_broken = self._fragmentation.get_num_broken_bonds()
        dohop = self._fragmentation.do_fmohop_fragmentation()
        if nbonds_broken > 0 and dohop:
            s += " $FMOHYB\n"
            basis_sets = self._fragmentation.get_qm_basis()
            nbas = len(basis_sets)
            nlayers = self._nlayers
            if nbas == 1:
                nlayers = nbas
            for ilayer in range(nlayers):
                basis = basis_sets[ilayer]
                path_to_hmo = os.path.join(self._directories['share'], "hmo/{0:s}".format(basis))
                with open(path_to_hmo, 'r') as f:
                    s += "".join(f.readlines())
            s += " $END\n"
        return s

    @staticmethod
    def write_gamess_system_group() -> str:
        return " $SYSTEM MWORDS=125 $END\n"

    @staticmethod
    def write_gamess_scf_group() -> str:
        return " $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"

    @staticmethod
    def write_gamess_gddi_group() -> str:
        return " $GDDI NGROUP=1 $END\n"

    def write_fmo_fmoprp_group(self) -> str:
        return " $FMOPRP NPRINT=9 NGUESS=%i $END\n" % self.calculate_fmo_nguess()

    def calculate_fmo_nguess(self) -> int:
        nguess = 2
        if self._nlayers > 1 and len(self._active_atoms) != 0:
            nguess += 128  # this is needed for multilayer formulation, but not FD
        return nguess

    def write_gamess_contrl_group(self) -> str:
        """ Returns the $CONTRL group

            if a geometry optimzation is requested this method
            also returns the $STATPT group
        """

        nbonds_broken = self._fragmentation.get_num_broken_bonds()
        dohop = self._fragmentation.do_fmohop_fragmentation()
        localize = ""
        if nbonds_broken > 0 and not dohop:
            localize = " LOCAL=BOYS"
        base = " $CONTRL NPRINT=-5 ISPHER=1%s\n         RUNTYP=%s\n $END\n"
        statpt = " $STATPT OPTTOL=5.0e-4 NSTEP=2000\n%s\n $END\n"
        if len(self._active_fragments) == 0 and self._active_atoms_distance <= 0.0:
            return base % (localize, "ENERGY")
        else:
            base_final = base % (localize, "OPTIMIZE")
            active_string = "      IACTAT(1)=%s" % list_of_ranges_to_string(list_to_ranges(self._active_atoms),
                                                                            maxlength=40,
                                                                            line_format="%5s",
                                                                            item_format="%s,",
                                                                            tuple_format="%s,%s,",
                                                                            terminator_format=None)
            statpt_final = statpt % active_string
            return "{0:s}{1:s}".format(base_final, statpt_final)

    def _get_active_atoms_from_fragments(self) -> List[int]:
        atoms = []
        fragments = self._fragmentation.get_fragments()
        for idx in self._active_fragments:
            atoms.extend(fragments[idx-1])
        return sorted(atoms)

    def _get_active_atoms_from_distance(self):
        atoms = []
        central_atoms = self._fragmentation.get_fragments()[self._central_fragment - 1]
        atoms.extend(central_atoms)
        all_atoms = range(1, len(self._fragmentation.get_ob_atoms()) + 1)
        for atom_idx in central_atoms:
            for atom_jdx in all_atoms:
                if atom_jdx in central_atoms or atom_jdx in atoms:
                    continue
                distance = self.get_distance_between_atoms(atom_idx, atom_jdx)
                if distance < self._active_atoms_distance:
                    atoms.append(atom_jdx)
                    continue
        atoms = remove_duplicates(atoms)
        return sorted(atoms)

    def write_gamess_basis_group(self):
        """ Returns a $BASIS group for GAMESS input.

            If there are multiple basis sets defined in the configuration
            file and there are multiple layers basis must be specified
            differently in the $DATA group. In that case this subroutine
            returns nothing

            This method will print warnings if there are discrepancies
            between number of layers and basis sets provided.
        """
        basis_sets = self._fragmentation.get_qm_basis()

        nbas = len(basis_sets)
        nlayers = self._nlayers

        if nbas != nlayers:
            if nbas > 1:
                raise ValueError("Error: FragIt [GAMESS-FMO] Number of basis sets ({}) and layers ({}) do not agree.".format(nbas, nlayers))

        if nbas > 1 and nlayers > 1:
            return ""

        if nbas == 0:
            print("Warning: FragIt [GAMESS-FMO] no basis set defined. Defaulting to 3-21G")
            basis_sets = ['3-21G']

        basis = basis_sets[0]

        return " $BASIS {0:s} $END\n".format(GAMESS_BASIS_GROUP[basis])

    def write_fmo_fmobnd_group(self):
        broken_bonds = self._fragmentation.get_explicitly_break_atom_pairs()
        if not isinstance(broken_bonds, list):
            raise TypeError
        if len(broken_bonds) == 0:
            return "\n"

        return " $FMOBND{0:s}\n $END\n".format(self.get_bond_group_data(broken_bonds))

    def get_bond_group_data(self, bonds) -> str:
        return "".join(["{0:s}".format(self.write_broken_bonds(bond)) for bond in bonds])

    def write_broken_bonds(self, bond_atoms: Tuple[int, int]) -> str:
        output_string = "\n{0:>10s}{1:10d}".format("-" + str(bond_atoms[0]), bond_atoms[1])
        dohop = self._fragmentation.do_fmohop_fragmentation()
        basis_sets = self._fragmentation.get_qm_basis()
        nbas = len(basis_sets)
        if dohop:
            for i in range(self._nlayers):
                if nbas == 1:
                    output_string += " {0:>s}".format(basis_sets[0])
                else:
                    output_string += " {0:>s}".format(basis_sets[i])
        return output_string

    def write_gamess_data_group(self):
        """ Returns the $DATA group

            This group usually only contains simple information as the
            basis set is provided through the $BASIS group. However, for
            multilayer calculations the basis set has to be specified here
            when multiple basis sets are requested.
        """
        return " $DATA\n%s\nc1\n%s $END\n" % (self._title, self.get_basis_set())

    def get_basis_set(self) -> str:
        return "".join([self.get_basis_set_for_atoms_in_layer(ilayer) for ilayer in range(1, self._nlayers + 1)])

    def get_basis_set_for_atoms_in_layer(self, layer: int) -> str:
        atom_numbers = remove_duplicates([atom.GetAtomicNum() for atom in self._fragmentation.get_ob_atoms()])
        atom_numbers.sort()
        atoms: List[str] = [Z2LABEL[atom_number] for atom_number in atom_numbers]
        return "".join([self.write_basis_set_for_atom(layer, atom) for atom in atoms])

    def write_basis_set_for_atom(self, layer: int, atom: str) -> str:
        """ Formats the basis set for a single atom

            the most common case is to just define the atom,
            the layer and the nuclear charge and the rest is
            handled by the $BASIS group.

            In multilayer cases where different basis sets are requested for each
            layer the basis set must be specified here.
        """
        s = "{0:s}-{1:d} {2:d}\n".format(atom, layer, LABEL2Z[atom])
        basis_sets = self._fragmentation.get_qm_basis()
        nbas = len(basis_sets)
        nlayers = self._nlayers

        if nbas > 1 and nlayers > 1:
            basis = basis_sets[layer - 1]
            try:
                basis_data = GAMESS_DATA_BASIS[basis]
            except KeyError:
                exit("Error: Basis set '{}' not defined. Aborting.".format(basis))
            try:
                _ = basis_data[atom]
            except KeyError:
                exit("Error: Basis set '{}' not defined for atom '{}'. Aborting.".format(basis, atom))

            s += "  {0:s}\n\n".format(basis_data[atom])

        return s

    def write_fmo_fmoxyz_group(self):
        s = " $FMOXYZ\n{0:s}\n $END\n"
        xyzstring = self.write_fmoxyz_atoms()[:-1]
        if len(xyzstring) == 0:
            s = ""
        return s.format(xyzstring)

    def write_fmoxyz_atoms(self) -> str:
        fragment_atoms = self._fragmentation.get_ob_atoms()
        fragments = self._fragmentation.get_fragments()
        fmo_fragments = []

        for i, fragment in enumerate(fragments):
            if i in self._water_fragments:
                continue
            fmo_fragments.extend(fragment)
        atoms = [fragment_atoms[i-1] for i in sorted(fmo_fragments)]
        return "".join([self.write_fmoxyz_atom(i, atom) for i, atom in enumerate(atoms, start=1)])

    def write_fmoxyz_atom(self, index: int, atom: openbabel.OBAtom) -> str:
        atom_label = "{0:7d}".format(index)
        if self._fragmentation.has_atom_names():
            names = self._fragmentation.get_atom_names()
            atom_label = "%7s" % (names[index - 1])
        return "%7s%7s%17f%13f%13f\n" % (atom_label,
                                         Z2LABEL[atom.GetAtomicNum()], atom.GetX(), atom.GetY(), atom.GetZ())

    def write_fmo_fmo_group(self) -> str:
        fmo_string = " $FMO\n%s%s%s\n%s\n%s\n%s\n%s\n%s\n%s\n $END\n"
        return fmo_string % (
            self.write_fmo_number_of_fragments(),
            self.write_fmo_defaults(),
            self.write_fmo_nlayers(),
            self.write_fmo_mplevl(),
            self.write_fmo_charges(),
            self.write_fmo_fragment_names(),
            self.write_fmo_indat(),
            self.write_fmo_layer(),
            self.write_fmo_active_fragment()
        )

    def write_fmo_number_of_fragments(self) -> str:
        """ Returns the number of FMO fragments in the calculations

            INFO: This number is modified if EFP waters are included
        """
        nfrag = len(self._fragmentation.get_fragments())
        nefpwaters = len(self._water_fragments)
        nfrag -= nefpwaters
        return "      NFRAG={0:d}\n".format(nfrag)

    def write_fmo_defaults(self) -> str:
        """ Returns FMO defaults """
        s = "      {0:s}\n".format("NBODY=2")
        nbonds_broken = self._fragmentation.get_num_broken_bonds()
        dohop = self._fragmentation.do_fmohop_fragmentation()
        if nbonds_broken > 0 and not dohop:
            s += "      {0:s}\n".format("RAFO(1)=1,1,1")
        s += "      {0:s}\n".format("RESDIM=2.0")
        s += "      {0:s}\n".format("RCORSD=2.0")
        return s

    def write_fmo_charges(self) -> str:
        return "      ICHARG(1)={0:s}".format(self._format_charges())

    def _format_charges(self) -> str:
        fragment_charges = self._fragmentation.get_fragment_charges()
        charges = []
        for i, charge in enumerate(fragment_charges):
            if i in self._water_fragments:
                continue
            charges.append(charge)
        output_2d_list = list_to_2d(charges, 10, "%3i")
        return list_2d_to_str(output_2d_list, ",", ",\n                 ")

    def write_fmo_fragment_names(self) -> str:
        return "      FRGNAM(1)={0:s}".format(self._format_fragment_names())

    def _format_fragment_names(self) -> str:
        """ Formats fragment names in FMO

            This group can technically be left out but at the moment it isn't so
        """
        names = self._fragmentation.get_fragment_names()
        fragnames = []
        fragment_index = 1
        for i, name in enumerate(names):
            if i in self._water_fragments:
                continue
            fragnames.append(" {0:5>s}{1:03d}".format(name, fragment_index))
            fragment_index += 1
        return list_2d_to_str(list_to_2d(fragnames, 5), ', ', ",\n                 ")

    def write_fmo_mplevl(self) -> str:
        return "      MPLEVL(1)=%s" % (
            list_2d_to_str(
                    list_to_2d([0 for _ in range(self._nlayers)], 10, '%i'),
                    ',', ",\n")
        )

    def write_fmo_indat(self) -> str:
        """ Returns the indices of fragments in an FMO calculation

            There is a sanity check going on here, that if the indices
            are not continous, i.e. 1, 2, 3, ... -> N the list will be
            rebuilt.
        """
        indat_base_string = "      INDAT(1)=0\n{0:s}"

        fragments = self._fragmentation.get_fragments()
        indices = []
        for i, fragment in enumerate(fragments):
            if i in self._water_fragments:
                continue
            indices.append(fragment)

        # we must check that the indices list is continous
        chklist = flatten(indices)
        chkval1 = sum(chklist)

        # the value of [sum_n=1^N n is 0.5*N*(N+1)] if it is continous.
        total = len(chklist)
        chkval2 = int(total*(total+1)/2)
        if chkval1 != chkval2:
            print("Warning: FragIt [GAMESS-FMO] Re-sequencing fragment indices.")
            # if we end up here, indices is not a continous series
            # which they must be for it to make sense in FMO.
            # So now we make a new index list with proper indices.
            new_indices = []
            i_start = 1
            for index in indices:
                i_end = i_start + len(index)
                new_indices.append(list(range(i_start, i_end)))
                i_start = i_end

            # copy new indices to indices list
            indices = new_indices[:]

        indat = "".join([list_of_ranges_to_string(list_to_ranges(frag)) for frag in indices])
        return indat_base_string.format(indat)

    def write_fmo_nlayers(self) -> str:
        return "      NLAYER=%i" % self._nlayers

    def write_fmo_active_fragment(self) -> str:
        if self._nlayers == 1:
            return ""
        if len(self._active_fragments) == 0:
            return ""
        return "      MODFD=1\n      IACTFG(1)=%s" % (list_2d_to_str(list_to_2d([i + 1 for i in self._active_fragments], 5, '%i'), ',', ',\n       '))

    def write_fmo_layer(self):
        if self._nlayers == 1 or self._central_fragment < 1:
            return ""
        layers = self._fragment_layers
        output_2d_list = list_to_2d(layers, 10, '%i')
        return "      LAYER(1)=%s" % list_2d_to_str(output_2d_list, ',', ",\n               ")
        
    def get_fragment_distances_vector(self, other_fragment: List[int]) -> List[float]:
        return [self.get_fragment_distance_to_fragment(fragment, other_fragment) for fragment in self._fragmentation.get_fragments()]

    def get_fragment_distance_to_fragment(self, fragment: List[int], other_fragment: List[int]) -> float:
        r_max = 1e30
        for atom_idx in fragment:
            for atom_jdx in other_fragment:
                r = self.get_distance_between_atoms(atom_idx, atom_jdx)
                if r < r_max:
                    r_max = r

        return r_max

    def get_distance_between_atoms(self, iat: int, jat: int) -> float:
        def get_ob_atom_vector(fragmentation, index: int) -> np.ndarray:
            atom = fragmentation.get_ob_atom(index)
            return np.array([atom.GetX(), atom.GetY(), atom.GetZ()])
        """ Returns the distance between two atoms """
        ivec = get_ob_atom_vector(self._fragmentation, iat)
        jvec = get_ob_atom_vector(self._fragmentation, jat)
        atom_vector = jvec - ivec
        r2: float = np.dot(atom_vector, atom_vector)
        return r2**0.5

    def get_layers_from_distances(self, distances: List[float]) -> List[int]:
        nfrags = len(self._fragmentation.get_fragments())
        fragment_layers = np.array([1 for i in range(nfrags)], dtype=int)
        layer = 2
        for distance in self._boundaries:
            selection = np.where(np.asarray(distances) < distance)
            for index in selection[0]:
                fragment_layers[index] = layer
            layer += 1
        return list(fragment_layers)

    def write_fmo_fmoefp_group(self) -> str:
        s = ""
        nefpwaters = len(self._water_fragments)
        if nefpwaters > 0:
            s += " $FMOEFP\n"
            s += "   NLEVEL=1\n"  # default to NLEVEL = 1 which minimizes EFP for each fragment and dimer.
            s += " $END\n"
        return s

    def write_gamess_efrag_group(self) -> str:
        s = ""
        nefpwaters = len(self._water_fragments)
        fragments = self._fragmentation.get_fragments()
        fragment_atoms = self._fragmentation.get_ob_atoms()
        if nefpwaters > 0:
            s += " $EFRAG\n"
            s += "COORD=CART\n"
            for water_fragment in self._water_fragments:
                s += "FRAGNAME=H2ORHF\n"
                for idx, atom_index in enumerate(fragments[water_fragment], start=1):
                    atom = fragment_atoms[atom_index-1]
                    symbol = Z2LABEL[atom.GetAtomicNum()]
                    s += "{0:s}{1:d}{2:17.6f}{3:13.6f}{4:13.6f}\n".format(symbol, idx, atom.GetX(), atom.GetY(), atom.GetZ())
            s += " $END\n"
        return s
