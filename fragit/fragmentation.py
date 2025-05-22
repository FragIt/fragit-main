"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2023 Casper Steinmann
"""
import sys
from typing import List, Tuple

from fragit.fragit_exceptions import OBNotFoundException

try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")

from fragit.util import remove_duplicates, flatten, difference, uniqifyListOfLists, is_tuple_values_in_either_list, LABEL2Z, Z2LABEL
from fragit.config import FragItConfig, FragItDataBase


class Fragmentation(FragItConfig):

    def __init__(self, mol: openbabel.OBMol, **kwargs):
        """ Readies a fragmentation object for an OpenBabel molecule

            Arguments:
            mol -- molecule

            Keyword Arguments:
            defaults -- FragItConfig object defaults
        """
        conffile = kwargs.get("conffile", None)
        defaults = kwargs.pop("defaults", FragItDataBase)
        FragItConfig.__init__(self, defaults=defaults, filename=conffile, **kwargs)
        self.mol = mol
        self.pat = openbabel.OBSmartsPattern()
        self._atom_names: List[str] = []
        self._residue_names: List[str] = []
        self._fragment_names: List[str] = []
        self._fragment_charges: List[int] = []
        self._fragments: List[List[int]] = []
        self._backbone_atoms: List[int] = []
        self._water_molecules: List[int] = []
        self._mergeable_atoms: List[int] = []
        self._atoms: List[openbabel.OBAtom] = []
        self._fragment_charges_filename = None  # to read in fragment charges
        self._fix_atoms_and_charges()
        self._nbonds_broken: int = 0
        self._total_charge: float = 0.0

    def _fix_atoms_and_charges(self):
        """ Removes metal atoms to make charge calculation work

            OpenBabel has some issues with metal atoms when trying
            to obtain the total charge / fragment charges. The
            strategy is to remove the offending metal atoms from
            the molecule
        """
        _metalAtoms = self._remove_metal_atoms()

        self.formalCharges = [0.0 for _ in range(self.mol.NumAtoms())]

        # We can use either of the following:
        # None to just have zero charges
        # openbabel to guess the charges
        # read a file with fragment charges
        model = self.get_charge_model()
        if "read" in model:
            tmp, self._fragment_charges_filename = model.split()
        else:

            print("FragIt: The charge model in use is: {}".format(model))
            charge_model = openbabel.OBChargeModel.FindType(model)
            if charge_model is None:
                raise ValueError("Error: FragIt [FRAGMENTATION] could not use the charge model '{0:s}'".format(model))

            if charge_model.ComputeCharges(self.mol):
                self.formalCharges = list(charge_model.GetPartialCharges())
            else:
                print("Info: FragIt [FRAGMENTATION] fragment charges are not available.")

            # add back the metals, use the formal charges
            if len(_metalAtoms) > 0:
                print("Info: FragIt [FRAGMENTATION] appending metal atoms.")
                for atom in _metalAtoms:
                    if not self.mol.AddAtom(atom):
                        raise Exception("Error: FragIt [FRAGMENTATION] encountered an error when reinserting the metals.")
                    else:
                        self._atoms.append(atom)
                        self.formalCharges.append(atom.GetFormalCharge())

    def _remove_metal_atoms(self):
        _metalAtoms = []

        removing = True
        while removing:
            added = 0
            if self.mol.NumAtoms() == 0:
                break
            for i in range(1, self.mol.NumAtoms()+1):
                atom = self.mol.GetAtom(i)
# PX: not make Cl to Cl-
                if atom.GetAtomicNum() in [1, 6, 7, 8, 9, 12, 15, 16, 17]:
                    if atom not in self._atoms:
                        self._atoms.append(atom)
                    added += 1
                else:
                    print("Info: FragIt [FRAGMENTATION] Temporarily removing atom with Z={}".format(atom.GetAtomicNum()))
                    # the atoms are most-likely counter ions, so we give them an
                    # appropriate formal charge (of +/- 1)
                    atomic_charge = 0  # default
                    if atom.GetAtomicNum() in [11, 19]:  # Na+ and K+:
                        atomic_charge = 1
                    #elif atom.GetAtomicNum() in [9, 17]:  # F- and Cl-
                    #    atomic_charge = -1

                    new_atom = openbabel.OBAtom()
                    new_atom.Duplicate(atom)
                    new_atom.SetFormalCharge(atomic_charge)

                    _metalAtoms.append(new_atom)
                    self.mol.DeleteAtom(atom)
                    break

                if added == self.mol.NumAtoms():
                    removing = False
                    break

        if len(self.get_explicitly_break_atom_pairs()) > 0:
            print("\n WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING")
            print("    You have requested to manually fragment atom indices. However,")
            print("    because you _also_ have metal atoms in your input file, atom")
            print("    indices might have shifted around and you will likely see an")
            print("    error below about atoms not being connected with a bond.\n")
            print("    A possible fix is to move metal atoms to the end of your file.\n")
            #print(" WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING")

        return _metalAtoms

    def begin_fragmentation(self):
        """ Performs necessary actions before fragmentation can begin

            This includes identifying fragments which could
            potentially be of importance later. It also identifies
            which atoms are potentially excluded from participating
            in fragmentation
        """
        self.print_pre_fragmentation_information()
        self.identify_backbone_atoms()
        self.identify_water_molecules()
        self.identify_mergeable_atoms()
        self.identify_residues()
        if self.use_atom_names():
            self.name_atoms()
        self.set_protected_atoms()

    def identify_backbone_atoms(self):
        """ Identifies protein-backbone atoms """
        pattern = "N([*])C([H])C(=O)"
        self.pat.Init(pattern)
        self.pat.Match(self.mol)
        self._backbone_atoms = remove_duplicates(flatten(self.pat.GetUMapList()))

    def identify_water_molecules(self):
        """ Identifies all water molecules in the system """
        pattern = "[OH2]"
        # maybe O([H])[H] will give us all the atoms, actually
        self.pat.Init(pattern)
        self.pat.Match(self.mol)
        self._water_molecules = remove_duplicates(flatten(self.pat.GetUMapList()))

    def identify_mergeable_atoms(self):
        """ Identifies fragmentation points that are to be merged with neighbour fragments

            The fragmentation patterns are "greedy" in the sense that
            they match _all_ points in a protein. Some fragments (like glycine)
            are too small to make physically sense in a calculation and is thus
            merged into neighbours by _not_ fragmenting the points.

            Note: This has 
        """
        patterns = self.get_merge_patterns()
        for pattern in patterns:
            if len(patterns[pattern]) == 0:
                continue
            value = self._get_atoms_to_protect(patterns[pattern])
            value.sort()
            self._mergeable_atoms.extend(value)

    def do_fragment_merging(self):
        fragments_to_merge = self.get_fragments_to_merge()
        if len(fragments_to_merge) == 0:
            return
        fragments = self.get_fragments()
        fragments_to_merge.reverse()
        for fragment_id in fragments_to_merge:
            #PX : for edge case,
            # when the first fragment is glycine,do not merge.
            if fragment_id == 0:
                continue
            previous_fragment = fragment_id-1
            ifrag = fragments.pop(fragment_id)
            jfrag = fragments[previous_fragment].extend(ifrag)

        self._fragments = fragments[:]
        self._clean_merged_bonds()

    def do_fragment_combination(self):
        fragments_to_combine = self.get_combine_fragments()
        if len(fragments_to_combine) == 0:
            return
        fragments_to_combine.reverse()
        fragments = self.get_fragments()
        combined_fragment = flatten([fragments.pop(index - 1) for index in fragments_to_combine])
        combined_fragment.sort()
        fragments.append(combined_fragment)
        self._fragments = fragments[:]
        self._clean_merged_bonds()

    def get_fragments_to_merge(self):
        fragments = self.get_fragments()
        fragments_to_merge = []
        for i, fragment in enumerate(fragments):
            for sid in self._mergeable_atoms:
                if sid in fragment and i not in fragments_to_merge:
                    fragments_to_merge.append(i)
        return fragments_to_merge

    def do_fragmentation(self):
        """ Performs the actual fragmentation based on the
            actions performed in beginFragmentation
        """
        self.break_bonds()
        self.build_fragments()

    def finish_fragmentation(self):
        self.determine_fragment_charges()
        self.name_fragments()

    def get_fragments(self) -> List[List[int]]:
        """
            Note:
            Probably should be a copy (instead of a reference)
            but there are some dependencies that in the QM/MM
            module that fails such a check.
        """
        return self._fragments

    def get_fragment_names(self) -> List[str]:
        return self._fragment_names

    def get_ob_atoms(self) -> List[openbabel.OBAtom]:
        return self._atoms

    def set_protected_atoms(self):
        self.apply_smart_protect_patterns()

    def apply_smart_protect_patterns(self):
        patterns = self.get_protect_patterns()
        for protectpattern in patterns.keys():
            pattern = patterns[protectpattern]
            if len(pattern) == 0:
                continue
            self.add_explicitly_protected_atoms(self._get_atoms_to_protect(pattern))

    def _get_atoms_to_protect(self, pattern: str) -> List[int]:
        self.pat.Init(pattern)
        self.pat.Match(self.mol)
        return flatten(self.pat.GetUMapList())

    def identify_residues(self):
        if len(self._residue_names) > 0:
            return self._residue_names
        patterns = {"WATER": "[H]O[H]",
                    "NH3+": "N[H3]",
                    "AMINO": "C(=O)NC",
                    "SUGAR": "C1C(CO)OC(O)C(O)C1(O)"
                    }
        result = []
        for residue in patterns.keys():
            pattern = patterns[residue]
            self.pat.Init(pattern)
            self.pat.Match(self.mol)
            for atoms in self.pat.GetUMapList():
                result.append({residue: atoms})

        self._residue_names = result
        return result

    def is_bond_protected(self, atom_pair: Tuple[int, int]) -> bool:
        """ Returns whether a bond between two atoms is protected """
        protected_atoms = self.get_explicitly_protected_atoms()
        for bond_atom in atom_pair:
            if bond_atom in protected_atoms:
                return True
        return False

    def break_bonds(self):
        """ Searches for and breaks bonds in the molecule """
        self._search_fragmentation_atom_pairs()
        self._delete_ob_mol_bonds()

    def _delete_ob_mol_bonds(self):
        """ Deletes all bonds in the OBMol instance where fragmentation points are found """
        for atom_pair in self.get_explicitly_break_atom_pairs():
            if not self.is_valid_explicit_bond(atom_pair):
                continue

            if self.is_bond_protected(atom_pair):
                continue

            self._delete_ob_mol_bond(atom_pair)

    def _delete_ob_mol_bond(self, atom_pair: Tuple[int, int]):
        """ Deletes the bond in the OBMol instance between a pair of atoms

            with version 2.5 of the openbabel API there has been a change in
            the way bonds are guesstimated and this requires us to increase
            the implicit hydrogen count on each atom by the bond order of the
            newly broken bond manually.

            Arguments:
            pair -- a tuple (i, j) of the atoms with the bond between them
        """
        bond = self.mol.GetBond(atom_pair[0], atom_pair[1])
        bond_order = bond.GetBondOrder()
        try:
            atom1 = self.get_ob_atom(atom_pair[0])
            atom2 = self.get_ob_atom(atom_pair[1])
            impl_h_count1 = atom1.GetImplicitHCount()
            impl_h_count2 = atom2.GetImplicitHCount()
        except AttributeError:  # not openbabel 2.5
            pass
        else:
            atom1.SetImplicitHCount(impl_h_count1 + bond_order)
            atom2.SetImplicitHCount(impl_h_count2 + bond_order)

        self.mol.DeleteBond(bond)

    def _search_fragmentation_atom_pairs(self):
        """ Searches for places to break bonds via SMARTS

            When bonds are found, they are added to the list of
            bonds that must be exclusively broken (along with user
            defined places)
        """
        break_patterns = self.get_break_patterns()
        for BOND_TYPE in break_patterns.keys():
            pattern = break_patterns[BOND_TYPE]
            if self.get_verbose():
                print("\n      Searching for '{0:s}' bonds.".format(BOND_TYPE))
                print("      Pattern: '{0:s}'".format(pattern))

            if len(pattern) == 0:
                continue

            # find occurrences of the fragmentation pattern
            self.pat.Init(pattern)
            self.pat.Match(self.mol)
            matches = self.pat.GetUMapList()
            self._nbonds_broken += len(matches)
            if self.get_verbose():
                print("      Found {0:d} matching bonds.".format(len(matches)))
            for pair in matches:
                self.add_fragmentation_atom_pair(pair)

    def add_fragmentation_atom_pair(self, atom_pair: Tuple[int, int]):
        """ Adds an atom pair to the list of fragmentation points

            A test is made to see if the atom pair is protected from
            fragmentation
        """
        if self.is_bond_protected(atom_pair):
            if self.get_verbose():
                s_pair = "({0[0]:d},{0[1]:d})".format(atom_pair)
                print("      bond pair {0:>15s} will not be broken.".format(s_pair))
            return
        self.add_explicitly_break_atom_pairs([atom_pair])

    def is_valid_explicit_bond(self, pair: Tuple[int, int]) -> bool:
        """ Returns whether a pair of atoms contains a bond between them

            Arguments:
            pair -- a tuple (i, j) of the atoms with the bond between them

            Returns:
            True or False depending on if there is a bond between the atoms

            Raises:
            ValueError if either i == j or of there is no bond between the atoms
        """
        if pair[0] == pair[1]:
            raise ValueError("Error: FragIt [FRAGMENTATION] Fragment pair '{0:s}' must be two different atoms.".format(str(pair)))
        if self.mol.GetBond(pair[0], pair[1]) is None:
            raise ValueError("Error: FragIt [FRAGMENTATION] Fragment pair '{0:s}' must be connected with a bond.".format(str(pair)))
        return True

    def build_fragments(self):
        """ Builds all fragments from the fragmentation information

            The process is done in three steps:
              * obtain unique fragments
              * find fragments that could not be grouped
              * check that the fragments are sane
        """
        self.get_unique_fragments()
        self.find_remaining_fragments()
        self.do_sanity_check_on_fragments()

    def get_unique_fragments(self):
        result = list()
        for pair in self.get_explicitly_break_atom_pairs():
            result.append(self.get_atoms_in_same_fragment(pair[1], 0))
            result.append(self.get_atoms_in_same_fragment(pair[0], 0))
        result = uniqifyListOfLists(result)
        self._fragments = sorted(result)

    def find_remaining_fragments(self):
        remaining_atoms = difference(list(range(1, self.mol.NumAtoms() + 1)), flatten(self._fragments))
        while len(remaining_atoms) > 0:
            newfrag = self.get_atoms_in_same_fragment(remaining_atoms[0])
            remaining_atoms = difference(remaining_atoms, newfrag)
            self._fragments.append(newfrag)

    def do_sanity_check_on_fragments(self):
        """ Performs some level of sanity check on the generated fragments. """
        for i, fragment in enumerate(self._fragments, start=1):
            if len(fragment) > self.get_maximum_fragment_size():
                raise ValueError("Error: FragIt [FRAGMENTATION] Size of fragment {0:d} is {1:d}. Maximum allowed is {2:d}.".format(i, len(fragment), self.get_maximum_fragment_size()))
            if len(fragment) == 0:
                raise ValueError("Error: FragIt [FRAGMENTATION] Fragment {0:d} is empty.".format(i))

    def _clean_merged_bonds(self):
        broken_bonds = self.get_explicitly_break_atom_pairs()
        fragments = self.get_fragments()
        for bond in broken_bonds:
            for fragment in fragments:
                if bond[0] in fragment and bond[1] in fragment:
                    self.pop_explicitly_break_atom_pairs(bond)

    def do_fragment_grouping(self):
        if len(self._fragments) == 0:
            raise ValueError("You must fragment the molecule first.")
        group_size = self.get_fragment_group_count()
        new_fragments: List[List[int]] = []
        last_fragment = None
        group_count = 0
        temporary_fragment = []
        for fragment in self._fragments:
            group_count += 1
            if len(temporary_fragment) + len(fragment) > self.get_maximum_fragment_size():
                # max size is reached
                is_joinable = False
            else:
                is_joinable = self.is_fragment_joinable(fragment, last_fragment)
            if group_count < group_size and is_joinable:
                # joined to previous fragment.
                temporary_fragment += fragment
                last_fragment = fragment
                continue
            if is_joinable:
                # The group is now big enough
                temporary_fragment += fragment
                new_fragments.append(sorted(temporary_fragment))
                # reset
                last_fragment = None
                temporary_fragment = []
                group_count = 0
            else:
                if len(temporary_fragment) != 0:
                    new_fragments.append(sorted(temporary_fragment))
                last_fragment = fragment
                group_count = 1
                temporary_fragment = fragment
        if temporary_fragment:
            new_fragments.append(sorted(temporary_fragment))

        self._fragments = new_fragments

    def is_fragment_joinable(self, frag1: List[int], frag2: List[int]) -> bool:
        if frag2 is None:
            return True
        for p in self.get_explicitly_break_atom_pairs():
            if is_tuple_values_in_either_list(p, frag1, frag2):
                self.pop_explicitly_break_atom_pairs(p)
                return True
        return False

    def get_fragment_charges(self):
        """ Returns the fragment charges """
        return self._fragment_charges

    def determine_fragment_charges(self):
        """ Computes and stores fragment charges """
        if self._fragment_charges_filename is not None:
            self._fragment_charges = [0 for _ in self._fragments]
            print("Reading fragment charges from {}".format(self._fragment_charges_filename))
            with open(self._fragment_charges_filename, 'r') as charges_in:
                for line in charges_in:
                    tokens = list(map(int, line.split()))
                    self._fragment_charges[tokens[0]-1] = tokens[1]
        else:
            self._fragment_charges = [self.get_integer_fragment_charge(fragment) for fragment in self._fragments]
            self._total_charge = sum(self._fragment_charges)
            self.validate_total_charge()

    def get_integer_fragment_charge(self, fragment: List[int]) -> int:
        charge = self.get_sum_of_atomic_charges_in_fragment(fragment)
        return int(round(charge, 0))

    def get_sum_of_atomic_charges_in_fragment(self, fragment: List[int]) -> float:
        """ Computes the sum of atomic charges in a fragment

            Arguments:
            fragment -- the fragment indices used to compute a fragment charge

            Returns:
            charge -- the fragment charge

            Raises:
            IndexError if the charges have not been computed
        """
        charge = 0.0
        try:
            charge = sum([self.formalCharges[atom_idx-1] for atom_idx in fragment])
        except IndexError:
            print("Error: FragIt [FRAGMENTATION] found that fragment %s has invalid charges." % fragment)
        return charge

    def validate_total_charge(self):
        total_charge2 = sum(self.formalCharges)
        total_charge2 = int(round(total_charge2, 0))
        if self._total_charge != total_charge2:
            s = "Error: FragIt [FRAGMENTATION] cannot determine the charges of fragments correctly.\n"
            s += "       This is likely a problem in your structure, fragmentation\n."
            s += "       patterns or OpenBabel. Or a combination of the above.\n"
            s += "       Total charge = {0:f}. Sum of fragment charges = {1:d}".format(self._total_charge, total_charge2)
            print(s)
            raise ValueError("Error: FragIt [FRAGMENTATION] Total charge = {0:f}. Sum of fragment charges = {1:d}".format(self._total_charge, total_charge2))

    def get_atoms_in_same_fragment(self, a1: int, a2: int = 0):
        """ The heart of FragIt.

            This method determines fragments in a system.
            It works by finding all possible atoms between
            two end-points in the molecular graph by calling
            the FindChildren method of OBMol.

            Note: There is some dark art going on here since
                  a2 is apparently always zero.
        """
        if not isinstance(a1, int) or not isinstance(a2, int):
            raise ValueError
        if a2 != 0:
            raise ValueError
        tmp = openbabel.vectorInt()
        self.mol.FindChildren(tmp, a2, a1)

        fragment = [value for value in tmp] + [a1]
        return sorted(fragment)

    def get_ob_atom(self, atom_index: int) -> openbabel.OBAtom:
        """ Retrieves an OpenBabel atom from the fragmentation object

            :param int atom_index: The atom index of the atom to retrieve.
            :returns: an openbabel atom
        """
        if not isinstance(atom_index, int):
            raise ValueError
        if atom_index < 1 or atom_index > self.mol.NumAtoms():
            raise IndexError("Error: FragIt[FRAGMENTATION] Index '{0:d}' out of range [{1:d},{2:d}]".format(atom_index, 1, self.mol.NumAtoms()))
        return self.mol.GetAtom(atom_index)

    def name_fragments(self) -> List[str]:
        names = list()
        for fragment in self.get_fragments():
            names.append(self.name_fragment(fragment))
        self._fragment_names = names
        return names

    def name_fragment(self, atoms: List[int]) -> str:
        charge_lbls = ["", "+", "-"]

        if len(atoms) == 0:
            raise ValueError("Error: FragIt [FRAGMENTATION] Cannot name empty fragments. Aborting.")

        if len(atoms) == 1:
            atom = self.mol.GetAtom(atoms[0])
            charge_lbl = charge_lbls[atom.GetFormalCharge()]
            element = Z2LABEL[atom.GetAtomicNum()]
            return "{0:s}{1:s}".format(element, charge_lbl)
        else:
            for residue in openbabel.OBResidueIter(self.mol):
                for atom in openbabel.OBResidueAtomIter(residue):
                    if atom.GetIdx() in atoms:
                        return str(residue.GetName())
        return "None"

    def get_backbone_atoms(self):
        return self._backbone_atoms

    def get_water_molecules(self):
        return self._water_molecules

    def name_atoms(self):
        """Attempts to name atoms """
        atoms_no_name = list(range(0, self.mol.NumAtoms()))

        # first try to name atoms according to biological
        # function, i.e. from a PDB file.

        # OpenBabel has problems with some elements, so we need to work around
        # those here.
        for residue in openbabel.OBResidueIter(self.mol):
            if residue.GetNumAtoms() > 0:
                for atom in openbabel.OBResidueAtomIter(residue):
                    atoms_no_name.remove(atom.GetId())
                    self._atom_names.append(residue.GetAtomID(atom))
            else:
                pass

        # if there are items left in "atoms_no_name" try to name them
        # print atoms_no_name
        if len(atoms_no_name) > 0:
            if self.get_verbose():
                print("Info: FragIt [FRAGMENTATION] will now try to name remaining {0:3d} atoms".format(len(atoms_no_name)))

            for i, index in enumerate(atoms_no_name):
                atom = self.mol.GetAtom(index+1)
                self._atom_names.append(atom.GetType())
                if self.get_verbose():
                    print("   atom {0:4d} is forced to have name '{1:s}'".format(index, atom.GetType()))

    def get_atom_names(self):
        return self._atom_names

    def has_atom_names(self):
        return len(self._atom_names) > 0

    def print_pre_fragmentation_information(self):
        if not self.get_verbose():
            return

        print("Info: FragIt [FRAGMENTATION] will break the following bonds defined in config file:")

        if len(self.get_explicitly_break_atom_pairs()) > 0:
            print("\n      Bonds explicitly defined between pairs of atoms:")
        for pair in self.get_explicitly_break_atom_pairs():
            print("      {0[0]:>4d} <-> {0[1]:d}".format(pair))

        print("")
        for pair in self.get_explicitly_break_atom_pairs():
            try:
                self.is_valid_explicit_bond(pair)
            except ValueError as e:
                print(e)
                sys.exit()  # we bomb here for now 

    def get_num_broken_bonds(self):
        return self._nbonds_broken
