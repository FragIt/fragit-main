"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2017 Casper Steinmann
"""
import os
import sys

try:
    import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")
import numpy

from .util import *
from .config import FragItConfig, FragItDataBase


class Fragmentation(FragItConfig):

    def __init__(self, mol, **kwargs):
        conffile = kwargs.get('conffile', None)
        defaults = kwargs.pop('defaults', FragItDataBase)
        FragItConfig.__init__(self, defaults=defaults, filename=conffile, **kwargs)
        self.mol     = mol
        self.obc     = openbabel.OBConversion()
        self.pat     = openbabel.OBSmartsPattern()
        self._atom_names = []
        self._residue_names = []
        self._fragment_names = []
        self._fragment_charges = []
        self._fragments = []
        self._backbone_atoms = []
        self._mergeable_atoms = []
        self._atoms = []
        self._fixAtomsAndCharges()
        self._elements = openbabel.OBElementTable()
        self._nbonds_broken = 0


    def _removeMetalAtoms(self):
        _metalAtoms = []

        removing = True
        while removing:
            added = 0
            if self.mol.NumAtoms() == 0:
                break
            for i in range(1, self.mol.NumAtoms()+1):
                atom = self.mol.GetAtom(i)
                if atom.GetAtomicNum() in [1,6,7,8,12,15,16]:
                    if atom not in self._atoms: self._atoms.append(atom)
                    added += 1
                else:
                    print("Info: FragIt [FRAGMENTATION] Temporarily removing atom with Z={}".format(atom.GetAtomicNum()))
                    # the atoms are most-likely counter ions, so we give them an
                    # appropriate formal charge (of +/- 1)
                    atomic_charge = 0 # default
                    if atom.GetAtomicNum() in [11, 19]: # Na+ and K+:
                        atomic_charge = 1
                    elif atom.GetAtomicNum() in [9, 17]: # F- and Cl-
                        atomic_charge = -1

                    new_atom = openbabel.OBAtom()
                    new_atom.Duplicate(atom)
                    new_atom.SetFormalCharge(atomic_charge)

                    _metalAtoms.append(new_atom)
                    self.mol.DeleteAtom(atom)
                    break

                if added == self.mol.NumAtoms():
                    removing = False
                    break

        if len(self.getExplicitlyBreakAtomPairs()) > 0:
            print("\n WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING")
            print("    You have requested to manually fragment atom indices. However,")
            print("    because you _also_ have metal atoms in your input file, atom")
            print("    indices might have shifted around and you will likely see an")
            print("    error below about atoms not being connected with a bond.\n")
            print("    A possible fix is to move metal atoms to the end of your file.\n")
            #print(" WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING")

        return _metalAtoms


    def _fixAtomsAndCharges(self):
        """Removes unwanted atoms to make the charge calculation
           work. Be sure to re-insert the atoms once it is done
        """
        _metalAtoms = self._removeMetalAtoms()

        # now lets do the charges (without the metals)
        model = self.getChargeModel()
        charge_model = openbabel.OBChargeModel.FindType(model)
        if charge_model is None:
            raise ValueError("Error: FragIt [FRAGMENTATION] could not use the charge model '{0:s}'".format(model))

        self.formalCharges = [0.0 for i in range(self.mol.NumAtoms())]
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


    def beginFragmentation(self):
        """ Performs the nescessary actions before fragmentation

            This includes identifying fragments which could
            potentially be of importance. It also identifies
            which atoms are potentially excluded from parti-
            cipating in fragmentation
        """
        self.printPreFragmentationInformation()
        self.identifyBackboneAtoms()
        self.identifyWaterMolecules()
        self.identifyMergeableAtoms()
        self.identifyResidues()
        if self.useAtomNames():
            self.nameAtoms()
        self.setProtectedAtoms()


    def identifyBackboneAtoms(self):
        pattern = "N([*])C([H])C(=O)"
        self.pat.Init(pattern)
        self.pat.Match(self.mol)
        self._backbone_atoms = Uniqify(self._listMatches(self.pat.GetUMapList()))


    def identifyWaterMolecules(self):
        """ Identifies all water molecules in the system """
        pattern = "[OH2]"
        # maybe O([H])[H] will give us all the atoms, actually
        self.pat.Init(pattern)
        self.pat.Match(self.mol)
        self._watermolecules = Uniqify(self._listMatches(self.pat.GetUMapList()))


    def identifyMergeableAtoms(self):
        patterns = self.getMergePatterns()
        for pattern in patterns:
            if len(patterns[pattern]) == 0: continue
            value = self._getAtomsToProtect(patterns[pattern])
            value.sort()
            self._mergeable_atoms.extend(value)


    def doFragmentMerging(self):
        fragments_to_merge = self.getFragmentsToMerge()
        if len(fragments_to_merge) == 0: return
        fragments = self.getFragments()
        fragments_to_merge.reverse()
        for fragment_id in fragments_to_merge:
            previous_fragment = fragment_id-1
            ifrag = fragments.pop(fragment_id)
            jfrag = fragments[previous_fragment].extend(ifrag)

        self._fragments = fragments[:]
        self._CleanMergedBonds()


    def doFragmentCombination(self):
        fragments_to_combine = self.getCombineFragments()
        if len(fragments_to_combine) == 0: return
        fragments_to_combine.reverse()
        fragments = self.getFragments()
        combined_fragment = ravel2D([fragments.pop(id-1) for id in fragments_to_combine])
        combined_fragment.sort()
        fragments.append(combined_fragment)
        self._fragments = fragments[:]
        self._CleanMergedBonds()


    def getFragmentsToMerge(self):
        fragments = self.getFragments()
        fragments_to_merge = []
        for i,fragment in enumerate(fragments):
            for sid in self._mergeable_atoms:
                if sid in fragment and i not in fragments_to_merge:
                    fragments_to_merge.append(i)
        return fragments_to_merge


    def doFragmentation(self):
        """ Performas the actual fragmentation based on the
            actions performed in beginFragmentation
        """
        self.breakBonds()
        self.determineFragments()


    def finishFragmentation(self):
        self.determineFragmentCharges()
        self.nameFragments()


    def getFragments(self):
        return self._fragments


    def getFragmentNames(self):
        return self._fragment_names


    def getAtoms(self):
        return self._atoms


    def setActiveFragments(self, active_fragments):
        if not isinstance(active_fragments, list): raise TypeError
        self.active_fragments = active_fragments


    def setProtectedAtoms(self):
        self.applySmartProtectPatterns()


    def clearProtectionPatterns(self):
        self.protected_atoms = list()


    def applySmartProtectPatterns(self):
        patterns = self.getProtectPatterns()
        for protectpattern in patterns.keys():
            pattern = patterns[protectpattern]
            if len(pattern) == 0: continue
            self.addExplicitlyProtectedAtoms(self._getAtomsToProtect(pattern))


    def _getAtomsToProtect(self,pattern):
        self.pat.Init(pattern)
        self.pat.Match(self.mol)
        return self._listMatches(self.pat.GetUMapList())


    def _listMatches(self,matches):
        results = []
        for match in matches:
            results.extend(match)
        return results


    def identifyResidues(self):
        if (len(self._residue_names) > 0):
            return self._residue_names
        patterns =     {
                "WATER"    : "[H]O[H]",
                "NH3+"  :"N[H3]",
                "AMINO"    : "C(=O)NC",
                "SUGAR" : "C1C(CO)OC(O)C(O)C1(O)"
                }
        result = []
        for residue in patterns.keys():
            pattern = patterns[residue]
            self.pat.Init(pattern)
            self.pat.Match( self.mol )
            for atoms in self.pat.GetUMapList():
                result.append({residue : atoms})

        self._residue_names = result
        return result


    def isBondProtected(self,bond_pair):
        protected_atoms = self.getExplicitlyProtectedAtoms()
        for bond_atom in bond_pair:
            if bond_atom in protected_atoms:
                return True
        return False


    def breakBonds(self):
        self._BreakPatternBonds()
        self._DeleteOBMolBonds()


    def _DeleteOBMolBonds(self):
        for pair in self.getExplicitlyBreakAtomPairs():
            if not self.isValidExplicitBond(pair):
                continue

            if self.isBondProtected(pair):
                continue

            self._DeleteOBMolBond(pair)


    def _DeleteOBMolBond(self,pair):
        bond = self.mol.GetBond(pair[0],pair[1])
        self.mol.DeleteBond(bond)


    def _BreakPatternBonds(self):
        breakPatterns = self.getBreakPatterns()
        for bondType in breakPatterns.keys():
            pattern = breakPatterns[bondType]
            if self.getVerbose():
                print("\n      Searching for '{0:s}' bonds.".format(bondType))
                print("      Pattern: '{0:s}'".format(pattern))

            if len(pattern) == 0:
                continue

            # find occurrences of the fragmentation pattern
            self.pat.Init(pattern)
            self.pat.Match(self.mol)
            matches = self.pat.GetUMapList()
            self._nbonds_broken += len(matches)
            if self.getVerbose():
                print("      Found {0:d} matching bonds.".format(len(matches)))
            for p in matches:
                self.realBondBreaker(bondType, p)


    def realBondBreaker(self, bondtype, bond_pair):
        if self.isBondProtected(bond_pair):
            if self.getVerbose():
                s_pair = "({0[0]:d},{0[1]:d})".format(bond_pair)
                print("      bond pair {0:>15s} will not be broken.".format(s_pair))
            return
        self.addExplicitlyBreakAtomPairs([bond_pair])


    def isValidExplicitBond(self, pair):
        if pair[0] == pair[1]:
            raise ValueError("Error: FragIt [FRAGMENTATION] Fragment pair '{0:s}' must be two different atoms.".format(str(pair)))
        if self.mol.GetBond(pair[0],pair[1]) is None:
            raise ValueError("Error: FragIt [FRAGMENTATION] Fragment pair '{0:s}' must be connected with a bond.".format(str(pair)))
        return True


    def determineFragments(self):
        """ Determines the individual fragments by gathering up atoms
            that belong to each fragment.
        """
        self.getUniqueFragments()
        self.findRemainingFragments()
        self.doSanityCheckOnFragments()


    def getUniqueFragments(self):
        result = list()
        for pair in self.getExplicitlyBreakAtomPairs():
            result.append(self.getAtomsInSameFragment(pair[1], 0))
            result.append(self.getAtomsInSameFragment(pair[0], 0))
        result = uniqifyListOfLists(result)
        self._fragments = sorted(result)


    def findRemainingFragments(self):
        remainingAtoms = listDiff(range(1, self.mol.NumAtoms() + 1), ravel2D(self._fragments))
        while len(remainingAtoms) > 0:
            newfrag = self.getAtomsInSameFragment(remainingAtoms[0])
            remainingAtoms = listDiff(remainingAtoms, newfrag)
            self._fragments.append(newfrag)


    def doSanityCheckOnFragments(self):
        """ Performs some level of sanity check on the generated fragments. """
        for i, fragment in enumerate(self._fragments, start=1):
            if len(fragment) > self.getMaximumFragmentSize():
                raise ValueError("Error: FragIt [FRAGMENTATION] Size of fragment {0:d} is {1:d}. Maximum allowed is {2:d}.".format(i, len(fragment), self.getMaximumFragmentSize()))
            if len(fragment) == 0:
                raise ValueError("Error: FragIt [FRAGMENTATION] Fragment {0:d} is empty.".format(i))


    def _CleanMergedBonds(self):
        broken_bonds = self.getExplicitlyBreakAtomPairs()
        fragments = self.getFragments()
        for bond in broken_bonds:
            for fragment in fragments:
                if bond[0] in fragment and bond[1] in fragment:
                    self.popExplicitlyBreakAtomPairs(bond)


    def doFragmentGrouping(self):
        if len(self._fragments) == 0: raise ValueError("You must fragment the molecule first.")
        groupSize = self.getFragmentGroupCount()
        newfragments = []
        lastFragment = None
        grpcount = 0
        tmp = []
        for fragment in self._fragments:
            grpcount += 1
            if len(tmp) + len(fragment) > self.getMaximumFragmentSize():
                #max size is reached
                isJoinable = False
            else:
                isJoinable     = self.isJoinable(fragment, lastFragment)
            if (grpcount < groupSize and isJoinable):
                #joined to previous fragment.
                tmp += fragment
                lastFragment =  fragment
                continue
            if (isJoinable):
                #The group is now big enough
                tmp += fragment
                newfragments.append(sorted(tmp))
                #reset
                lastFragment = None
                tmp = []
                grpcount = 0
            else:
                if (len(tmp) != 0):
                    newfragments.append(sorted(tmp))
                lastFragment     = fragment
                grpcount     = 1
                tmp         = fragment
                #raise ValueError("Maximum fragment size too small, unable to continue")
        if (tmp != []):
            newfragments.append(sorted(tmp))

        self._fragments = newfragments


    def isJoinable(self,frag1,frag2):
        if (frag2 == None): return True
        for p in self.getExplicitlyBreakAtomPairs():
            if tupleValuesInEitherList(p,frag1,frag2):
                self.popExplicitlyBreakAtomPairs(p)
                return True
        return False


    def getFragmentCharges(self):
        return self._fragment_charges


    def determineFragmentCharges(self):
        self._fragment_charges = [self.getIntegerFragmentCharge(fragment) for fragment in self._fragments]
        self.total_charge = sum(self._fragment_charges)
        self.validateTotalCharge()


    def getIntegerFragmentCharge(self, fragment):
        charge = self.getSumOfAtomicChargesInFragment(fragment)
        return int(round(charge,0))


    def getSumOfAtomicChargesInFragment(self,fragment):
        charge = 0.0
        try:
            charge = sum([self.formalCharges[atom_idx-1] for atom_idx in fragment])
        except IndexError:
            print("Error: FragIt [FRAGMENTATION] found that fragment %s has invalid charges." % (fragment))
        return charge


    def validateTotalCharge(self):
        total_charge2 = sum(self.formalCharges)
        total_charge2 = int(round(total_charge2,0))
        if (self.total_charge != total_charge2):
            s  = "Error: FragIt [FRAGMENTATION] cannot determine the charges of fragments correctly.\n"
            s += "       This is likely a problem in your structure, fragmentation\n."
            s += "       patterns or OpenBabel. Or a combination of the above.\n"
            s += "       Total charge = {0:d}. Sum of fragment charges = {1:d}".format(self.total_charge, total_charge2)
            print(s)
            raise ValueError("Error: FragIt [FRAGMENTATION] Total charge = {0:d}. Sum of fragment charges = {1:d}".format(self.total_charge, total_charge2))


    def getAtomsInSameFragment(self, a1, a2 = 0):
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
        fragment = sorted(toList(tmp) + [a1])
        return fragment


    def getOBAtom(self, atom_index):
        if not isinstance(atom_index, int): raise ValueError
        if atom_index < 1 or atom_index > self.mol.NumAtoms():
            raise ValueError("Error: FragIt[FRAGMENTATION] Index '{0:d}' out of range [{1:d},{2:d}]".format(atom_index,1,self.mol.NumAtoms()))
        return self.mol.GetAtom(atom_index)


    def nameFragments(self):
        names = list()
        for fragment in self.getFragments():
            names.append( self.tryNameFragment(fragment) )
        self._fragment_names = names
        return names


    def tryNameFragment(self, atoms):
        matched_atoms     = 0
        frag_name     = False
        residues = self.identifyResidues()
        charge_lbls = ["", "+", "-"]

        if len(atoms) == 0:
            raise ValueError("Error: FragIt [FRAGMENTATION] Cannot name empty fragments. Aborting.")

        if len(atoms) == 1:
            atom = self.mol.GetAtom(atoms[0])
            charge_lbl = charge_lbls[atom.GetFormalCharge()]
            element = self._elements.GetSymbol(atom.GetAtomicNum())
            return "{0:s}{1:s}".format(element, charge_lbl)
        else:
            for residue in openbabel.OBResidueIter( self.mol ):
                for atom in openbabel.OBResidueAtomIter( residue ):
                    if atom.GetIdx() in atoms:
                        return residue.GetName()
        return "None"


    def getBackboneAtoms(self):
        return self._backbone_atoms


    def getWaterMolecules(self):
        return self._watermolecules

    def nameAtoms(self):
        """attempt to name atoms """
        has_residues = False
        residue_has_atoms = False
        atoms_no_name = list(range(0, self.mol.NumAtoms()))

        # first try to name atoms according to biological
        # function, i.e. from a PDB file.

        # OpenBabel has problems with some elements so we need to work around
        # those here.
        for residue in openbabel.OBResidueIter(self.mol):
            has_residues = True
            if residue.GetNumAtoms() > 0:
                residue_has_atoms = True
                for atom in openbabel.OBResidueAtomIter( residue ):
                    atoms_no_name.remove(atom.GetId())
                    self._atom_names.append( residue.GetAtomID( atom ) )
            else:
                pass

        # if there are items left in "atoms_no_name" try to name them
        #print atoms_no_name
        if len(atoms_no_name) > 0:
            if self.getVerbose():
                print("Info: FragIt [FRAGMENTATION] will now try to name remaining {0:3d} atoms".format(len(atoms_no_name)))

            for i, id in enumerate(atoms_no_name):
                atom = self.mol.GetAtom(id+1)
                self._atom_names.append(atom.GetType())
                if self.getVerbose():
                    print("   atom {0:4d} is forced to have name '{1:s}'".format(id, atom.GetType()))


    def getAtomNames(self):
        return self._atom_names


    def hasAtomNames(self):
        return len(self._atom_names) > 0


    def printFragment(self, i_frag, atoms):
        print("{0:d}".format(len(atoms)))
        print("Fragment {0:4d} has {1:4d} atoms".format(i_frag+1, len(atoms)))
        for iat, atom in enumerate(atoms):
            print("{0:4d}{1:5d}".format(iat, atom))


    def printPreFragmentationInformation(self):
        if not self.getVerbose():
            return

        print("Info: FragIt [FRAGMENTATION] will break the following bonds defined in config file:")

        if len(self.getExplicitlyBreakAtomPairs()) > 0:
            print("\n      Bonds explicitly defined between pairs of atoms:")
        for pair in self.getExplicitlyBreakAtomPairs():
            print("      {0[0]:>4d} <-> {0[1]:d}".format(pair))

        print("")
        for pair in self.getExplicitlyBreakAtomPairs():
            try:
                self.isValidExplicitBond(pair)
            except ValueError as e:
                print(e)
                sys.exit()  # we bomb here for now 


    def getNumBrokenBonds(self):
        return self._nbonds_broken
