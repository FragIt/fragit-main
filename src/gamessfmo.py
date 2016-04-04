"""
**********************************************************************
gamessfmo.py - the GAMESS FMO output writer

Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2016 Casper Steinmann

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
from numpy import sqrt, dot, where, array
from writer import Standard
from util import WriteStringToFile

from util import file_extension,is_list,listTo2D,join2D,is_int
from util import listToRanges,listOfRangesToString,Uniqify,ravel2D
from util import deepLength

class GamessFMO(Standard):
    def __init__(self, fragmentation):
        Standard.__init__(self,fragmentation)

        # initialize layered stuff
    def setup(self):
        self._setupLayeredInformation()
        self._setupActiveFragmentsInformation()
        self._validateMultiLayerInformation()
        if self._do_pymol: self._dump_pymol()
        if self._do_jmol: self._dump_jmol()

    def _setupLayeredInformation(self):
        self._fragment_layers = self._getFragmentLayersFromFragment()

    def _getFragmentLayersFromFragment(self):
        fragments = self._fragmentation.getFragments()
        if self._central_fragment == 0: return array([1 for i in fragments])
        other_fragment = fragments[self._central_fragment-1]
        distances = self._getFragmentDistancesVector(other_fragment)
        layers = self._getLayersFromDistances(distances)
        return layers

    def _dump_pymol(self):
        from pymol import PymolTemplate
        pt = PymolTemplate(self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _dump_jmol(self):
        from jmol import JmolTemplate
        pt = JmolTemplate(self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _setTemplateData(self, template):
        template.setFragmentsData(self._fragmentation.getFragments())
        template.setBufferData(self._fragment_layers)
        template.setActiveData(self._active_atoms)
        template.setBackboneData(self._fragmentation.getBackboneAtoms())
        template.setPairData(self._fragmentation.getExplicitlyBreakAtomPairs())

    def _writeTemplateFile(self, template):
        template.override()
        template.write()

    def _setupActiveFragmentsInformation(self):
        active_atoms = self._getActiveAtomsFromFragments()
        if self._active_atoms_distance > 0.0:
            active_atoms = self._getActiveAtomsFromDistance()
            if self._verbose:
                print("Info: FragIt [GAMESS-FMO] found {0:d} atoms which should be active".format(len(active_atoms)))

        self._active_atoms = active_atoms[:]
        atoms = self._active_atoms[:]
        frags = self._fragmentation.getFragments()
        fragment_layers = self._fragment_layers[:]
        if self._central_fragment == 0:
            self._active_frags = []
            return

        # 1) If there are active atoms, we must find the associated fragments
        # 2) The associated fragments are labelled active
        # 3) The active fragments have their atoms made flexible
        # we now have region A
        if len(self._active_atoms) > 0:
            active_frags = [] #self._active_fragments[:]
            active_frags.append(self._central_fragment -1) # central must also be active
            for atom in atoms:
                ifrg = self._getFragmentFromAtom(atom)
                active_frags.extend([ifrg])
            active_frags = Uniqify(active_frags)
            active_frags = sorted(active_frags)
            self._active_fragments = active_frags[:]

            # promote active fragments to layer 2
            for active_fragment_id in active_frags:
                self._fragment_layers[active_fragment_id] = 2
    
            # add active fragment atoms to list of active atoms
            fragments = self._fragmentation.getFragments()
            for frag in active_frags:
                atoms.extend(fragments[frag])
            atoms = Uniqify(atoms)
            atoms = sorted(atoms)
            if self._verbose and len(atoms) != len(self._active_atoms):
                print("Info: FragIt [GAMESS-FMO] active region is now {0:d} atoms ({1:d} fragments) ".format(len(active_atoms), len(active_frags)))

        # Optionally freeze backbone atoms in the active region
        if self._freeze_backbone:
            for item in self._fragmentation.getBackboneAtoms():
                if item in atoms:
                    atoms.remove(item)
                    continue

        atoms = Uniqify(atoms)
        atoms = sorted(atoms)
        self._active_atoms = atoms[:]

    def _getFragmentFromAtom(self, atom):
        for i, fragment in enumerate(self._fragmentation.getFragments()):
            if atom in fragment:
                return i

    def _validateMultiLayerInformation(self):
        """ Validates multilayer (and FD) information and attemps
            to fix any issues should be be present.

            This method makes sure that the buffer region, b, around the active region A is
            does not have close contacts between A and F.
        """
        active_atoms = self._active_atoms[:]
        active_fragments = self._active_fragments[:]
        fragment_layers = self._fragment_layers[:]
        fragments = self._fragmentation.getFragments()
        for atom in active_atoms:  # this is for specifying something lame on the command line
            if not self._isAtomInActiveLayer(atom):
                #print "removing atom %i" % atom
                self._active_atoms.remove(atom)

        # Here we make sure that fragments in A are not physically close
        # to fragments in F.
        # We do this by extending B (which includes A) with a buffer region, b,
        # with the buffer distance that a user wants.
        # Technically we promote the fragments in b to layer 2
        if(len(active_atoms) > 0 and len(active_fragments) > 0):
            for active_fragment_index in active_fragments:
                active_fragment = fragments[active_fragment_index]
                distance_vector = self._getFragmentDistancesVector(active_fragment)
                selection = where( (0.1 < distance_vector) & (distance_vector < self._buffer_maximum_distance))
                self._fragment_layers[selection] = 2

            if self._verbose:
                frags = list()
                for i,k in enumerate(self._fragment_layers):
                    if k == 2:
                        frags.append(i)

                atms = list()
                for frag in frags:
                    atms.extend(fragments[frag])

                print("Info: FragIt [GAMESS-FMO] region B is {0:d} atoms ({1:d} fragments)".format(len(atms), len(frags)))

        self._active_atoms = active_atoms[:]

    def _isAtomInActiveLayer(self, atom):
        frags = self._fragmentation.getFragments()
        active_fragments = self._active_fragments[:]
        fragment_layers = self._fragment_layers[:]
        for i,fragment in enumerate(frags):
            if (atom in fragment) and (fragment_layers[i] == 2):
                return True
        return False

    def writeFile(self, filename):
        outStringTemplate = "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n"
        outString = outStringTemplate % (    self.SYSTEMgroup(),self.GDDIgroup(),self.SCFgroup(),
                            self.CONTRLgroup(),self.BASISgroup(),self.FMOPRPgroup(),
                            self.FMOgroup(),self.FMOBNDgroup(),self.DATAgroup(),
                            self.FMOXYZgroup())
        WriteStringToFile(filename, outString)

    def SYSTEMgroup(self):
        return " $SYSTEM MWORDS=125 $END"

    def SCFgroup(self):
        return " $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END"

    def GDDIgroup(self):
        return " $GDDI NGROUP=1 $END"

    def FMOPRPgroup(self):
        return self._get_FMOPRP_basestring() % self._calculateOrbitalGuess()

    def _get_FMOPRP_basestring(self):
        return " $FMOPRP NPRINT=9 NGUESS=%i $END"

    def _calculateOrbitalGuess(self):
        nguess = 2 # project orbitals out of huckel guess
        if self._nlayers > 1 and len(self._active_atoms) != 0:
            nguess += 128 # this is needed for multilayer formulation, but not FD
        return nguess

    def CONTRLgroup(self):
        localize = " LOCAL=BOYS"
        if len(self._fragmentation.getExplicitlyBreakAtomPairs()) == 0:
            localize = ""
        base = " $CONTRL NPRINT=-5 ISPHER=1%s\n         RUNTYP=%s\n $END"
        statpt = " $STATPT OPTTOL=5.0e-4 NSTEP=2000\n%s\n $END"
        if(len(self._active_fragments) == 0 and self._active_atoms_distance <= 0.0):
            return base % (localize, "ENERGY")
        else:
            active_string = self._getActiveAtomsString(self._active_atoms)
            base_final = base % (localize, "OPTIMIZE")
            statpt_final = statpt % active_string
            final = "%s\n%s" % (base_final, statpt_final)
            return final

    def _getActiveAtomsFromFragments(self):
        atoms = []
        fragments = self._fragmentation.getFragments()
        for idx in self._active_fragments:
            atoms.extend( fragments[idx-1] )
        return sorted(atoms)

    def _getActiveAtomsFromDistance(self):
        atoms = []
        central_atoms = self._fragmentation.getFragments()[self._central_fragment-1]
        atoms.extend(central_atoms)
        all_atoms = range(1,len(self._fragmentation.getAtoms())+1)
        for atom_idx in central_atoms:
            for atom_jdx in all_atoms:
                if atom_jdx in central_atoms: continue
                if atom_jdx in atoms: continue
                R = self._getDistanceBetweenAtoms(atom_idx,atom_jdx)
                if R < self._active_atoms_distance:
                    atoms.append(atom_jdx)
                    continue
        atoms = Uniqify(atoms)
        return sorted(atoms)

    def _getActiveAtomsString(self, atoms):
        active = "      IACTAT(1)=%s"
        #active_atoms ="".join([listOfRangesToString(listToRanges(frag)) for frag in atoms])
        active_atoms = listOfRangesToString(listToRanges(self._active_atoms),
                        maxlength=40,
                        line_format="%5s",
                        item_format="%s,",
                        tuple_format="%s,%s,",
                        terminator_format=None)
        return active % active_atoms

    def BASISgroup(self):
        default_basis = " $BASIS GBASIS=N21 NGAUSS=3 $END"
        #basis = self._fragmentation.getBasisSet()
        #if len(basis) == 0:
        #    basis = default_basis
        return default_basis

    def FMOBNDgroup(self):
        broken_bonds = self._fragmentation.getExplicitlyBreakAtomPairs()
        if not is_list(broken_bonds):
            raise TypeError
        if len(broken_bonds) == 0:
            return "\n"

        return " $FMOBND%s\n $END" % self._getBondGroupData(broken_bonds)

    def _getBondGroupData(self, bonds):
        return "".join(["%s" % self._formatBrokenBond(bond) for bond in bonds])

    def _formatBrokenBond(self,bond_atoms):
        return "\n%10s%10s" % ("-"+str(bond_atoms[0]),bond_atoms[1])

    def DATAgroup(self):
        return " $DATA\n%s\nc1\n%s $END" % (self._title, self._getBasisSet())

    def _getBasisSet(self):
        return "".join([self._getBasisAtoms(ilayer) for ilayer in range(1, self._nlayers+1)])

    def _getBasisAtoms(self, ilayer):
        atom_numbers = Uniqify([atom.GetAtomicNum() for atom in self._fragmentation.getAtoms()])
        atom_numbers.sort()
        atoms = [self._elements.GetSymbol(atom_number) for atom_number in atom_numbers]
        return "".join([self._formatSingleAtomBasis(ilayer,atom) for atom in atoms])

    def _formatSingleAtomBasis(self,ilayer,atom):
        return "%s-%i %i\n" % (atom,ilayer,self._elements.GetAtomicNum(atom))

    ## the $FMOXYZ group
    def FMOXYZgroup(self):
        return " $FMOXYZ\n%s\n $END" % self._formatAtoms()[:-1]

    def _formatAtoms(self):
        return "".join([self._formatSingleAtom(i+1,atom) for i,atom in enumerate(self._fragmentation.getAtoms())])

    def _formatSingleAtom(self, index, atom):
        strindex = "%7i" % (index)
        if self._fragmentation.hasAtomNames():
             names = self._fragmentation.getAtomNames()
             strindex = "%7s" % (names[index-1])
        return "%7s%7s%17f%13f%13f\n" % (strindex,
             self._elements.GetSymbol(atom.GetAtomicNum()), atom.GetX(), atom.GetY(), atom.GetZ())

    ## The $FMO group
    def FMOgroup(self):
        fmoString = " $FMO\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n $END"
        fmo = fmoString % (    self._getFMODefaults(),self._getFMONLayer(),self._getFMONFrag(),
                  self._getFMOICharg(),self._getFMOFrgnam(),self._getFMOMplevl(),
                  self._getFMOIndat(),self._getFMOLayer(),self._getFMOActfg())
        return fmo

    def _getFMODefaults(self):
        return "      %s\n      %s\n      %s" % ("NBODY=2","RAFO(1)=1,1,1","RESDIM=2.0 RCORSD=2.0")

    def _getFMONFrag(self):
        return "      NFRAG=%i" % len(self._fragmentation.getFragments())

    def _getFMOICharg(self):
        return "      ICHARG(1)=%s" % (self._formatCharges())

    def _formatCharges(self):
        list2D = listTo2D(self._fragmentation.getFragmentCharges(), 10, '%3i')
        return join2D(list2D, ',', ",\n                 ")

    def _getFMOFrgnam(self):
        return "      FRGNAM(1)=%s" % (self._formatFragmentNames())

    def _formatFragmentNames(self):
        fragnames = list()
        names = self._fragmentation.getFragmentNames()
        for i in range(1, len(names) + 1):
            s = "%5s%03i" % (names[i-1], i)
            fragnames.append(s)
        return join2D(listTo2D(fragnames, 5), ', ', ",\n                 ")

    def _getFMOMplevl(self):
        return "      MPLEVL(1)=%s" % (join2D(listTo2D([0 for i in range(self._nlayers)],10,'%i'),',',",\n"))

    def _getFMOIndat(self):
        frags = self._fragmentation.getFragments()
        indat_string = "      INDAT(1)=0\n%s"
        indat ="".join([listOfRangesToString(listToRanges(frag)) for frag in frags])
        return indat_string % indat

    def _getFMONLayer(self):
        return "      NLAYER=%i" % self._nlayers

    def _getFMOActfg(self):
        if self._nlayers == 1: return ""
        if len(self._active_fragments) == 0: return ""
        return "      MODFD=1\n      IACTFG(1)=%s" % (join2D(listTo2D([i+1 for i in self._active_fragments],5,'%i'),',',',\n       '))

    def _getFMOLayer(self):
        if self._nlayers == 1 or self._central_fragment < 1: return ""
        layers = self._fragment_layers
        list2D = listTo2D(layers, 10, '%i')
        return "      LAYER(1)=%s" % join2D(list2D, ',', ",\n               ")
        
    def _getFragmentDistancesVector(self, other_fragment):
        return array([self._getFragmentDistanceToFragment(fragment, other_fragment) for fragment in self._fragmentation.getFragments()])

    def _getFragmentDistanceToFragment(self, fragment, other_fragment):
        r_max = 1e30
        for atom_idx in fragment:
            for atom_jdx in other_fragment:
                r = self._getDistanceBetweenAtoms(atom_idx, atom_jdx)
                if r < r_max:
                    r_max = r

        return r_max

    def _getDistanceBetweenAtoms(self, iat, jat):
        ivec = self._getAtomVector(iat)
        jvec = self._getAtomVector(jat)
        atom_vector = jvec - ivec
        return sqrt( dot( atom_vector, atom_vector ) )

    def _getAtomVector(self, index):
        atom = self._fragmentation.getOBAtom(index)
        return array([atom.GetX(), atom.GetY(), atom.GetZ()])

    def _getLayersFromDistances(self, distances):
        nfrags = len(self._fragmentation.getFragments())
        fragment_layers = array([1 for i in range(nfrags)])
        layer = 2
        for distance in self._boundaries:
            fragment_layers[where(distances<distance)] = layer
            layer +=1
        return fragment_layers

