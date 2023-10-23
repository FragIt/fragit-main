"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2017 Casper Steinmann
"""
import os

from numpy import sqrt, dot, where, array

from .writer import Standard
from .util import WriteStringToFile
from .util import file_extension,listTo2D,join2D
from .util import listToRanges,listOfRangesToString,Uniqify,ravel2D
from .util import deepLength, Z2LABEL, LABEL2Z

# $BASIS set data depending on basis set
GAMESS_BASIS_GROUP = dict()
GAMESS_BASIS_GROUP['STO-3G']      =  "GBASIS=STO NGAUSS=3"
GAMESS_BASIS_GROUP['3-21G']       =  "GBASIS=N21 NGAUSS=3"
GAMESS_BASIS_GROUP['6-31G']       =  "GBASIS=N31 NGAUSS=6"
GAMESS_BASIS_GROUP['6-31G*']      =  "GBASIS=N31 NGAUSS=6 NDFUNC=1"
GAMESS_BASIS_GROUP['6-31G(d)']    =  "GBASIS=N31 NGAUSS=6 NDFUNC=1"
GAMESS_BASIS_GROUP['6-31+G*']     =  "GBASIS=N31 NGAUSS=6 NDFUNC=1 DIFFSP=.T."
GAMESS_BASIS_GROUP['6-31+G(d)']   =  "GBASIS=N31 NGAUSS=6 NDFUNC=1 DIFFSP=.T."
GAMESS_BASIS_GROUP['cc-pVDZ']     =  "GBASIS=CCD"
GAMESS_BASIS_GROUP['cc-pVTZ']     =  "GBASIS=CCT"
GAMESS_BASIS_GROUP['aug-cc-pVDZ'] =  "GBASIS=ACCD"
GAMESS_BASIS_GROUP['aug-cc-pVTZ'] =  "GBASIS=ACCT"

# basis set data for atoms in $DATA group
GAMESS_DATA_BASIS = dict()
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
        Standard.__init__(self,fragmentation, directories)

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
        self._water_fragments = []
        layers = array([1 for i in fragments])

        if self._fragmentation.doFMOEFPWaters():
            #1 - figure out which fragments these oxygens belong to.
            #    If this is not a multilayer calculation, ALL
            #    waters are promoted (demoted?) to EFP waters
            water_oxygen = self._fragmentation.getWaterMolecules()
            self._water_fragments = sorted([self._getFragmentFromAtom(i) for i in water_oxygen])


        if self._central_fragment == 0:
            return layers
        other_fragment = fragments[self._central_fragment-1]
        distances = self._getFragmentDistancesVector(other_fragment)
        layers = self._getLayersFromDistances(distances)

        if self._fragmentation.doFMOEFPWaters():
            #3 - if it is a multilayer calculation, all waters in
            #    the layer specified by getFMOEFPWatersFromLayer
            #    are promoted (demoted?) to EFP waters and multilayer
            #    is switched off.
            #    self._fragment_layers[active_fragment_id] = 2
            water_fragments = []
            for idx, water_fragment in enumerate(sorted(self._water_fragments)):
                if layers[water_fragment] == self._fragmentation.getFMOEFPWatersFromLayer():
                    water_fragments.append(water_fragment)
            self._water_fragments = water_fragments[:]
            layers = array([1 for i in fragments])
            self._nlayers = 1
        return layers

    def _dump_pymol(self):
        from .pymol import PymolTemplate
        pt = PymolTemplate(self._directories, self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _dump_jmol(self):
        from .jmol import JmolTemplate
        pt = JmolTemplate(self._directories, self._input_filename, self._output_filename)
        self._setTemplateData(pt)
        self._writeTemplateFile(pt)

    def _setTemplateData(self, template):
        template.setFragmentsData(self._fragmentation.getFragments())
        template.setBufferData(self._fragment_layers)
        template.setActiveData(self._active_atoms)
        template.setBackboneData(self._fragmentation.getBackboneAtoms())
        template.setPairData(self._fragmentation.getExplicitlyBreakAtomPairs())
        template.setFragmentCharges(self._fragmentation.getFragmentCharges())

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
        outStringTemplate = "%s%s%s%s%s%s%s%s%s%s%s%s%s"
        outString = outStringTemplate % (self.SYSTEMgroup(),
                                         self.GDDIgroup(),
                                         self.SCFgroup(),
                                         self.CONTRLgroup(),
                                         self.BASISgroup(),
                                         self.FMOPRPgroup(),
                                         self.FMOgroup(),
                                         self.FMOBNDgroup(), 
                                         self.DATAgroup(),
                                         self.FMOHYBgroup(), 
                                         self.FMOXYZgroup(),
                                         self.FMOEFPgroup(),
                                         self.EFRAGgroup())
        WriteStringToFile(filename, outString)

    def FMOHYBgroup(self):
        """ Generates the FMOHYB input group """
        s = ""
        nbonds_broken = self._fragmentation.getNumBrokenBonds()
        dohop = self._fragmentation.doFMOHOPFragmentation()
        if nbonds_broken > 0 and dohop:
            s += " $FMOHYB\n"
            basis_sets = self._fragmentation.getQMBasis()
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

    def SYSTEMgroup(self):
        return " $SYSTEM MWORDS=125 $END\n"

    def SCFgroup(self):
        return " $SCF CONV=1E-7 DIRSCF=.T. NPUNCH=0 DIIS=.F. SOSCF=.T. $END\n"

    def GDDIgroup(self):
        return " $GDDI NGROUP=1 $END\n"

    def FMOPRPgroup(self):
        return self._get_FMOPRP_basestring() % self._calculateOrbitalGuess()

    def _get_FMOPRP_basestring(self):
        return " $FMOPRP NPRINT=9 NGUESS=%i $END\n"

    def _calculateOrbitalGuess(self):
        nguess = 2 # project orbitals out of huckel guess
        if self._nlayers > 1 and len(self._active_atoms) != 0:
            nguess += 128 # this is needed for multilayer formulation, but not FD
        return nguess

    def CONTRLgroup(self):
        """ Returns the $CONTRL group

            if a geometry optimzation is requested this method
            also returns the $STATPT group
        """

        nbonds_broken = self._fragmentation.getNumBrokenBonds()
        dohop = self._fragmentation.doFMOHOPFragmentation()
        localize = ""
        if nbonds_broken > 0 and not dohop:
            localize = " LOCAL=BOYS"
        base = " $CONTRL NPRINT=-5 ISPHER=1%s\n         RUNTYP=%s\n $END\n"
        statpt = " $STATPT OPTTOL=5.0e-4 NSTEP=2000\n%s\n $END\n"
        if(len(self._active_fragments) == 0 and self._active_atoms_distance <= 0.0):
            return base % (localize, "ENERGY")
        else:
            base_final = base % (localize, "OPTIMIZE")
            active_string = self._getActiveAtomsString(self._active_atoms)
            statpt_final = statpt % active_string
            return "{0:s}{1:s}".format(base_final, statpt_final)

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
        """ Returns a $BASIS group for GAMESS input.

            If there are multiple basis sets defined in the configuration
            file and there are multiple layers basis must be specified
            differently in the $DATA group. In that case this subroutine
            returns nothing

            This method will print warnings if there are discrepancies
            between number of layers and basis sets provided.
        """
        basis_sets = self._fragmentation.getQMBasis()

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

    def FMOBNDgroup(self):
        broken_bonds = self._fragmentation.getExplicitlyBreakAtomPairs()
        if not isinstance(broken_bonds, list):
            raise TypeError
        if len(broken_bonds) == 0:
            return "\n"

        return " $FMOBND{0:s}\n $END\n".format(self._getBondGroupData(broken_bonds))

    def _getBondGroupData(self, bonds):
        return "".join(["{0:s}".format(self._formatBrokenBond(bond)) for bond in bonds])

    def _formatBrokenBond(self,bond_atoms):
        s = "\n{0:>10s}{1:10d}".format("-"+str(bond_atoms[0]),bond_atoms[1])
        dohop = self._fragmentation.doFMOHOPFragmentation()
        basis_sets = self._fragmentation.getQMBasis()
        nbas = len(basis_sets)
        if dohop:
           for i in range(self._nlayers):
               if nbas == 1:
                   s += " {0:>s}".format(basis_sets[0])
               else:
                   s += " {0:>s}".format(basis_sets[i])
        return s

    def DATAgroup(self):
        """ Returns the $DATA group

            This group usually only contains simple information as the
            basis set is provided through the $BASIS group. However, for
            multilayer calculations the basis set has to be specified here
            when multiple basis sets are requested.
        """
        return " $DATA\n%s\nc1\n%s $END\n" % (self._title, self._getBasisSet())

    def _getBasisSet(self):
        return "".join([self._getBasisAtoms(ilayer) for ilayer in range(1, self._nlayers+1)])

    def _getBasisAtoms(self, ilayer):
        atom_numbers = Uniqify([atom.GetAtomicNum() for atom in self._fragmentation.getAtoms()])
        atom_numbers.sort()
        atoms = [Z2LABEL[atom_number] for atom_number in atom_numbers]
        return "".join([self._formatSingleFMOXYZAtomBasis(ilayer,atom) for atom in atoms])

    def _formatSingleFMOXYZAtomBasis(self,ilayer,atom):
        """ Formats the basis set for a single atom

            the most common case is to just define the atom,
            the layer and the nuclear charge and the rest is
            handled by the $BASIS group.

            In multilayer cases where different basis sets are requested for each
            layer the basis set must be specified here.
        """
        s = "{0:s}-{1:d} {2:d}\n".format(atom,ilayer,LABEL2Z[atom])
        basis_sets = self._fragmentation.getQMBasis()
        nbas = len(basis_sets)
        nlayers = self._nlayers

        if nbas > 1 and nlayers > 1:
            basis = basis_sets[ilayer-1]
            try:
                basis_data = GAMESS_DATA_BASIS[basis]
            except KeyError:
                exit("Error: Basis set '{}' not defined. Aborting.".format(basis))
            try:
                atom_basis_data = basis_data[atom]
            except KeyError:
                exit("Error: Basis set '{}' not defined for atom '{}'. Aborting.".format(basis, atom))

            s += "  {0:s}\n\n".format(basis_data[atom])

        return s

    ## the $FMOXYZ group
    def FMOXYZgroup(self):
        s = " $FMOXYZ\n{0:s}\n $END\n"
        xyzstring = self._formatFMOXYZAtoms()[:-1]
        if len(xyzstring) == 0:
            s = ""
        return s.format(xyzstring)

    def _formatFMOXYZAtoms(self):
        fragment_atoms = self._fragmentation.getAtoms()
        fragments = self._fragmentation.getFragments()
        fmo_fragments = []

        for i, fragment in enumerate(fragments):
            if i in self._water_fragments:
                continue
            fmo_fragments.extend(fragment)
        atoms = [fragment_atoms[i-1] for i in sorted(fmo_fragments)]
        return "".join([self._formatSingleFMOXYZAtom(i,atom) for i, atom in enumerate(atoms, start=1)])

    def _formatSingleFMOXYZAtom(self, index, atom):
        strindex = "{0:7d}".format(index)
        if self._fragmentation.hasAtomNames():
             names = self._fragmentation.getAtomNames()
             strindex = "%7s" % (names[index-1])
        return "%7s%7s%17f%13f%13f\n" % (strindex,
             Z2LABEL[atom.GetAtomicNum()], atom.GetX(), atom.GetY(), atom.GetZ())

    ## The $FMO group
    def FMOgroup(self):
        fmoString = " $FMO\n%s%s%s\n%s\n%s\n%s\n%s\n%s\n%s\n $END\n"
        fmo = fmoString % ( self._getFMONFrag(), self._getFMODefaults(),self._getFMONLayer(),self._getFMOMplevl(),
                  self._getFMOICharg(), self._getFMOFrgnam(),
                  self._getFMOIndat(), self._getFMOLayer(), self._getFMOActfg())
        return fmo

    def _getFMONFrag(self):
        """ Returns the number of FMO fragments in the calculations

            INFO: This number is modified if EFP waters are included
        """
        nfrag = len(self._fragmentation.getFragments())
        nefpwaters = len(self._water_fragments)
        nfrag -= nefpwaters
        return "      NFRAG={0:d}\n".format(nfrag)

    def _getFMODefaults(self):
        """ Returns FMO defaults """
        s = "      {0:s}\n".format("NBODY=2")
        nbonds_broken = self._fragmentation.getNumBrokenBonds()
        dohop = self._fragmentation.doFMOHOPFragmentation()
        if nbonds_broken > 0 and not dohop:
            s += "      {0:s}\n".format("RAFO(1)=1,1,1")
        s += "      {0:s}\n".format("RESDIM=2.0")
        s += "      {0:s}\n".format("RCORSD=2.0")
        return s

    def _getFMOICharg(self):
        return "      ICHARG(1)={0:s}".format(self._formatCharges())

    def _formatCharges(self):
        fragment_charges = self._fragmentation.getFragmentCharges()
        charges = []
        for i, charge in enumerate(fragment_charges):
            if i in self._water_fragments:
                continue
            charges.append(charge)
        list2D = listTo2D(charges, 10, '%3i')
        return join2D(list2D, ',', ",\n                 ")

    def _getFMOFrgnam(self):
        return "      FRGNAM(1)={0:s}".format(self._formatFragmentNames())

    def _formatFragmentNames(self):
        """ Formats fragment names in FMO

            This group can technically be left out but at the moment it isn't so
        """
        names = self._fragmentation.getFragmentNames()
        fragnames = []
        fragment_index = 1
        for i, name in enumerate(names):
            if i in self._water_fragments:
                continue
            fragnames.append(" {0:5>s}{1:03d}".format(name, fragment_index))
            fragment_index += 1
        return join2D(listTo2D(fragnames, 5), ', ', ",\n                 ")

    def _getFMOMplevl(self):
        return "      MPLEVL(1)=%s" % (join2D(listTo2D([0 for i in range(self._nlayers)],10,'%i'),',',",\n"))

    def _getFMOIndat(self):
        """ Returns the indices of fragments in an FMO calculation

            There is a sanity check going on here, that if the indices
            are not continous, i.e. 1, 2, 3, ... -> N the list will be
            rebuilt.
        """
        indat_base_string = "      INDAT(1)=0\n{0:s}"

        fragments = self._fragmentation.getFragments()
        indices = []
        for i, fragment in enumerate(fragments):
            if i in self._water_fragments:
                continue
            indices.append(fragment)

        # we must check that the indices list is continous
        chklist = ravel2D(indices)
        chkval1 = sum(chklist)

        # the value of [sum_n=1^N n is 0.5*N*(N+1)] if it is continous.
        N = len(chklist)
        chkval2 = int(N*(N+1)/2)
        if chkval1 != chkval2:
            print("Warning: FragIt [GAMESS-FMO] Re-sequencing fragment indices.")
            # if we end up here, indices is not a continous series
            # which is must be for it to make sense in FMO.
            # So now we make a new index list with proper indices.
            new_indices = []
            i_start = 1
            for index in indices:
                i_end = i_start + len(index)
                new_indices.append(list(range(i_start, i_end)))
                i_start = i_end

            # copy new indices to indices list
            indices = new_indices[:]

        indat ="".join([listOfRangesToString(listToRanges(frag)) for frag in indices])
        return indat_base_string.format(indat)

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
        """ Returns the distance between two atoms """
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

    def FMOEFPgroup(self):
        s = ""
        nefpwaters = len(self._water_fragments)
        if nefpwaters > 0:
            s += " $FMOEFP\n"
            s += "   NLEVEL=1\n" # default to NLEVEL = 1 which minimizes EFP for each fragment and dimer.
            s += " $END\n"
        return s

    def EFRAGgroup(self):
        s = ""
        nefpwaters = len(self._water_fragments)
        fragments = self._fragmentation.getFragments()
        fragment_atoms = self._fragmentation.getAtoms()
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
