"""
Copyright (C) 2013-2017 Casper Steinmann
"""
from .fragit_exceptions import OBNotFoundException
try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")
import numpy

from .util import calculate_hydrogen_position
from .util import ravel2D, listDiff, Uniqify
from .util import getOBAtomVector


class QMMM(object):
    """ Actions related to performing additional refinement specific to QM/MM codes
    """

    def __init__(self, fragmentation, qmlist):
        if not isinstance(qmlist, list):
            raise TypeError("Expected a list of fragments to be included in the QM region")

        if len(qmlist) == 0:
            raise ValueError("List of fragments to be included in QM region is empty")

        self._fragmentation = fragmentation
        self._qmfrags = map(int, qmlist)
        self._fd = FragmentDistances(fragmentation)
        retract = lambda L: [l-1 for l in L]
        self._qmfrags = retract(self._qmfrags)
        self._qmfrags.sort()
        self._qmfrags.reverse()
        self._qm_charges = fragmentation.getFragmentCharges()

    def pop_qm_fragment(self):
        """ Remove the qm fragments from the fragmentation. Adds hydrogens to both the
            qm-fragment that is returned and to the neighbouring fragments if bonds were cut.
        """

        # at this point, fragments have been generated
        fragments = self._fragmentation.getFragments()

        # investigate which fragments should be included in the QM-region
        qmfrags = self._qmfrags[:]
        distance = self._fragmentation.getActiveAtomsDistance()
        qmfrags = Uniqify(self._add_fragments_to_QM(qmfrags))
        qmfrags.sort()
        qmfrags.reverse()

        qm_region_charge = 0

        fragments_for_qm_no_hydrogens = []

        # instead of removing the qm-fragment, we rather adopt a quite novel approach
        # in which we signal to the user of the API that the fragment is not to be used
        # further by setting its original atom numbers to -1
        for idx in qmfrags:
            qm_region_charge += self._qm_charges[idx]
            old_fragment = fragments.pop(idx) 
            fragments.insert(idx, [-1 for i in old_fragment])
            fragments_for_qm_no_hydrogens.insert(0,old_fragment[:])

        # for simplicity, let us just squash the qm fragments into one big fragment
        fragments_for_qm_no_hydrogens = ravel2D(fragments_for_qm_no_hydrogens)

        breaks = self._fragmentation.getExplicitlyBreakAtomPairs()
        if len(breaks) == 0:
            return (fragments_for_qm_no_hydrogens, qm_region_charge)

        # below here: add hydrogens to qm-fragment and to the rest of the capped structure
        fragment_for_qm = fragments_for_qm_no_hydrogens[:]

        # first, fix the QM-fragment, removing any bond-breaks that reside
        # inside (or bordering) the qm-fragment. if breaks are bordering,
        # add appropriate hydrogen atoms.
        lenqmfrags = len(fragments_for_qm_no_hydrogens)
        remove_breaks = []
        for ibreak, bbreak in enumerate(breaks):
            lendiff = len(listDiff(fragments_for_qm_no_hydrogens, list(bbreak)))
            difflen = lenqmfrags - lendiff

            # the break is not present in the qm-region, leave it
            if difflen == 0: continue
            # mark the ibreak'th item for removal. no hydrogens to be added
            if difflen == 2: remove_breaks.append(ibreak)
            # the break is bordering the qm-region and the mm-region
            if difflen == 1:

                for iibreak in bbreak:
                    # fix the qm-fragment first
                    if iibreak in fragments_for_qm_no_hydrogens:
                        new_atoms = self.satisfyValency(fragments_for_qm_no_hydrogens, iibreak, bbreak)
                        if len(new_atoms) > 0:
                            print("Info: FragIt adds", len(new_atoms), "atom(s) to the QM fragment.")
                            self._fragmentation._atom_names.extend(['  H '] * len(new_atoms))
                            fragment_for_qm.extend(new_atoms)

                    # then fix the fragments themselves
                    # INFO/WARNING: this is a lists of lists thing. BE CAREFULL
                    for ifragment, fragment in enumerate(fragments):
                        if iibreak in fragment:
                            new_atoms = self.satisfyValency(fragment, iibreak, bbreak)
                            print("Info: FragIt adds", len(new_atoms), "atom(s) to MM fragment", ifragment+1)
                            self._fragmentation._atom_names.extend(['  H '] * len(new_atoms))
                            fragments[ifragment].extend(new_atoms)

                # also mark the ibreak'th item for removal
                remove_breaks.append(ibreak)

        for ibreak in remove_breaks:
            self._fragmentation.popExplicitlyBreakAtomPairs(breaks[ibreak])

        return (fragment_for_qm, qm_region_charge)

    def satisfyValency(self, fragment, idx_heavy, bbreak):
        """ Satisfies the valency of a fragment that had a bond cut by
            adding a hydrogen atom along the same vector as that bond

            Arguments:
            fragment -- the fragment which valency needs to be satisfied
            iheavy -- index of the heavy atom on which the hydrogen is to be satisfied
            bbreak -- the indices of the atoms that were broken

            Returns a new fragment with all atoms in the correct place.
        """

        assert idx_heavy in bbreak

        # heavy is the heavy atom that wants a hydrogen
        # the following section of code figures out which
        # part of the bond we are investigating
        heavy = self._fragmentation.getOBAtom(idx_heavy)
        idx_light = 0
        if bbreak.index(idx_heavy) == 0: idx_light = 1
        light = self._fragmentation.getOBAtom(bbreak[idx_light])

        # now we add the new hydrogen atom
        try:
            rval = heavy.GetValence()
        except AttributeError: # OpenBabel3 has a different name
            rval = heavy.GetExplicitDegree()

        try:
            ival = heavy.GetImplicitValence()
        except AttributeError: # openbabel3 API changes
            ival = rval + heavy.GetImplicitHCount()

        new_atoms = []
        if ival != rval:
            if self._fragmentation.mol.AddHydrogens(heavy):
                for nbatom in openbabel.OBAtomAtomIter(heavy):
                    if nbatom.GetIdx() not in fragment:
                        x,y,z = calculate_hydrogen_position(heavy, light)
                        nbatom.SetVector(x,y,z)
                        new_atoms.append(nbatom.GetIdx())
        return new_atoms

    def _add_fragments_to_QM(self, qmfragments):
        """ Adds nearby fragments to the QM region if they fulfill
            certain requirements:

                hydrogen binds to QM region
                covalently bound to QM region

            Arguments:
            qmfragments -- list of QM fragments

            Returns:
            updated list of QM fragments with possible neighbours included.
        """
        qmfrags = qmfragments[:]
        for qmfragment in qmfragments:
            qmfrags.extend( self._fd.getHydrogenBoundFragments(qmfragment) )
            qmfrags.extend( self._fd.getCovalentlyBoundFragments(qmfragment) )
            qmfrags.extend( self._fd.getFragmentsWithinDistanceFrom(qmfragment) )
        return qmfrags


class FragmentDistances(object):
    def __init__(self, fragmentation):
        self._fragmentation = fragmentation
        self._fragments = fragmentation.getFragments()
        self._fragment_donors = [self._donors_from_fragment(i) for i in range(len(self._fragments))]
        self._fragment_acceptors = [self._acceptors_from_fragment(i) for i in range(len(self._fragments))]

        #print len(self._fragment_acceptors), len(self._fragments)

    def getHydrogenBoundFragments(self, idx):
        """ Obtains all hydrogen bound fragments to the idx'th fragment

            Arguments:
            idx -- the fragment index to which hydrogen bonds are to be found

            Returns:
            A list of fragments hydrogen bound to the idx'th fragment.
        """
        hb_fragments = []
        donors = []
        acceptors = []

        if self._fragmentation.doQMMMHydrogenBondDonors():
            donors = self._fragment_donors[idx]
        if self._fragmentation.doQMMMHydrogenBondAcceptors():
            acceptors = self._fragment_acceptors[idx]

        if len(donors) > 0 or len(acceptors) > 0:
            print("Info: FragIt will include all hydrogen bound molecules in the QM region.")

        # first we will find any donor (current fragment) -> acceptor (whole system) pairs
        isDonor = False
        for ifg, _acceptors in enumerate(self._fragment_acceptors):
            if ifg == idx: continue
            for A in _acceptors:
                if isDonor: break
                for H,D in donors:
                    isDonor = self._isHydrogenBond(D, H, A)
                    if isDonor:
                        hb_fragments.append( ifg )
                        break
            isDonor = False

        # then we will find acceptors (current fragment) -> donor (whole system)
        isAcceptor = False
        for ifg, _donors in enumerate(self._fragment_donors):
            if ifg == idx: continue
            for H,D in _donors:
                if isAcceptor: break
                for A in acceptors:
                    isAcceptor = self._isHydrogenBond(D, H, A)
                    if isAcceptor:
                        hb_fragments.append( ifg )
                        break
            isAcceptor = False

        return hb_fragments

    def getCovalentlyBoundFragments(self, idx):
        """ Returns list of fragments that are covalently connected to the idx fragment

            Arguments:
            idx -- the fragment index to use when looking for other covalently bound fragments

            Returns:
            A list of fragments covalently bound to the idx'th fragment.
        """
        atomids = self._fragments[idx]
        bonds = self._fragmentation.getExplicitlyBreakAtomPairs()

        other_fragment_atomids = []
        other_fragments = []

        if not self._fragmentation.doQMMMIncludeCovalent():
            return other_fragments

        print("Info: FragIt will include all fragments covalently bound to the QM region.")

        # let us find the nearby fragments that are covalently connected
        # currently, this only works with nearest neighbours
        for (l,r) in bonds:
            if l in atomids:
                other_fragment_atomids.append(r)
            if r in atomids:
                other_fragment_atomids.append(l)

        # lets find the associated fragments
        for ifg,atoms in enumerate(self._fragments):
            if len(atoms) > len(listDiff(atoms, other_fragment_atomids)):
                other_fragments.append(ifg)

        return other_fragments


    def _isHydrogenBond(self, D, H, A):
        """ angle > 110 is according to:
            http://pac.iupac.org/publications/pac/pdf/2011/pdf/8308x1637.pdf

            Arguments:
            D -- donor
            H -- hydrogen
            A -- acceptor
        """
        RAH = A.GetDistance(H)
        RAD = A.GetDistance(D)
        if RAH < self._fragmentation.getHBondDistanceMin() and RAD < self._fragmentation.getHBondDistanceMax():
            ANG = D.GetAngle(H,A) # in degrees
            if ANG > self._fragmentation.getHBondAngle():
                return True
        return False

    def _donors_from_fragment(self, idx):
        """ Searches a for a hydrogen-bond donor (aha a hydrogen). we consider only
            X --- H-N
            X --- H-O
            as potential donors. the atoms will be returned as a list of tuples (H, ?)
        """
        obatoms = [self._fragmentation.getOBAtom(i) for i in self._fragments[idx]]
        # lets locate all the hydrogen atoms connected to a nitrogen or oxygen
        donors = []
        for obatom in obatoms:
            if obatom.GetAtomicNum() == 1 and obatom.IsHbondDonorH():
                for otheratom in obatoms:
                    if obatom == otheratom: continue
                    if obatom.IsConnected(otheratom) and otheratom.IsHbondDonor():
                        donors.append((obatom, otheratom))
                        break

        return donors


    def _acceptors_from_fragment(self, idx):
        """ Searches a for a hydrogen-bond acceptor (aha a carboxyl oxygen). we consider only
            X --- O=C
            as potential donors. the atoms will be returned
        """
        obatoms = [self._fragmentation.getOBAtom(i) for i in self._fragments[idx]]
        acceptors = []
        for obatom in obatoms:
            if obatom.IsHbondAcceptor():
                acceptors.append(obatom)
        return acceptors


    def _is_doner(self, atom):
        """ Returns whether or not an atom is a hydrogen bond donor
        
            Arguments:
            atom -- OpenBabel atom
        """
        return atom.IsNitrogen() or atom.IsOxygen()


    def getFragmentsWithinDistanceFrom(self, idx):
        """ Returns all fragments that have at least an atom with a distance R
            of the fragment idx

            Arguments:
            idx -- the fragment index to use when looking for nearby fragments

            Returns:
            A list of fragments having at least an atom within a certain distance
        """
        other_fragments = []
        if not self._fragmentation.doQMMMIncludeAllWithin():
            return other_fragments

        print("Info: FragIt will include all fragments within R = {0:5.2f} angstrom in the QM region.".format(self._fragmentation.getQMMMIncludeAllWithinDistance()))

        vectors = self._get_vectors_from_fragment(idx)

        for i_frag, other_fragment in enumerate(self._fragments):
            if i_frag == idx:
                continue

            other_vectors = self._get_vectors_from_fragment(i_frag)
            R2 = self._get_min_distances2(vectors, other_vectors)
            R = numpy.sqrt(R2)
            if R < self._fragmentation.getQMMMIncludeAllWithinDistance():
                print("Info: FragIt includes fragment {0:5d} (R = {1:6.2f}) in the QM region.".format(i_frag+1, R))
                other_fragments.append(i_frag)

        return other_fragments


    def _get_vectors_from_fragment(self, idx):
        this_atoms = self._fragments[idx]
        this_vectors = []
        for atom_index in this_atoms:
            atom = self._fragmentation.getOBAtom(atom_index)
            this_vectors.append(getOBAtomVector(atom))

        return this_vectors

    def _get_min_distances2(self, d1s, d2s):
        """ Returns the minimum distance squared (hence the 2) between
            two sets of coordinates d1s and d2s

            Arguments:
            d1s -- the first set of coordinates
            d2s -- the second set of coordinates
        """
        d2s = numpy.array(d2s)
        R_min = 1.0e30
        for i, r1 in enumerate(d1s):
            dr = d2s-r1
            for r in dr:
                R_min = min(R_min, r.dot(r))

        return R_min
