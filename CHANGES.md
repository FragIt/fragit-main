FragIt v1.9.0 Release Notes
===========================

This is a much-needed major revision to remove bugs
and squash problems related to the old-age of FragIt

Updates since  v1.8.0

  * FragIt can now read fragment charges from an
    external file. This is specified in the input
    configuration files in

    [fragmentation]
    chargemodel = read test.dat

    where the MMFF94 keyword is replaced with
    read <filename>. Each line in the supplied
    file gives two values. The first value is the
    fragment index and the second value is the
    fragment charge. The following example shows
    that the charge of the first fragment  is -1
    and the fifth fragment has charge +1. Fragments
    not listed are assumed to be neutral.
```text
    1 -1
    5  1
```
   
  * Includes some changes necessary for openbabel
    cooperation with especially version 3.

  * Fixed a few errors in the jmol script and template.

  * Switch to GitHub actions for build system testing.

  * Updated the GAMESS to HMO script to python3

Developer Changes

 * Made substantial changes to align FragIt with how
python packages should be organized

 * code updates and improvements all around

FragIt v1.8.0 Release Notes
===========================

This major update brings several quality of life
improvements in FragIt.

Updates since  v1.7.2

  * Now full (automated) testing of multiple python
    versions (2.7, 3.5 and 3.6 as we speak) coupled
    with multiple (master and 2.4.1 release) versions
    of openbabel

  * includes preliminary fixes for the upcoming openbabel
    2.5 release (no idea when it is due) specifically
    regarding OBElementsTable and implicit valence info
    which was relocated in the openbabel API.

  * It is now possible to show the charges of fragments
    in PyMol. The command

    ColorFragments charge

    in the PyMol terminal will color the fragments with
    red and blue for negative and positive charges of the
    fragments, respectively.

  * FragIt now dumps it's input file in standard output
    so one can actually see the options being used. This
    is mainly to be inline with other programs and to
    make sure users understand what options were set.

  * FragIt now correctly warns people using explicit
    fragmentation points in combination with metal atoms
    that atoms are shifted around. This procedure has
    not been completely automated.

  * A lot of interal API documentation stuff happened.


FragIt v1.7.2 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.7.1

  * Mixing metals/ions with explicit bond cutting
    can (and is still) problematic. For now, print
    a warning with possible fixes so that the user
    not only sees a "bond not found" error. (#26)
                   -- Thanks to Morten S Nørby
                      for reporting this issue.


FragIt v1.7.1 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.7.0

  * Fixes a python3 issue with the pymol scripts and
    template system.

  * PyMol scripts are now named .py instead of .pymol
    to signify that they are different from the .pml
    scripts. PyMol 1.8.6 enforces this rule.

  * the calculation of fragment labels has changed to
    a more efficient algorithm. Per default, fragment
    labels are hidden.

  In addition, some behind the scene updates were made
  to the template system


FragIt v1.7.0 Release Notes
===========================

This major update brings several changes to default keywords
and to the command line interface for FragIt.

Updates since v1.6.1

  * FragIt is now compatible with python 3.5 after
    minor modifications to the source code mostly
    regarding imports and print statements.

  * Travis-CI builds are failing for python 3.5 for
    some compilation/install of openbabel problem
    so testing is done manually for now. This is
    also why it is compatible but not entirely
    supported for now.

  * Templates are now located appropriately in a
    'share' directory in the installation directory

  * The PyMol template received some colorscheme
    upgrades for multilayer and FMO/FD (FDD) runs.

  * The GAMESS-FMO writer now supports the hybrid
    orbital projection (HOP) scheme for fragmentation
    as well as the AFO scheme that has always been
    present. HOP orbitals are located in the
    'hmo' folder in the share directory. Limited to
    support of only SP3 C-C bonds but that should
    cover 99 % of all current use cases.

    HOP is disabled per default. You can enable it
    by specifying

    [fmo]
    dohop = True

    in the configuration file.

  * The GAMESS-FMO writer now has improved
    support for multilayer calculations as well.
    Through a new section [qm] in the configuration
    it is possible to specify the basis set for an
    FMO calculation.

    In the case of multilayer runs FragIt
    supports the multilayer notation from FMO
    as well so one specifies basis sets on a
    per layer basis separated by colons:

    [qm]
    basis = STO-3G:3-21G
    
    In multilayer calculations the appropriate
    basis set will be added to the $DATA group
    and leave out the $BASIS group.

    For regular FMO runs just specifying basis
    sets as

    [qm]
    basis = 6-31G(d)

    will suffice.

  * The GAMESS-FMO writer now has support to
    generate FMO/EFP input files. The feature
    is enabled by

    [fmo]
    efpwaters = 1

    in the configuration file. This will, if
    no multilayer options are specified, gen-
    erate an FMO input file with *all* water
    molecules demoted from FMO fragments to
    EFP fragments.

    This new option also supports a slightly
    more involved use case scenario where
    all water molecules some distance R from
    a central molecule, such as a ligand in
    a protein, are demoted to EFP water. The
    options are

    [fmo]
    efpwaters = 1

    [output]
    centralfragment = 1
    boundaries = 3.0

    which will convert all water molecules not
    within a distance R = 3.0 angstrom from 
    fragment 1 to EFP water.

    Please note that this will *disable* the
    multilayer calculation as that option
    seems a little exuberant.


FragIt v1.6.1 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.6.0

  * Fixes a default value in the base configuration
    for FMO input files such that atom names are
    ignored

FragIt v1.6.0 Release Notes
===========================

This major update brings several changes to default keywords
and to the command line interface for FragIt

Updates since v1.5.1

  * Removed --make-config option from FragIt commandline
    tool. Instead, a fragit-conf tool is supplied which
    allows for the generation of different types of
    configuration files for many types of fragmentation jobs.

  * Added fragit-conf tool to generate configuration files
    for FragIt.
    See fragit-conf -h for details.

  * Added multiple default configurations for FragIt.
    FragIt (through fragit-conf) now supports:
    - fragment molecular orbital (FMO) specific fragmentation
      patterns.
    - polarizable embedding (PE) specific fragmentation
      patterns. This is specifically made for use through
      the polarizable embedding assistant script (PEAS).
    The default configuration settings is still FMO.
    See fragit-conf -h.

  * Default values were changed for some parameters for the
    fragmentation group. The new values are:
    - maxfragsize = 100
    - useatomnames = False  - ONLY FMO.

  * The setup script is now more self-contained.

FragIt v1.5.1 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.5.0

  * Fixed issue with FMO writer which dumped errornous
    keyword in FMO input file even if no bonds were
    broken. (#15)  -- Thanks to Anders S. Christensen
                      for reporting this issue.

FragIt v1.5.0 Release Notes
===========================

This major update brings changes to defaults in keywords and
a major overhaul to the information printed.

The changes to the defaults is mostly to make sure that users
of the API get something sensible back per default without having
to make too many changes before they get started.

Updates since v1.4.4

  * Default values were changed for some parameters for the QM/MM
    group. The new values are:
    - includehbonddonors = False
    - includehbondacceptors = False
    - includecovalent = False

  * Default values were changed for a single parameter for the OUTPUT
    group. The new value is:
    - useatomnames = True

  * A new option to enable printout (especially handy when using fragit
    as an API from other programs) in the OUTPUT group:
    - verbose = False

    when set to true, FragIt will print information on its doings.

  * FragIt can now correctly inform users if it fails because OpenBabel
    is missing. Previously this was not the case. Issue 13.

FragIt v1.4.4 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.4.3

  * Include within a distance in QM/MM.
  * Adds atom names for added hydrogens.
  * Correct naming of atoms added to cap.
  * Cleanup in the code. Removing legacy space / tabs etc.

FragIt v1.4.3 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.4.2

  * Interfaced with Travis CI for continous integration.
  * Build status on README.
  * Moved tests out of the 'src' directory and reconfigured
    tests such that nosetests is now supported.

FragIt v1.4.2 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.4.1

  * FragIt now attemps to set the formal charges of
    metal and counter ions.
  * Naming of atoms should be greatly improved.

FragIt v1.4.1 Release Notes
===========================

This is a bug-fix and minor improvements release.

Updates since  v1.4.0

  * Fixed a bug which could make FragIt enter an infinite
    loop if the files used were only counter ions or metals

FragIt v1.4.0 Release Notes
===========================

Updates since v1.3.7

  * FragIt now correctly reads the configuration file early
    such that settings are used (instead of defaults).

FragIt v1.3.7 Release Notes
===========================

This is a bug-fix and minor improvements release.

  * Fixes an error where naming of atoms were not strings but None.
  * Better naming of atoms and a proper warning when naming fails.

FragIt v1.3.6 Release Notes
===========================

This is a bug-fix and minor improvements release.

  * Fixes error when a QM region not connected with covalent bonds
    were extracted.

FragIt v1.3.5 Release Notes
===========================

This is a bug-fix and minor improvements release.

  * Fixes error in the way QM/MM regions were extracted.

FragIt v1.3.4 Release Notes
===========================

This is a bug-fix and minor improvements release.

  * Adds option to not include covalent neighbours to QM fragments.
    It is not always wanted to include covalently bound neighbours.

FragIt v1.3.3 Release Notes
===========================

This is a bug-fix and minor improvements release.

 * P-H hydrogen bond distance added. This is used when using
   FragIt to cap DNA with the MFCC capping procedure.

 * Adds magnesium ions to ignore metals list.

 * Attempt to add residue names in a less fancy way. The old naming
   scheme gave "none" too often. Now, attempt to use info from
   PDB files if that exists

 * Fixed color space for fragments.

 * Expose QMMM features through config files. Options were added for
   inclusion based upon hydrogen bonds.

FragIt v1.3.2 Release Notes
===========================

This is a bug-fix and minor improvements release.

 * Adds support for atomnames read from PDB files for instance.

 * Boolean options are read as strings. There was a problem with
   options not being parsed correctly.

FragIt v1.3.1 Release Notes
===========================

This is a bug-fix and minor improvements release.

 * Phosphorus atoms are not to be excluded as a metal atom. This concerns mostly DNA.

FragIt v1.3.0 Release Notes
===========================

Updates since v1.2.1

 * stabilization of the API for 3rd party uses

 * API has support for a QM-region in the XYZ-MFCC
   writer. This QM-region can be expanded by including
   covalently bound or hydrogen bound neighbours.

   This support is slated to be included more generally
   in a future fragit 1.4 release series.

 * A visually more pleasing color scheme. We are no
   longer living in an 8-bit world.

 * FragIt now exists more cleanly if it cannot find itself.

FragIt v1.2.1 Release Notes
===========================

Updates since v1.2.0

 * fixed naming scheme in xyz and xyz-mfcc writers
   to be similar

FragIt v1.2.0 Release Notes
===========================

Updates since v1.0.4

 * MFCC style fragmentation is beginning to be work
   however, it will take some time before this matures

 * Slight extension to the API to allow for better
   interoperability with external software wishing
   to use FragIt.

 * Squash fragments

 * Found a reliable workflow - including bug fixes
   to treat metal ions. There are some constraints
   still.

FragIt v1.0.4 Release Notes
===========================

Updates since v1.0.3

 * specify charge models via the command
   line interface (--charge-model=) or
   configuration files.

FragIt v1.0.3 Release Notes
===========================

Updates since v1.0.2

 * Renamed GAMESS writer GAMESS-FMO due
   to ambiguity.

FragIt v1.0.2 Release Notes
===========================

Updates since v1.0.1

 * Most of the options in the command line
   version of FragIt now uses the underlying
   FragItConfig object for defaults.

FragIt v1.0.0 Release Notes
===========================

 * The code is now distributed as a library

 * Rename option for distance based layering to
     --output-active-distance

 * Add option to disable protection patterns
   from command line

 * Add fragmentation of DNA

 * Add code and patterns to merge glycine to
   neighbor fragments. Other fragments can
   also be merged via the new mergepatterns
   group in the configuration files

 * Add possibility to generate configuration
   file without needing an input file. In the
   same instance, file-specific fragmentation
   info has been removed from the configuration
   files.

 * Add fix so custom patterns work consistently
   in all cases.

 * Add fix to make explicit pairs of atoms always
   work
