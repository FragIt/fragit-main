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

Updates since v1.2.0

 * fixed naming scheme in xyz and xyz-mfcc writers
   to be similar

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

Updates since v1.0.3

 * specify charge models via the command
   line interface (--charge-model=) or
   configuration files.

Updates since v1.0.2
--------------------

 * Renamed GAMESS writer GAMESS-FMO due
   to ambiguity.

Updates since v1.0.1
--------------------

 * Most of the options in the command line
   version of FragIt now uses the underlying
   FragItConfig object for defaults.

Updates since v0.9
------------------

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


Fixes since v0.9
------------------

 * Add fix so custom patterns work consistently
   in all cases.

 * Add fix to make explicit pairs of atoms always
   work
