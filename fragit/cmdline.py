"""
Copyright (C) 2011-2016 Casper Steinmann
"""
import sys
import argparse

from fragit.util import file_to_mol, file_basename
from fragit.fragmentation import Fragmentation
from fragit.outputformats import get_writer_and_extension, supported_output_formats
from fragit.strings import version_str, doc_str


def main(directories, argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # load defaults so we can use them below
    from .config import FragItConfig, ConfigSettings
    cfg = FragItConfig(defaults=ConfigSettings['BARE'])

    parser = argparse.ArgumentParser(description=doc_str)

    parser.add_argument("-o", "--output", dest="outputfile", type=str, default="", metavar="FILENAME")
    parser.add_argument("--version", action="version", version="%(prog)s {0:s}".format(version_str))
    parser.add_argument("input", type=str, metavar="INPUTFILE")


    configuration_group = parser.add_argument_group('Configuration File Options', description="""
It is possible to use configuration files with FragIt.
This removes some issues with having many options defined on the commandline.
To use FragIt as an API all options must also be provided in the form of a configuration file.
The simplest way to use configuration files is to first generate it using the fragit-conf program.
Then a configuration file is use by specifying the --use-config flag when running fragit.""")
    configuration_group.add_argument("--use-config", dest="useconfigfile", type=str, default="", metavar="FILENAME",
                    help="Specify configuration file to use. This will ignore other command line parameters.")

    general_group = parser.add_argument_group("Fragmentation Options", description="""
Controls how to fragmentat the provided input molecule.
From the command line the default fragmentation patterns will be used.
It is recommended to use configuration file options as outlined above
for increased flexibility.
""")
    general_group.add_argument("-m", "--maxfragsize", dest="maxFragmentSize", type=int, default=cfg.get_maximum_fragment_size(), metavar="MAX FRAGMENT SIZE",
                               help="The maximum fragment size allowed [default: %(default)s]")
    general_group.add_argument("-g", "--groupcount", dest="groupcount", type=int, default=cfg.get_fragment_group_count(), metavar="GROUP n FRAGMENTS INTO ONE",
                               help="Specify number of consecutive fragments to combine into a single fragment [default: %(default)s]")
    general_group.add_argument("--disable-protection", dest="disable_protection", action="store_true", default=False,
                help="Specify this flag to disable the use protection patterns.")
    general_group.add_argument("--merge-glycine", action="store_true", dest="merge_glycine", default=False,
                help="Merge a glycine to the neighbor fragment when fragmenting proteins.")
    general_group.add_argument("--charge-model", dest="charge_model", default=cfg.get_charge_model(),
                               help="Charge model to use [%(default)s]")
    general_group.add_argument("--combine-fragments", dest="combinefragments", type=str, default="",metavar="LIST OF INTEGERS",
                       help="Combines several fragments into one.")

    output = parser.add_argument_group("Output Options", description="""
Specifies which output format (and options of that format) to make when fragit is finished with the fragmentation.
Some options only apply for some formats.
The XYZ-MFCC format performes the molecular fractionation with conjugate caps procedure such that the produced fragments
are capped (when nescessary) and concaps are created (when nescessary).
""")
    output.add_argument("--output-format", dest="format", type=str, default=cfg.get_writer(),
                        help="Output format [%(default)s]", choices=supported_output_formats().keys())
    output.add_argument("--output-boundaries", dest="boundaries", type=str, default="",metavar="list of floats",
                help="Specifies boundaries for multiple layers. Must be used with --central-fragment option")
    output.add_argument("--output-central-fragment", dest="central_fragment", type=int, default=cfg.get_central_fragment_id(), metavar="integer",
                        help="Specifies the fragment to use as the central one. Used in combination with --output-boundaries to make layered inputs")
    output.add_argument("--output-active-distance", dest="active_atoms_distance", type=float, default=cfg.get_active_atoms_distance(), metavar="float",
                        help="Atoms within this distance from --output-central-fragment will be active. Use with --output-buffer-distance to add buffer region between active and frozen parts. [default: %(default)s]")
    output.add_argument("--output-buffer-distance", dest="maximum_buffer_distance", type=float, default=cfg.get_buffer_distance(), metavar="float",
                        help="Maximum distance in angstrom from active fragments from which to include nearby fragments as buffers. This option adds and extends to --output-boundaries. [default: %(default)s]")
    output.add_argument("--output-freeze-backbone", dest="freeze_backbone", action="store_true", default=cfg.get_freeze_backbone(),
                        help="Option to freeze the backbone of the active region.")
    output.add_argument("--output-jmol-script", dest="output_jmol_script", action="store_true", default=cfg.get_write_jmol_script(),
                        help="Write a complimentary jmol script for visualization.")
    output.add_argument("--output-pymol-script", dest="output_pymol_script", action="store_true", default=cfg.get_write_pymol_script(),
                        help="Write a complimentary pymol script for visualization.")

    args = parser.parse_args()

    infile = args.input

    molecule = file_to_mol(infile)
    conffile = None
    if len(args.useconfigfile) > 0:
        conffile = args.useconfigfile
    fragmentation = Fragmentation(molecule, conffile=conffile)


    # if there is a config file, read it and ignore other command line options
    if len(args.useconfigfile) > 0:
        (writer, output_extension) = get_writer_and_extension(fragmentation.get_output_format())
    else:
        fragmentation.set_charge_model(args.charge_model)
        fragmentation.set_maximum_fragment_size(args.maxFragmentSize)
        fragmentation.set_output_format(args.format)
        if args.groupcount > 1: fragmentation.set_fragment_group_count(args.groupcount)

        (writer, output_extension) = get_writer_and_extension(args.format)

    outfile = "{0:s}{1:s}".format(file_basename(infile), output_extension)
    if len(args.outputfile) > 0:
        outfile = args.outputfile


    # do the fragmentation procedure
    fragmentation.set_combine_fragments(args.combinefragments)
    if args.disable_protection:
        fragmentation.clear_protect_patterns()
    if args.merge_glycine:
        fragmentation.enable_merge_glycine_pattern()
    fragmentation.begin_fragmentation()
    fragmentation.do_fragmentation()
    fragmentation.do_fragment_merging()
    fragmentation.do_fragment_combination()
    if fragmentation.get_fragment_group_count() > 1:
        fragmentation.do_fragment_grouping()
    fragmentation.finish_fragmentation()

    # write to file
    out = writer(fragmentation, directories)

    # set options from command line
    boundaries = args.boundaries
    central_fragment = args.central_fragment
    active_atoms_distance = args.active_atoms_distance
    maximum_buffer_distance = args.maximum_buffer_distance
    freeze_backbone = args.freeze_backbone
    output_pymol_script = args.output_pymol_script
    output_jmol_script = args.output_jmol_script

    # set options from config file
    if len(args.useconfigfile) > 0:
        boundaries = fragmentation.get_boundaries()
        central_fragment = fragmentation.get_central_fragment_id()
        output_pymol_script = fragmentation.get_write_pymol_script()
        output_jmol_script = fragmentation.get_write_jmol_script()
        freeze_backbone = fragmentation.get_freeze_backbone()
        maximum_buffer_distance = fragmentation.get_buffer_distance()
        active_atoms_distance = fragmentation.get_active_atoms_distance()

    # set the options
    out.set_boundaries_from_string(boundaries)
    out.set_central_fragment_id(central_fragment)
    out.set_active_atoms_distance(active_atoms_distance)
    out.set_buffer_max_distance(maximum_buffer_distance)
    if freeze_backbone: out.set_freeze_backbone()
    if output_pymol_script: out.set_pymol_output(infile, outfile)
    if output_jmol_script: out.set_jmol_output(infile, outfile)

    out.setup()
    out.write_file(outfile)

    # write configuration file
    #if len(args.makeconfigfile) > 0:
    #    fragmentation.setBoundaries(boundaries)
    #    fragmentation.writeConfigurationToFile(args.makeconfigfile)
