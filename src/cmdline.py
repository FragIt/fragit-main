"""
Copyright (C) 2011-2016 Casper Steinmann
"""
import sys
import argparse

from .util import fileToMol, file_basename
from .fragmentation import Fragmentation
from .outputformats import get_writer_and_extension, supported_output_formats
from .strings import version_str, doc_str


def main(directories, argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # load defaults so we can use them below
    from .config import FragItConfig
    cfg = FragItConfig()

    #parser = argparse.ArgumentParser(version=version_str,
    #                                 description=doc_str)
    #                                 #usage=strings.usage)
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
    general_group.add_argument("-m", "--maxfragsize", dest="maxFragmentSize", type=int, default=cfg.getMaximumFragmentSize(),metavar="MAX FRAGMENT SIZE",
                help="The maximum fragment size allowed [default: %(default)s]")
    general_group.add_argument("-g", "--groupcount", dest="groupcount", type=int, default=cfg.getFragmentGroupCount(),metavar="GROUP n FRAGMENTS INTO ONE",
                help="Specify number of consecutive fragments to combine into a single fragment [default: %(default)s]")
    general_group.add_argument("--disable-protection", dest="disable_protection", action="store_true", default=False,
                help="Specify this flag to disable the use protection patterns.")
    general_group.add_argument("--merge-glycine", action="store_true", dest="merge_glycine", default=False,
                help="Merge a glycine to the neighbor fragment when fragmenting proteins.")
    general_group.add_argument("--charge-model", dest="charge_model", default=cfg.getChargeModel(),
                      help="Charge model to use [%(default)s]")
    general_group.add_argument("--combine-fragments", dest="combinefragments", type=str, default="",metavar="LIST OF INTEGERS",
                       help="Combines several fragments into one.")

    output = parser.add_argument_group("Output Options", description="""
Specifies which output format (and options of that format) to make when fragit is finished with the fragmentation.
Some options only apply for some formats.
The XYZ-MFCC format performes the molecular fractionation with conjugate caps procedure such that the produced fragments
are capped (when nescessary) and concaps are created (when nescessary).
""")
    output.add_argument("--output-format", dest="format", type=str, default=cfg.getWriter(),
                help="Output format [%(default)s]", choices=supported_output_formats().keys())
    output.add_argument("--output-boundaries", dest="boundaries", type=str, default="",metavar="list of floats",
                help="Specifies boundaries for multiple layers. Must be used with --central-fragment option")
    output.add_argument("--output-central-fragment", dest="central_fragment", type=int, default=cfg.getCentralFragmentID(), metavar="integer",
                help="Specifies the fragment to use as the central one. Used in combination with --output-boundaries to make layered inputs")
    output.add_argument("--output-active-distance", dest="active_atoms_distance", type=float, default=cfg.getActiveAtomsDistance(), metavar="float",
                help="Atoms within this distance from --output-central-fragment will be active. Use with --output-buffer-distance to add buffer region between active and frozen parts. [default: %(default)s]")
    output.add_argument("--output-buffer-distance", dest="maximum_buffer_distance", type=float, default=cfg.getBufferDistance(), metavar="float",
                help="Maximum distance in angstrom from active fragments from which to include nearby fragments as buffers. This option adds and extends to --output-boundaries. [default: %(default)s]")
    output.add_argument("--output-freeze-backbone", dest="freeze_backbone", action="store_true", default=cfg.getFreezeBackbone(),
                help="Option to freeze the backbone of the active region.")
    output.add_argument("--output-jmol-script", dest="output_jmol_script", action="store_true", default=cfg.getWriteJmolScript(),
                help="Write a complimentary jmol script for visualization.")
    output.add_argument("--output-pymol-script", dest="output_pymol_script", action="store_true", default=cfg.getWritePymolScript(),
                help="Write a complimentary pymol script for visualization.")

    args = parser.parse_args()

    infile = args.input

    molecule = fileToMol(infile)
    conffile = None
    if len(args.useconfigfile) > 0:
        conffile = args.useconfigfile
    fragmentation = Fragmentation(molecule, conffile=conffile)


    # if there is a config file, read it and ignore other command line options
    if len(args.useconfigfile) > 0:
        (writer, output_extension) = get_writer_and_extension(fragmentation.getOutputFormat())
    else:
        fragmentation.setChargeModel(args.charge_model)
        fragmentation.setMaximumFragmentSize(args.maxFragmentSize)
        fragmentation.setOutputFormat(args.format)
        if args.groupcount > 1: fragmentation.setFragmentGroupCount(args.groupcount)

        (writer, output_extension) = get_writer_and_extension(args.format)

    outfile = "%s%s" % (file_basename(infile), output_extension)
    if len(args.outputfile) > 0:
        outfile = args.outputfile


    # do the fragmentation procedure
    fragmentation.setCombineFragments(args.combinefragments)
    if args.disable_protection:
        fragmentation.clearProtectPatterns()
    if args.merge_glycine:
        fragmentation.enableMergeGlycinePattern()
    fragmentation.beginFragmentation()
    fragmentation.doFragmentation()
    fragmentation.doFragmentMerging()
    fragmentation.doFragmentCombination()
    if fragmentation.getFragmentGroupCount() > 1:
        fragmentation.doFragmentGrouping()
    fragmentation.finishFragmentation()

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
        boundaries = fragmentation.getBoundaries()
        central_fragment = fragmentation.getCentralFragmentID()
        output_pymol_script = fragmentation.getWritePymolScript()
        output_jmol_script = fragmentation.getWriteJmolScript()
        freeze_backbone = fragmentation.getFreezeBackbone()
        maximum_buffer_distance = fragmentation.getBufferDistance()
        active_atoms_distance = fragmentation.getActiveAtomsDistance()

    # set the options
    out.setBoundariesFromString(boundaries)
    out.setCentralFragmentID(central_fragment)
    out.setActiveAtomsDistance(active_atoms_distance)
    out.setBufferMaxDistance(maximum_buffer_distance)
    if freeze_backbone: out.setFreezeBackbone()
    if output_pymol_script: out.setPymolOutput(infile,outfile)
    if output_jmol_script: out.setJmolOutput(infile,outfile)

    out.setup()
    out.writeFile(outfile)

    # write configuration file
    #if len(args.makeconfigfile) > 0:
    #    fragmentation.setBoundaries(boundaries)
    #    fragmentation.writeConfigurationToFile(args.makeconfigfile)
