"""
**********************************************************************
cmdline.py

Copyright (C) 2011-2012 Casper Steinmann

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
# python 2.6
import sys
from optparse import OptionParser, OptionGroup

from util import fileToMol, file_basename
from fragmentation import Fragmentation

from outputformats import get_writer_and_extension


import strings

def main(argv=None):
	if argv is None:
		argv = sys.argv[1:]

	parser = OptionParser(usage=strings.usage,
				description=strings.description,
				version=strings.version_str)
	parser.add_option("-o", "--output", dest="outputfile", type=str, default="", metavar="filename")

	configuration = OptionGroup(parser, "Configuration")
	general = OptionGroup(parser, "Fragmentation")
        output = OptionGroup(parser, "Output")

	configuration.add_option("--use-config", dest="useconfigfile", type=str, default="", metavar="filename",
					help="Specify configuration file to use. This will ignore other command line parameters.")
        configuration.add_option("--make-config", dest="makeconfigfile", type=str, default="", metavar="filename",
					help="Specify a filename to use as a configuration file. Use command line options to modify defaults. It is possible to use this command without specifying an input file to generate a clean configuration file.")

	general.add_option("-m", "--maxfragsize", dest="maxFragmentSize", type=int, default=50,metavar="integer",
				help="The maximum fragment size allowed [default: %default]")
	general.add_option("-g", "--groupcount", dest="groupcount", type=int, default=1,metavar="integer",
				help="Specify number of consecutive fragments to combine into a single fragment [default: %default]")
	general.add_option("--disable-protection", dest="disable_protection", action="store_true", default=False,
				help="Specify this flag to disable the use protection patterns.")
	general.add_option("--merge-glycine", action="store_true", dest="merge_glycine", default=False,
				help="Merge a glycine to the neighbor fragment when fragmenting proteins.")
	output.add_option("--output-format", dest="format", type=str, default="GAMESS",
				help="Output format [%default]")
	output.add_option("--output-boundaries", dest="boundaries", type=str, default="",metavar="list of floats",
				help="Specifies boundaries for multiple layers. Must be used with --central-fragment option")
	output.add_option("--output-central-fragment", dest="central_fragment", type=int, default=0, metavar="integer",
				help="Specifies the fragment to use as the central one. Used in combination with --output-boundaries to make layered inputs")
	output.add_option("--output-active-distance", dest="active_atoms_distance", type=float, default=0.0, metavar="float",
				help="Atoms within this distance from --output-central-fragment will be active. Use with --output-buffer-distance to add buffer region between active and frozen parts. [default: %default]")
	output.add_option("--output-buffer-distance", dest="maximum_buffer_distance", type=float, default=3.0, metavar="float",
				help="Maximum distance in angstrom from active fragments from which to include nearby fragments as buffers. This option adds and extends to --output-boundaries. [default: %default]")
	output.add_option("--output-freeze-backbone", dest="freeze_backbone", action="store_true", default=False,
				help="Option to freeze the backbone of the active region.")
	output.add_option("--output-jmol-script", dest="output_jmol_script", action="store_true", default=False,
				help="Write a complimentary jmol script for visualization.")
	output.add_option("--output-pymol-script", dest="output_pymol_script", action="store_true", default=False,
				help="Write a complimentary pymol script for visualization.")

	parser.add_option_group(configuration)
	parser.add_option_group(general)
	parser.add_option_group(output)
	(options, args) = parser.parse_args(argv)

	if len(args) == 0 and len(options.makeconfigfile) > 0:
		from config import FragItConfig
		cfg = FragItConfig()
		cfg.writeConfigurationToFile(options.makeconfigfile)
		sys.exit()

	if len(args) != 1:
		parser.print_help()
		sys.exit()

	infile = args[0]

	molecule = fileToMol(infile)
	fragmentation = Fragmentation(molecule)


	# if there is a config file, read it and ignore other command line options
        if len(options.useconfigfile) > 0:
		fragmentation.readConfigurationFromFile(options.useconfigfile)
		(writer, output_extension) = get_writer_and_extension(fragmentation.getOutputFormat())
	else:
		fragmentation.setMaximumFragmentSize(options.maxFragmentSize)
		if options.groupcount > 1: fragmentation.setFragmentGroupCount(options.groupcount)

		(writer, output_extension) = get_writer_and_extension(options.format)

	outfile = "%s%s" % (file_basename(infile), output_extension)
	if len(options.outputfile) > 0:
		outfile = options.outputfile


	# do the fragmentation procedure
	if options.disable_protection:
		fragmentation.clearProtectPatterns()
	if options.merge_glycine:
		fragmentation.enableMergeGlycinePattern()
	fragmentation.beginFragmentation()
	fragmentation.doFragmentation()
	fragmentation.doFragmentMerging()
	if fragmentation.getFragmentGroupCount() > 1:
		fragmentation.doFragmentGrouping()
	fragmentation.finishFragmentation()

	# write to file
	out = writer(fragmentation)

	# set options from command line
	boundaries = options.boundaries
	central_fragment = options.central_fragment
	active_atoms_distance = options.active_atoms_distance
	maximum_buffer_distance = options.maximum_buffer_distance
	freeze_backbone = options.freeze_backbone
	output_pymol_script = options.output_pymol_script
	output_jmol_script = options.output_jmol_script

	# set options from config file
	if len(options.useconfigfile) > 0:
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
	if len(options.makeconfigfile) > 0:
		fragmentation.setBoundaries(boundaries)
		fragmentation.writeConfigurationToFile(options.makeconfigfile)
