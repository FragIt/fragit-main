"""
**********************************************************************
outputformats.py

Copyright (C) 2011-2012 Casper Steinmann

This file is part of the FragIt project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
"""

import gamess

def get_writer_and_extension(theformat):
	formats = supported_output_formats()
	extensions = supported_output_fileexts()
	if not formats.has_key(theformat):
		raise ValueError("The format you requested is not available")
	return (formats[theformat],extensions[theformat])

## Returns ALL supported output formats
def supported_output_formats():
	formats = dict()
	formats['GAMESS'] = gamess.Gamess
	return formats

## Returns ALL supported output formats
def supported_output_fileexts():
	formats = dict()
	formats['GAMESS'] = ".inp"
	return formats

