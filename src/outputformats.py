"""
**********************************************************************
outputformats.py

Copyright (C) 2011-2013 Casper Steinmann

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

import gamessfmo
import xyzmfcc

def get_writer_and_extension(theformat):
	formats = supported_output_formats()
	extensions = supported_output_fileexts()
	if not formats.has_key(theformat):
		raise ValueError("The format you requested is not available")
	return (formats[theformat],extensions[theformat])

## Returns ALL supported output formats
def supported_output_formats():
	formats = dict()
	formats['GAMESS-FMO'] = gamessfmo.GamessFMO
	formats['XYZ-MFCC'] = xyzmfcc.XYZMFCC
	return formats

## Returns ALL supported output formats
def supported_output_fileexts():
	formats = dict()
	formats['GAMESS-FMO'] = ".inp"
	formats['XYZ-MFCC'] = ".xyz"
	return formats

