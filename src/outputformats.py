"""
Copyright (C) 2011-2016 Casper Steinmann
"""

import gamessfmo
import xyzmfcc
import xyz

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
    formats['XYZ'] = xyz.XYZ
    return formats

## Returns ALL supported output formats
def supported_output_fileexts():
    formats = dict()
    formats['GAMESS-FMO'] = ".inp"
    formats['XYZ-MFCC'] = ".xyz"
    formats['XYZ'] = ".xyz"
    return formats

