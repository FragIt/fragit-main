"""
Copyright (C) 2011-2016 Casper Steinmann
"""

from .gamessfmo import GamessFMO
from .xyzmfcc import XYZMFCC
from .xyz import XYZ

def get_writer_and_extension(theformat):
    formats = supported_output_formats()
    extensions = supported_output_fileexts()
    if theformat not in formats:
        raise ValueError("The format you requested is not available")
    return (formats[theformat],extensions[theformat])

## Returns ALL supported output formats
def supported_output_formats():
    formats = dict()
    formats['GAMESS-FMO'] = GamessFMO
    formats['XYZ-MFCC'] = XYZMFCC
    formats['XYZ'] = XYZ
    return formats

## Returns ALL supported output formats
def supported_output_fileexts():
    formats = dict()
    formats['GAMESS-FMO'] = ".inp"
    formats['XYZ-MFCC'] = ".xyz"
    formats['XYZ'] = ".xyz"
    return formats

