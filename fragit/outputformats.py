"""
Copyright (C) 2011-2016 Casper Steinmann
"""

from fragit.gamessfmo import GamessFMO
from fragit.xyzmfcc import XYZMFCC
from fragit.xyz import XYZ
from fragit.writer import Standard

from typing import Dict, Tuple, Type


def get_writer_and_extension(theformat: str) -> Tuple[Type[Standard], str]:
    formats = supported_output_formats()
    extensions = supported_output_fileexts()
    if theformat not in formats:
        raise ValueError("The format you requested is not available")
    return formats[theformat], extensions[theformat]


def supported_output_formats() -> Dict[str, Type[Standard]]:
    formats: Dict[str, Type[Standard]] = {
        "GAMESS-FMO": GamessFMO,
        "XYZ-MFCC": XYZMFCC,
        "XYZ": XYZ
    }
    return formats


def supported_output_fileexts():
    formats = dict()
    formats['GAMESS-FMO'] = ".inp"
    formats['XYZ-MFCC'] = ".xyz"
    formats['XYZ'] = ".xyz"
    return formats

