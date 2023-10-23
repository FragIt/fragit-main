"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2016 Casper Steinmann
"""
import os
import string

from .fragit_exceptions import OBNotFoundException
try:
    import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")

import numpy

# converts nuclear charge to atom label
Z2LABEL = {
 1: 'H',                                                             2: 'He',
 3: 'Li',  4: 'Be',  5: 'B',   6: 'C',   7: 'N',  8: 'O',  9: 'F',  10: 'Ne',
11: 'NA', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar'
}

# converts an atomic label to a nuclear charge
LABEL2Z = {}
for key in Z2LABEL:
    LABEL2Z[Z2LABEL[key]] = key


def file_exists(filename):
    if not isinstance(filename, str):
        raise TypeError("filename provided must be a string.")
    try:
        f = open(filename,'r')
    except IOError:
        return False
    f.close()
    return True

def file_extension(path_to_file):
    (filename,extension) = getFilenameAndExtension(path_to_file)
    return extension

def file_basename(path_to_file):
    (filename,extension) = getFilenameAndExtension(path_to_file)
    return filename

def getFilenameAndExtension(path_to_file):
    if not isinstance(path_to_file, str):
        raise TypeError
    basename = os.path.split(path_to_file)[1]
    return os.path.splitext(basename)


def Uniqify(thelist):
    invalid_types = [int, float, bool, str]
    is_invalid = False
    for invalid_type in invalid_types:
        is_invalid = is_invalid or isinstance(thelist, invalid_type)

    if is_invalid:
        raise TypeError
    unique_set = set(thelist)
    return list(unique_set)

def uniqifyListOfLists(thelist):
    result = list()
    keys = list()
    
    for sublist in thelist:
        if (sublist[0] in keys): continue
        result.append(sublist)
        keys.append(sublist[0])
        
    return result

def ravel2D(input):
    raveled_list = list()
    for i in input:
        for j in i:
            raveled_list.append(j)
    return raveled_list

def deepLength(input):
    count = 0
    for i in input:
        for j in i:
            count += 1
    return count

def listDiff(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    return list(set1 - set2)

def listTo2D(list1D, sublength, elmFormat = None):
    list2D     = list()
    tmplist = list()
    
    for element in list1D:
        if (elmFormat is None):    tmplist.append(element)
        else:            tmplist.append(elmFormat % element)
        if (len(tmplist) == sublength):
            list2D.append(tmplist)
            tmplist = list()
    
    # add the rest
    if (len(tmplist) != 0):
        list2D.append(tmplist)

    # the unusual empty list
    if (len(list2D) == 0):
        list2D = [[]]

    return list2D

def join2D(list2D, item_divisor, list_divisor):
    if not isStringList(list2D):
        raise TypeError 

    tmplist = list()
    for x in list2D:
        tmplist.append(item_divisor.join(x))

    return list_divisor.join(tmplist)
    
def joinIntList(glue, intlist):
    if not isIntegerList(intlist):
        raise TypeError
    result = ""
    for i in intlist:
        i = int(i)
        if (result != ""): result += glue
        result += str(i)
    return result

def intlistToString(intlist):
    return joinIntList(",", intlist)


def valuelist_from_string(string, op, seperator=","):
    if not isinstance(string, str):
        raise TypeError("Argument 'string' to valuelist_from_string is not of type str")

    if not isinstance(op, type):
        raise TypeError("Argument 'op' to valuelist_from_string is not of type int or float")

    if len(string) == 0:
        return []

    return [op(value) for value in string.split(seperator)]

def intlistFromString(string, seperator=","):
    """ Returns a list of integers from a string

        Arguments:
        string -- the string to convert to list of integers
        seperator -- the seperator to use when converting the string to a list

        Returns list of integers
    """
    return valuelist_from_string(string, int, seperator)

def floatlistFromString(string, seperator=","):
    """ Returns a list of floats from a string

        Arguments:
        string -- the string to convert to list of floats
        seperator -- the seperator to use when converting the string to a list

        Returns list of floats 
    """
    return valuelist_from_string(string, float, seperator)


def isStringList( string_list ):
    if not isinstance(string_list, list):
        raise TypeError

    for item in string_list:
        if isinstance(item, list):
            if not isStringList(item):
                return False
        elif not isinstance(item, str):
            return False

    return True

def isIntegerList(string_list):
    if not isinstance(string_list, list) and not isinstance(string_list, tuple):
        raise TypeError

    for item in string_list:
        if isinstance(item, list):
            if not isIntegerList(item):
                return False
        elif not isinstance(item, int):
            return False

    return True

def WriteStringToFile(filename, string):
    if not isinstance(filename, str):
        raise TypeError
    if not isinstance(string, str):
        raise TypeError
    with open(filename, 'w') as f:
        f.write(string)


def ReadStringFromFile(filename):
    if not isinstance(filename, str):
        raise TypeError
    f = open(filename,'r')
    string = f.read()
    f.close()
    return string

def ReadStringListFromFile(filename):
    string = ReadStringFromFile(filename)
    string = string.replace("\r","") # windows line endings
    return string.split("\n")


def tupleValuesInEitherList(the_tuple,list1,list2):
    if not isinstance(the_tuple, tuple):
        raise ValueError
    in_lhs = the_tuple[0] in list1 and the_tuple[1] in list2
    in_rhs = the_tuple[1] in list1 and the_tuple[0] in list2
    return in_lhs or in_rhs

def listToRanges(the_list):
    result = []
    if len(the_list) == 0: return result
    first = None
    last = None
    for value in the_list:
        if first is None:
            first = value
            last = value
            continue

        if last+1 != value:
            RangeAction(result,value,first,last)
            first = value
        last = value
    RangeAction(result,value,first,last)
    return result

def RangeAction(result,value,first,last):
    if first == last:
        result.append(last)
    elif first+1 == last:
        result.extend([first,last])
    else:
        result.append((first,last))

def listOfRangesToString(atomList, maxlength = 72, line_format = "%10s", item_format="%7s", tuple_format="%7s%7s",
            terminator_format="%7i\n"):
    result = ''
    line = line_format % ''
    for i in atomList:
        if isinstance(i, int):
            tmp = item_format % i
        else:
            tmp = tuple_format % (i[0], ( '-' + str(i[1])))
        if (len(line+tmp+"\n") >= maxlength):
            result += line + "\n"
            line = line_format % '' + tmp
        else:
            line += tmp
    if( terminator_format is not None ):
        terminator = terminator_format % 0
        if (len(line+terminator) >= maxlength):
            result += line + "\n"
            line    = line_format % ''
        result += line + terminator
    else:
        result += line[:-1]
    return result

def fileToMol(filename):
    file_format = OBGetFormatFromFilename(filename)
    mol = OBMoleculeFromFilenameAndFormat(filename, file_format)
    OBCheckMoleculeConsistency(mol)
    return mol

def OBGetFormatFromFilename(filename):
    return file_extension(filename)[1:]

def OBMoleculeFromFilenameAndFormat(filename, file_format):
    obc = openbabel.OBConversion()
    obc.SetInFormat(file_format)
    mol = openbabel.OBMol()
    obc.ReadFile(mol, filename)
    return mol

def OBCheckMoleculeConsistency(molecule):
    if molecule.NumAtoms() < 1:
        raise ValueError("Molecule has no atoms.")

def getOBAtomVector(atom):
    return numpy.array([atom.GetX(), atom.GetY(), atom.GetZ()])

def calculate_hydrogen_position(heavy, light):
    """ Positions a hydrogen atom in the "correct" position between two points
    """
    table = {6: 1.09, 7: 1.01, 8: 0.96, 16: 1.35, 15: 1.42}
    alpha = table[heavy.GetAtomicNum()]
    p1 = numpy.array([heavy.GetX(), heavy.GetY(), heavy.GetZ()])
    p2 = numpy.array([light.GetX(), light.GetY(), light.GetZ()])
    n = numpy.linalg.norm(p2-p1)
    return p1 + alpha/n * (p2 - p1)

def shares_elements(a, b):
    """Returns True if lists (sets) a and b shares elements. Otherwise false.
    """
    sa = set(a)
    sb = set(b)
    return len(sa & sb) > 0

def directories(from_file):
    """ Sets up directories needed internally in FragIt.
        This pertains to especially the share directory
        that holds the template files.

        Arguments:
        ----------
        from_file -- the basename to use to extract the files

        Returns:
        --------
        dictionary with paths. The following keys are available:
        path -- the base path for the executable
        bin -- the path for the binary
        share -- the path of the share directory

    """
    abs_path = os.path.abspath(from_file)
    bin_path = os.path.dirname(abs_path)
    path = os.path.dirname(bin_path)
    share_path = os.path.join(path, 'share')
    return {'path': path, 'bin': bin_path, 'share': share_path}


def substitute_file(from_file, to_file, substitutions):
    """ Substitute contents in from_file with substitutions and
        output to to_file using string.Template class

        Arguments:
        ----------
        from_file -- template file to load
        to_file -- substituted file
        substitutions -- dictionary of substitutions.
    """
    with open(from_file, "r") as f_in:
        source = string.Template(f_in.read())

        with open(to_file, "w") as f_out:
            outcome = source.safe_substitute(substitutions)
            f_out.write(outcome)
