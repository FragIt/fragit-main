"""
Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2016 Casper Steinmann
"""
import os
import string
from typing import List, Tuple, Any, Union, Callable, Optional

from fragit.fragit_exceptions import OBNotFoundException
try:
    from openbabel import openbabel
except ImportError:
    raise OBNotFoundException("OpenBabel not found. Please install OpenBabel to use FragIt.")

import numpy as np

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


def file_extension(path_to_file: str) -> str:
    (filename, extension) = get_filename_and_extension(path_to_file)
    return extension


def file_basename(path_to_file: str) -> str:
    (filename, extension) = get_filename_and_extension(path_to_file)
    return filename


def get_filename_and_extension(path_to_file: str) -> Tuple[str, str]:
    if not isinstance(path_to_file, str):
        raise TypeError
    basename = os.path.split(path_to_file)[1]
    return os.path.splitext(basename)


def remove_duplicates(input_list: Union[List, Tuple]) -> List:
    """ Removes any duplicates """
    invalid_types = [int, float, bool, str]
    is_invalid = False
    for invalid_type in invalid_types:
        is_invalid = is_invalid or isinstance(input_list, invalid_type)

    if is_invalid:
        raise TypeError
    unique_set = set(input_list)
    return list(unique_set)


def uniqifyListOfLists(thelist: List[List[int]]) -> List[List[int]]:
    result = list()
    keys = list()
    
    for sublist in thelist:
        if sublist[0] in keys:
            continue
        result.append(sublist)
        keys.append(sublist[0])
        
    return result


def flatten(input_list: List[List[Any]]) -> List[Any]:
    """ Flattens a list of lists into a list """
    flattened_list = list()
    for i in input_list:
        for j in i:
            flattened_list.append(j)
    return flattened_list


def difference(first: List[int], second: List[int]) -> List[int]:
    set1 = set(first)
    set2 = set(second)
    return list(set1 - set2)


def list_to_2d(input_list: List, length: int, element_format: Optional[str] = None) -> List[List]:
    output_list: List[List] = list()
    temp_list: List = list()
    
    for element in input_list:
        if element_format is None:
            temp_list.append(element)
        else:
            temp_list.append(element_format % element)
        if len(temp_list) == length:
            output_list.append(temp_list)
            temp_list = list()
    
    # add the rest
    if len(temp_list) != 0:
        output_list.append(temp_list)

    # the unusual empty list
    if len(output_list) == 0:
        output_list = [[]]

    return output_list


def list_2d_to_str(input_list: List[List[str]], element_divisor: str, list_divisor: str) -> str:
    if not is_string_list(input_list):
        raise TypeError 

    output_list = list()
    x: List[str]
    for x in input_list:
        output_list.append(element_divisor.join(x))

    return list_divisor.join(output_list)


def list_from_string(input_string: str, op: Callable, seperator: str = ","):
    if not isinstance(input_string, str):
        raise TypeError("Argument 'string' to valuelist_from_string is not of type str")

    if not isinstance(op, type):
        raise TypeError("Argument 'op' to valuelist_from_string is not of type int or float")

    if len(input_string) == 0:
        return []

    return [op(value) for value in input_string.split(seperator)]


def int_list_from_string(input_string: str, seperator:str = ","):
    """ Returns a list of integers from a string

        Arguments:
        string -- the string to convert to list of integers
        seperator -- the seperator to use when converting the string to a list

        Returns list of integers
    """
    return list_from_string(input_string, int, seperator)


def float_list_from_string(input_string: str, seperator: str = ","):
    """ Returns a list of floats from a string

        Arguments:
        string -- the string to convert to list of floats
        seperator -- the seperator to use when converting the string to a list

        Returns list of floats 
    """
    return list_from_string(input_string, float, seperator)


def is_string_list(string_list: Union[List[str], List[List[str]]]) -> bool:
    if not isinstance(string_list, list):
        raise TypeError

    for item in string_list:
        if isinstance(item, list):
            if not is_string_list(item):
                return False
        elif not isinstance(item, str):
            return False

    return True


def is_integer_list(string_list: Union[List[int], List[List[int]]]) -> bool:
    if not isinstance(string_list, list) and not isinstance(string_list, tuple):
        raise TypeError

    for item in string_list:
        if isinstance(item, list):
            if not is_integer_list(item):
                return False
        elif not isinstance(item, int):
            return False

    return True


def write_string_to_file(filename: str, value: str):
    if not isinstance(filename, str):
        raise TypeError
    if not isinstance(value, str):
        raise TypeError
    with open(filename, 'w') as f:
        f.write(value)


def read_string_from_file(filename):
    if not isinstance(filename, str):
        raise TypeError
    f = open(filename,'r')
    string = f.read()
    f.close()
    return string


def read_string_list_from_file(filename):
    string = read_string_from_file(filename)
    string = string.replace("\r","") # windows line endings
    return string.split("\n")


def is_tuple_values_in_either_list(the_tuple: Tuple[int, int], list1: List[int], list2: List[int]) -> bool:
    if not isinstance(the_tuple, tuple):
        raise ValueError
    in_lhs = the_tuple[0] in list1 and the_tuple[1] in list2
    in_rhs = the_tuple[1] in list1 and the_tuple[0] in list2
    return in_lhs or in_rhs


def list_to_ranges(the_list: List[int]) -> List[Union[int, Tuple[int, int]]]:
    """ Converts a list of numbers into possible ranges defined as tuples

    example: [1,2,3,4] -> [(1,4)]

    :param the_list:
    :return:
    """
    def expand_range(output: List[Union[int, Tuple[int, int]]], first: int, last: int):
        if first == last:
            output.append(last)
        elif first + 1 == last:
            output.extend([first, last])
        else:
            output.append((first, last))

    result: List[Union[int, Tuple[int, int]]] = []

    if len(the_list) == 0:
        return result

    first_pass = True
    first: int = 0
    last: int = 0
    for value in the_list:
        if first_pass:
            first_pass = False
            first = value
            last = value
            continue

        if last+1 != value:
            expand_range(result, first, last)
            first = value
        last = value
    expand_range(result, first, last)
    return result


def list_of_ranges_to_string(atom_list: List[Union[int, Tuple[int, int]]],
                             maxlength: int = 72,
                             line_format: str = "%10s",
                             item_format: str = "%7s",
                             tuple_format: str = "%7s%7s",
                             terminator_format= "%7i\n"
                             ):
    result = ""
    line = line_format % ""
    for atom_index in atom_list:
        # atom indices are either direct indices
        if isinstance(atom_index, int):
            tmp = item_format % atom_index
        else:
            # or they are tuples requiring a special format
            tmp = tuple_format % (atom_index[0], ('-' + str(atom_index[1])))
        if len(line+tmp+"\n") >= maxlength:
            result += line + "\n"
            line = line_format % "" + tmp
        else:
            line += tmp
    if terminator_format is not None:
        terminator = terminator_format % 0
        if len(line+terminator) >= maxlength:
            result += line + "\n"
            line = line_format % ""
        result += line + terminator
    else:
        result += line[:-1]
    return result


def file_to_mol(filename: str) -> openbabel.OBMol:
    """ Converts a filename to an openbabel molecule """
    file_format = file_extension(filename)[1:]
    obc = openbabel.OBConversion()
    obc.SetInFormat(file_format)
    mol = openbabel.OBMol()
    obc.ReadFile(mol, filename)
    OBCheckMoleculeConsistency(mol)
    return mol


def OBCheckMoleculeConsistency(molecule: openbabel.OBMol):
    if molecule.NumAtoms() < 1:
        raise ValueError("Molecule has no atoms.")


def calculate_hydrogen_position(heavy: openbabel.OBAtom, light: openbabel.OBAtom) -> np.ndarray:
    """ Positions a hydrogen atom in the "correct" position between two points
    """
    table = {6: 1.09, 7: 1.01, 8: 0.96, 16: 1.35, 15: 1.42}
    alpha: float = table[heavy.GetAtomicNum()]
    p1 = np.array([heavy.GetX(), heavy.GetY(), heavy.GetZ()])
    p2 = np.array([light.GetX(), light.GetY(), light.GetZ()])
    n = np.linalg.norm(p2-p1)
    return np.asarray(p1 + alpha/n * (p2 - p1))


def shares_elements(a: List, b: List) -> bool:
    """Returns True if lists (sets) a and b shares elements. Otherwise, false.
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
