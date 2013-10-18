"""
**********************************************************************
util.py

Copyright (C) 2010-2011 Mikael W. Ibsen
Some portions Copyright (C) 2011-2012 Casper Steinmann

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
import os
import openbabel
import numpy

def isEqual(a, b):
	return a == b and type(a) == type(b)

def file_exists(filename):
	if not is_string(filename):
		raise TypeError
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
	if not is_string(path_to_file):
		raise TypeError
	basename = os.path.split(path_to_file)[1]
	return os.path.splitext(basename)

def TupleToStringTypeRepresentation(tuple_of_vars):
	if not is_tuple( tuple_of_vars ):
	    raise TypeError
	if (len(tuple_of_vars) == 0): return ""

	tuple_string = toString(tuple_of_vars[0])
	for i in range(1, len(tuple_of_vars)):
		tuple_string += ", " + toString(tuple_of_vars[i])
		
	return tuple_string

def toString(msg):
	if is_int(msg)			:	msg = "Int: " 			+ str(msg)
	if is_float(msg)		:	msg = 'Float: '			+ str(msg)
	if is_list(msg)			:	msg = 'List: '			+ str(msg)
	if msg == None			:	msg = 'None: '			+ str(msg)
	if not is_string(msg)	:	msg = str(type(msg)) + ': '	+ str(msg)
	return msg

def is_string(input):
	return type(input) == type('a')

def is_bool(input):
	return type(input) == type(True)

def is_int(input):
	return type(input) == type(1)

def is_float(input):
	return type(input) == type(1.0)

def is_list(input):
	return type(input) == type([])

def is_tuple(input):
	return type(input) == type((1,2))

def is_dictionary(input):
	return type(input) == type(dict())

def toList(input):
	return [p for p in input]

def maximum_value(list_of_values):
	if not is_list( list_of_values ):
		raise TypeError
	if len(list_of_values) == 0: return 0
	max_value = max(list_of_values)
	min_value = min(list_of_values)
	if abs(min_value) > abs(max_value):
		return min_value
	return max_value

def Uniqify(thelist):
	if is_int(thelist) or is_float(thelist) or is_bool(thelist) or is_string(thelist):
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
	if not is_list(list1) or not is_list(list2): raise TypeError
	if len(list2) > len(list1): raise ValueError("The first argument should be the list with the most items.")

	difference_list = list()
	for i in list1:
		if i not in list2:
			difference_list.append(i)
	return difference_list

def lenOfLists(lol):
	result = list()
	for i in lol:
		result.append(len(i))
	return result

def listTo2D(list1D, sublength, elmFormat = None):
	list2D 	= list()
	tmplist = list()
	
	for element in list1D:
		if (elmFormat is None):	tmplist.append(element)
		else:			tmplist.append(elmFormat % element)
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
	
	
def intlistFromString(string):
	if not is_string(string):
		raise TypeError
	if (string == ""): return []
	intlist = string.split(",")
	result = list()
	for i in intlist:
		i = int(i)
		result.append(i)

	return result

def floatlistFromString(string):
	if not is_string(string):
		raise TypeError
	if(string == ""): return []
	floatlist = string.split(",")
	result = list()
	for f in floatlist:
		f = float(f)
		result.append(f)
	return result



def listOfDoubleIntTupleToString(diTuple, sep = ";"):
	assert is_list(diTuple)

	tmp = list()

	for (a1,a2) in diTuple:
		assert is_int(a1)
		assert is_int(a2)
		tmp.append( "(" + str(a1) + "," + str(a2) + ")")

	return ";".join(tmp)


def listOfDoubleIntTupleFromString(string, sep = ";"):
	if not is_string(string):
		raise TypeError
	result = list()

	if (string == ""): return list()

	for pair in string.split(";"):
		pair = pair[1:-1]
		ints = pair.split(",")
		result_tuple = (int(ints[0]), int(ints[1]))
		result.append(result_tuple)
	return result

def isStringList( string_list ):
	if not is_list( string_list ):
		raise TypeError

	for item in string_list:
		if is_list( item ):
			if not isStringList( item ):
				return False
		elif not is_string( item ):
			return False

	return True

def isIntegerList(string_list):
	if not is_list(string_list) and not is_tuple(string_list):
		raise TypeError

	for item in string_list:
		if is_list(item):
			if not isIntegerList(item):
				return False
		elif not is_int(item):
			return False

	return True

def WriteStringToFile(filename, string):
	if not is_string(string):
		raise TypeError
	f = open(filename, 'w')
	f.write(string)
	f.close()

def WriteStringListToFile(filename,thelist):
	if not isStringList(thelist):
		raise TypeError
	f = open(filename,'w')
	for item in thelist:
		f.write("%s\n" % item)
	f.close()

def ReadStringFromFile(filename):
	f = open(filename,'r')
	string = f.read()
	f.close()
	return string

def ReadStringListFromFile(filename):
	string = ReadStringFromFile(filename)
	string = string.replace("\r","") # windows line endings
	return string.split("\n")


def tupleValuesInEitherList(the_tuple,list1,list2):
	if not is_tuple(the_tuple): raise ValueError
	lhs = the_tuple[0] in list1 and the_tuple[1] in list2
	rhs = the_tuple[1] in list1 and the_tuple[0] in list2
	return lhs or rhs

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
		if (is_int(i)):
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

def calculate_hydrogen_position(heavy, light):
	""" Positions a hydrogen atom in the "correct" position between two points
	"""
	table = {6: 1.09, 7: 1.01, 8: 0.96, 16: 1.35}
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
