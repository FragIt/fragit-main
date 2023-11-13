"""
Copyright (C) 2011-2017 Casper Steinmann
"""
import os
import unittest
from typing import Callable

import fragit.util as util


class TestUtilModule(unittest.TestCase):

    def setUp(self):
        self.btrue = True
        self.bfalse = False

        self.ione = 1
        self.fone = 1.0
        self.isix = 6

        self.lsimple = list(range(self.ione, self.isix))
        self.tsimple = tuple(self.lsimple)
        self.ssimple = "lorem ipsum dolor."
        self.dsimple = {"a": "value", "b": "myvalue"}

    @staticmethod
    def delete_file(filename):
        try:
            f = open(filename)
        except IOError:
            return
        else:
            f.close()
            os.remove(filename)

    def type_error_base_tests(self, func: Callable):
        self.assertRaises(TypeError, func, self.ione)
        self.assertRaises(TypeError, func, self.fone)
        self.assertRaises(TypeError, func, self.btrue)
        self.assertRaises(TypeError, func, self.tsimple)
        self.assertRaises(TypeError, func, self.lsimple)

    def test_file_extension(self):
        """ file extension """
        self.type_error_base_tests(util.file_extension)

        teststring1 = "file.txt"
        self.assertEqual(util.file_extension(teststring1), ".txt")

        teststring2 = "/home/user/file2.ext"
        self.assertEqual(util.file_extension(teststring2), ".ext")

    def test_file_basename(self):
        self.type_error_base_tests(util.file_basename)

        teststring1 = "file.txt"
        self.assertEqual(util.file_basename(teststring1), "file")

        teststring2 = "/home/user/file2.ext"
        self.assertEqual(util.file_basename(teststring2), "file2")

    def test_Uniqify(self):
        """ unique items in a list """
        self.assertRaises(TypeError, util.remove_duplicates, self.ione)
        self.assertRaises(TypeError, util.remove_duplicates, self.fone)
        self.assertRaises(TypeError, util.remove_duplicates, self.btrue)
        self.assertRaises(TypeError, util.remove_duplicates, self.ssimple)
        self.assertEqual(util.remove_duplicates(self.lsimple), self.lsimple)
        self.assertEqual(util.remove_duplicates([1, 2, 2, 3, 4, 5, 5]), self.lsimple)

    def test_uniqifyListOfLists(self):
        """ unique list of lists """
        test_array = [[1, 2, 3, 4],
                      [5, 2, 3, 4],
                      [1, 2, 3, 4]]
        self.assertEqual(util.uniqifyListOfLists(test_array), [[1, 2, 3, 4], [5, 2, 3, 4]])

    def test_ravel2DEmpty(self):
        test_array = []
        self.assertEqual(util.flatten(test_array), [])

    def test_ravel2D(self):
        test_array = [[1, 2, 3, 4], [1, 2, 3, 4]]
        self.assertEqual(util.flatten(test_array), [1, 2, 3, 4, 1, 2, 3, 4])

    def test_ravel2DOddSize(self):
        test_array = [[1, 2, 3, 4], [1, 2, 3, 4, 5]]
        self.assertEqual(util.flatten(test_array), [1, 2, 3, 4, 1, 2, 3, 4, 5])

    def test_ravel2DEmptySublists(self):
        test_array = [[1, 2, 3, 4], [], [1, 2, 3, 4]]
        self.assertEqual(util.flatten(test_array), [1, 2, 3, 4, 1, 2, 3, 4])

    def test_listDiff(self):
        test_array1 = list(range(10))
        test_array2 = list(range(3, 7))
        self.assertEqual(util.difference(test_array1, test_array2), [0, 1, 2, 7, 8, 9])

    def test_listTo2D(self):
        # check empty, it should return the correct dimension still
        test_array = []
        final_array = [[]]
        self.assertEqual(util.list_to_2d(test_array, 5), final_array)

        test_array = list(range(20))
        final_array = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14], [15, 16, 17, 18, 19]]
        self.assertEqual(util.list_to_2d(test_array, 5), final_array)

        # test the odd sized conversion too
        test_array = list(range(19))
        final_array = [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11, 12, 13, 14], [15, 16, 17, 18]]
        self.assertEqual(util.list_to_2d(test_array, 5), final_array)

    def test_join2D(self):
        test_array = [["a", "b", "c"], ["ab", "cd", "ef"]]
        self.assertEqual(util.list_2d_to_str(test_array, "|", "--"), 'a|b|c--ab|cd|ef')

        test_array = [["1", "2", 4], ["1", "2", "5"]]
        self.assertRaises(TypeError, util.list_2d_to_str, test_array, "|", "--")

    def test_intlistFromString(self):
        self.type_error_base_tests(util.int_list_from_string)
        self.assertEqual(util.int_list_from_string(""), [])
        self.assertEqual(util.int_list_from_string("1,2,3"), [1, 2, 3])
        self.assertEqual(util.int_list_from_string("-1,2,3"), [-1, 2, 3])

    def test_floatlistToString(self):
        self.type_error_base_tests(util.float_list_from_string)
        self.assertEqual(util.float_list_from_string(""), [])
        self.assertEqual(util.float_list_from_string("1.0,2.0,3.0"), [1.0, 2.0, 3.0])
        self.assertEqual(util.float_list_from_string("-1.0,2.0,3.0"), [-1.0, 2.0, 3.0])

    def test_isStringList(self):
        self.assertRaises(TypeError, util.is_string_list, self.ione)
        self.assertRaises(TypeError, util.is_string_list, self.fone)
        self.assertRaises(TypeError, util.is_string_list, self.btrue)
        self.assertRaises(TypeError, util.is_string_list, self.tsimple)
        test_array = ["a", "2", "cv"]
        self.assertEqual(util.is_string_list(test_array), True)
        test_array = [["1", "f"], ["3", "5"]]
        self.assertEqual(util.is_string_list(test_array), True)
        test_array = [["1", "f"], ["3", 5]]
        self.assertEqual(util.is_string_list(test_array), False)

    def test_isIntegerList(self):
        self.assertRaises(TypeError, util.is_integer_list, self.ione)
        self.assertRaises(TypeError, util.is_integer_list, self.fone)
        self.assertRaises(TypeError, util.is_integer_list, self.btrue)
        self.assertRaises(TypeError, util.is_integer_list, self.ssimple)

        test_array = [1, 2, 3]
        self.assertEqual(util.is_integer_list(test_array), True)

        test_array = [[2, 4], [5, 1]]
        self.assertEqual(util.is_integer_list(test_array), True)

        test_array = [["1", 2], [4, 5]]
        self.assertEqual(util.is_integer_list(test_array), False)

        test_array = [[1.0, 2], [4, 5]]
        self.assertEqual(util.is_integer_list(test_array), False)

    def test_WriteStringToFile(self):
        test_string = "Hello World"
        test_filename = "test.dat"
        util.write_string_to_file(test_filename, test_string)
        f = open(test_filename, "r")
        read_string = f.read()
        f.close()
        self.assertEqual(read_string, test_string)

        # only strings can be written
        self.assertRaises(TypeError, util.write_string_to_file, test_filename, self.ione)
        self.assertRaises(TypeError, util.write_string_to_file, test_filename, self.fone)
        self.assertRaises(TypeError, util.write_string_to_file, test_filename, self.btrue)
        self.assertRaises(TypeError, util.write_string_to_file, test_filename, self.tsimple)
        self.assertRaises(TypeError, util.write_string_to_file, test_filename, self.lsimple)

        # filename must be a string
        self.assertRaises(TypeError, util.write_string_to_file, self.ione, test_string)
        self.assertRaises(TypeError, util.write_string_to_file, self.fone, test_string)
        self.assertRaises(TypeError, util.write_string_to_file, self.btrue, test_string)
        self.assertRaises(TypeError, util.write_string_to_file, self.tsimple, test_string)
        self.assertRaises(TypeError, util.write_string_to_file, self.lsimple, test_string)

        self.delete_file(test_filename)

    def test_ReadStringFromFile(self):
        test_string = "Hello World"
        test_filename = "test.dat"
        with open(test_filename, 'w') as f:
            f.write(test_string)

        read_string = util.read_string_from_file(test_filename)
        self.assertEqual(read_string, test_string)

        # filename must be a string
        self.assertRaises(TypeError, util.read_string_from_file, self.ione)
        self.assertRaises(TypeError, util.read_string_from_file, self.fone)
        self.assertRaises(TypeError, util.read_string_from_file, self.btrue)
        self.assertRaises(TypeError, util.read_string_from_file, self.tsimple)
        self.assertRaises(TypeError, util.read_string_from_file, self.lsimple)
        self.assertRaises(TypeError, util.read_string_from_file, self.dsimple)

        self.delete_file(test_filename)
        
    def test_ReadStringListFromFile(self):
        test_string = "Hello World\nWelcome Home"
        test_list = ["Hello World", "Welcome Home"]
        test_filename = "test.dat"
        f = open(test_filename, "w")
        f.write(test_string)
        f.close()
        read_list = util.read_string_list_from_file(test_filename)
        self.assertEqual(read_list, test_list)

    def test_tupleValuesInEitherList(self):
        list1 = [1, 2, 3, 4]
        list2 = [5, 6, 7, 8]
        self.assertEqual(util.is_tuple_values_in_either_list((1, 5), list1, list2), True)
        self.assertEqual(util.is_tuple_values_in_either_list((5, 1), list1, list2), True)
        self.assertEqual(util.is_tuple_values_in_either_list((1, 9), list1, list2), False)
        self.assertEqual(util.is_tuple_values_in_either_list((6, 9), list1, list2), False)

    def test_listToRanges(self):
        list1 = list(range(1, 10))
        self.assertEqual(util.list_to_ranges(list1), [(1, 9)])
        list2 = [1, 2] + list(range(5, 10))
        self.assertEqual(util.list_to_ranges(list2), [1, 2, (5, 9)])
        list3 = [1, 2, 3] + list(range(5, 10))
        self.assertEqual(util.list_to_ranges(list3), [(1, 3), (5, 9)])
        list4 = [1, 2, 4] + list(range(7, 10))
        self.assertEqual(util.list_to_ranges(list4), [1, 2, 4, (7, 9)])
        list5 = [2, 3] + list(range(5, 9)) + [12, 13] + list(range(15, 20)) + [31]
        self.assertEqual(util.list_to_ranges(list5), [2, 3, (5, 8), 12, 13, (15, 19), 31])
        list6 = [2, 3] + list(range(5, 9)) + [12, 13] + list(range(15, 20)) + [31, 32]
        self.assertEqual(util.list_to_ranges(list6), [2, 3, (5, 8), 12, 13, (15, 19), 31, 32])


def suite():
    s = unittest.TestSuite()
    s.addTest(unittest.makeSuite(TestUtilModule))
    return s


if __name__ == '__main__':
    unittest.main()
