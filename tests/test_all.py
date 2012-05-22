import unittest

import test_util
import test_Fragmentation
import test_OutputFormats

suite = unittest.TestSuite()
suite.addTest(test_util.suite())
suite.addTest(test_Fragmentation.suite())
suite.addTest(test_OutputFormats.suite())

unittest.TextTestRunner().run(suite)
