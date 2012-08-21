import unittest

import test_util
import test_fragmentation
import test_outputformats
import test_standardwriter
import test_gamessoutput

if __name__ == '__main__':
  suite = unittest.TestSuite()
  suite.addTest(test_util.suite())
  suite.addTest(test_fragmentation.suite())
  suite.addTest(test_outputformats.suite())
  suite.addTest(test_standardwriter.suite())
  suite.addTest(test_gamessoutput.suite())
  unittest.TextTestRunner().run(suite)
