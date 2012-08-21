import unittest

import test_util
import test_Fragmentation
import test_OutputFormats
import test_StandardWriter
import test_GamessOutput

if __name__ == '__main__':
  suite = unittest.TestSuite()
  suite.addTest(test_util.suite())
  suite.addTest(test_Fragmentation.suite())
  suite.addTest(test_OutputFormats.suite())
  suite.addTest(test_StandardWriter.suite())
  suite.addTest(test_GamessOutput.suite())
  unittest.TextTestRunner().run(suite)
