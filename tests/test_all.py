import unittest

from .test_util import suite as tu_suite
from .test_fragmentation import suite as tf_suite
from .test_outputformats import suite as to_suite
from .test_writers import suite as tw_suite
from .test_gamessfmooutput import suite as tg_suite
from .test_qmmm_interface import suite as qm_suite
from .test_mfcc_interface import suite as mf_suite

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(tu_suite)
    suite.addTest(tf_suite)
    suite.addTest(to_suite)
    suite.addTest(tw_suite)
    suite.addTest(tg_suite)
    suite.addTest(qm_suite)
    suite.addTest(mf_suite)
    unittest.TextTestRunner().run(suite)
