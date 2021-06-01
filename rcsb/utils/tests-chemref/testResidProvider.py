##
# File:    ResidProviderTests.py
# Author:  J. Westbrook
# Date:    18-Mar-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Test utilities to read the RESID resource file and build loadable documents and identifier
    correspondences.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.ResidProvider import ResidProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ResidProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testReadResidData(self):
        residP = ResidProvider(cachePath=self.__workPath, useCache=False)
        ok = residP.testCache()
        self.assertTrue(ok)
        version = residP.getVersion()
        self.assertTrue(version is not None)
        #


def readResidData():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ResidProviderTests("testReadResidData"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readResidData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
