##
# File:    CODProviderTests.py
# Author:  J. Westbrook
# Date:    8-Feb-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing Crystallographic Open Database (COD) molecule data.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.chemref.CODProvider import CODProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class CODProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testFetchCODSmiles(self):
        cP = CODProvider(cachePath=self.__cachePath, useCache=False)
        ok = cP.testCache()
        self.assertTrue(ok)

    def testXFetchCif(self):
        cP = CODProvider(cachePath=self.__cachePath, useCache=True)
        ok = cP.testCache()
        self.assertTrue(ok)
        for codId in ["2002916", "2002925", "2002926"]:
            ok = cP.fetchCif(codId)
            self.assertTrue(ok)


def fetchCODSmiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CODProviderTests("testFetchCODSmiles"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchCODSmiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
