##
# File:    RcsbLigandScoreProviderTests.py
# Author:  J. Westbrook
# Date:    8-Feb-2021
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities managing supporting data for RCSB ligand scoring.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.chemref.RcsbLigandScoreProvider import RcsbLigandScoreProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class RcsbLigandScoreProviderTests(unittest.TestCase):
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

    def testFetchScoreFiles(self):
        logger.info("Testing reload of ligand score files with empty cache (expect failure, since data needs to be pre-built or restored from BL)")
        rlscP = RcsbLigandScoreProvider(cachePath=self.__cachePath, useCache=False)
        ok = rlscP.testCache()
        self.assertFalse(ok)
        #
        logger.info("Testing reload of ligand score files from GitHub fallback (only do this for sake of tests, since official source should be BL)")
        rlscP.reload(useCache=False, useFallback=True)
        ok = rlscP.testCache()
        self.assertTrue(ok)
        #
        meanD, stdD, loadingD = rlscP.getParameterStatistics()
        for ky, val in meanD.items():
            logger.info("%-20s Mean %.4f stddev %.4f loading %.4f", ky, val, stdD[ky], loadingD[ky])
        #
        fitRank = rlscP.getFitScoreRanking(-1.9)
        geoRank = rlscP.getGeometryScoreRanking(-1.5)
        self.assertGreaterEqual(fitRank, 0.90)
        self.assertGreaterEqual(geoRank, 0.90)
        logger.info("fitRank (-1.9) %.4f  geoRank (-1.5) %.4f", fitRank, geoRank)


def fetchLigandScoreFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(RcsbLigandScoreProviderTests("testFetchScoreFiles"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = fetchLigandScoreFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
