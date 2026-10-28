##
# File:    BirdProviderTests.py
# Author:  J. Westbrook
# Date:    1-Jun-2021
# Version: 0.001
#
# Update:
# 28-Oct-2026 dwp  Add URL selection (current and beta archive layouts) and URL availability tests
##
"""
Tests for utilities to read and process BIRD dictionary definitions.

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

from importlib.metadata import version as get_package_version
from rcsb.utils.chemref.BirdProvider import BirdProvider
from rcsb.utils.io.FileUtil import FileUtil

__version__ = get_package_version("rcsb.utils.chemref")

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class BirdProviderTests(unittest.TestCase):
    def setUp(self):
        self.__startTime = time.time()
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__cachePath = os.path.join(HERE, "test-output")
        #
        logger.debug("Running tests on version %s", __version__)
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testBuildBirdCache(self):
        """Test build Bird definition cache"""
        bP = BirdProvider(cachePath=self.__cachePath, useCache=False, molLimit=None)
        ok = bP.testCache(minCount=1000)
        self.assertTrue(ok)
        bD = bP.getBirdD()
        for ky in ["releaseStatus", "representAs", "chemCompId", "class"]:
            self.assertTrue(ky in bD["PRD_000288"])
        self.assertEqual(bP.getReleaseStatus("PRD_000288"), "REL")
        self.assertEqual(bP.getRepresentation("PRD_000288"), "single molecule")
        self.assertEqual(bP.getChemCompId("PRD_000288"), "0GJ")
        self.assertEqual(bP.getClassType("PRD_000288"), "Inhibitor")
        bP = BirdProvider(cachePath=self.__cachePath, useCache=True, molLimit=None)
        ok = bP.testCache(minCount=1000)
        self.assertTrue(ok)

    def testBuildUrlTarget(self):
        """Test URL selection for the current (.../pub) and beta archive layouts"""
        curUrl = "https://files.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz"
        betaUrl = "https://files-beta.wwpdb.org/pub/wwpdb/refdata/bird/prd/prd-all.cif.gz"
        for basePdbRepoUrl, urlTarget in [
            (None, curUrl),
            ("", curUrl),
            ("https://files.wwpdb.org/pub", curUrl),
            ("https://files.wwpdb.org/pub/", curUrl),
            ("http://mirror.example.org/pub", "http://mirror.example.org/pub/pdb/data/bird/prd/prd-all.cif.gz"),
            ("https://files-beta.wwpdb.org", betaUrl),
            ("https://files-beta.wwpdb.org/", betaUrl),
            ("https://files-beta.wwpdb.org/pub", betaUrl),
            ("http://mirror-beta.example.org", "http://mirror-beta.example.org/pub/wwpdb/refdata/bird/prd/prd-all.cif.gz"),
            ("http://mirror-beta.example.org/pub/", "http://mirror-beta.example.org/pub/wwpdb/refdata/bird/prd/prd-all.cif.gz"),
            ("https://files.wwpdb.org", "https://files.wwpdb.org/pub/wwpdb/refdata/bird/prd/prd-all.cif.gz"),
        ]:
            self.assertEqual(BirdProvider.buildUrlTarget(basePdbRepoUrl), urlTarget, msg=repr(basePdbRepoUrl))

    def testUrlTargetsExist(self):
        """Test that the source file is available on both the current (files.wwpdb.org) and beta (files-beta.wwpdb.org) archives"""
        fU = FileUtil()
        for basePdbRepoUrl in ["https://files.wwpdb.org/pub", "https://files-beta.wwpdb.org"]:
            urlTarget = BirdProvider.buildUrlTarget(basePdbRepoUrl)
            logger.info("Checking %s", urlTarget)
            self.assertTrue(fU.exists(urlTarget), msg=urlTarget)


def buildCacheFiles():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(BirdProviderTests("testBuildBirdCache"))
    suiteSelect.addTest(BirdProviderTests("testBuildUrlTarget"))
    suiteSelect.addTest(BirdProviderTests("testUrlTargetsExist"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = buildCacheFiles()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
