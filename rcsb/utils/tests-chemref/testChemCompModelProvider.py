##
# File:    ChemCompModelProviderTests.py
# Author:  J. Westbrook
# Date:    1-Nov-2018
# Version: 0.001
#
# Update:
# 28-Oct-2026 dwp  Add URL selection (current and beta archive layouts) and URL availability tests
#
##
"""
Various utilities for extracting data from CCDC correspondence data files.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.ChemCompModelProvider import ChemCompModelProvider
from rcsb.utils.io.FileUtil import FileUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ChemCompModelProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testReadChemCompModelRef(self):
        ccm = ChemCompModelProvider(cachePath=self.__workPath, useCache=False)
        ok = ccm.testCache()
        self.assertTrue(ok)
        rD = ccm.getMapping()
        self.assertGreater(len(rD), 1000)
        aD = ccm.getMapping()
        self.assertGreater(len(aD), 1000)
        #
        #
        logger.info("Model match length %d", len(rD))
        for ccId in rD:
            if len(rD[ccId]) > 1:
                logger.debug("Match length for %r %d", ccId, len(rD[ccId]))

    def testBuildUrlTarget(self):
        """Test URL selection for the current (.../pub) and beta archive layouts"""
        curUrl = "https://files.wwpdb.org/pub/pdb/data/component-models/complete/chem_comp_model.cif.gz"
        betaUrl = "https://files-beta.wwpdb.org/pub/wwpdb/pdb/data/component-models/complete/chem_comp_model.cif.gz"
        for basePdbRepoUrl, urlTarget in [
            (None, curUrl),
            ("", curUrl),
            ("https://files.wwpdb.org/pub", curUrl),
            ("https://files.wwpdb.org/pub/", curUrl),
            ("http://mirror.example.org/pub", "http://mirror.example.org/pub/pdb/data/component-models/complete/chem_comp_model.cif.gz"),
            ("https://files-beta.wwpdb.org", betaUrl),
            ("https://files-beta.wwpdb.org/", betaUrl),
            ("https://files-beta.wwpdb.org/pub", betaUrl),
            ("http://mirror-beta.example.org", "http://mirror-beta.example.org/pub/wwpdb/pdb/data/component-models/complete/chem_comp_model.cif.gz"),
            ("http://mirror-beta.example.org/pub/", "http://mirror-beta.example.org/pub/wwpdb/pdb/data/component-models/complete/chem_comp_model.cif.gz"),
            ("https://files.wwpdb.org", "https://files.wwpdb.org/pub/wwpdb/pdb/data/component-models/complete/chem_comp_model.cif.gz"),
        ]:
            self.assertEqual(ChemCompModelProvider.buildUrlTarget(basePdbRepoUrl), urlTarget, msg=repr(basePdbRepoUrl))

    def testUrlTargetsExist(self):
        """Test that the source file is available on both the current (files.wwpdb.org) and beta (files-beta.wwpdb.org) archives"""
        fU = FileUtil()
        for basePdbRepoUrl in ["https://files.wwpdb.org/pub", "https://files-beta.wwpdb.org"]:
            urlTarget = ChemCompModelProvider.buildUrlTarget(basePdbRepoUrl)
            logger.info("Checking %s", urlTarget)
            self.assertTrue(fU.exists(urlTarget), msg=urlTarget)


def readChemCompModelInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompModelProviderTests("testReadChemCompModelRef"))
    suiteSelect.addTest(ChemCompModelProviderTests("testBuildUrlTarget"))
    suiteSelect.addTest(ChemCompModelProviderTests("testUrlTargetsExist"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readChemCompModelInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
