##
# File:    ChemCompModelProviderTests.py
# Author:  J. Westbrook
# Date:    1-Nov-2018
# Version: 0.001
#
# Update:
#
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


def readChemCompModelInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompModelProviderTests("testReadChemCompModelRef"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readChemCompModelInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
