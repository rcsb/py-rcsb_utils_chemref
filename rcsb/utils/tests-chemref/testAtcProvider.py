##
# File:    AtcProviderTests.py
# Author:  J. Westbrook
# Date:    1-Nov-2018
# Version: 0.001
#
# Update:
#
#
##
"""
Various utilities for extracting data from the ATC hierarchy.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.AtcProvider import AtcProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class AtcProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testReadAtcInfo(self):
        atcP = AtcProvider(cachePath=self.__cachePath, useCache=False)
        version = atcP.getVersion()
        self.assertTrue(version is not None)
        tnL = atcP.getTreeNodeList()
        logger.info("length of tree list %d", len(tnL))
        self.assertGreaterEqual(len(tnL), 6000)
        idLin = atcP.getIdLineage("B05DB")
        logger.info("ID Lineage %r", idLin)
        atcName = atcP.getAtcName("B05DB")
        logger.info("ATC Name %r", atcName)
        nameLin = atcP.getNameLineage("B05DB")
        logger.info("Name Lineage %r", nameLin)
        # for tn in tnL:
        #    logger.info(">>>node %r", tn)


def readAtcInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(AtcProviderTests("testReadAtcInfo"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readAtcInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
