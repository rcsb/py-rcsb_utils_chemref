##
# File:    PubChemProviderTests.py
# Author:  J. Westbrook
# Date:    27-Jul-2018
#
# Update:
#
#
##
"""
Tests for accessors for managing PubChem extracted annotations.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.PubChemProvider import PubChemProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PubChemProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        self.__pcAnnotD = self.__mU.doImport(os.path.join(self.__dataPath, "pubchem_mapped_annotations.json"), fmt="json")

    def tearDown(self):
        pass

    def testPubChemAnnotBootstrap(self):
        pcP = PubChemProvider(cachePath=self.__cachePath, useCache=False)
        ok = pcP.load(self.__pcAnnotD["identifiers"], "identifiers", fmt="json", indent=3)
        self.assertTrue(ok)
        riD = pcP.getIdentifiers()
        logger.info("riD (%d)", len(riD))
        ok = pcP.testCache(minCount=23)
        self.assertTrue(ok)
        pcP = PubChemProvider(cachePath=self.__cachePath, useCache=True)
        ok = pcP.testCache(minCount=23)
        self.assertTrue(ok)


def pubchemAnnotSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PubChemProviderTests("testPubChemAnnotBootstrap"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = pubchemAnnotSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
