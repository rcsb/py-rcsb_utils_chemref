##
# File:    DrugCentralProviderTests.py
# Author:  J. Westbrook
# Date:    11-Nov-2020
#
# Update:
#
#
##
"""
Tests for accessors for DrugCentral small molecule data.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.DrugCentralProvider import DrugCentralProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class DrugCentralProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)

    def tearDown(self):
        pass

    def testFetchDrugCentralData(self):
        try:
            dcP = DrugCentralProvider(cachePath=self.__cachePath, useCache=False)
            ok = dcP.testCache()
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def drugCentralDataSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DrugCentralProviderTests("testFetchDrugCentralData"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = drugCentralDataSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
