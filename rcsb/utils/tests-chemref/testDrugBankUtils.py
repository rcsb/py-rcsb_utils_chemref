##
# File:    DrugBankUtilsTests.py
# Author:  J. Westbrook
# Date:    17-Oct-2018
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities extracting data from DrugBank repository exported data files.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest

from rcsb.utils.chemref.DrugBankUtils import DrugBankUtils
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class DrugBankUtilsTests(unittest.TestCase):
    def setUp(self):
        mockTopPath = os.path.join(TOPDIR, "mock-data")
        configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info"
        cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=mockTopPath)
        self.__user = cfgOb.getSecret("DRUGBANK_AUTH_USERNAME", sectionName=configName)
        self.__pw = cfgOb.getSecret("DRUGBANK_AUTH_PASSWORD", sectionName=configName)
        self.__workPath = os.path.join(HERE, "test-output")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testReadDrugBankInfo(self):
        dbu = DrugBankUtils(dirPath=self.__workPath, useCache=False, clearCache=True, username=self.__user, password=self.__pw)
        dbMapD = dbu.getMapping()
        logger.info("Mapping length %d", len(dbMapD))
        self.assertGreaterEqual(len(dbMapD["id_map"]), 5850)
        dbDocL = dbu.getDocuments()
        self.assertGreaterEqual(len(dbDocL), 5850)

    @unittest.skip("Long test")
    def testReReadDrugBankInfo(self):
        dbu = DrugBankUtils(dirPath=self.__workPath, useCache=False, clearCache=True, username=self.__user, password=self.__pw)
        dbMapD = dbu.getMapping()
        logger.info("Mapping length %d", len(dbMapD))
        self.assertGreaterEqual(len(dbMapD["id_map"]), 5850)
        dbDocL = dbu.getDocuments()
        self.assertGreaterEqual(len(dbDocL), 5850)
        #
        dbu = DrugBankUtils(dirPath=self.__workPath, useCache=True, clearCache=False)
        dbMapD = dbu.getMapping()
        logger.info("Mapping length %d", len(dbMapD))
        self.assertGreaterEqual(len(dbMapD["id_map"]), 5850)
        dbDocL = dbu.getDocuments()
        self.assertGreaterEqual(len(dbDocL), 5850)


def readDrugBankInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DrugBankUtilsTests("testReadDrugBankInfo"))
    suiteSelect.addTest(DrugBankUtilsTests("testReReadDrugBankInfo"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readDrugBankInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
