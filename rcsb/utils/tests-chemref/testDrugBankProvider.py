##
# File:    DrugBankProviderTests.py
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

from rcsb.utils.chemref.DrugBankProvider import DrugBankProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class DrugBankProviderTests(unittest.TestCase):
    def setUp(self):
        configPath = os.path.join(HERE, "test-data", "drugbank-config-example.yml")
        configName = "site_info_configuration"
        cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
        self.__user = cfgOb.get("_DRUGBANK_AUTH_USERNAME", sectionName=configName)
        self.__pw = cfgOb.get("_DRUGBANK_AUTH_PASSWORD", sectionName=configName)
        self.__cachePath = os.path.join(HERE, "test-output")
        #
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)\n", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testReadDrugBankInfo(self):
        dbu = DrugBankProvider(cachePath=self.__cachePath, useCache=False, username=self.__user, password=self.__pw)
        dbMapD = dbu.getMapping()
        version = dbu.getVersion()
        logger.info("Drugbank %r mapping length %d", version, len(dbMapD))
        self.assertGreaterEqual(len(dbMapD["id_map"]), 5850)
        dbDocL = dbu.getDocuments()
        self.assertGreaterEqual(len(dbDocL), 5850)

    @unittest.skip("Long test")
    def testReReadDrugBankInfo(self):
        dbu = DrugBankProvider(cachePath=self.__cachePath, useCache=False, username=self.__user, password=self.__pw)
        dbMapD = dbu.getMapping()
        logger.info("Mapping length %d", len(dbMapD))
        self.assertGreaterEqual(len(dbMapD["id_map"]), 5850)
        dbDocL = dbu.getDocuments()
        self.assertGreaterEqual(len(dbDocL), 5850)
        atcD = dbMapD["db_atc_map"]
        logger.info("atcD length %d", len(atcD))
        self.assertGreaterEqual(len(atcD), 3100)
        logger.info("atcD %r", atcD)
        #
        dbu = DrugBankProvider(cachePath=self.__cachePath, useCache=True)
        dbMapD = dbu.getMapping()
        logger.info("Mapping length %d", len(dbMapD))
        self.assertGreaterEqual(len(dbMapD["id_map"]), 5850)
        dbDocL = dbu.getDocuments()
        self.assertGreaterEqual(len(dbDocL), 5850)


def readDrugBankInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DrugBankProviderTests("testReadDrugBankInfo"))
    suiteSelect.addTest(DrugBankProviderTests("testReReadDrugBankInfo"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readDrugBankInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
