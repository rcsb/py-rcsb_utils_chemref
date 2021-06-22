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
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testReadAbbrevDrugBankInfo(self):
        urlTarget = os.path.join(HERE, "test-data", "full_database.zip")
        logger.info("Loading abbreviated Drugbank file %s", urlTarget)
        dbu = DrugBankProvider(urlTarget=urlTarget, cachePath=self.__cachePath, useCache=False, username=self.__user, password=self.__pw)
        dbMapD = dbu.getMapping()
        version = dbu.getVersion()
        logger.info("Drugbank %r mapping length %d", version, len(dbMapD))
        self.assertGreaterEqual(len(dbMapD["id_map"]), 340)
        dbDocL = dbu.getDocuments()
        self.assertGreaterEqual(len(dbDocL), 340)

    def testReadAbbrevDrugBankFeature(self):
        urlTarget = os.path.join(HERE, "test-data", "full_database.zip")
        logger.info("Loading abbreviated Drugbank file %s", urlTarget)
        dbu = DrugBankProvider(urlTarget=urlTarget, cachePath=self.__cachePath, useCache=False, username=self.__user, password=self.__pw)
        dbId = "DB00114"
        tS = dbu.getFeature(dbId, "smiles")
        logger.debug("%s SMILES is %r", dbId, tS)
        self.assertEqual(tS, "CC1=NC=C(COP(O)(O)=O)C(C=O)=C1O")
        fnL = ["name", "description", "indication", "pharmacodynamics", "mechanism-of-action"]
        for fn in fnL:
            fS = dbu.getFeature("DB00114", fn)
            logger.debug("%s %s is %r", dbId, fn, fS)
        # logger.info("%r %s", dbId, pprint.pformat(dbu.getD(dbId)))

    @unittest.skip("Long test")
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
