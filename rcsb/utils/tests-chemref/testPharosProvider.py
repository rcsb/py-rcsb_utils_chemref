##
# File:    PharosProviderTests.py
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

__docformat__ = "google.en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.PharosProvider import PharosProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PharosProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        self.__phL = self.__mU.doImport(os.path.join(self.__dataPath, "Pharos-chembl-id-20200730.list"), fmt="list")
        self.__stashUrl = None
        self.__stashRemotePath = os.path.join(self.__cachePath, "stash-remote")

    def tearDown(self):
        pass

    def testPharosAnnotBootstrap(self):
        minCount = 266600
        pcP = PharosProvider(cachePath=self.__cachePath, useCache=False)
        ok = pcP.load(self.__phL, "identifiers", fmt="json", indent=3)
        self.assertTrue(ok)
        riD = pcP.getIdentifiers()
        logger.info("riD (%d)", len(riD))
        ok = pcP.testCache(minCount=minCount)
        self.assertTrue(ok)
        ok = pcP.toStash(self.__stashUrl, self.__stashRemotePath)
        self.assertTrue(ok)
        ok = pcP.fromStash(self.__stashUrl, self.__stashRemotePath)
        self.assertTrue(ok)
        ok = pcP.reload()
        self.assertTrue(ok)
        ok = pcP.testCache(minCount=minCount)
        self.assertTrue(ok)
        #
        pcP = PharosProvider(cachePath=self.__cachePath, useCache=True)
        ok = pcP.testCache(minCount=minCount)
        self.assertTrue(ok)

    @unittest.skip("Private test")
    def testStashRemote(self):
        minCount = 266600
        configPath = os.path.join(self.__dataPath, "stash-config-example.yml")
        configName = "site_info_configuration"
        cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName)
        userName = cfgOb.get("_STASH_AUTH_USERNAME", sectionName=configName)
        password = cfgOb.get("_STASH_AUTH_PASSWORD", sectionName=configName)
        basePath = cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=configName)
        url = cfgOb.get("STASH_SERVER_URL", sectionName=configName)
        urlFallBack = cfgOb.get("STASH_SERVER_FALLBACK_URL", sectionName=configName)
        #
        pcP = PharosProvider(cachePath=self.__cachePath, useCache=True)
        ok = pcP.testCache(minCount=minCount)
        self.assertTrue(ok)
        ok = pcP.toStash(url, basePath, userName=userName, password=password)
        self.assertTrue(ok)
        pcP.toStash(urlFallBack, basePath, userName=userName, password=password)
        self.assertTrue(ok)
        #
        ok = pcP.fromStash(url, basePath, userName=userName, password=password)
        self.assertTrue(ok)
        ok = pcP.fromStash(urlFallBack, basePath, userName=userName, password=password)
        self.assertTrue(ok)
        ok = pcP.reload()
        self.assertTrue(ok)
        #


def pharosAnnotSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PharosProviderTests("testPharosAnnotBootstrap"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = pharosAnnotSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
