##
# File:    ChemRefDataPrepUtilsTests.py
# Author:  J. Westbrook
# Date:    4-Dec-2018
# Version: 0.001
#
# Updates:
#
##
"""
Tests for loading external chemical reference annotations.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"


import logging
import os
import time
import unittest

from rcsb.utils.chemref.ChemRefDataPrep import ChemRefDataPrep
from rcsb.utils.config.ConfigUtil import ConfigUtil

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))


class ChemRefDataPrepTests(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super(ChemRefDataPrepTests, self).__init__(methodName)
        self.__verbose = True

    def setUp(self):
        #
        #
        mockTopPath = os.path.join(TOPDIR, 'mock-data')
        logger.info("mockTopPath %s" % mockTopPath)
        configPath = os.path.join(mockTopPath, 'config', 'dbload-setup-example.yml')
        configName = 'site_info'
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=mockTopPath)
        # self.__cfgOb.dump()
        self.__resourceName = "MONGO_DB"
        self.__readBackCheck = True
        self.__numProc = 2
        self.__chunkSize = 10
        self.__documentLimit = 1000

        self.__exampleFilePath = os.path.join(mockTopPath, 'DrugBank', 'full_example.xml')
        self.__fullDrugBankFilePath = os.path.join(mockTopPath, 'DrugBank', 'full_database.xml.gz')
        #
        # self.__drugBankMappingFile = os.path.join(HERE, 'test-output', 'drugbank_pdb_mapping.json')
        # self.__drugBankDescriptorFile = os.path.join(HERE, 'test-output', 'drugbank_inchikey_mapping.json')
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s" % (self.id(),
                                            time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)\n" % (self.id(),
                                                              time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                              endTime - self.__startTime))

    def testGetDrugBankDataExample(self):
        extResource = "DrugBank"
        crl = ChemRefDataPrep(self.__cfgOb)
        dL = crl.fetchDocuments(extResource, self.__exampleFilePath)
        self.assertEqual(len(dL), 2)
        #
        oL = []
        for d in dL:
            oD = crl.buildDocument(extResource, d)
            logger.debug("loadable document %r" % (list(oD.items())))
            oL.append(oD)
        self.assertEqual(len(oL), 2)

    def testGetDrugBankDataFull(self):
        extResource = "DrugBank"
        oL = []
        crl = ChemRefDataPrep(self.__cfgOb)
        dL = crl.fetchDocuments(extResource, self.__fullDrugBankFilePath)
        for d in dL:
            oL.append(crl.buildDocument(extResource, d))
        #
        logger.debug("Processed document length %d" % len(oL))
        self.assertGreater(len(oL), 11000)

    def testGetDrugBankLoadableDocs(self):
        extResource = "DrugBank"
        crl = ChemRefDataPrep(self.__cfgOb)
        dL = crl.fetchDocuments(extResource, self.__exampleFilePath)
        self.assertEqual(len(dL), 2)
        exIdD = {}
        for d in dL:
            exIdD[d['drugbank_id']] = True
        dList = crl.getDocuments(extResource, exIdD)
        logger.info("document list length %d" % len(dList))


def accessionMappingSuite():
    suiteSelect = unittest.TestSuite()
    if True:
        suiteSelect.addTest(ChemRefDataPrepTests("testGetDrugBankDataExample"))
        suiteSelect.addTest(ChemRefDataPrepTests("testGetDrugBankDataFull"))
        suiteSelect.addTest(ChemRefDataPrepTests("testGetDrugBankLoadableDocs"))
    return suiteSelect


if __name__ == '__main__':
    #
    if (True):
        mySuite = accessionMappingSuite()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
