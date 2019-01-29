##
# File:    ChemCompDataPrepUtilsTests.py
# Author:  J. Westbrook
# Date:    29-Jan-2019
# Version: 0.001
#
# Updates:
#
##
"""
Tests for processing chemical componet reference definitions.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"


import logging
import os
import time
import unittest

from rcsb.utils.chemref.ChemCompDataPrep import ChemCompDataPrep
from rcsb.utils.config.ConfigUtil import ConfigUtil

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))


class ChemCompDataPrepTests(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super(ChemCompDataPrepTests, self).__init__(methodName)
        self.__verbose = True

    def setUp(self):
        #
        #
        mockTopPath = os.path.join(TOPDIR, 'mock-data')
        logger.info("mockTopPath %s" % mockTopPath)
        configPath = os.path.join(mockTopPath, 'config', 'dbload-setup-example.yml')
        configName = 'site_info'
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=mockTopPath)
        #
        self.__pathJsonContext = os.path.join(HERE, 'test-output', 'pdbx-contexts.json')

        self.__extendedChemCompFileNames = ['NAG-rev.cif']
        self.__workPath = os.path.join(HERE, 'test-output')
        self.__dataPath = os.path.join(HERE, 'test-data')
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s" % (self.id(),
                                            time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

    def tearDown(self):

        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)\n" % (self.id(),
                                                              time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                              endTime - self.__startTime))

    def testProcessContext(self):
        ccdp = ChemCompDataPrep(self.__cfgOb)
        catD, atD = ccdp.getContextInfo()
        logger.debug("category index %r" % catD.items())
        logger.debug("attribute index %r" % atD.items())
        self.assertTrue('WWPDB_LOCAL' in catD)
        self.assertTrue('WWPDB_LOCAL' in atD)

    def testExportContext(self):
        ccdp = ChemCompDataPrep(self.__cfgOb)
        ok = ccdp.exportContexts(self.__pathJsonContext)
        self.assertTrue(ok)

    def testFilterContext(self):
        ccdp = ChemCompDataPrep(self.__cfgOb)
        ok = ccdp.exportContexts(self.__pathJsonContext)
        self.assertTrue(ok)
        contextD = ccdp.importContexts(self.__pathJsonContext)
        for fn in self.__extendedChemCompFileNames:
            inpPath = os.path.join(self.__dataPath, fn)
            outPath = os.path.join(self.__workPath, fn)
            ccdp.filterByContext(inpPath, outPath, contextD)


def updateProcessingSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompDataPrepTests("testProcessContext"))
    suiteSelect.addTest(ChemCompDataPrepTests("testExportContext"))
    suiteSelect.addTest(ChemCompDataPrepTests("testFilterContext"))
    return suiteSelect


if __name__ == '__main__':
    #
    if (True):
        mySuite = updateProcessingSuite()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
