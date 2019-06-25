##
# File:    ChemCompModelUtilsTests.py
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

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import json
import logging
import os
import unittest

from rcsb.utils.chemref.ChemCompModelUtils import ChemCompModelUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ChemCompModelUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        self.__modelFile = os.path.join(self.__dirPath, "chem_comp_models", "chem_comp_model.cif.gz")
        #
        self.__modelMappingFile = os.path.join(HERE, "test-output", "ccdc_pdb_mapping.json")

    def tearDown(self):
        pass

    def testReadChemCompModelRef(self):
        ccm = ChemCompModelUtils()
        rD = ccm.getMapping(self.__modelFile)
        #
        logger.info("Model match length %d", len(rD))
        for ccId in rD:
            if len(rD[ccId]) > 1:
                logger.info("Match length for %r %d", ccId, len(rD[ccId]))

        self.__serializeJson(self.__modelMappingFile, rD)
        #
        #

    def __serializeJson(self, filePath, oD):
        with open(filePath, "w") as outfile:
            json.dump(oD, outfile, indent=3)


def readChemCompModelInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompModelUtilsTests("testReadChemCompModelRef"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readChemCompModelInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
