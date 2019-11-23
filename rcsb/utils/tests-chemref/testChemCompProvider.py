##
# File:    ChemCompProviderTests.py
# Author:  J. Westbrook
# Date:    1-Nov-2018
# Version: 0.001
#
# Update:
#
#
##
"""
Various utilities for extracting data from chemical component definitions.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.ChemCompProvider import ChemCompProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ChemCompProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testReadChemComplRef(self):
        ccP = ChemCompProvider(dirPath=os.path.join(self.__workPath, "chem_ref"), useCache=False)
        ok = ccP.testCache()
        self.assertTrue(ok)
        #
        rD = ccP.getAbbridged()
        self.assertGreaterEqual(len(rD), 30000)
        #
        nsL = []
        for ccId in ccP.getComponentIds():
            tS = ccP.getParentComponent(ccId)
            if tS:
                nsL.append(tS)
        logger.info("nsL %r", len(nsL))
        self.assertGreaterEqual(len(nsL), 1550)
        #
        ccP = ChemCompProvider(dirPath=os.path.join(self.__workPath, "chem_ref"), useCache=True)
        ok = ccP.testCache()
        self.assertTrue(ok)


def readChemCompInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompProviderTests("testReadChemCompRef"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readChemCompInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
