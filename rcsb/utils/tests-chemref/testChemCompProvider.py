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

__docformat__ = "google.en"
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

    def testReadChemCompRef(self):
        ccP = ChemCompProvider(cachePath=self.__workPath, useCache=False)
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
        ccP = ChemCompProvider(cachePath=self.__workPath, useCache=True)
        ok = ccP.testCache()
        self.assertTrue(ok)
        numAtoms = ccP.getAtomCount("ATP")
        numChiralAtoms = ccP.getAtomCountChiral("ATP")
        numHeavyAtoms = ccP.getAtomCountHeavy("ATP")
        logger.debug("%r %r %r", numAtoms, numChiralAtoms, numHeavyAtoms)
        self.assertEqual(numAtoms, 47)
        self.assertEqual(numChiralAtoms, 6)
        self.assertEqual(numHeavyAtoms, 31)

        ccIdL = ccP.getComponentIds()
        iCount = 0
        fwCount = 0
        for ccId in ccIdL:
            if ccP.getAtomCountHeavy(ccId) < 1:
                logger.info("%s has no heavy atoms", ccId)
                iCount += 1
            if not ccP.getFormulaWeight(ccId):
                logger.info("%s has no formula weight", ccId)
                fwCount += 1
        self.assertGreaterEqual(3, iCount)
        self.assertGreaterEqual(3, fwCount)


def readChemCompInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompProviderTests("testReadChemCompRef"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readChemCompInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
