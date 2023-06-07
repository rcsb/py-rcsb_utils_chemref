##
# File:    ChemCompProviderTests.py
# Author:  J. Westbrook
# Date:    1-Nov-2018
# Version: 0.001
#
# Updates:
# 28-Mar-2022 bv Fix pylint issue with "Iterated dict modified inside for loop body" in testChemCompProvider
#  7-Jun-2022 aae Include bond count in chem comp data and fix typo
#
##
"""
Various utilities for extracting data from chemical component definitions.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import datetime
import logging
import os
import unittest
import copy

from rcsb.utils.chemref.ChemCompProvider import ChemCompProvider
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.TimeUtil import TimeUtil


HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class ChemCompProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__dataPath = os.path.join(HERE, "test-data")

    def tearDown(self):
        pass

    @unittest.skip("maintenance test")
    def testProcessReleaseData(self):
        """Test process raw chemical component release data extracted from released and obsolete entries"""
        try:
            fpOut = os.path.join(self.__workPath, "chem-comp-release-data-summary.json")
            fpRel = os.path.join(self.__dataPath, "chem-comp-release-data.json.gz")
            fpObs = os.path.join(self.__dataPath, "pdbx_obsolete-chem-comp-release-data.json")
            tU = TimeUtil()
            mU = MarshalUtil(workPath=self.__workPath)
            relTupL = mU.doImport(fpRel, fmt="json")
            obsTupL = mU.doImport(fpObs, fmt="json")
            logger.info("rel length (%d) obs length (%d)", len(relTupL), len(obsTupL))
            cD = {}
            for entryId, tS, ccIdL in obsTupL + relTupL:
                if len(tS) < 4:
                    continue
                dtObj = tU.getDateTimeObj(tS)
                for ccId in ccIdL:
                    if ccId not in cD:
                        cD[ccId] = (entryId, dtObj)
                    elif dtObj < cD[ccId][1]:
                        cD[ccId] = (entryId, dtObj)
            #
            logger.info("ALA %r %s", cD["ALA"][0], cD["ALA"][1])
            logger.info("cD (%d)", len(cD))
            cloneD = copy.deepcopy(cD)
            for ccId in cloneD:
                cD[ccId] = (cloneD[ccId][0], cloneD[ccId][1].strftime("%Y-%m-%d"))

            tS = datetime.datetime.now().isoformat()
            vS = datetime.datetime.now().strftime("%Y-%m-%d")
            ok = mU.doExport(fpOut, {"version": vS, "created": tS, "release_data": cD}, fmt="json", indent=3)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

    def testReadChemCompRef(self):
        ccP = ChemCompProvider(cachePath=self.__workPath, useCache=False)
        ok = ccP.testCache()
        self.assertTrue(ok)
        #
        rD = ccP.getAbridged()
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
        numBonds = ccP.getBondCount("ATP")
        logger.debug("%r %r %r %r", numAtoms, numChiralAtoms, numHeavyAtoms, numBonds)
        self.assertEqual(numAtoms, 47)
        self.assertEqual(numChiralAtoms, 6)
        self.assertEqual(numHeavyAtoms, 31)
        self.assertEqual(numBonds, 49)

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
        #
        tS = ccP.getReleaseDate("ALA")
        self.assertEqual(tS, "1973-05-03")
        fEntryId = ccP.getEntryOfFirstRelease("ALA")
        self.assertEqual(fEntryId, "2pti")


def readChemCompInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompProviderTests("testReadChemCompRef"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readChemCompInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
