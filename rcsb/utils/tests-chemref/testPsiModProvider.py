##
# File:    PsiModProviderTests.py
# Author:  J. Westbrook
# Date:    10-Dec-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for various utilities for extracting data from PSI-MOD obo export files

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.PsiModProvider import PsiModProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PsiModProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testReloadPsiMod(self):
        """Test load from source"""
        try:
            psimodP = PsiModProvider(cachePath=self.__workPath, useCache=False)
            ok = psimodP.testCache()
            self.assertTrue(ok)
            version = psimodP.getVersion()
            self.assertTrue(version is not None)
            #
            rL = psimodP.getRootNodes()
            logger.info("root nodes %d %r", len(rL), rL)
            self.assertEqual(len(rL), 1)
            #
            modId = "MOD:00000"
            nL = psimodP.getSuccessors(modId)
            self.assertEqual(len(nL), 0)
            nL = psimodP.getPredecessors(modId)
            self.assertIsNotNone(nL)
            logger.info("root children %r", nL)
            ok = psimodP.exists(modId)
            self.assertTrue(ok)
            nm = psimodP.getName(modId)
            logger.info("name %r", nm)
            self.assertEqual(nm, "protein modification")
            #
            doTree = True
            if doTree:
                modIdL = [("MOD:00032", 1), ("MOD:00538", 1), ("MOD:01156", 1), ("MOD:01157", 1)]
                for modId, numParents in modIdL:
                    nm = psimodP.getName(modId)
                    self.assertIsNotNone(nm)
                    nL = psimodP.getAdjacentParents(modId)
                    self.assertIsNotNone(nL)
                    nL = psimodP.getPredecessors(modId)
                    self.assertIsNotNone(nL)
                    nL = psimodP.getSuccessors(modId)
                    self.assertIsNotNone(nL)
                    linL = psimodP.getDescendants(modId, includeSelf=True)
                    self.assertEqual(len(linL), numParents + 1)
                    logger.debug("%a Lineage(%d) %r", modId, len(linL), linL)

                trL = psimodP.exportTreeNodeList([tup[0] for tup in modIdL])
                logger.debug("trL %r", trL)
                logger.info("Length of tree node list %d", len(trL))
                self.assertGreaterEqual(len(trL), 5)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def readPsiMod():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PsiModProviderTests("testReloadPsiMod"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = readPsiMod()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
