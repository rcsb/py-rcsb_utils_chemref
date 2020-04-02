##
# File:    PubChemUtilsTests.py
# Author:  J. Westbrook
# Date:    30-Mar-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Test access to PubChem fetch and search utilities.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.PubChemUtils import PubChemUtils
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class PubChemUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testFetchCompound(self):
        pcU = PubChemUtils()
        retStatus, jD = pcU.fetchCompound(identifier="2244", identifierType="cid")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        retStatus, jD = pcU.fetchCompound(identifier="2-acetyloxybenzoic acid", identifierType="name")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        #

    def testFetchDrugCompound(self):
        pcU = PubChemUtils()
        retStatus, vD = pcU.fetchCompound(identifier="123631", identifierType="cid", searchType="view")
        self.assertTrue(retStatus)
        mU = MarshalUtil()
        ok = mU.doExport(os.path.join(self.__workPath, "123631-pubchem-view.json"), vD, fmt="json", indent=3)
        ok = pcU.traversePubChemCompoundView(vD)
        self.assertTrue(ok)
        pcD = pcU.parsePubChemCompoundView(vD)
        logger.info("pcD %r", pcD)

    def testSearchCompound(self):
        pcU = PubChemUtils()
        retStatus, jD = pcU.fetchCompound(identifier="CC(=O)OC1=CC=CC=C1C(=O)O", identifierType="smiles", searchType="fastidentity")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        retStatus, jD = pcU.fetchCompound(identifier="InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)", identifierType="inchi", searchType="fastidentity")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        retStatus, jD = pcU.fetchCompound(identifier="BSYNRYMUTXBXSQ-UHFFFAOYSA-N", identifierType="inchikey")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)

    def testFetchCompoundXrefs(self):
        pcU = PubChemUtils()
        retStatus, jD = pcU.fetchCompound(identifier="2244", identifierType="cid", searchType="xrefs")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)

    def testFetchCompoundView(self):
        pcU = PubChemUtils()
        retStatus, jD = pcU.fetchCompound(identifier="2244", identifierType="cid", searchType="view")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)


def fetchPubChemData():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompound"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = fetchPubChemData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
