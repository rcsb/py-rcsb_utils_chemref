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

from rcsb.utils.chemref.PubChemUtils import PubChemUtils, ChemicalIdentifier
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PubChemUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")

    def tearDown(self):
        pass

    def testFetchCompound(self):
        pcU = PubChemUtils()
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="2244", identifierType="cid")
        retStatus, jD = pcU.fetch(chemId, returnType="record")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="2-acetyloxybenzoic acid", identifierType="name")
        retStatus, jD = pcU.fetch(chemId, returnType="record")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        #

    def testFetchDrugCompound(self):
        pcU = PubChemUtils()
        chemId = ChemicalIdentifier(idCode="test", identifierType="cid", identifier="123631")
        retStatus, vD = pcU.fetch(chemId, returnType="view")
        self.assertTrue(retStatus)
        mU = MarshalUtil()
        ok = mU.doExport(os.path.join(self.__workPath, "123631-pubchem-view.json"), vD, fmt="json", indent=3)
        # ok = pcU.traversePubChemCompoundView(vD)
        self.assertTrue(ok)
        pcD = pcU.parsePubChemCompoundView(vD)
        logger.debug("pcD %r", pcD)

    def __containsCid(self, jD, cid):
        try:
            if "PC_Compounds" in jD:
                for rD in jD["PC_Compounds"]:
                    ok = rD["id"]["id"]["cid"] = cid
                    if ok:
                        return ok
            return False
        except Exception:
            pass
        return False

    def testSearchCompound(self):
        pcU = PubChemUtils()
        cId = 2244
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="CC(=O)OC1=CC=CC=C1C(=O)O", identifierType="smiles", identifierSource="3d-model")
        retStatus, jD = pcU.fetch(chemId, searchType="fastidentity", returnType="record")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        ok = self.__containsCid(jD, cId)
        self.assertTrue(ok)

        #
        chemId = ChemicalIdentifier(
            idCode="aspirin|test", identifier="InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)", identifierType="inchi", identifierSource="3d-model"
        )
        retStatus, jD = pcU.fetch(chemId, searchType="fastidentity")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        ok = self.__containsCid(jD, cId)
        self.assertTrue(ok)
        #
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="BSYNRYMUTXBXSQ-UHFFFAOYSA-N", identifierType="inchikey", identifierSource="3d-model")
        retStatus, jD = pcU.fetch(chemId, searchType="lookup", returnType="record")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)
        ok = self.__containsCid(jD, cId)
        self.assertTrue(ok)

    def testFetchCompoundOtherReturnTypes(self):
        pcU = PubChemUtils()
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="2244", identifierType="cid")
        for returnType, kys in [("classification", "Hierarchies"), ("property", "PropertyTable"), ("xrefs", "InformationList"), ("synonyms", "InformationList")]:
            retStatus, jD = pcU.fetch(chemId, returnType=returnType)
            logger.debug("returnType %r keys %s\n", returnType, jD.keys())
            logger.debug(">>> jD = %r\n", jD)
            self.assertTrue(retStatus)
            self.assertEqual(kys, list(jD.keys())[0])

    def testFetchCompoundView(self):
        pcU = PubChemUtils()
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="2244", identifierType="cid")
        retStatus, jD = pcU.fetch(chemId, returnType="view")
        self.assertTrue(retStatus)
        logger.debug("jD = %r", jD)


def fetchPubChemData():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompound"))
    suiteSelect.addTest(PubChemUtilsTests("testFetchDrugCompound"))
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompoundView"))
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompoundOtherReturnTypes"))
    suiteSelect.addTest(PubChemUtilsTests("testSearchCompound"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = fetchPubChemData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
