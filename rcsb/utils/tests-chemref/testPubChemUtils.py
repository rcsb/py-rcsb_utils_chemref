##
# File:    PubChemUtilsTests.py
# Author:  J. Westbrook
# Date:    1-May-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Test access to PubChem fetch and search utilities.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.chemref.PubChemUtils import PubChemUtils, ChemicalIdentifier

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class PubChemUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output", "PubChem")

    def tearDown(self):
        pass

    def testFetchCompoundRecord(self):
        cId = "2244"
        cName = "2-acetyloxybenzoic acid"
        #
        pcU = PubChemUtils()
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier=cId, identifierType="cid")
        rawPath = os.path.join(self.__workPath, "%s-pubchem-record-raw.json" % cId)
        extractedPath = os.path.join(self.__workPath, "%s-pubchem-record-extracted.json" % cId)
        retStatus, rDL = pcU.fetch(chemId, returnType="record", storeRawResponsePath=rawPath, storeResponsePath=extractedPath)
        self.assertTrue(retStatus)
        ok = self.__containsCid(rDL, cId)
        self.assertTrue(ok)
        #
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier=cName, identifierType="name")
        retStatus, rDL = pcU.fetch(chemId, returnType="record")
        self.assertTrue(retStatus)
        ok = self.__containsCid(rDL, cId)
        self.assertTrue(ok)
        #

    def __containsCid(self, rDL, cid):
        try:
            for rD in rDL:
                if "cid" in rD and rD["cid"] == cid:
                    return True
            return False
        except Exception:
            pass
        return False

    def testSearchCompoundRecord(self):
        pcU = PubChemUtils()
        cId = "2244"
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="CC(=O)OC1=CC=CC=C1C(=O)O", identifierType="smiles", identifierSource="3d-model")
        retStatus, rDL = pcU.fetch(chemId, searchType="fastidentity", returnType="record")
        self.assertTrue(retStatus)
        ok = self.__containsCid(rDL, cId)
        self.assertTrue(ok)
        #
        chemId = ChemicalIdentifier(
            idCode="aspirin|test", identifier="InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)", identifierType="inchi", identifierSource="3d-model"
        )
        retStatus, rDL = pcU.fetch(chemId, searchType="fastidentity")
        self.assertTrue(retStatus)
        ok = self.__containsCid(rDL, cId)
        self.assertTrue(ok)
        #
        chemId = ChemicalIdentifier(idCode="aspirin|test", identifier="BSYNRYMUTXBXSQ-UHFFFAOYSA-N", identifierType="inchikey", identifierSource="3d-model")
        retStatus, rDL = pcU.fetch(chemId, searchType="lookup", returnType="record")
        self.assertTrue(retStatus)
        ok = self.__containsCid(rDL, cId)
        self.assertTrue(ok)
        #
        chemId = ChemicalIdentifier(idCode="Dummy", identifier="BSYNRYMUTXBXSQ-UHFFFAOYSA-Q", identifierType="inchikey", identifierSource="3d-model")
        retStatus, rDL = pcU.fetch(chemId, searchType="lookup", returnType="record")
        self.assertFalse(retStatus)
        self.assertFalse(rDL)

    # @unittest.skip("Skipping until api is more reliable")
    def testFetchCompoundAltReturnTypes(self):
        pcU = PubChemUtils()
        cIdList = ["123631", "2244"]
        for cId in cIdList:
            chemId = ChemicalIdentifier(idCode=cId, identifier=cId, identifierType="cid")
            for returnType, _ in [("classification", "Hierarchies"), ("property", "PropertyTable"), ("xrefs", "InformationList"), ("synonyms", "InformationList")]:
                rawResponsePath = os.path.join(self.__workPath, "%s-pubchem-%s-raw.json" % (cId, returnType))
                extractedResponsePath = os.path.join(self.__workPath, "%s-pubchem-%s-extract.json" % (cId, returnType))
                retStatus, rDL = pcU.fetch(chemId, returnType=returnType, storeRawResponsePath=rawResponsePath, storeResponsePath=extractedResponsePath)
                self.assertTrue(retStatus)
                ok = self.__containsCid(rDL, cId)
                self.assertTrue(ok)

    def testFetchCompoundView(self):
        cIdList = ["2244", "123631"]
        for cId in cIdList:
            pcU = PubChemUtils()
            chemId = ChemicalIdentifier(idCode="test", identifierType="cid", identifier=cId)
            rawResponsePath = os.path.join(self.__workPath, "%s-pubchem-view-raw.json" % cId)
            extractedResponsePath = os.path.join(self.__workPath, "%s-pubchem-view-extracted.json" % cId)
            retStatus, vL = pcU.fetch(chemId, returnType="view", storeRawResponsePath=rawResponsePath, storeResponsePath=extractedResponsePath)
            self.assertTrue(retStatus)
            self.assertGreaterEqual(len(vL), 1)

    # @unittest.skip("Skipping until api is more reliable")
    def testFetchCompoundExtTable(self):
        try:
            cIdList = ["2244", "123631"]
            extTable = cId = None
            for cId in cIdList:
                pcU = PubChemUtils()
                chemId = ChemicalIdentifier(idCode="test", identifierType="cid", identifier=cId)
                for extTable in ["dgidb", "pathway", "fdaorangebook", "clinicaltrials", "bioactivity"]:
                    rawResponsePath = os.path.join(self.__workPath, "%s-pubchem-%s-raw.json" % (cId, extTable))
                    extractedResponsePath = os.path.join(self.__workPath, "%s-pubchem-%s-extracted.json" % (cId, extTable))
                    retStatus, vL = pcU.fetch(chemId, returnType=extTable, storeRawResponsePath=rawResponsePath, storeResponsePath=extractedResponsePath)
                    self.assertTrue(retStatus)
                    self.assertGreaterEqual(len(vL), 1)
        except Exception as e:
            logger.exception("Failing with cId %r exTable %r %s", cId, extTable, str(e))
            self.fail()

    def testAssemble(self):
        cIdList = ["2244", "123631"]
        for cId in cIdList:
            pcU = PubChemUtils()
            chemId = ChemicalIdentifier(idCode=cId, identifierType="cid", identifier=cId)
            retStatus, retDL = pcU.assemble(chemId, exportPath=os.path.join(self.__workPath, "PubChem"))
            self.assertTrue(retStatus)
            self.assertTrue("record" in retDL[0]["data"])
            self.assertTrue("dgidb" in retDL[0]["data"])


def fetchPubChemData():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompoundRecord"))
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompoundView"))
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompoundAltReturnTypes"))
    suiteSelect.addTest(PubChemUtilsTests("testFetchCompoundExtTable"))
    suiteSelect.addTest(PubChemUtilsTests("testSearchCompoundRecord"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = fetchPubChemData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
