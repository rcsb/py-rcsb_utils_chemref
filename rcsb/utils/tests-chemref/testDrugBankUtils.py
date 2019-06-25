##
# File:    DrugBankUtilsTests.py
# Author:  J. Westbrook
# Date:    17-Oct-2018
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for utilities extracting data from DrugBank repository exported data files.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import json
import logging
import os
import unittest

from rcsb.utils.chemref.DrugBankUtils import DrugBankUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


class DrugBankUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        self.__exampleFile = os.path.join(self.__dirPath, "DrugBank", "full_example.xml")
        self.__fullDrugBankFile = os.path.join(self.__dirPath, "DrugBank", "full_database.xml.gz")
        #
        self.__drugBankMappingFile = os.path.join(HERE, "test-output", "drugbank_pdb_mapping.json")
        self.__drugBankDescriptorFile = os.path.join(HERE, "test-output", "drugbank_inchikey_mapping.json")

    def tearDown(self):
        pass

    def testReadDrugBankExample(self):
        dbu = DrugBankUtils()
        rL = dbu.read(self.__exampleFile)
        #
        logger.info("DrugBank example length %d", len(rL))
        logger.debug("DrugBank example keys %r", rL[0].keys())
        logger.debug("DrugBank example aliases %r", rL[0]["aliases"])
        self.assertTrue(len(rL), 2)

    def testMatchPdbDrugBank(self):
        dbu = DrugBankUtils()
        dL = dbu.read(self.__fullDrugBankFile)
        logger.info("DrugBank full record length %d", len(dL))
        dbD = {}
        mD = {}
        for dD in dL:
            dbId = dD["drugbank_id"]
            pdbIds = ""
            if "external_identifiers" in dD:
                for exD in dD["external_identifiers"]:
                    if exD["resource"] == "PDB":
                        logger.debug("dbId %s pdbids %r ccids %r", dbId, pdbIds, exD["identifier"])
                        if exD["identifier"] not in mD:
                            mD[exD["identifier"]] = []
                        mD[exD["identifier"]] = {"drugbank_id": dbId, "aliases": list(dD["aliases"])}
                        #
                        if "atc_codes" in dD and dD["atc_codes"]:
                            mD[exD["identifier"]]["atc_codes"] = dD["atc_codes"]

                        if "target_interactions" in dD:
                            for tid in dD["target_interactions"]:
                                tD = {}
                                tD["type"] = tid["category"]
                                tD["name"] = tid["name"]
                                tD["organism"] = tid["organism"]
                                if tid["actions"]:
                                    tD["actions"] = tid["actions"]
                                if tid["known_action"]:
                                    tD["known_action"] = tid["known_action"]
                                if "uniprot_ids" in tid:
                                    tD["uniprot_ids"] = tid["uniprot_ids"]
                                #
                                if "target_interactions" not in mD[exD["identifier"]]:
                                    mD[exD["identifier"]]["target_interactions"] = []
                                mD[exD["identifier"]]["target_interactions"].append(tD)
        logger.info("Match length is %d", len(mD))
        self.assertGreater(len(mD), 5000)
        dbD["id_map"] = mD
        #
        inD = {}
        for dD in dL:
            dbId = dD["drugbank_id"]
            if "inchikey" in dD and dD["inchikey"] and len(dD["inchikey"]) > 13:
                if dD["inchikey"] not in inD:
                    inD[dD["inchikey"]] = []
                inD[dD["inchikey"]].append({"drugbank_id": dbId, "inchikey": dD["inchikey"], "name": dD["name"]})
        #
        logger.info("Drugbank InChIKey dictionary length %d", len(inD))
        #
        self.assertGreater(len(inD), 9000)
        dbD["inchikey_map"] = inD
        self.__serializeJson(self.__drugBankMappingFile, dbD)
        #
        #

    def __serializeJson(self, filePath, oD):
        with open(filePath, "w") as outfile:
            json.dump(oD, outfile, indent=3)


def readDrugBankInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DrugBankUtilsTests("testReadDrugBankExample"))
    suiteSelect.addTest(DrugBankUtilsTests("testMatchPdbDrugBank"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = readDrugBankInfo()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
