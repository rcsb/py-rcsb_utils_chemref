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
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), 'rcsb', 'mock-data')
        self.__exampleFile = os.path.join(self.__dirPath, 'DrugBank', 'full_example.xml.gz')
        self.__fullDrugBankFile = os.path.join(self.__dirPath, 'DrugBank', 'full_database.xml.gz')
        #
        self.__fullDrugBankMappingFile = os.path.join(HERE, 'test-output', 'drugbank_pdb_mapping.json')

    def tearDown(self):
        pass

    def testReadDrugBankExample(self):
        dbu = DrugBankUtils()
        rL = dbu.read(self.__exampleFile)
        #
        logger.info("DrugBank example length %d" % len(rL))
        logger.info("DrugBank example keys %r" % rL[0].keys())
        logger.info("DrugBank example aliases %r" % rL[0]['aliases'])

    def testReadDrugBankFull(self):
        dbu = DrugBankUtils()
        dL = dbu.read(self.__fullDrugBankFile)
        logger.info("DrugBank full record length %d" % len(dL))
        mD = {}
        for d in dL:
            dbId = d['drugbank_id']
            pdbIds = ''
            # if 'pdb_entries' in d:
            #    pdbIds = d['pdb_entries']
            if 'external_identifiers' in d:
                for exD in d['external_identifiers']:
                    if exD['resource'] == 'PDB':
                        logger.debug("dbId %s pdbids %r ccids %r" % (dbId, pdbIds, exD['identifier']))
                        if exD['identifier'] not in mD:
                            mD[exD['identifier']] = []
                        mD[exD['identifier']] = {"drugbank_id": dbId, "aliases": list(d['aliases'])}

                        if 'target_interactions' in d:
                            for tid in d['target_interactions']:
                                tD = {}
                                tD['type'] = tid['category']
                                tD['name'] = tid['name']
                                tD['organism'] = tid['organism']
                                if len(tid['actions']):
                                    tD['actions'] = tid['actions']
                                if len(tid['known_action']):
                                    tD['known_action'] = tid['known_action']
                                if 'uniprot_ids' in tid:
                                    tD['uniprot_ids'] = tid['uniprot_ids']
                                #
                                if 'target_interactions' not in mD[exD['identifier']]:
                                    mD[exD['identifier']]['target_interactions'] = []
                                mD[exD['identifier']]['target_interactions'].append(tD)

        logger.info("Match length is %d" % len(mD))
        self.__serializeJson(self.__fullDrugBankMappingFile, mD)
        #
        #

    def __serializeJson(self, filePath, oD):
        with open(filePath, "w") as outfile:
            json.dump(oD, outfile, indent=3)


def readDrugBankInfo():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DrugBankUtilsTests("testReadDrugBankExample"))
    suiteSelect.addTest(DrugBankUtilsTests("testReadDrugBankFull"))
    return suiteSelect


if __name__ == '__main__':

    if True:
        mySuite = readDrugBankInfo()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
