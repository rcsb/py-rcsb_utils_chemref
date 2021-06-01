##
# File: PharosReadSqlDump.py
# Date: 30-Mar-2020
#
# Adhoc class to extract and export a list ChEMBL identifers in selected tables from a dump of the
# Pharos relational database (TCRDv6.7).
##

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time

from rcsb.utils.io.MarshalUtil import MarshalUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class PharosReadSqlDump(object):
    """Adhoc class to extract and export a list ChEMBL identifers in selected tables from a dump of the
    Pharos relational database (TCRDv6.7).
    """

    def __init__(self):
        #
        self.__workPath = "."
        self.__mU = MarshalUtil(workPath=self.__workPath)
        #
        self.__startTime = time.time()
        logger.debug("Starting at %s", time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        #

    def readActivity(self, inpFileName="cmpd_activity.tdd"):
        chemblIdS = set()
        pubChemIdS = set()
        rowDictL = self.__mU.doImport(inpFileName, fmt="tdd", rowFormat="dict")
        logger.info("row 0 %r", rowDictL[0])
        for rowDict in rowDictL:
            if rowDict["catype"].upper() == "CHEMBL":
                chemblIdS.add(rowDict["cmpd_id_in_src"])
            pubChemIdS.add(rowDict["cmpd_id_in_src"])
        return chemblIdS, pubChemIdS

    def readDrug(self, inpFileName="drug_activity.tdd"):
        chemblIdS = set()
        rowDictL = self.__mU.doImport(inpFileName, fmt="tdd", rowFormat="dict")
        # logger.info("row 0 %r", rowDictL[0])
        for rowDict in rowDictL:
            if rowDict["cmpd_chemblid"] and rowDict["cmpd_chemblid"].startswith("CHEM"):
                chemblIdS.add(rowDict["cmpd_chemblid"])
        #
        return chemblIdS

    def listChEmbl(self):
        chS, pcS = self.readActivity()
        logger.info("ChEMBL (%d) PubChem (%d)", len(chS), len(pcS))
        #
        chDS = self.readDrug()
        logger.info("ChEMBL (%d) ", len(chDS))
        #
        logger.info("ChEMBL (%d) ", len(chS.union(chDS)))
        chL = list(chS.union(chDS))
        self.__mU.doExport("./Pharos-chembl-id.list", chL, fmt="list")


if __name__ == "__main__":
    rp = PharosReadSqlDump()
    rp.listChEmbl()
