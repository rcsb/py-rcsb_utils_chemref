##
#  File:           DrugCentralProvider.py
#  Date:           11-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for DrugCentral small molecule data.

"""

import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class DrugCentralProvider:
    """Accessors for DrugCentral small molecule data."""

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirPath = os.path.join(cachePath, "CACHE", "DrugCentral")
        drugCentralUrl = "https://unmtid-shinyapps.net/download/DrugCentral/20200516/structures.smiles.tsv"
        # drugCentralUrl = "http://unmtid-shinyapps.net/download/structures.smiles.tsv"
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__retD = self.__reload(drugCentralUrl, self.__dirPath, useCache=useCache)
        #

    def testCache(self):
        logger.info("Length %d", len(self.__retD))
        return self.__retD and len(self.__retD) > 3900

    def getInChIKey(self, drugCentralId):
        rVal = None
        try:
            rVal = self.__retD[drugCentralId]["inchikey"]
        except Exception:
            pass
        return rVal

    def getSMILES(self, drugCentralId):
        rVal = None
        try:
            rVal = self.__retD[drugCentralId]["smiles"]
        except Exception:
            pass
        return rVal

    def getCAS(self, drugCentralId):
        rVal = None
        try:
            rVal = self.__retD[drugCentralId]["cas"]
        except Exception:
            pass
        return rVal

    def getINN(self, drugCentralId):
        rVal = None
        try:
            rVal = self.__retD[drugCentralId]["inn"]
        except Exception:
            pass
        return rVal

    def __reload(self, drugCentralUrl, dirPath, useCache=True, fmt="json"):
        retD = {}
        startTime = time.time()
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        drugCentralDataPath = os.path.join(dirPath, "drugcentral_data." + fmt)
        if useCache and self.__mU.exists(drugCentralDataPath):
            logger.info("useCache %r using %r", useCache, drugCentralDataPath)
            retD = self.__mU.doImport(drugCentralDataPath, fmt=fmt)
            ok = True
        else:
            #
            inpFileName = fU.getFileName(drugCentralUrl)
            inpPath = os.path.join(dirPath, inpFileName)
            logger.info("Fetching url %s path %s", drugCentralUrl, inpPath)
            ok = fU.get(drugCentralUrl, inpPath)
            logger.info("Completed fetches at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            if ok:
                retD = self.__drugCentralParse(inpPath)
                self.__mU.doExport(drugCentralDataPath, retD, fmt=fmt, indent=3)
        return retD

    def __drugCentralParse(self, filePath):
        #     Fields:   SMILES  InChI   InChIKey        ID      INN     CAS_RN
        retD = {}
        rowL = self.__mU.doImport(filePath, fmt="tdd", rowFormat="list")
        logger.info("DrugCentral data length %d", len(rowL))
        #
        for row in rowL[1:]:
            lenR = len(row)
            if lenR < 3:
                logger.debug("Bad row %r", row)
                continue
            casRn = row[5] if lenR > 5 else None
            inn = row[4] if lenR > 4 else None
            retD[row[3]] = {"smiles": row[0], "inchiKey": row[2], "cas": casRn, "inn": inn}
        return retD
