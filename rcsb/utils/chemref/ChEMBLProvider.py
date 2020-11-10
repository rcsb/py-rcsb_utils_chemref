##
#  File:           ChEMBLProvider.py
#  Date:           9-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for ChEMBL small molecule data.

"""

import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChEMBLProvider:
    """Accessors for ChEMBL small molecule data."""

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirPath = os.path.join(cachePath, "CACHE", "ChEMBL")

        chemblDbUrl = kwargs.get("ChEMBLDbUrl", "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__retD = self.__reload(chemblDbUrl, self.__dirPath, useCache=useCache)
        #

    def testCache(self):
        return self.__retD and len(self.__retD) > 10000

    def getInChIKey(self, chemblId):
        rVal = None
        try:
            rVal = self.__retD[chemblId]["inchikey"]
        except Exception:
            pass
        return rVal

    def getSMILES(self, chemblId):
        rVal = None
        try:
            rVal = self.__retD[chemblId]["smiles"]
        except Exception:
            pass
        return rVal

    def __reload(self, chemblDbUrl, dirPath, useCache=True, fmt="json"):
        startTime = time.time()
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        baseVersion = 27
        # ChEMBL current version 27,...
        # template:  chembl_27_chemreps.txt.gz
        #
        inpFileName = "chembl_" + str(baseVersion) + ".fa.gz"
        #

        chemblDataPath = os.path.join(dirPath, "chembl_data." + fmt)

        #
        if useCache and self.__mU.exists(chemblDataPath):
            logger.info("useCache %r using %r", useCache, chemblDataPath)
            retD = self.__mU.doImport(chemblDataPath, fmt=fmt)
            ok = True
        else:
            #
            for vers in range(baseVersion, baseVersion + 10):
                inpFileName = "chembl_" + str(vers) + "_chemreps.txt.gz"
                inpPath = os.path.join(dirPath, inpFileName)
                url = os.path.join(chemblDbUrl, inpFileName)
                ok = fU.get(url, inpPath)
                logger.info("Fetching url %s path %s", url, inpPath)
                if ok:
                    break
            #
            logger.info("Completed fetches at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            retD = self.__chemblParse(inpPath)
            self.__mU.doExport(chemblDataPath, retD, fmt=fmt, indent=3)
        return retD

    def __chemblParse(self, filePath):
        retD = {}
        rowL = self.__mU.doImport(filePath, fmt="tdd", rowFormat="list")
        logger.info("ChEMBL data length %d", len(rowL))
        for row in rowL:
            if len(row) < 4:
                logger.debug("Bad row %r", row)
                continue
            retD[row[0]] = {"smiles": row[1], "inchiKey": row[3]}
        return retD
