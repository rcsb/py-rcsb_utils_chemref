##
#  File:           CODProvider.py
#  Date:           8-Feb-2021 jdw
#
#  Updated:
#
##
"""
Accessors for Crystallographic Open Database (COD) molecule data.

"""

import logging
import os.path
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class CODProvider:
    """Accessors for COD molecule data."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "COD-molecules")
        self.__codCifDirPath = os.path.join(self.__dirPath, "CIF")
        self.__useCache = kwargs.get("useCache", True)
        codSmilesDumpUrl = kwargs.get("CODSmilesDumpUrl", "http://www.crystallography.net/cod/smi/allcod.smi")
        #
        self.__codCifTemplateUrl = "http://www.crystallography.net/cod/"
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        ok = self.__reload(self.__dirPath, codSmilesDumpUrl, self.__useCache)
        #
        logger.info("COD SMILES data status (%r)", ok)
        #

    def testCache(self, minCount=200000):
        _ = minCount
        return True

    def getSmilesPath(self):
        return os.path.join(self.__dirPath, "cod-molecules.smi")

    def __reload(self, dirPath, codSmilesDumpUrl, useCache):
        startTime = time.time()

        ok = False
        fU = FileUtil()
        codSmilesDumpFileName = "allcod.smi"
        codSmilesDumpPath = os.path.join(dirPath, codSmilesDumpFileName)
        #
        fU.mkdir(dirPath)
        codSmilesDataPath = self.getSmilesPath()
        #
        logger.info("useCache %r CODSmilesDumpPath %r", useCache, codSmilesDumpPath)
        if useCache and fU.exists(codSmilesDataPath):
            #
            ok = True
        else:
            logger.info("Fetching url %s path %s", codSmilesDumpUrl, codSmilesDumpPath)
            ok = fU.get(codSmilesDumpUrl, codSmilesDumpPath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            #
            numSmiles = self.__reformatCodSmilesData(codSmilesDumpPath, codSmilesDataPath)
            ok = numSmiles > 200000
        # ---
        return ok

    def __reformatCodSmilesData(self, inpFilePath, outFilePath):
        """Reformat COD SMILES data -

        Args:
            inpFilePath (str): input file path for COD molecule data
            outFilePath (str): reformatted file path for COD molecule data

        Returns:
            (int):  number of reformatted SMILES records
        """
        sCount = 0
        try:
            nl = "\n"
            with open(inpFilePath, "r", encoding="utf-8") as ifh, open(outFilePath, "w", encoding="utf-8") as ofh:
                for line in ifh.readlines():
                    ff = line[:-1].split("\t")
                    ff.reverse()
                    ofh.write("\t".join(ff) + nl)
                    sCount += 1
            #
            logger.info("Parsed COD molecule references (%d)", sCount)
        except Exception as e:
            logger.exception("Failing using %r with %s", inpFilePath, str(e))
        return sCount

    def getCodSmilesList(self):
        smiTupL = []
        inpFilePath = self.getSmilesPath()
        try:
            with open(inpFilePath, "r", encoding="utf-8") as ifh:
                for line in ifh.readlines():
                    ff = line[:-1].split("\t")
                    smiTupL.append((ff[0], ff[1]))
            return smiTupL
        except Exception as e:
            logger.exception("Failing for %r with %s", inpFilePath, str(e))

    def getCifPath(self, codId):
        return os.path.join(self.__codCifDirPath, codId[-1], codId + ".cif")

    def fetchCif(self, codId):
        ok = False
        try:
            startTime = time.time()
            fU = FileUtil()
            fU.mkdir(self.__codCifDirPath)
            codCifPath = self.getCifPath(codId)
            #
            logger.info("useCache %r codCifPath %r", self.__useCache, codCifPath)
            if self.__useCache and fU.exists(codCifPath):
                ok = True
            else:
                codCifUrl = os.path.join(self.__codCifTemplateUrl, codId + ".cif")
                logger.info("Fetching url %s path %s", codCifUrl, codCifPath)
                ok = fU.get(codCifUrl, codCifPath)
                logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        except Exception as e:
            logger.exception("Failing for %r with %s", codId, str(e))
            # ---
        return ok
