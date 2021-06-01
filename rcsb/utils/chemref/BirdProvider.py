##
# File:    BirdProvider.py
# Author:  J. Westbrook
# Date:    1-Jun-2021
#
# Updates:
#
##
"""
Utilities to read and serialize portions of the dictionary of PDBx/mmCIF BIRD definitions.
"""
__docformat__ = "google.en"
__author__ = "John Westbrook"
__email__ = "john.westbrook@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.IoUtil import getObjSize
from rcsb.utils.io.MarshalUtil import MarshalUtil

# from rcsb.utils.io.SingletonClass import SingletonClass

logger = logging.getLogger(__name__)


class BirdProvider(object):
    """Utilities to read and serialize portions of the dictionary of PDBx/mmCIF BIRD definitions."""

    def __init__(self, **kwargs):
        # Default source target locators
        self.__birdUrlTarget = kwargs.get("birdUrlTarget", None)
        self.__birdUrlTarget = self.__birdUrlTarget if self.__birdUrlTarget else "http://ftp.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz"
        #
        cachePath = kwargs.get("cachePath", ".")
        dirPath = os.path.join(cachePath, "bird")
        useCache = kwargs.get("useCache", True)
        molLimit = kwargs.get("molLimit", 0)
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__birdD = self.__reload(self.__birdUrlTarget, dirPath, useCache=useCache, molLimit=molLimit)

    def testCache(self, minCount=None, logSizes=False):
        if logSizes and self.__birdD:
            logger.info("birdMolD object size %.2f MB", getObjSize(self.__birdD) / 1000000.0)
        ok = self.__birdD and len(self.__birdD) >= minCount if minCount else self.__birdD is not None
        return ok

    def getBirdD(self):
        return self.__birdD

    def getRepresentation(self, prdId):
        try:
            return self.__birdD[prdId]["representAs"]
        except Exception as e:
            logger.debug("Get representation %r failing with %s", prdId, str(e))
        return None

    def getReleaseStatus(self, prdId):
        try:
            return self.__birdD[prdId]["releaseStatus"]
        except Exception as e:
            logger.debug("Get status %r failing with %s", prdId, str(e))
        return None

    def getChemCompId(self, prdId):
        try:
            return self.__birdD[prdId]["chemCompId"]
        except Exception as e:
            logger.debug("Get compId %r failing with %s", prdId, str(e))
        return None

    def getClassType(self, prdId):
        try:
            return self.__birdD[prdId]["class"]
        except Exception as e:
            logger.debug("Get class %r failing with %s", prdId, str(e))
        return None

    def __reload(self, birdUrlTarget, dirPath, useCache=False, molLimit=None):
        """Reload or create serialized data dictionary of BIRD molecule definitions.

        Args:
            birdUrlTarget (str): target url for bird dictionary resource file
            dirPath (str): path to the directory containing cache files
            useCache (bool):
            molLimit (int): maximum number of definitions to process

         Returns:
            (list): list of bird data containers
        """
        #
        startTime = time.time()
        birdDataFilePath = os.path.join(dirPath, "bird-definitions.json")
        _, fExt = os.path.splitext(birdDataFilePath)
        birdDataFormat = "json" if fExt == ".json" else "pickle"
        #
        if useCache and self.__mU.exists(birdDataFilePath):
            birdObjD = self.__mU.doImport(birdDataFilePath, fmt=birdDataFormat)
        else:
            # Source component data files ...
            birdFilePath = self.__fetchUrl(birdUrlTarget, dirPath, useCache=useCache)
            birdObjD = self.__readBirdDefinitions(birdFilePath, molLimit=molLimit)
            ok = self.__mU.doExport(birdDataFilePath, birdObjD, fmt=birdDataFormat)
            logger.info("Storing %d definitions (status=%r) path: %s ", len(birdObjD), ok, birdDataFilePath)
        #
        endTime = time.time()
        logger.info("Loaded/reloaded %d definitions (%.4f seconds)", len(birdObjD), endTime - startTime)
        return birdObjD

    def __fetchUrl(self, urlTarget, dirPath, useCache=False):
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        filePath = os.path.join(dirPath, fn)
        if not (useCache and fU.exists(filePath)):
            startTime = time.time()
            ok2 = fU.get(urlTarget, filePath)
            endTime = time.time()
            if ok2:
                logger.info("Fetched %s for resource file %s (status = %r) (%.4f seconds)", urlTarget, filePath, ok2, endTime - startTime)
            else:
                logger.error("Failing fetch for %s for resource file %s (status = %r) (%.4f seconds)", urlTarget, filePath, ok2, endTime - startTime)
        #
        return filePath

    def __readBirdDefinitions(self, birdFilePath, molLimit=None):
        birdD = {}
        try:
            startTime = time.time()
            logger.info("Reading definitions from %s", birdFilePath)
            dataContainerL = self.__mU.doImport(birdFilePath, fmt="mmcif")
            endTime = time.time()
            logger.info("Read %s with %d BIRD definitions (%.4f seconds)", birdFilePath, len(dataContainerL), endTime - startTime)
            # -------
            startTime = time.time()
            dataContainerL = dataContainerL[:molLimit] if molLimit else dataContainerL
            for dataContainer in dataContainerL:
                prdId = prdReleaseStatus = representAs = None
                if dataContainer.exists("pdbx_reference_molecule"):
                    prdObj = dataContainer.getObj("pdbx_reference_molecule")
                    prdId = prdObj.getValueOrDefault("prd_id", 0, defaultValue=None)
                    prdReleaseStatus = prdObj.getValueOrDefault("release_status", 0, defaultValue=None)
                    representAs = prdObj.getValueOrDefault("represent_as", 0, defaultValue=None)
                    chemCompId = prdObj.getValueOrDefault("chem_comp_id", 0, defaultValue=None)
                    classType = prdObj.getValueOrDefault("class", 0, defaultValue=None)
                if prdId:
                    birdD[prdId] = {"releaseStatus": prdReleaseStatus, "representAs": representAs, "chemCompId": chemCompId, "class": classType}
            endTime = time.time()
            logger.info("Processed %d definitions (%.4f seconds)", len(birdD), endTime - startTime)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return birdD
