##
#  File:           PubChemProvider.py
#  Date:           27-Jul-2020 jdw
#
#  Updated:
# 21-Jul-2021 jdw  Make this provider a subclass of StashableBase
#
##
"""
Accessors for PubChem mapped annotations.

"""

import datetime
import logging
import os.path
import time

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class PubChemProvider(StashableBase):
    """Accessors for PubChem managed annotations.

    dirPath -> CACHE/PubChem/

                             mapped_annotations/pubchem_mapped_annotations.json

                             stash/pubchem_mapped_annotations.tar.gz

    """

    def __init__(self, **kwargs):
        dirName = "PubChem-mapping"
        cachePath = kwargs.get("cachePath", ".")
        super(PubChemProvider, self).__init__(cachePath, [dirName])
        self.__dirPath = os.path.join(cachePath, dirName)
        #
        self.__version = datetime.datetime.now().strftime("%Y-%m-%d")
        useCache = kwargs.get("useCache", True)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__pcD = self.__reload(fmt="json", useCache=useCache)
        #

    def testCache(self, minCount=1):
        if minCount == 0 or minCount is None:
            return True
        if self.__pcD and ("identifiers" in self.__pcD) and len(self.__pcD["identifiers"]) >= minCount:
            logger.info("PubChem identifiers (%d)", len(self.__pcD["identifiers"]))
            return True
        return False

    def getIdentifiers(self):
        """Return a dictionary of related identifiers organized by chemical component/BIRD id.

        Returns:
            (dict): {ccId: {'idtype1': ids, 'idtype1': ids}, ... }
        """
        try:
            return self.__pcD["identifiers"] if self.__pcD["identifiers"] else {}
        except Exception as e:
            logger.error("Failing with %r", str(e))
        return {}

    def __getAnnotFilePath(self, fmt="json"):
        stashBaseFileName = "pubchem_mapped_annotations"
        fExt = ".json" if fmt == "json" else ".pic"
        fp = os.path.join(self.__dirPath, stashBaseFileName + fExt)
        return fp

    def load(self, pcObj, contentType, fmt="json", indent=0):
        """Load the input object of type contentType into the PubChem provider cache.

        Args:
            pcObj (dict): PubChem annotation object.
            contentType (str): PubChem annotation content type (identifiers|...)
            fmt (str, optional): format string (json|pickle). Defaults to 'json'.

        Returns:
            (bool): True for success for False otherwise
        """
        ok = False
        try:
            if contentType in ["identifiers"]:
                tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
                self.__pcD = {"version": self.__version, "created": tS, contentType: pcObj}
                #
                annotFilePath = self.__getAnnotFilePath(fmt=fmt)
                kwargs = {"indent": indent} if fmt == "json" else {}
                ok = self.__mU.doExport(annotFilePath, self.__pcD, fmt=fmt, **kwargs)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reload(self):
        """Reload from the current cache file."""
        ok = False
        try:
            self.__pcD = self.__reload(fmt="json", useCache=True)
            ok = self.__pcD is not None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reload(self, fmt="json", useCache=True):
        annotFilePath = self.__getAnnotFilePath(fmt=fmt)
        #
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        pcD = {"version": self.__version, "created": tS, "identifiers": {}}
        logger.info("useCache %r annotFilePath %r", useCache, annotFilePath)
        if useCache and self.__mU.exists(annotFilePath):
            pcD = self.__mU.doImport(annotFilePath, fmt=fmt)
        return pcD
