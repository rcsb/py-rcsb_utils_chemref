##
#  File:           PubChemProvider.py
#  Date:           27-Jul-2020 jdw
#
#  Updated:
#
##
"""
Accessors for PubChem extracted annotations.

"""

import logging
import os.path
import time

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashUtil import StashUtil


logger = logging.getLogger(__name__)


class PubChemProvider:
    """ Accessors for PubChem extracted annotations.

        dirPath -> CACHE/PubChem/

                                 mapped_annotations/pubchem_mapped_annotations.json

                                 stash/pubchem_mapped_annotations.tar.gz

    """

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirPath = os.path.join(cachePath, "PubChem")
        #
        #  - Configuration for stash services -
        #
        #    Local target directory name to be stashed.  (subdir of dirPath)
        self.__stashDir = "mapped_annotations"
        #
        #    Remote path to the directory containing the bundle file -
        # self.__stashRemoteDirPath = kwargs.get("stashRemoteDirPath", None)
        #
        #  Optional configuration for stash services -
        #
        # self.__stashUrl = kwargs.get("stashUrl", None)
        # self.__stashUserName = kwargs.get("stashUserName", None)
        # self.__stashPassword = kwargs.get("stashPassword", None)
        # self.__stashRemotePrefix = kwargs.get("stashRemotePrefix", None)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__pcD = self.__reload(fmt="json", useCache=useCache)
        #

    def testCache(self, minCount=1):
        if self.__pcD and minCount and ("identifiers" in self.__pcD) and len(self.__pcD["identifiers"]) >= minCount:
            return True
        return False

    def getIdentifiers(self):
        try:
            return self.__pcD["identifiers"]
        except Exception as e:
            logger.info("Failing with %r", str(e))
        return None

    def __getAnnotFilePath(self, fmt="json"):
        stashBaseFileName = "pubchem_mapped_annotations"
        fExt = ".json" if fmt == "json" else ".pic"
        fp = os.path.join(self.__dirPath, self.__stashDir, stashBaseFileName + fExt)
        return fp

    def load(self, pcObj, contentType, fmt="json", indent=0):
        """Load the input object into the PubChem provider cache.

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

    def __reload(self, fmt="json", useCache=True):
        annotFilePath = self.__getAnnotFilePath(fmt=fmt)
        #
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        pcD = {"version": self.__version, "created": tS, "identifiers": None}
        logger.info("useCache %r annotFilePath %r", useCache, annotFilePath)
        if useCache and self.__mU.exists(annotFilePath):
            pcD = self.__mU.doImport(annotFilePath, fmt=fmt)
        return pcD

    def toStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "pubchem_mapped_annotations")
            ok = stU.makeBundle(self.__dirPath, [self.__stashDir])
            if ok:
                ok = stU.storeBundle(url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.exception("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok

    def fromStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        try:
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "pubchem_mapped_annotations")
            ok = stU.fetchBundle(self.__dirPath, url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.exception("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok
