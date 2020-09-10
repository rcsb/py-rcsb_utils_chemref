##
#  File:           PharosProvider.py
#  Date:           30-Jul-2020 jdw
#
#  Updated:
#
##
"""
Accessors for Pharos ChEMBL compound assignments.

"""

import logging
import os.path
import time

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashUtil import StashUtil

logger = logging.getLogger(__name__)


class PharosProvider:
    """Accessors for Pharos comound assignments."""

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirPath = os.path.join(cachePath, "Pharos")
        #
        #  - Configuration for stash services -
        #
        #    Local target directory name to be stashed.  (subdir of dirPath)
        self.__stashDir = "mapped_annotations"
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__phD = self.__reload(fmt="json", useCache=useCache)
        #

    def testCache(self, minCount=0):
        if minCount == 0:
            return True
        if self.__phD and minCount and ("identifiers" in self.__phD) and len(self.__phD["identifiers"]) >= minCount:
            return True
        return False

    def getIdentifiers(self):
        """Return a dictionary of related identifiers organized by chemical component/BIRD id.

        Returns:
            (dict): {phId: True}
        """
        try:
            return self.__phD["identifiers"] if self.__phD["identifiers"] else {}
        except Exception as e:
            logger.error("Failing with %r", str(e))
        return {}

    def __getAnnotFilePath(self, fmt="json"):
        stashBaseFileName = "pharos_mapped_annotations"
        fExt = ".json" if fmt == "json" else ".pic"
        fp = os.path.join(self.__dirPath, self.__stashDir, stashBaseFileName + fExt)
        return fp

    def load(self, chemblIdList, contentType, fmt="json", indent=0):
        """Load the input object of type contentType into the Pharos provider cache.

        Args:
            chemblIdList (list): Pharos ChEMBL ID assignments
            contentType (str): Pharos annotation content type (identifiers|...)
            fmt (str, optional): format string (json|pickle). Defaults to 'json'.

        Returns:
            (bool): True for success for False otherwise
        """
        ok = False
        try:
            if contentType in ["identifiers"]:
                phD = {chemblId: True for chemblId in chemblIdList}
                tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
                self.__phD = {"version": self.__version, "created": tS, contentType: phD}
                #
                annotFilePath = self.__getAnnotFilePath(fmt=fmt)
                kwargs = {"indent": indent} if fmt == "json" else {}
                ok = self.__mU.doExport(annotFilePath, self.__phD, fmt=fmt, **kwargs)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reload(self):
        """Reload from the current cache file."""
        ok = False
        try:
            self.__phD = self.__reload(fmt="json", useCache=True)
            ok = self.__phD is not None
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

    def toStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Copy tar and gzipped bundled cache data to remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "pharos_mapped_annotations")
            ok = stU.makeBundle(self.__dirPath, [self.__stashDir])
            if ok:
                ok = stU.storeBundle(url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok

    def fromStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Restore local cache from a tar and gzipped bundle to fetched from a remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        try:
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "pharos_mapped_annotations")
            ok = stU.fetchBundle(self.__dirPath, url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok
