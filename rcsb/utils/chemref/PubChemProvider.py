##
#  File:           PubChemProvider.py
#  Date:           27-Jul-2020 jdw
#
#  Updated:
#
##
"""
Accessors for PubChem mapped annotations.

"""

import logging
import os.path
import time

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashUtil import StashUtil

logger = logging.getLogger(__name__)


class PubChemProvider:
    """ Accessors for PubChem managed annotations.

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
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__pcD = self.__reload(fmt="json", useCache=useCache)
        #

    def testCache(self, minCount=0):
        if minCount == 0:
            return True
        if self.__pcD and minCount and ("identifiers" in self.__pcD) and len(self.__pcD["identifiers"]) >= minCount:
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
        fp = os.path.join(self.__dirPath, self.__stashDir, stashBaseFileName + fExt)
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
        """ Reload from the current cache file.
        """
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
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "pubchem_mapped_annotations")
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
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "pubchem_mapped_annotations")
            ok = stU.fetchBundle(self.__dirPath, url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok
