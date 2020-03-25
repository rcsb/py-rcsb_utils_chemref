##
# File:    ChemCompProvider.py
# Author:  J. Westbrook
# Date:    22-Nov-2019
#
# Updates:
#
##
"""
Utilities to provide essential data items for chemical component definitions.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChemCompProvider(object):
    """Utilities to provide essential data items for chemical component definitions.
    """

    def __init__(self, **kwargs):
        urlTarget = kwargs.get("ccUrlTarget", "http://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz")
        # self.__birdUrlTarget = kwargs.get("birdUrlTarget", "ftp://ftp.wwpdb.org/pub/pdb/data/bird/prd/prdcc-all.cif.gz")
        #
        dirPath = os.path.join(kwargs.get("cachePath", "."), "chem_comp")
        useCache = kwargs.get("useCache", True)
        ccdFileName = kwargs.get("ccdFileName", "ccd_abbridged_definitions.json")
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__ccdD = self.__reload(urlTarget, dirPath, ccdFileName, useCache=useCache)

    def testCache(self):
        logger.info("Lengths map %d", len(self.__ccdD))
        if len(self.__ccdD) > 30000:
            return True
        return False

    def getParentComponent(self, ccId):
        try:
            tS = self.__ccdD[ccId]["mon_nstd_parent_comp_id"]
            return tS if tS and tS not in [".", "?"] else None
        except Exception:
            pass
        return None

    def getType(self, ccId):
        try:
            tS = self.__ccdD[ccId]["type"]
            return tS if tS and tS not in [".", "?"] else None
        except Exception:
            pass
        return None

    def getSubComponentList(self, ccId):
        try:
            tS = self.__ccdD[ccId]["pdbx_subcomponent_list"]
            return tS if tS and tS not in [".", "?"] else None
        except Exception:
            pass
        return None

    def getComponentIds(self):
        try:
            return list(self.__ccdD.keys())
        except Exception:
            pass
        return []

    def getAbbridged(self):
        return self.__ccdD

    def __reload(self, urlTarget, dirPath, ccdFileName, useCache=True):
        """Reload input mmCIF model mapping resource file and return a container list.

        Args:
            urlTarget (str): target url for resource file
            dirPath (str): path to the directory containing cache files
            ccdFileName (str): mapping file name
            useCache (bool, optional): flag to use cached files. Defaults to True.

        Returns:
            (dict): mapping dictionary (pdb_ccId -> CCDC details)
        """
        mD = {}
        #
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        filePath = os.path.join(dirPath, fn)
        ccdFilePath = os.path.join(dirPath, ccdFileName)
        self.__mU.mkdir(dirPath)
        #
        if not useCache:
            for fp in [filePath, ccdFilePath]:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and fU.exists(ccdFilePath):
            mD = self.__mU.doImport(ccdFilePath, fmt="json")
        else:
            ok = True
            if not (useCache and fU.exists(filePath)):
                logger.info("Fetching url %s for resource file %s", urlTarget, filePath)
                ok = fU.get(urlTarget, filePath)
            if ok:
                cL = self.__mU.doImport(filePath, fmt="mmcif")
                mD = self.__buildAbbridged(cL)
                ok = self.__mU.doExport(ccdFilePath, mD, fmt="json", indent=3)
        #
        return mD

    def __buildAbbridged(self, cL):
        """Return a dictionary of abbridged CCD info.
        """
        atNameList = [
            "id",
            "name",
            "type",
            "pdbx_type",
            "formula",
            "mon_nstd_parent_comp_id",
            "pdbx_synonyms",
            "pdbx_formal_charge",
            "pdbx_initial_date",
            "pdbx_modified_date",
            "pdbx_ambiguous_flag",
            "pdbx_release_status",
            "pdbx_replaced_by",
            "pdbx_replaces",
            "formula_weight",
            "one_letter_code",
            "three_letter_code",
        ]
        retD = {}
        for dataContainer in cL:
            tD = {}
            logger.debug("Processing model %r", dataContainer.getName())
            cObj = dataContainer.getObj("chem_comp")
            ccId = cObj.getValue("id", 0)
            for atName in atNameList:
                tD[atName] = cObj.getValueOrDefault(atName, 0, defaultValue=None)
            retD[ccId] = tD
        return retD
