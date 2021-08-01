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
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class ChemCompProvider(StashableBase):
    """Utilities to provide essential data items for chemical component definitions."""

    def __init__(self, **kwargs):
        dirName = "chem_comp"
        cachePath = kwargs.get("cachePath", ".")
        super(ChemCompProvider, self).__init__(cachePath, [dirName])
        #
        urlTarget = kwargs.get("ccUrlTarget", "http://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz")
        # self.__birdUrlTarget = kwargs.get("birdUrlTarget", "ftp://ftp.wwpdb.org/pub/pdb/data/bird/prd/prdcc-all.cif.gz")
        #
        dirPath = os.path.join(cachePath, dirName)
        useCache = kwargs.get("useCache", True)
        ccdFileName = kwargs.get("ccdFileName", "ccd_abbridged_definitions.json")
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        #
        relDataUrlTarget = kwargs.get("relDataUrlTarget", "https://github.com/rcsb/py-rcsb_exdb_assets/raw/development/fall_back/chem-comp-release-data-summary.json")
        self.__ccdRelD = self.__reloadRelData(relDataUrlTarget, dirPath, useCache=useCache)
        #
        self.__ccdD = self.__reload(urlTarget, dirPath, ccdFileName, useCache=useCache)
        #

    def testCache(self, minCount=30000):

        if self.__ccdD and self.__ccdRelD and len(self.__ccdD) > minCount and len(self.__ccdRelD) > minCount:
            logger.info("ChemComp details map length (%d) release data (%d)", len(self.__ccdD), len(self.__ccdRelD))
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

    def getAtomCountHeavy(self, ccId):
        try:
            return self.__ccdD[ccId]["atom_count_heavy"]
        except Exception:
            pass
        return 0

    def getAtomCount(self, ccId):
        try:
            return self.__ccdD[ccId]["atom_count"]
        except Exception:
            pass
        return 0

    def getAtomCountChiral(self, ccId):
        try:
            return self.__ccdD[ccId]["atom_count_chiral"]
        except Exception:
            pass
        return 0

    def getFormulaWeight(self, ccId):
        try:
            return float(self.__ccdD[ccId]["formula_weight"])
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

    def getReleaseDate(self, ccId):
        try:
            return self.__ccdRelD[ccId.upper()][1]
        except Exception:
            pass
        return None

    def getEntryOfFirstRelease(self, ccId):
        try:
            return self.__ccdRelD[ccId.upper()][0].lower()
        except Exception:
            pass
        return None

    def __reload(self, urlTarget, dirPath, ccdFileName, useCache=True):
        """Reload abbreviated chemical component details resource file.

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
        elif not useCache:
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

    def __reloadRelData(self, urlTarget, dirPath, useCache=True):
        """Reload chemical component release data details.

        Args:
            urlTarget (str): target url for resource file
            dirPath (str): path to the directory containing cache files
            ccdFileName (str): mapping file name
            useCache (bool, optional): flag to use cached files. Defaults to True.

        Returns:
            (dict): mapping dictionary (pdb_ccId -> CCDC details)
        """
        cD = {}
        #
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        filePath = os.path.join(dirPath, fn)
        self.__mU.mkdir(dirPath)
        #
        if not useCache:
            for fp in [filePath]:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and fU.exists(filePath):
            cD = self.__mU.doImport(filePath, fmt="json")
        elif not useCache:
            ok = True
            if not (useCache and fU.exists(filePath)):
                logger.info("Fetching url %s for resource file %s", urlTarget, filePath)
                ok = fU.get(urlTarget, filePath)
            if ok:
                cD = self.__mU.doImport(filePath, fmt="json")
        #
        return cD["release_data"] if "release_data" in cD else {}

    def __buildAbbridged(self, cL):
        """Return a dictionary of abbridged CCD info."""
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
            tD["formula_weight"] = float(tD["formula_weight"]) if tD["formula_weight"] and "formula_weight" in tD else None
            retD[ccId] = tD
            #
            #  - add some counts -
            #
            numAtoms = 0
            numAtomsHeavy = 0
            numAtomsChiral = 0
            try:
                cObj = dataContainer.getObj("chem_comp_atom")
                numAtoms = cObj.getRowCount()
                numAtomsHeavy = 0
                numAtomsChiral = 0
                for ii in range(numAtoms):
                    el = cObj.getValue("type_symbol", ii)
                    if el != "H":
                        numAtomsHeavy += 1
                    chFlag = cObj.getValue("pdbx_stereo_config", ii)
                    if chFlag != "N":
                        numAtomsChiral += 1
            except Exception:
                logger.warning("Missing chem_comp_atom category for %s", ccId)
                numAtoms = 0
                numAtomsHeavy = 0
                numAtomsChiral = 0
            #
            retD[ccId]["atom_count"] = numAtoms
            retD[ccId]["atom_count_chiral"] = numAtomsChiral
            retD[ccId]["atom_count_heavy"] = numAtomsHeavy
            #
        return retD
