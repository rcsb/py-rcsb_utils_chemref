##
#  File:           RcsbLigandScoreProvider.py
#  Date:           10-Feb-2021 jdw
#
#  Update2:
#  25-Jul-2022 dwp  Change location of fall back files to point to master (not development) branch on py-rcsb_exdb_assets
#  23-Feb-2026 dwp  Adjust provider to not use GitHub fallback files, and instead expect pre-built data from rcsb.workflow.refstats or BL
##
"""
Accessors for RCSB Ligand quality score supporting data.

"""

import bisect
import logging
import math
import os.path
import statistics
import time

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class RcsbLigandScoreProvider(StashableBase):
    """Accessors for RCSB Ligand quality score supporting data."""

    def __init__(self, **kwargs):
        self.__dirName = "rcsb-ligand-score"
        cachePath = kwargs.get("cachePath", ".")
        super().__init__(cachePath, [self.__dirName])
        self.__dirPath = os.path.join(cachePath, self.__dirName)
        #
        useCache = kwargs.get("useCache", True)
        useFallback = kwargs.get("useFallback", False)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__ligandScoreDL, self.__ligandExcludeD = self.__reload(self.__dirPath, useCache, useFallback=useFallback)
        #
        self.__meanD = {}
        self.__stdD = {}
        self.__loadingD = {}
        self.__geoScoreList = None
        self.__fitScoreList = None

    def getLigandScoreDataPath(self):
        """Return the path to final desired output file.

        Note that this is called by rcsb.workflow.refstats.LigandQualityReferenceGenerator,
        in order to support backup and restore functionalities to BL.
        """
        return os.path.join(self.__dirPath, "ligand_score_reference.csv")

    def testCache(self, minCount=200000):
        if self.__ligandScoreDL and self.__ligandExcludeD:
            logger.info("Ligand score (%d) exclude (%d)", len(self.__ligandScoreDL), len(self.__ligandExcludeD))
            if len(self.__ligandScoreDL) > minCount and len(self.__ligandExcludeD) > 0:
                return True
        return False

    def getLigandExcludeList(self):
        return list(self.__ligandExcludeD.keys())

    def isLigandExcluded(self, ccId):
        return ccId in self.__ligandExcludeD

    def getFitScoreRanking(self, score):
        try:
            if not self.__fitScoreList:
                self.__fitScoreList = sorted([float(tD["fit_pc1"]) for tD in self.__ligandScoreDL])
                logger.debug("Sorted model fit score (%d) range %.3f : %.3f", len(self.__fitScoreList), self.__fitScoreList[0], self.__fitScoreList[-1])
                logger.debug("Sorted model fit score (%d) range %.3f : %.3f", len(self.__fitScoreList), min(self.__fitScoreList), max(self.__fitScoreList))
            frac = bisect.bisect(self.__fitScoreList, score) / float(len(self.__fitScoreList) - 1)
            return 1.0 - frac
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return 0

    def getGeometryScoreRanking(self, score):
        try:
            if not self.__geoScoreList:
                self.__geoScoreList = sorted([float(tD["geo_pc1"]) for tD in self.__ligandScoreDL])
                logger.debug("Sorted model geometry score (%d) range %.3f : %.3f", len(self.__geoScoreList), self.__geoScoreList[0], self.__geoScoreList[-1])
                logger.debug("Sorted model geometry score (%d) range %.3f : %.3f", len(self.__geoScoreList), min(self.__geoScoreList), max(self.__geoScoreList))
            frac = bisect.bisect(self.__geoScoreList, score) / float(len(self.__geoScoreList) - 1)
            return 1.0 - frac
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return 0

    def reload(self, useCache=True, useFallback=False):
        self.__ligandScoreDL, self.__ligandExcludeD = self.__reload(self.__dirPath, useCache, useFallback=useFallback)

    def __reload(self, dirPath, useCache, useFallback=False):
        startTime = time.time()
        ligandScoreDL = []
        ligandExcludeD = {}
        #
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        ligandScoreFilePath = self.getLigandScoreDataPath()
        #
        if useCache and fU.exists(ligandScoreFilePath):
            ok = True
        #
        if not ok and useFallback:
            rcsbLigandScoreFallbackUrl = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/rcsb_ligand_score/ligand_score_reference.csv"
            logger.info("Fetching url %s path %s", rcsbLigandScoreFallbackUrl, ligandScoreFilePath)
            ok = fU.get(rcsbLigandScoreFallbackUrl, ligandScoreFilePath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        #
        if ok:
            ligandScoreDL = self.__mU.doImport(ligandScoreFilePath, fmt="csv", rowFormat="dict")
            ligandExcludeD = self.getLigandExcludeDict()
            logger.info("Reloaded ligand score list (%d) and exclude list (%d)", len(ligandScoreDL), len(ligandExcludeD))
            # ---
        logger.info("Completed reload (useCache %r) status %r at %s (%.4f seconds)", useCache, ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return ligandScoreDL, ligandExcludeD

    def getParameterStatistics(self):
        """Return the mean, standard deviation and paramter loadings for score model parameters

        Returns:
            (dict,dict,dict): mean, std. dev., and loading ({"rsr": v, "rscc": v, "mogul_bonds_rmsz": v, "mogul_angles_rmsz": v})
        """
        if not (self.__meanD and self.__stdD and self.__loadingD):
            self.__meanD, self.__stdD, self.__loadingD = self.__calcParameterStatistics()
        return self.__meanD, self.__stdD, self.__loadingD

    def __calcParameterStatistics(self):
        """Calculatefe the mean, standard deviation and paramter loadings for score model parameters

        Returns:
            (dict,dict,dict): mean, std. dev., and loading ({"rsr": v, "rscc": v, "mogul_bonds_rmsz": v, "mogul_angles_rmsz": v})
        """
        meanD = {}
        stdD = {}
        loadingD = {}
        try:
            #
            for ky in ["rsr", "rscc", "mogul_bonds_rmsz", "mogul_angles_rmsz"]:
                tL = [float(tD[ky]) for tD in self.__ligandScoreDL]
                meanD[ky] = statistics.mean(tL)
                stdD[ky] = statistics.stdev(tL)
                loadingD[ky] = math.sqrt(2.0) / 2.0 if ky != "rscc" else -math.sqrt(2.0) / 2.0
                #
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return meanD, stdD, loadingD

    def getLigandExcludeDict(self):
        ligExcludeL = [
            "N",
            "UNK",
            "NCO",
            "IRI",
            "SO4",
            "ZN",
            "MG",
            "CL",
            "CA",
            "NA",
            "PO4",
            "ACT",
            "MN",
            "K",
            "NI",
            "FE",
            "CU",
            "CD",
            "IOD",
            "FE2",
            "CO",
            "NO3",
            "HG",
            "FLC",
            "BR",
            "SCN",
            "CO3",
            "CAC",
            "BCT",
            "NH4",
            "CU1",
            "BA",
            "SR",
            "OH",
            "ALF",
            "NO2",
            "CS",
            "PT",
            "MLT",
            "OXL",
            "SO3",
            "VO4",
            "YB",
            "LI",
            "F",
            "RB",
            "OAA",
            "WO4",
            "PB",
            "SM",
            "PR",
            "YT3",
            "IUM",
            "RU",
            "TB",
            "PO3",
            "EMC",
            "3CO",
            "PD",
            "Y1",
            "OS",
            "AR",
            "LU",
            "IR",
            "EU",
            "EU3",
            "CR",
            "IR3",
            "2PO",
            "AUC",
            "LCP",
            "GA",
            "RE",
            "RH3",
            "3NI",
            "SE4",
            "PT4",
            "PBM",
            "AL",
            "D8U",
            "ER3",
            "RHD",
            "VN3",
            "RH",
            "SB",
            "TH",
            "4TI",
            "V",
            "BS3",
            "PTN",
            "ND4",
            "AM",
            "0BE",
            "TCN",
            "CF",
            "LCO",
            "CUL",
            "ZCM",
            "DY",
            "SFL",
            "PDV",
            "IN",
            "OS4",
            "4PU",
            "TA0",
            "YB2",
            "ZR",
            "GOL",
            "EDO",
            "PEG",
            "DMS",
            "ACE",
            "MPD",
            "MES",
            "TRS",
            "PG4",
            "PGE",
            "NH2",
            "FMT",
            "SF4",
            "EPE",
            "CIT",
            "BME",
            "ACY",
            "IMD",
            "1PE",
            "MLI",
            "FES",
            "UNX",
            "IPA",
            "MRD",
            "TLA",
            "UNL",
            "P6G",
            "POP",
            "F3S",
            "BTB",
            "EOH",
            "DTT",
            "NHE",
            "PYR",
            "BEF",
            "DIO",
            "MLA",
            "PGO",
            "2PE",
            "XE",
            "B3P",
            "PE4",
            "TAR",
            "PPV",
            "TAM",
            "PG0",
            "O",
            "CUA",
            "0QE",
            "P33",
            "AU",
            "AF3",
            "AG",
            "ARS",
            "BO3",
            "ICS",
            "CF0",
            "TBR",
            "AU3",
            "TAS",
            "BO4",
            "2T8",
            "AST",
            "ART",
            "BF2",
            "6BP",
            "D6N",
            "BF4",
            "ICE",
            "8AR",
            "0KA",
            "RMO",
            "HG2",
            "82N",
            "ICH",
            "ICZ",
            "ICG",
            "8P8",
            "NOB",
            "202",
            "HOH",
            "DOD",
            "1PG",
            "15P",
            "12P",
            "PG6",
            "PE5",
            "PE8",
            "7PE",
            "PE3",
            "PG5",
            "ETE",
            "XPE",
            "PEU",
            "P15",
            "7PG",
            "P4K",
            "33O",
            "9FO",
        ]
        #
        ligandExcludeD = {lig: True for lig in ligExcludeL}
        return ligandExcludeD
