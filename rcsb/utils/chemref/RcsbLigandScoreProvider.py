##
#  File:           RcsbLigandScoreProvider.py
#  Date:           10-Feb-2021 jdw
#
#  Updated:
#
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
        dirName = "rcsb-ligand-score"
        self.__cachePath = kwargs.get("cachePath", ".")
        super(RcsbLigandScoreProvider, self).__init__(self.__cachePath, [dirName])
        #
        self.__dirPath = os.path.join(self.__cachePath, dirName)
        self.__useCache = kwargs.get("useCache", True)
        rcsbLigandScoreUrl = kwargs.get("rcsbLigandScoreUrl", "https://github.com/rcsb/py-rcsb_exdb_assets/raw/development/fall_back/rcsb_ligand_score/ligand_score_reference.csv")
        rcsbLigandExcludeUrl = kwargs.get("rcsbLigandExcludeUrl", "https://github.com/rcsb/py-rcsb_exdb_assets/raw/development/fall_back/rcsb_ligand_score/ligand_score_exclude.list")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__ligandScoreDL, self.__ligandExcludeD = self.__reload(self.__dirPath, rcsbLigandScoreUrl, rcsbLigandExcludeUrl, self.__useCache)
        #
        self.__meanD = {}
        self.__stdD = {}
        self.__loadingD = {}
        self.__geoScoreList = None
        self.__fitScoreList = None

    def testCache(self):
        if self.__ligandScoreDL and self.__ligandExcludeD:
            logger.info("Ligand score (%d) exclude (%d)", len(self.__ligandScoreDL), len(self.__ligandExcludeD))
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

    def __reload(self, dirPath, rcsbLigandScoreUrl, rcsbLigandExcludeUrl, useCache):
        startTime = time.time()
        ligandScoreDL = []
        ligandExcludeD = {}
        #
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        fn = os.path.basename(rcsbLigandScoreUrl)
        ligandScoreFilePath = os.path.join(dirPath, fn)
        fn = os.path.basename(rcsbLigandExcludeUrl)
        ligandExcludeFilePath = os.path.join(dirPath, fn)
        #
        if useCache and fU.exists(ligandScoreFilePath) and fU.exists(ligandExcludeFilePath):
            ok = True
        elif not useCache:
            logger.info("Fetching url %s path %s", rcsbLigandScoreUrl, ligandScoreFilePath)
            ok = fU.get(rcsbLigandScoreUrl, ligandScoreFilePath)
            logger.info("Fetching url %s path %s", rcsbLigandExcludeUrl, ligandExcludeFilePath)
            ok = fU.get(rcsbLigandExcludeUrl, ligandExcludeFilePath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            #
        if ok:
            ligandScoreDL = self.__mU.doImport(ligandScoreFilePath, fmt="csv", rowFormat="dict")
            ligExcludeL = self.__mU.doImport(ligandExcludeFilePath, fmt="list")
            ligandExcludeD = {lig: True for lig in ligExcludeL}
            # ---
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
