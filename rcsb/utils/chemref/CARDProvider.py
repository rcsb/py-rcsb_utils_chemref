##
#  File:           CARDProvider.py
#  Date:           28-Nov-2020 jdw
#
#  Updated:
#
##
"""
Accessors for CARD molecule data.

"""

import logging
import os.path
import time

import networkx
import obonet

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class CARDProvider:
    """Accessors for CARD molecule data."""

    def __init__(self, **kwargs):
        #
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__dirPath = os.path.join(self.__cachePath, "CARD-molecules")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__oD, self.__version = self.__reload(self.__dirPath, **kwargs)
        #

    def testCache(self, minCount=350):
        if self.__oD and len(self.__oD) >= minCount:
            logger.info("CARD mapping length (%d)", len(self.__oD))
            return True
        else:
            return False

    def __reload(self, dirPath, **kwargs):
        oD = None
        startTime = time.time()
        useCache = kwargs.get("useCache", True)
        #
        cardDumpUrl = kwargs.get("CARDDumpUrl", "https://card.mcmaster.ca/latest/ontology")
        ok = False
        fU = FileUtil()
        cardDumpFileName = "card-ontology.tar.bz2"
        cardDumpPath = os.path.join(dirPath, cardDumpFileName)
        cardDumpDirPath = os.path.join(dirPath, "dump")
        #
        fU.mkdir(dirPath)
        cardDataPath = os.path.join(dirPath, "card-molecules.json")
        #
        logger.info("useCache %r CARDDumpPath %r", useCache, cardDumpPath)
        if useCache and self.__mU.exists(cardDataPath):
            oD = self.__mU.doImport(cardDataPath, fmt="json")
        else:
            logger.info("Fetching url %s path %s", cardDumpUrl, cardDumpPath)
            ok = fU.get(cardDumpUrl, cardDumpPath)
            fU.mkdir(cardDumpDirPath)
            fU.uncompress(cardDumpPath, outputDir=cardDumpDirPath)
            fU.unbundleTarfile(os.path.join(cardDumpDirPath, cardDumpFileName[:-4]), dirPath=cardDumpDirPath)
            logger.info("Completed fetch (%r) at %s (%.4f seconds)", ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            #
            oD = self.__parseCardData(os.path.join(cardDumpDirPath, "aro.obo"))
            ok = self.__mU.doExport(cardDataPath, oD, fmt="json", indent=3)
        # ---
        return oD

    def __parseCardData(self, filePath):
        """Parse CARD ontology file describing antibiotic molecules.

        Args:
            filePath (str): output file path for CARD molecule data

        Returns:
            (dict, version): dictionary of select CARD data, version string


                Note: https://card.mcmaster.ca/aro/0000003
        """
        oD = {}
        version = None
        try:
            cardGraph = obonet.read_obo(filePath)
            logger.info("CARD graph nodes (%d) edges (%d) directed %r", len(cardGraph), cardGraph.number_of_edges(), networkx.is_directed_acyclic_graph(cardGraph))
            for ky, gD in cardGraph.nodes(data=True):
                # logger.info("ky %r gD %r", ky, gD.keys())
                if "xref" in gD:
                    # oD["name"] = gD["name"]
                    # oD["pubChemCId"] = gD["xref"]
                    oD[ky] = {"name": gD["name"], "pubChemCId": gD["xref"], "description": gD["def"]}
            logger.info("Parsed CARD molecule references (%d)", len(oD))
        except Exception as e:
            logger.exception("Failing using %r with %s", filePath, str(e))
        return oD, version
