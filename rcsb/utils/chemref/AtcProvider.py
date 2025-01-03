##
#  File:           AtcProvider.py
#  Date:           3-Apr-2019 jdw
#
#  Updates:
#   5-Jun-2021 jdw  Update ATC data source and fallback
#  20-Jul-2021 jdw  Make this provider a subclass of StashableBase
#   3-Jan-2022 dwp  Update data loading methods to address recent changes in source NCBO ATC files
#  25-Jul-2022 dwp  Revert last change - source NCBO ATC files were updated again to restore previous format
#  27-Jul-2022 dwp  Adapt to possible naming schemes for class IDs (may be identified with "ATC" or "UATC")
#  22-Oct-2024 dwp  Add standardization method for ATC data to handle potential stringified arrays
##
"""
  Extract ATC term descriptions from NCBO ATC flat files.

"""

import collections
import logging
import os.path
import sys
import ast

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase


logger = logging.getLogger(__name__)


class AtcProvider(StashableBase):
    """Extract term descriptions and ATC classifications from ATC flat files."""

    def __init__(self, **kwargs):
        self.__dirName = "atc"
        self.__cachePath = kwargs.get("cachePath", ".")
        super(AtcProvider, self).__init__(self.__cachePath, [self.__dirName])
        atcDirPath = os.path.join(self.__cachePath, self.__dirName)
        useCache = kwargs.get("useCache", True)

        urlTarget = "https://data.bioontology.org/ontologies/ATC/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv"
        urlTargetFallback = "https://github.com/rcsb/py-rcsb_exdb_assets/raw/master/fall_back/ATC.csv.gz"
        self.__version = kwargs.get("AtcVersion", "2023")
        #
        self.__mU = MarshalUtil(workPath=atcDirPath)
        self.__atcD = self.__reload(urlTarget, urlTargetFallback, atcDirPath, useCache=useCache)
        #

    def getVersion(self):
        return self.__version

    def __reload(self, urlTarget, urlTargetFallback, atcDirPath, useCache=True):
        pyVersion = sys.version_info[0]
        atcFilePath = os.path.join(atcDirPath, "atc-py%s.pic" % str(pyVersion))
        #
        atcD = {}
        if useCache and self.__mU.exists(atcFilePath):
            atcD = self.__mU.doImport(atcFilePath, fmt="pickle")
            # logger.debug("ATC name length %d parent length %d assignments %d", len(atcD["names"]), len(atcD["parents"]), len(atcD["assignments"]))
            # nD = atcD["names"]
            # pD = atcD["parents"]

        elif not useCache:
            fn = "ATC.csv.gz"
            fp = os.path.join(atcDirPath, fn)
            logger.debug("Fetch ATC term descriptions from source %s", fp)
            fileU = FileUtil(workPath=atcDirPath)
            fileU.mkdir(atcDirPath)
            try:
                ok = fileU.get(urlTarget, fp)
                logger.info("ATC fetch status is %r", ok)
                if not ok:
                    ok = fileU.get(urlTargetFallback, fp)
                    logger.info("ATC fallback fetch status is %r", ok)
                atcL = self.__mU.doImport(fp, fmt="csv", rowFormat="dict", uncomment=False)
                #
                columnL = list(atcL[0].keys()) if atcL else []
                nD = self.__extractNames(atcL)
                pD = self.__extractHierarchy(atcL)
                atcD = {"names": nD, "parents": pD}
                ok = self.__mU.doExport(atcFilePath, atcD, fmt="pickle")
                logger.info("ATC cache status %r data length %d columns %r names %d parents %d", ok, len(atcL), columnL, len(nD), len(pD))
            except Exception as e:
                logger.exception("Failing for %r with %s", fp, str(e))
            #
        return atcD

    def getAtcVersion(self):
        return self.__version

    def testCache(self):
        if self.__atcD and "names" in self.__atcD and "parents" in self.__atcD and (len(self.__atcD["names"]) > 6100) and (len(self.__atcD["parents"]) > 6100):
            logger.info("ATC names %d parents %d", len(self.__atcD["names"]), len(self.__atcD["parents"]))
            return True
        return False

    def getAtcName(self, atcId):
        try:
            return self.__atcD["names"][atcId]["name"]
        except Exception:
            logger.info("Undefined ATC  %r", atcId)
        return None

    def getIdLineage(self, atcId):
        pList = []
        try:
            pList.append(atcId)
            pt = self.__atcD["parents"][atcId]
            while (pt is not None) and (pt != 0):
                pList.append(pt)
                pt = self.__atcD["parents"][pt]
        except Exception as e:
            logger.exception("Failing for %r with %s", atcId, str(e))
        #
        pList.reverse()
        return pList

    def getNameLineage(self, atcId):
        try:
            return [self.getAtcName(cId) for cId in self.getIdLineage(atcId)]
        except Exception as e:
            logger.exception("Failing for %r with %s", atcId, str(e))
        return None

    def getTreeNodeList(self, filterD=None):
        return self.__exportTreeNodeList(self.__atcD["names"], self.__atcD["parents"], filterD=filterD)

    def __extractNames(self, atcL):
        """
        Example public data set from BioPortal -
        Class ID,Preferred Label,Synonyms,Definitions,Obsolete,CUI,Semantic Types,Parents,ATC LEVEL,Is Drug Class,Semantic type UMLS property

        http://purl.bioontology.org/ontology/ATC/A03AX13,silicones,,,false,C0037114,http://purl.bioontology.org/ontology/STY/T109|http://purl.bioontology.org/ontology/STY/T122,http://purl.bioontology.org/ontology/UATC/A03AX,5,,http://purl.bioontology.org/ontology/STY/T109|http://purl.bioontology.org/ontology/STY/T122
        http://purl.bioontology.org/ontology/ATC/J01DB07,cefatrizine,,,false,C0007545,http://purl.bioontology.org/ontology/STY/T195|http://purl.bioontology.org/ontology/STY/T109,http://purl.bioontology.org/ontology/UATC/J01DB,5,,http://purl.bioontology.org/ontology/STY/T195|http://purl.bioontology.org/ontology/STY/T109

        Occasionally changes from "ATC" to "UATC", so capture these too:
        http://purl.bioontology.org/ontology/UATC/A03AX13,silicones,,,false,C0037114,http://purl.bioontology.org/ontology/STY/T109|http://purl.bioontology.org/ontology/STY/T122,http://purl.bioontology.org/ontology/UATC/A03AX,5,,http://purl.bioontology.org/ontology/STY/T109|http://purl.bioontology.org/ontology/STY/T122
        http://purl.bioontology.org/ontology/UATC/J01DB07,cefatrizine,,,false,C0007545,http://purl.bioontology.org/ontology/STY/T195|http://purl.bioontology.org/ontology/STY/T109,http://purl.bioontology.org/ontology/UATC/J01DB,5,,http://purl.bioontology.org/ontology/STY/T195|http://purl.bioontology.org/ontology/STY/T109
        """
        nD = {}
        nsATC = "http://purl.bioontology.org/ontology/ATC/"
        nsUATC = "http://purl.bioontology.org/ontology/UATC/"
        for rD in atcL:
            if nsATC in rD["Class ID"]:
                ns = nsATC
            elif nsUATC in rD["Class ID"]:
                ns = nsUATC
            else:
                continue
            idCode = rD["Class ID"].replace(ns, "")
            name = self.__standardizePreferredLabel(rD["Preferred Label"])
            synonyms = rD["Synonyms"]
            definition = rD["Definitions"]
            level = int(rD["ATC LEVEL"]) if rD["ATC LEVEL"].isdigit() else None
            nD[idCode] = {"name": name, "synonyms": synonyms, "definition": definition, "level": level}
        logger.info("Length of name dictionary %d", len(nD))
        return nD

    def __standardizePreferredLabel(self, label):
        """Temporary method to standardize the "Preferred Label" value of an ATC datum,
        due to recent unintended changes to source ATC data in which the column is appearing
        as a stringified array instead of a simple string. E.g.,:

        Before:
            http://purl.bioontology.org/ontology/ATC/B05DB,Hypertonic solutions,...<remainder_omitted>

        Now (as of Version 2024AA):
            http://purl.bioontology.org/ontology/ATC/B05DB,"[""Hypertonic solutions""]",...<remainder_omitted>

        Once this is resolved at the source data level, can stop using this method.

        Args:
            label (str): raw label

        Returns:
            label: processed label
        """
        try:
            # Attempt to parse the label as a list
            parsedLabel = ast.literal_eval(label)
            # If it's a list, return the first element (assuming it's a single-item array)
            if isinstance(parsedLabel, list) and len(parsedLabel) == 1:
                return parsedLabel[0]
            else:
                logger.info("parsed label %r", parsedLabel)
        except (ValueError, SyntaxError):
            # If not a list or parsing fails, return the original label
            return label
        return label

    def __extractHierarchy(self, atcL):
        """ """
        pD = {}
        nsATC = "http://purl.bioontology.org/ontology/ATC/"
        nsUATC = "http://purl.bioontology.org/ontology/UATC/"
        for rD in atcL:
            if nsATC in rD["Class ID"]:
                ns = nsATC
            elif nsUATC in rD["Class ID"]:
                ns = nsUATC
            else:
                continue
            idCode = rD["Class ID"].replace(ns, "")
            pIdCode = rD["Parents"].replace(ns, "")
            pIdCode = None if "http://www.w3.org/2002/07/owl#Thing" in pIdCode else pIdCode
            pD[idCode] = pIdCode
        #
        logger.info("Length of parent dictionary %d", len(pD))
        # for ky, vl in pD.items():
        #    logger.info(" %r  %r", ky, vl)
        return pD

    def __exportTreeNodeList(self, nD, pD, filterD=None):
        """Create node list from the ATC  parent and name/description dictionaries."""
        #
        # pL = [0]
        pL = []
        logger.info("nD %d pD %d", len(nD), len(pD))
        # create child dictionary
        cD = {}
        for ctId, ptId in pD.items():
            if not ptId:
                pL.append(ctId)
                continue
            cD.setdefault(ptId, []).append(ctId)

        #
        logger.debug("cD %d", len(cD))
        #
        idL = []
        for rootId in sorted(pL):
            visited = set([rootId])
            queue = collections.deque(visited)
            while queue:
                tId = queue.popleft()
                idL.append(tId)
                if tId not in cD:
                    # logger.warning("No children for Atc tId %r" % tId)
                    continue
                for childId in cD[tId]:
                    if childId not in visited:
                        queue.append(childId)
                        visited.add(childId)
        #
        uD = {}
        dL = []
        for tId in idL:
            if filterD and tId not in filterD:
                continue
            # id plus parents
            idlL = self.getIdLineage(tId)
            for iDepth, idl in enumerate(idlL, 0):
                if idl in uD:
                    continue
                uD[idl] = True
                displayName = nD[idl]["name"] if idl in nD else None
                ptId = pD[idl] if idl in pD else None
                #
                # dD = {"id": tId, "name": displayName, "lineage": lL, "parents": [ptId], "depth": len(lL)}
                #
                if ptId:
                    dD = {"id": idl, "name": displayName, "parents": [ptId], "depth": iDepth}
                else:
                    dD = {"id": idl, "name": displayName, "depth": iDepth}
                #
                dL.append(dD)

        return dL
