##
# File:    ResidProvider.py
# Author:  J. Westbrook
# Date:    18-Mar-2020
#
# Updates:
#
##
"""
Utilities to read the RESID resource file and build loadable documents and identifier
    correspondences.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.utils.chemref.ResidReader import ResidReader
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ResidProvider(object):
    """Utilities to read the RESID resource file and build loadable documents and identifier
    correspondences.
    """

    def __init__(self, **kwargs):
        #
        urlTarget = kwargs.get("residUrlTarget", "ftp://ftp.pir.georgetown.edu/pir_databases/other_databases/resid/RESIDUES.XML")
        #
        self.__dirPath = kwargs.get("dirPath", ".")
        useCache = kwargs.get("useCache", True)
        #
        dirPath = kwargs.get("dirPath", ".")
        useCache = kwargs.get("useCache", True)
        residFileName = kwargs.get("residFileName", "resid_correspondences_definitions.json")
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__residD = self.__reload(urlTarget, dirPath, residFileName, useCache=useCache)

    def testCache(self):
        logger.info("RESID map length %d", len(self.__residD))
        if len(self.__residD) > 200:
            return True
        return False

    def getComponentIdsOLD(self):
        try:
            return list(self.__residD.keys())
        except Exception:
            pass
        return []

    def getMapping(self):
        return self.__residD

    def __reload(self, urlTarget, dirPath, residFileName, useCache=True):
        """Reload input mmCIF model mapping resource file and return a container list.

        Args:
            urlTarget (str): target url for resource file
            dirPath (str): path to the directory containing cache files
            residFileName (str): mapping file name
            useCache (bool, optional): flag to use cached files. Defaults to True.

        Returns:
            (dict): mapping dictionary (pdb_ccId -> RESID details)
        """
        mD = {}
        #
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        filePath = os.path.join(dirPath, fn)
        residFilePath = os.path.join(dirPath, residFileName)
        self.__mU.mkdir(dirPath)
        #
        if not useCache:
            for fp in [filePath, residFilePath]:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and fU.exists(residFilePath):
            mD = self.__mU.doImport(residFilePath, fmt="json")
        else:
            ok = True
            if not (useCache and fU.exists(filePath)):
                logger.info("Fetching url %s for resource file %s", urlTarget, filePath)
                ok = fU.get(urlTarget, filePath)
            if ok:
                logger.info("Reading %r", filePath)
                xTree = self.__mU.doImport(filePath, fmt="xml")
                rr = ResidReader()
                dbObjL = rr.read(xTree)
                mD = self.__buildResidMapping(dbObjL)
                ok = self.__mU.doExport(residFilePath, mD, fmt="json", indent=3)
        #
        return mD

    def __getPdbCrossRef(self, xRefList):
        ccId = None
        for xRef in xRefList:
            if xRef.startswith("PDBHET"):
                ccId = xRef.split(":")[1]
                break
        return ccId

    def __buildResidMapping(self, dbObjL):
        """Return a dictionary of RESID to chemical component mapping.
        doc = {
            "residCode": entryElement.findtext("{ns}Header/{ns}Code".format(ns=self.__ns)),
            "names": [name.text for name in entryElement.findall("{ns}Names/{ns}Name".format(ns=self.__ns))],
            "nameXrefs": [xref.text for xref in entryElement.findall("{ns}Names/{ns}Xref".format(ns=self.__ns))],
            "ontRefs": [xref.text for xref in entryElement.findall("{ns}SequenceCode/{ns}Xref".format(ns=self.__ns))],
            "genEnzymes": [xref.text for xref in entryElement.findall("{ns}GeneratingEnzyme/{ns}EnzymeName".format(ns=self.__ns))],
            "features": [feature.text for feature in entryElement.findall("{ns}Features/{ns}Feature".format(ns=self.__ns))],
        }
        """
        retD = {}
        logger.info("Processing (%d) RESID entries", len(dbObjL))
        for dbObj in dbObjL:
            ccId = self.__getPdbCrossRef(dbObj["nameXrefs"])
            if not ccId:
                continue
            #
            residCode = dbObj["residCode"]
            residName = dbObj["names"][0]
            related = []
            for xRef in dbObj["nameXrefs"]:
                db = xRef.split(":")[0]
                if db == "PDBHET":
                    continue
                dbCode = xRef.split(":")[1]
                #
                related.append({"resourceName": db, "resourceCode": dbCode})
            #
            modResL = []
            for feature in dbObj["features"]:
                if feature.startswith("MOD_RES"):
                    modResL.append(feature[len("MOD_RES") :])
            genEnzymes = dbObj["genEnzymes"]
            #
            ontRefs = dbObj["ontRefs"]
            retD.setdefault(ccId, []).append(
                {"residCode": residCode, "residName": residName, "relatedResource": related, "modRes": modResL, "genEnzymes": genEnzymes, "ontRefs": ontRefs}
            )

        return retD
