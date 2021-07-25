##
# File: DrugBankProvider.py
# Author:  J. Westbrook
# Date:  8-Aug-2019
#
# Build loadable correspondence data from DrugBank.
#
# Update:
# 20-Jul-2021 jdw Make this provider a subclass of StashableBase
##

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"


import logging
import os
import time

from rcsb.utils.chemref.DrugBankReader import DrugBankReader
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class DrugBankProvider(StashableBase):
    """Utilities to read the DrugBank resource file and build loadable documents and identifier
    correspondences.
    """

    def __init__(self, **kwargs):
        self.__dirName = "DrugBank"
        self.__cachePath = kwargs.get("cachePath")
        super(DrugBankProvider, self).__init__(self.__cachePath, [self.__dirName])
        #
        self.__dbMapD, self.__dbObjL = self.__reload(**kwargs)
        self.__dbD = None
        #

    def testCache(self):
        try:
            if self.__dbMapD and self.__dbObjL and "id_map" in self.__dbMapD and (len(self.__dbMapD["id_map"]) > 340) and (len(self.__dbObjL) > 995):
                logger.info("DrugBank lengths map %d dbObjL %d", len(self.__dbMapD["id_map"]), len(self.__dbObjL))
                return True
        except Exception as e:
            logger.debug("Failing with %s", str(e))

        return False

    def getMapping(self):
        return self.__dbMapD

    def getVersion(self):
        return self.__dbMapD["version"] if self.__dbMapD and "version" in self.__dbMapD else None

    def getDocuments(self, mapD=None):
        """Get documents subject to corresponding identifiers in the input mapping dictionary.

        Args:
            mapD ([type], optional): mapping dictionary {drugbank_id: ccId}. Defaults to None.

        Returns:
            (list): loadable DrugBank documents
        """
        if mapD:
            return self.__buildDocuments(self.__dbObjL, mapD)
        else:
            dbIdD = {}
            if "id_map" in self.__dbMapD:
                for ccId, dD in self.__dbMapD["id_map"].items():
                    if "drugbank_id" in dD:
                        dbIdD.setdefault(dD["drugbank_id"], []).append(ccId)
            return self.__buildDocuments(self.__dbObjL, dbIdD)

    def getFeature(self, drugbankId, featureName):
        """For the input DrugBank code return the value of the input feature.

        Args:
            dbId (str): DrugBank identifier code
            featureName (str): DrugBank feature name

        Returns:
            any: feature value
        """
        try:
            # --
            if not self.__dbD:
                self.__dbD = {}
                for dD in self.__dbObjL:
                    self.__dbD[dD["drugbank_id"]] = dD
            # --
            return self.__dbD[drugbankId][featureName]
        except Exception as e:
            logger.error("Failing for %r %r with %s", drugbankId, featureName, str(e))
        return None

    def getD(self, drugbankId):
        try:
            # --
            if not self.__dbD:
                self.__dbD = {}
                for dD in self.__dbObjL:
                    self.__dbD[dD["drugbank_id"]] = dD
            # --
            return self.__dbD[drugbankId]
        except Exception as e:
            logger.error("Failing for %r with %s", drugbankId, str(e))
            return None

    def getDrugbankIds(self):
        try:
            # --
            if not self.__dbD:
                self.__dbD = {}
                for dD in self.__dbObjL:
                    self.__dbD[dD["drugbank_id"]] = dD
            # --
            return list(self.__dbD.keys())
        except Exception as e:
            logger.error("Failing with %s", str(e))
            return None

    def __reload(self, **kwargs):
        """Reload DrugBank mapping data and optionally supporting repository data file.

        Args:
            urlTarget (str): target url for resource file
            dirPath (str): path to the directory containing cache files
            useCache (bool, optional): flag to use cached files. Defaults to True.
            useDownload (bool, optional): flag to use a downloaded DrugBank resource file. Default to True.

        Returns:
            (dict, list): identifiers mapping dictionary (pdb_ccId -> DrugBank details), DrugBank object list
        """
        startTime = time.time()
        logger.info("Starting db reload at %s", time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        dirPath = os.path.join(self.__cachePath, self.__dirName)
        useCache = kwargs.get("useCache", True)
        mappingFilePath = os.path.join(dirPath, "drugbank_pdb_mapping.json")
        docListFilePath = os.path.join(dirPath, "drugbank_documents.pic")
        #
        dbMapD = {}
        dbObjL = []
        fU = FileUtil()
        mU = MarshalUtil(workPath=dirPath)
        mU.mkdir(dirPath)
        #
        if useCache and fU.exists(mappingFilePath) and fU.exists(docListFilePath):
            logger.debug("Using cached %r", mappingFilePath)
            dbMapD = mU.doImport(mappingFilePath, fmt="json")
            dbObjL = mU.doImport(docListFilePath, fmt="pickle")
            # done all cached -
            endTime = time.time()
            logger.info(
                "Completed cache recovery (%d/%d) at %s (%.4f seconds)",
                len(dbObjL),
                len(dbMapD["id_map"]),
                time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                endTime - startTime,
            )
            return dbMapD, dbObjL
        #
        elif useCache:
            return dbMapD, dbObjL
        #
        # Rebuild cache file from source
        urlTarget = kwargs.get("urlTarget", "https://go.drugbank.com/releases/latest/downloads/all-full-database")
        filePath = os.path.join(dirPath, "full database.xml")
        useDownload = kwargs.get("useDownload", True)
        username = kwargs.get("username", None)
        password = kwargs.get("password", None)
        for fp in [filePath, mappingFilePath, docListFilePath]:
            try:
                os.remove(fp)
            except Exception:
                pass
        ok = fU.exists(filePath)
        if not ok:
            if not username or not password:
                logger.warning("Missing credentials for DrugBank file download...")
            zipFilePath = os.path.join(dirPath, "full_database.zip")
            if useDownload and fU.exists(zipFilePath):
                logger.info("Using existing downloaded file %r", zipFilePath)
            else:
                logger.info("Fetching url %s to resource file %s", urlTarget, filePath)
                ok = fU.get(urlTarget, zipFilePath, username=username, password=password)
                endTime = time.time()
                logger.info("Completed db fetch at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
            #
            fp = fU.uncompress(zipFilePath, outputDir=dirPath)
            ok = fp.endswith("full database.xml")
            endTime = time.time()
            logger.info("Completed unzip at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)

        if ok:
            logger.debug("Reading %r", filePath)
            xTree = mU.doImport(filePath, fmt="xml")
            endTime = time.time()
            logger.info("Completed xml read at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
            dbr = DrugBankReader()
            version, dbObjL = dbr.read(xTree)
            endTime = time.time()
            logger.info("Completed parsing (%d) (%r) at %s (%.4f seconds)", len(dbObjL), version, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)

            dbMapD = self.__buildMapping(dbObjL)
            dbMapD["version"] = version
            ok = mU.doExport(mappingFilePath, dbMapD, fmt="json", indent=3, enforceAscii=False)
            ok = mU.doExport(docListFilePath, dbObjL, fmt="pickle")
            if ok:
                fU.remove(zipFilePath)
                fU.remove(filePath)
            endTime = time.time()
            logger.info(
                "Completed db %d/%d processing at %s (%.4f seconds)", len(dbObjL), len(dbMapD["id_map"]), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime
            )
        else:
            logger.error("Drugbank resource file missing %r", fp)
        #
        return dbMapD, dbObjL

    def __buildDocuments(self, dbObjL, dbIdD=None):
        """Build loadable documents subject to corresponding identifiers in the input mapping dictionary.

        Args:
            dbObjL (list): list of extracted DrugBank objects
            dbIdD (dict, optional): dictionary of DrugBank identifier mapping. Defaults to None.

        Returns:
            (list): DrugBank loadable documents
        """
        dbIdD = dbIdD if dbIdD else {}
        dbDocL = []
        for dbObj in dbObjL:
            if "drugbank_id" in dbObj:
                if dbIdD and dbObj["drugbank_id"] not in dbIdD:
                    continue
                lD = self.__buildDocument(dbObj)
                dbDocL.append(lD)
            else:
                logger.debug("dbObj.keys() %r", list(dbObj.keys()))
        return dbDocL
        #

    def __buildDocument(self, dbObj):
        """Construct local loadable object from input DrugBank extracted data object
        conforming to the following category organization:

         _drugbank_info.drugbank_id
         _drugbank_info.name
         _drugbank_info.description
         _drugbank_info.synonyms
         _drugbank_info.brand_names
         _drugbank_info.affected_organisms
         _drugbank_info.indication
         _drugbank_info.pharmacology
         _drugbank_info.mechanism_of_action
         _drugbank_info.cas_number
         _drugbank_info.drug_categories
         _drugbank_info.drug_groups
         _drugbank_info.atc_codes


         _drugbank_target.ordinal
         _drugbank_target.name
         _drugbank_target.interaction_type
         _drugbank_target.target_actions
         _drugbank_target.organism_common_name
         _drugbank_target.reference_database_name
         _drugbank_target.reference_database_accession_code
         _drugbank_target.seq_one_letter_code

        """
        oD = {}
        # oD["_drugbank_id"] = dbObj["drugbank_id"]
        oD["drugbank_container_identifiers"] = {"drugbank_id": dbObj["drugbank_id"]}

        dbiD = {}
        textKeys = [
            ("drugbank_id", "drugbank_id"),
            ("name", "name"),
            ("description", "description"),
            ("indication", "indication"),
            ("pharmacology", "pharmacology"),
            ("mechanism_of_action", "mechanism_of_action"),
            ("cas_number", "cas_number"),
        ]
        listKeys = [
            ("drug_categories", "drug_categories"),
            ("groups", "drug_groups"),
            ("aliases", "synonyms"),
            ("products", "brand_names"),
            ("affected_organisms", "affected_organisms"),
            ("atc_codes", "atc_codes"),
        ]
        # For category drugbank_info
        for textKey, docKey in textKeys:
            if textKey in dbObj and dbObj[textKey]:
                dbiD[docKey] = dbObj[textKey].replace("\r", "").replace("\n", " ")
        #
        for listKey, docKey in listKeys:
            if listKey in dbObj and dbObj[listKey]:
                dbiD[docKey] = list(dbObj[listKey])
            else:
                logger.debug("MISSING KEY %s", listKey)
        #
        oD["drugbank_info"] = dbiD
        #
        # For category drugbank_target -
        #
        tL = []
        if "target_interactions" in dbObj:
            for ii, tid in enumerate(dbObj["target_interactions"], 1):
                tD = {}
                tD["ordinal"] = ii
                tD["interaction_type"] = tid["category"]
                tD["name"] = tid["name"]
                if "organism" in tid:
                    tD["organism_common_name"] = tid["organism"]
                if "actions" in tid and tid["actions"]:
                    tD["target_actions"] = tid["actions"]
                if "amino-acid-sequence" in tid and tid["amino-acid-sequence"] and tid["amino-acid-sequence"]:
                    tD["seq_one_letter_code"] = "".join(tid["amino-acid-sequence"].split("\n")[1:])
                if "uniprot_ids" in tid and tid["uniprot_ids"]:
                    tD["reference_database_name"] = "UniProt"
                    tD["reference_database_accession_code"] = tid["uniprot_ids"]
                #
                tL.append(tD)
        #
        if tL:
            oD["drugbank_target"] = tL
        #
        return oD

    def __buildMapping(self, dbObjL):
        """Build PDB->DrugBank identifier correspondences.

        Args:
            dbObjL (list): list of extracted DrugBank data objects

        Returns:
            (dict): dictionary DrugBank identifier correspondences
        """
        logger.info("DrugBank full record length %d", len(dbObjL))
        dbMapD = {}
        mD = {}
        for dD in dbObjL:
            dbId = dD["drugbank_id"]
            pdbIds = ""
            if "external_identifiers" in dD:
                for exD in dD["external_identifiers"]:
                    if exD["resource"] == "PDB":
                        logger.debug("dbId %s pdbids %r ccids %r", dbId, pdbIds, exD["identifier"])
                        if exD["identifier"] not in mD:
                            mD[exD["identifier"]] = []
                        mD[exD["identifier"]] = {"drugbank_id": dbId, "aliases": list(dD["aliases"])}
                        #
                        if "atc_codes" in dD and dD["atc_codes"]:
                            mD[exD["identifier"]]["atc_codes"] = dD["atc_codes"]
                        #
                        if "products" in dD and dD["products"]:
                            mD[exD["identifier"]]["brand_names"] = dD["products"]
                        #
                        if "target_interactions" in dD:
                            for tid in dD["target_interactions"]:
                                tD = {}
                                tD["type"] = tid["category"]
                                tD["name"] = tid["name"]
                                tD["organism"] = tid["organism"]
                                if tid["actions"]:
                                    tD["actions"] = tid["actions"]
                                if tid["known_action"]:
                                    tD["known_action"] = tid["known_action"]
                                if "uniprot_ids" in tid:
                                    tD["uniprot_ids"] = tid["uniprot_ids"]
                                #
                                if "target_interactions" not in mD[exD["identifier"]]:
                                    mD[exD["identifier"]]["target_interactions"] = []
                                mD[exD["identifier"]]["target_interactions"].append(tD)
        logger.info("Match length is %d", len(mD))
        dbMapD["id_map"] = mD
        #
        inD = {}
        atcD = {}
        for dD in dbObjL:
            dbId = dD["drugbank_id"]
            if "inchikey" in dD and dD["inchikey"] and len(dD["inchikey"]) > 13:
                if dD["inchikey"] not in inD:
                    inD[dD["inchikey"]] = []
                inD[dD["inchikey"]].append({"drugbank_id": dbId, "inchikey": dD["inchikey"], "name": dD["name"]})
            #
            if "atc_codes" in dD and dD["atc_codes"]:
                atcD[dbId] = dD["atc_codes"]
        logger.info("Drugbank InChIKey dictionary length %d", len(inD))
        logger.info("Drugbank ATC  dictionary length %d", len(atcD))
        #
        dbMapD["inchikey_map"] = inD
        dbMapD["db_atc_map"] = atcD
        return dbMapD
