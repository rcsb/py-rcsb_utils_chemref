##
# File: DrugBankUtils.py
# Author:  J. Westbrook
# Date:  8-Aug-2019
#
# Build loadable correspondence data from DrugBank.
#
# Update:
#
##

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"


import logging
import os

from rcsb.utils.chemref.DrugBankReader import DrugBankReader
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

# from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class DrugBankUtils(object):
    """ Utilities to read the DrugBank resource file and build loadable documents and identifier
    correspondences.
    """

    def __init__(self, **kwargs):
        #
        urlTarget = kwargs.get("urlTarget", "https://www.drugbank.ca/releases/latest/downloads/all-full-database")
        dirPath = kwargs.get("dirPath", ".")
        useCache = kwargs.get("useCache", True)
        clearCache = kwargs.get("clearCache", False)
        username = kwargs.get("username", None)
        password = kwargs.get("password", None)

        mappingFileName = kwargs.get("mappingFileName", "drugbank_pdb_mapping.json")
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__dbMapD, self.__dbObjL = self.__reload(urlTarget, dirPath, mappingFileName, useCache=useCache, clearCache=clearCache, username=username, password=password)
        #

    def getMapping(self):
        return self.__dbMapD

    def getDocuments(self, mapD=None):
        """[summary]

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
                        dbIdD[dD["drugbank_id"]] = ccId
            return self.__buildDocuments(self.__dbObjL, dbIdD)

    def __reload(self, urlTarget, dirPath, mappingFileName, useCache=True, clearCache=False, username=None, password=None):
        """ Reload DrugBank mapping data and optionally supporting repository data file.

        Args:
            urlTarget (str): target url for resource file
            dirPath (str): path to the directory containing cache files
            mappingFileName (str): mapping file name
            useCache (bool, optional): flag to use cached files. Defaults to True.
            clearCache (bool, optional): flag to clear any cached files. Defaults to False.

        Returns:
            (dict, list): identifiers mapping dictionary (pdb_ccId -> DrugBank details), DrugBank object list
        """
        dbMapD = {}
        dbObjL = []
        fU = FileUtil()
        mappingFilePath = os.path.join(dirPath, mappingFileName)
        filePath = os.path.join(dirPath, "full database.xml")
        #
        if clearCache:
            try:
                os.remove(filePath)
                os.remove(mappingFilePath)
            except Exception:
                pass
        #
        if useCache and fU.exists(mappingFilePath):
            logger.debug("Using cached %r", mappingFilePath)
            dbMapD = self.__mU.doImport(mappingFilePath, fmt="json")

        ok = True
        if not fU.exists(filePath):
            if not username or not password:
                logger.warning("Missing credentials for DrugBank file download...")
            logger.debug("Fetching url %s to resource file %s", urlTarget, filePath)
            zipFilePath = os.path.join(dirPath, "full_database.zip")
            ok = fU.get(urlTarget, zipFilePath, username=username, password=password)
            fp = fU.uncompress(zipFilePath, outputDir=dirPath)
            ok = fp.endswith("full database.xml")
        if ok:
            logger.debug("Reading %r", filePath)
            xTree = self.__mU.doImport(filePath, fmt="xml")
            dbr = DrugBankReader()
            dbObjL = dbr.read(xTree)
            dbMapD = self.__buildMapping(dbObjL)
            ok = self.__mU.doExport(mappingFilePath, dbMapD, fmt="json", indent=3)
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
        oD["_drugbank_id"] = dbObj["drugbank_id"]
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
        logger.debug("DrugBank full record length %d", len(dbObjL))
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
        for dD in dbObjL:
            dbId = dD["drugbank_id"]
            if "inchikey" in dD and dD["inchikey"] and len(dD["inchikey"]) > 13:
                if dD["inchikey"] not in inD:
                    inD[dD["inchikey"]] = []
                inD[dD["inchikey"]].append({"drugbank_id": dbId, "inchikey": dD["inchikey"], "name": dD["name"]})
        #
        logger.info("Drugbank InChIKey dictionary length %d", len(inD))
        #
        dbMapD["inchikey_map"] = inD
        return dbMapD
