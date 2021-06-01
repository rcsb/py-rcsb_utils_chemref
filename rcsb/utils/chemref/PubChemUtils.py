##
# File:    PubChemUtils.py
# Date:    1-May-2020
#
# Updates:
#
#   28-May-2020 jdw use alternative url fetch library.
##
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

# pylint: disable=too-many-lines

import logging
import os
import time
from collections import namedtuple, defaultdict

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.UrlRequestUtil import UrlRequestUtil

try:
    from urllib.parse import quote
except ImportError:
    from urllib2 import quote

logger = logging.getLogger(__name__)

ExternalAnnotationFields = ("featureType", "name", "value", "description", "dataType", "units", "provReferenceId", "provSourceName")
ExternalAnnotation = namedtuple("ExternalAnnotation", ExternalAnnotationFields, defaults=(None,) * len(ExternalAnnotationFields))

ChemicalIdentifierFields = ("idCode", "identifierSource", "identifierType", "identifier", "indexName")
ChemicalIdentifier = namedtuple("ChemicalIdentifier", ChemicalIdentifierFields, defaults=(None,) * len(ChemicalIdentifierFields))


class PubChemUtils(object):
    """Manage search and fetch queries for PubChem compound, substance and assay data and related annotations."""

    def __init__(self, **kwargs):
        self.__verbose = kwargs.get("verbose", False)
        self.__delaySeconds = float(kwargs.get("delaySeconds", "0.15"))
        self.__urlPrimary = kwargs.get("urlPrimary", "https://pubchem.ncbi.nlm.nih.gov")
        #

    def assemble(self, chemicalIdentifier, **kwargs):
        """Build a PubChem exchange data object.

        Args:
            chemicalIdentifier (namedtuple): ChemicalIdentifier(identifierSource, identifierType, identifier)
            exportPath (str, optional): path to export intermediate fetch results (default: None)
            exportIntermediates (bool, optional): flag controlling export of intermediate results (default: False)
            matchIdOnly (bool, optional): flag to stop after performing record lookup to obtain matching identifiers
            contentTypes (list): PubChem content types to fetch (default: ["view", "classification", "property", "xrefs", \
                                                                           "synonyms", "dgidb", "pathway", \
                                                                           "fdaorangebook", "clinicaltrials", "bioactivity"])


        Return:
             (bool, dict): status, PubChem exchange data object

        """
        exportPath = kwargs.get("exportPath", None)
        exportIntermediates = kwargs.get("exportIntermediates", False)
        matchIdOnly = kwargs.get("matchIdOnly", False)
        contentTypes = kwargs.get("contentTypes", ["view", "classification", "property", "xrefs", "synonyms", "dgidb", "pathway", "fdaorangebook", "clinicaltrials", "bioactivity"])

        retStatus = False
        retD = {}
        assemDL = []
        # -- First identifier the PubChem compound CID for the input chemical identifier.
        if chemicalIdentifier.identifierType in ["smiles", "inchi"]:
            searchType = "fastidentity"
        elif chemicalIdentifier.identifierType in ["cid", "name", "inchikey"]:
            searchType = "lookup"
        else:
            pass
        #
        retStatus, rDL = self.fetch(chemicalIdentifier, searchType=searchType, returnType="record", delaySeconds=self.__delaySeconds)
        if not retStatus:
            return retStatus, assemDL
        #
        # PubChem compound CIDs
        cidList = []
        for rD in rDL:
            cid = rD["cid"] if "cid" in rD else None
            if cid:
                cidList.append(cid)
                retD.setdefault(cid, {}).update({"record": rD})
        #
        logger.info(
            "%s searchType %r idType %r idSource %r finds PubChem compounds (%d) %r",
            chemicalIdentifier.idCode,
            searchType,
            chemicalIdentifier.identifierType,
            chemicalIdentifier.identifierSource,
            len(cidList),
            cidList,
        )
        if not cidList:
            return retStatus, assemDL

        if matchIdOnly:
            for cid in retD:
                assemDL.append({"cid": cid, "data": retD[cid]})
            return retStatus, assemDL
        #
        # --  Build data object for each CID
        #
        for cid in cidList:
            chemId = ChemicalIdentifier(idCode=cid, identifier=cid, identifierType="cid")
            for returnType in contentTypes:
                rawResponsePath, extractedResponsePath = None, None
                if exportPath and exportIntermediates:
                    rawResponsePath = os.path.join(exportPath, "%s-pubchem-%s-raw.json" % (cid, returnType))
                    extractedResponsePath = os.path.join(exportPath, "%s-pubchem-%s-extract.json" % (cid, returnType))
                tStatus, rDL = self.fetch(chemId, returnType=returnType, storeRawResponsePath=rawResponsePath, storeResponsePath=extractedResponsePath, delaySeconds=self.__delaySeconds)
                # retStatus = tStatus and retStatus
                if tStatus and rDL:
                    if returnType in ["view", "classification", "property", "xrefs", "synonyms"] and len(rDL[0]) > 1:
                        retD[cid][returnType] = rDL[0]
                    else:
                        retD[cid][returnType] = rDL
        #
        for cid in retD:
            assemDL.append({"cid": cid, "data": retD[cid]})
        #
        if exportPath:
            resultPath = os.path.join(exportPath, "%s-pubchem-extracted-all.json" % cid)
            mU = MarshalUtil()
            mU.doExport(resultPath, assemDL, fmt="json", indent=3)

        return retStatus, assemDL

    def fetch(self, chemicalIdentifier, searchType="lookup", returnType="record", storeResponsePath=None, storeRawResponsePath=None, delaySeconds=0.15):
        """Fetch PubChem content (returnType) for the input chemical identifier using the input searchType.

        Args:
            chemicalIdentifer (ChemicalIdentifier Tuple): chemical identifier tuple
            searchType (string): search type (lookup|fastsimilarity|fastidentity)
            returnType (string): object type to return (default: record)

        Return:
            (bool, list): status, return object list (selected items extracted from each returned record type)

        """
        if self.__verbose:
            logger.info(
                "fetching identifier %r nameSpace %r searchType %r returnType %r delay %r",
                chemicalIdentifier.identifier,
                chemicalIdentifier.identifierType,
                searchType,
                returnType,
                delaySeconds,
            )
        #
        if delaySeconds:
            time.sleep(delaySeconds)
        #
        ok = False
        retL = None
        if returnType in ["dgidb", "pathway", "fdaorangebook", "clinicaltrials", "bioactivity"] and chemicalIdentifier.identifierType == "cid":
            response, retCode = self.__doSgdRequest(chemicalIdentifier.identifier, returnType=returnType)
        elif returnType == "view":
            response, retCode = self.__doPugViewRequest(chemicalIdentifier.identifier, nameSpace=chemicalIdentifier.identifierType, domain="compound")
        else:
            response, retCode = self.__doPugRequest(
                chemicalIdentifier.identifier, nameSpace=chemicalIdentifier.identifierType, domain="compound", searchType=searchType, returnType=returnType
            )
        if storeRawResponsePath and response:
            mU = MarshalUtil()
            mU.doExport(storeRawResponsePath, response, fmt="json", indent=3)

        ok = retCode in [200] and response and len(response) > 0
        #
        if ok and returnType == "record":
            retL = self.__parsePubChemRecord(response)
        elif ok and returnType == "xrefs":
            retL = self.__parsePubChemXrefs(response)
        elif ok and returnType == "property":
            retL = self.__parsePubChemProperties(response)
        elif ok and returnType == "classification" and chemicalIdentifier.identifierType == "cid":
            retL = self.__parsePubChemClassifications(chemicalIdentifier.identifier, response)
        elif ok and returnType == "synonyms":
            retL = self.__parsePubChemSynonyms(response)
        elif ok and returnType == "view":
            retL = [self.__parsePubChemCompoundView(response)]
        elif ok and returnType in ["dgidb", "pathway", "fdaorangebook", "clinicaltrials", "bioactivity"]:
            retL = self.__parseSdqResponse(returnType, response)
        #
        if storeResponsePath and retL:
            mU = MarshalUtil()
            mU.doExport(storeResponsePath, retL, fmt="json", indent=3)
        #
        return ok, retL

    def __doPugRequest(self, identifier, nameSpace="cid", domain="compound", searchType="lookup", returnType="record"):
        """Wrapper for PubChem PUG API requests

        https://pubchem.ncbi.nlm.nih.gov/rest/pug/<input specification>/<operation specification>/[<output specification>][?<operation_options>]

        Input:
                <input specification> = <domain>/<namespace>/<identifiers>
                <domain> = substance | compound | assay | <other inputs>
                compound domain <namespace> = cid | name | smiles | inchi | sdf | inchikey | formula | <structure search> | <xref> | listkey | <fast search>
        Operations:
            compound domain <operation specification> = record | <compound property> | synonyms | sids | cids | aids | assaysummary | classification | <xrefs> | description | conformers
            <compound property> = property / [comma-separated list of property tags]

            substance domain <operation specification> = record | synonyms | sids | cids | aids | assaysummary | classification | <xrefs> | description
            <xrefs> = xrefs / [comma-separated list of xrefs tags]

        Examples:
            https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/smiles/C1=NC2=C(N1)C(=O)N=C(N2)N/JSON
            https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/C1=NC2=C(N1)C(=O)N=C(N2)N/JSON?Threshold=99
            https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/MolecularFormula,MolecularWeight/CSV
            with “cid=1,2,3,4,5” in the POST body, and  “Content-Type: application/x-www-form-urlencoded”

        """
        #
        httpCodesCatch = [404]
        requestType = "POST"
        outputType = "JSON"
        sslCert = "enabled"
        timeOutSeconds = 10
        retCode, ret = None, None
        try:
            baseUrl = self.__urlPrimary
            pD = {}
            ureq = UrlRequestUtil()
            if domain == "compound":
                if nameSpace in ["cid", "name", "inchikey"] and returnType in ["record"] and searchType in ["lookup"] and requestType == "GET":
                    uId = quote(identifier.encode("utf8"))
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, uId, outputType])
                    ret, retCode = ureq.getUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["record"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
                #
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["xrefs"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, "RegistryID,RN", outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["property"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, "MolecularFormula,XLogP,TPSA,Volume3D", outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["synonyms"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
                #
                elif nameSpace in ["cid"] and returnType in ["classification"] and searchType in ["lookup"]:
                    # Needs to be specifically targeted on a particular compound ...
                    uId = quote(identifier.encode("utf8"))
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, uId, returnType, outputType])
                    # pD = {"classification_type": "simple"}
                    pD = {}
                    ret, retCode = ureq.getUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
                #
                elif nameSpace in ["cid"] and returnType in ["classification"] and searchType in ["lookup"] and requestType == "POST":
                    # Needs to be specifically targeted on a particular compound ...
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, outputType])
                    # This is a long request return server codes may be observed 500
                    # pD = {nameSpace: identifier, "classification_type": "simple"}
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
                #
                elif nameSpace in ["smiles", "inchi"] and returnType == "record" and searchType in ["fastidentity", "fastsimilarity_2d"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, searchType, nameSpace, outputType])
                    if searchType == "fastsimilarity_2d":
                        pD = {nameSpace: identifier, "Threshold": 95}
                    elif searchType == "fastidentity":
                        identityType = "same_stereo_isotope"
                        if identityType in [
                            "same_connectivity",
                            "same_tautomer",
                            "same_stereo",
                            "same_isotope",
                            "same_stereo_isotope",
                            "nonconflicting_stereo",
                            "same_isotope_nonconflicting_stereo",
                        ]:
                            pD = {nameSpace: identifier, "identity_type": identityType}
                        else:
                            pD = {nameSpace: identifier}
                    ret, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
            elif domain in ["substance", "assay"]:
                logger.error("Fetch not implemented for domain %s", domain)
            #

        except Exception as e:
            logger.error("Failing identifier %r returnType %r nameSpace %r with (retCode %r) %s", identifier, returnType, nameSpace, retCode, str(e))

        return ret, retCode

    def __doPugViewRequest(self, identifier, nameSpace="cid", domain="compound"):
        """Wrapper for PubChem PUG_VIEW API requests

        Example:
            https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/1234/JSON

        """
        requestType = "GET"
        outputType = "JSON"
        retCode, ret = None, None
        httpCodesCatch = [404]
        sslCert = "enabled"
        timeOutSeconds = 10
        try:
            baseUrl = self.__urlPrimary
            #
            pD = {}
            ureq = UrlRequestUtil()
            if nameSpace in ["cid"] and requestType == "GET":
                uId = quote(identifier.encode("utf8"))
                endPoint = "/".join(["rest", "pug_view", "data", domain, uId, outputType])
                ret, retCode = ureq.getUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
            elif nameSpace in ["cid"] and requestType == "POST":
                endPoint = "/".join(["rest", "pug_view", "data", domain, outputType])
                pD = {nameSpace: identifier}
                ret, retCode = ureq.postUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
        except Exception as e:
            logger.error("Failing for identifier %r returnType view with (retCode %r) %s", identifier, retCode, str(e))
        #
        return ret, retCode

    def __doSgdRequest(self, identifier, returnType="dgidb"):
        """Wrapper for PubChem SGD API requests

            Example:

        /sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={"download":"*","collection":"dgidb","where":{"ands":[{"cid":"2244"}]},
                          "order":["relevancescore,desc"],"start":1,"limit":10000000,"downloadfilename":"CID_2244_dgidb"}
        /sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={"download":"*","collection":"pathway","where":{"ands":[{"cid":"2244"},{"core":"1"}]},
                          "order":["name,asc"],"start":1,"limit":10000000,"downloadfilename":"CID_2244_pathway"}

        """
        nameSpace = "cid"
        requestType = "GET"
        retCode, ret = None, None
        httpCodesCatch = [404]
        sslCert = "enabled"
        timeOutSeconds = 10
        try:
            baseUrl = self.__urlPrimary
            pD = {}
            ureq = UrlRequestUtil()
            if nameSpace in ["cid"] and requestType == "GET" and returnType in ["dgidb", "pathway", "fdaorangebook", "clinicaltrials", "bioactivity"]:
                # uId = quote(identifier.encode("utf8"))
                endPoint = "/".join(["sdq", "sdqagent.cgi"])
                pD = {"infmt": "json", "outfmt": "json", "query": '{"select":"*","collection":"%s","where":{"ands":{"cid":"%s"}}}' % (returnType, identifier)}
                ret, retCode = ureq.getUnWrapped(baseUrl, endPoint, pD, httpCodesCatch=httpCodesCatch, returnContentType="JSON", sslCert=sslCert, timeOut=timeOutSeconds)
        except Exception as e:
            logger.error("Failing identifier %r return type %r with (retCode %r) %s", identifier, returnType, retCode, str(e))
        #
        logger.debug("identifier %s returnType %r [%r] result not empty %r", identifier, returnType, retCode, ret is not None)
        #  Display the column names returned from the various SDQ collections -
        # logger.debug("[%r] %s %s", retCode, ret["SDQOutputSet"][0]["collection"], ret["SDQOutputSet"][0]["rows"][0].keys())
        return ret, retCode

    def __parseSdqResponse(self, returnType, vD):
        #
        # Preliminary mapping to capture what is currently deemed useful from these collections.
        #
        cMapD = {
            "dgidb": [
                ("pmids", None),
                ("gid", None),
                ("srcid", None),
                ("geneid", "PubChemGeneId"),
                ("genename", "PubChemGeneName"),
                ("geneclaimname", None),
                ("interactionclaimsource", "provSource"),
                ("interactiontypes", "interactionType"),
                ("sid", "PubChemSid"),
                ("cid", "PubChemCid"),
                ("drugname", "drugName"),
                ("drugclaimname", None),
                ("drugclaimprimaryname", None),
                ("drugchemblid", None),
                ("dois", None),
            ],
            "pathway": [
                ("taxid", "ncbiTaxId"),
                ("pathwayid", "PubChemPathwayId"),
                ("srcid", None),
                ("core", None),
                ("cids", "PubChemCids"),
                ("geneids", "PubChemGeneIds"),
                ("ecs", "enzymeClassifications"),
                ("name", "pathwayName"),
                ("protacxns", None),
                ("taxname", "ncbiTaxonomy"),
                ("pwtype", "pathwayType"),
                ("category", "pathwayCategory"),
                ("url", "pathwayLinkoutUrl"),
                ("source", "provSource"),
                ("extid", "pathwayId"),
                ("externalid", None),
            ],
            "fdaorangebook": [
                ("marketingstatus", "marketing"),
                ("gid", None),
                ("srcid", None),
                ("cid", "PubChemCid"),
                ("activeingredient", "activeIngredient"),
                ("tradename", "tradeName"),
                ("dosageformroute", "formRoute"),
                ("strength", "strength"),
                ("applicationnumber", "FDA Application Number"),
                ("uspatent", "USPatent"),
                ("approvaldate", "FDA Approval Date"),
                ("applicant", "applicant"),
                ("url", "FDAApplicationLinkoutUrl"),
            ],
            "clinicaltrials": [
                ("updatedate", "updateDate"),
                ("gid", None),
                ("srcid", None),
                ("ctid", "clinicalTrialId"),
                ("cids", "PubChemCids"),
                ("link", "clinicalTrialsLinkoutUrl"),
                ("title", "title"),
                ("conditions", "targetConditions"),
                ("interventions", "interventions"),
                ("phase", "trialPhase"),
                ("status", "trialStatus"),
                ("diseaseids", None),
            ],
            "bioactivity": [
                ("cid", "PubChemCid"),
                ("sid", "PubChemSid"),
                ("aid", "PubChemAssayId"),
                ("geneid", "PubChemGeneId"),
                ("taxid", "ncbiTaxonomyId"),
                ("baid", "bioAssayId"),
                ("mid", None),
                ("rnai", None),
                ("hasdrc", None),
                ("aidsrcname", "assaySource"),
                ("ecs", None),
                ("repacxn", None),
                ("aidname", "assayName"),
                ("aidtype", "assayType"),
                ("acvalue", "assayValueMicroMolar"),
                ("aidmdate", "modificationDate"),
                ("activity", "activity"),
                ("protacxn", None),
                ("targetname", "targetName"),
                ("targeturl", "targetUrl"),
            ],
        }
        rL = []
        try:
            mTupL = cMapD[returnType] if returnType in cMapD else []
            for rowD in vD["SDQOutputSet"][0]["rows"]:
                rD = {}
                for srcA, destA in mTupL:
                    if destA and srcA in rowD:
                        rD[destA] = rowD[srcA]
                rL.append(rD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return rL

    def __traversePubChemCompoundView(self, vD):
        """
        "Record": {
                "RecordType": "CID",
                "RecordNumber": 123631,
                "RecordTitle": "Gefitinib",
                "Section": []

        """
        try:
            qD = vD["Record"] if "Record" in vD else {}
            #
            logger.info("keys %s", qD.keys())
            logger.info("Record.RecordType %r", self.__getKeyValues(vD, ["Record.RecordType"]))
            logger.info("Record.RecordNumber %r", self.__getKeyValues(vD, ["Record.RecordNumber"]))
            logger.info("Record.RecordTitle %r", self.__getKeyValues(vD, ["Record.RecordTitle"]))

            for tD in qD["Section"]:
                logger.info("TOCHeading %r Description %r ...", tD["TOCHeading"], tD["Description"][:30])
                if "Section" in tD:
                    for ttD in tD["Section"]:
                        # logger.debug("ttD.keys() %r", ttD.keys())
                        if "Description" in ttD:
                            logger.info("  >> TOCHeading %r Description %r ...", ttD["TOCHeading"], ttD["Description"][:30])
                        else:
                            logger.info("  >> TOCHeading %r Description None ...", ttD["TOCHeading"])
                        if "Section" in ttD:
                            for tttD in ttD["Section"]:
                                if "Description" in tttD:
                                    logger.info("        >> TOCHeading %r Description %r ...", tttD["TOCHeading"], tttD["Description"][:30])
                                else:
                                    logger.info("        >> TOCHeading %r Description  None ...", tttD["TOCHeading"])
                                if "Section" in tttD:
                                    for ttttD in tttD["Section"]:
                                        if "Description" in ttttD:
                                            logger.info("           >> TOCHeading %r Description %r ...", ttttD["TOCHeading"], ttttD["Description"][:30])
                                        else:
                                            logger.info("           >> TOCHeading %r Description  None ...", ttttD["TOCHeading"])
                                        if "Section" in ttttD:
                                            for tttttD in ttttD["Section"]:
                                                if "Description" in tttttD:
                                                    logger.info("              >> TOCHeading %r Description %r ...", tttttD["TOCHeading"], tttttD["Description"][:30])
                                                else:
                                                    logger.info("              >> TOCHeading %r Description  None ...", tttttD["TOCHeading"])

            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return False

    def __parsePubChemCompoundView(self, vD):
        """
        "Record": {
                "RecordType": "CID",
                "RecordNumber": 123631,
                "RecordTitle": "Gefitinib",
                "Section": []

             "Information": [
                       {
                          "ReferenceNumber": 55,
                          "Reference": [
                             "Computed by InChI 1.0.5 (PubChem release 2019.06.18)"
                          ],
                          "Value": {
                             "StringWithMarkup": [
                                {
                                   "String": "InChI=1S/C22H24ClFN4O3/c1-29-20-13-19-16...",
            'Refererence': [
                {
                    "ReferenceNumber": 158,
                    "SourceName": "The National Institute for Occupational Safety and Health (NIOSH)",
                    "SourceID": "npgd0010",
                    "Name": "Acetylsalicylic acid",
                    "Description": "The NIOSH Pocket Guide to Chemical Hazards is intended as a chemicals/Read more: https://www.cdc.gov/niosh/npg/",
                    "URL": "https://www.cdc.gov/niosh/npg/npgd0010.html",
                    "LicenseNote": "The information provided using CDC Web site egulations.",
                    "LicenseURL": "https://www.cdc.gov/Other/disclaimer.html",
                    "ANID": 2266310
                },
                ]                                  }
        """
        pcD = {}
        refD = {}
        mainSectionNameD = {
            "Names and Identifiers": {
                "Computed Descriptors": {"ky": "descriptors", "subSections": []},
                "Other Identifiers": {"ky": "identifiers", "subSections": []},
                "Synonyms": {"ky": "synonyms", "subSections": ["MeSH Entry Terms"]},
            },
            "Chemical and Physical Properties": {
                "Computed Properties": {"ky": "properties", "subSections": []},
                "Experimental Properties": {"ky": "properties", "subSections": []},
            },
            "Drug and Medication Information": {"Drug Indication": {"ky": "drugInfo", "subSections": []}, "Drug Classes": {"ky": "drugInfo", "subSections": []}},
            "Pharmacology and Biochemistry": {
                "Pharmacology": {"ky": "pharmaInfo", "subSections": []},
                "MeSH Pharmacological Classification": {"ky": "pharmaInfo", "subSections": []},
                "FDA Pharmacological Classification": {"ky": "pharmaInfo", "subSections": []},
                "ATC Code": {"ky": "pharmaInfo", "subSections": []},
            },
            # note that classifications are suppressed in the view and must be obtained by other means.
            # "Classification": {"Ontologies": {"ky": "ontologies", "subSections": []}},
        }
        #
        try:
            qD = vD["Record"] if "Record" in vD else {}
            logger.debug("keys %s", qD.keys())
            logger.debug("Record.RecordType %r", self.__getKeyValues(vD, ["Record.RecordType"]))
            logger.debug("Record.RecordNumber %r", self.__getKeyValues(vD, ["Record.RecordNumber"]))
            logger.debug("Record.RecordTitle %r", self.__getKeyValues(vD, ["Record.RecordTitle"]))
            #
            for ref in qD["Reference"]:
                sourceName = ref["SourceName"] if "SourceName" in ref else None
                sourceId = ref["SourceId"] if "SourceId" in ref else None
                name = ref["Name"] if "Name" in ref else None
                description = ref["Description"] if "Description" in ref else None
                url = ref["URL"] if "URL" in ref else None
                licenseInfo = ref["LicenseNote"] if "LicsenceNote" in ref else None
                licenseUrl = ref["LicenseURL"] if "LicenseURL" in ref else None
                refD[ref["ReferenceNumber"]] = {
                    "sourceName": sourceName,
                    "sourceId": sourceId,
                    "name": name,
                    "description": description,
                    "url": url,
                    "licenseInfo": licenseInfo,
                    "licenseUrl": licenseUrl,
                }
            #
            for sectionD in qD["Section"]:
                # Only process main sections identified in mainSectionNameD
                if "TOCHeading" in sectionD and sectionD["TOCHeading"] in mainSectionNameD:
                    sectionName = sectionD["TOCHeading"]
                    sectionMapD = mainSectionNameD[sectionName]
                    logger.debug("Section %r Description %r ...", sectionName, sectionD["Description"][:30])
                    #
                    if "Section" in sectionD:
                        for subSectionD in sectionD["Section"]:
                            subSectionName = subSectionD["TOCHeading"]
                            logger.debug("  >> subSection %r", subSectionName)

                            pcKy = sectionMapD[subSectionName]["ky"] if subSectionD["TOCHeading"] in sectionMapD else None
                            inclSubSections = sectionMapD[subSectionName]["subSections"] if subSectionName in sectionMapD else []
                            # ---
                            if pcKy and "Section" in subSectionD:
                                for subSubSectionD in subSectionD["Section"]:
                                    ky = subSubSectionD["TOCHeading"]
                                    if inclSubSections and ky not in inclSubSections:
                                        logger.debug("        >> SKIPPED subsubSection %r (%r) %r", ky, pcKy, inclSubSections)
                                        continue
                                    logger.debug("        >> subsubSection %r (%r)", ky, pcKy)
                                    # ----
                                    pcD = self.__parseViewInfoList(pcKy, ky, subSubSectionD["Information"], pcD)
                                    #
                            elif pcKy and "Information" in subSectionD:
                                pcD = self.__parseViewInfoList(pcKy, subSectionD["TOCHeading"], subSectionD["Information"], pcD)

                        # ---
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        rD = {}
        for ky in pcD:
            for ky2, pTupS in pcD[ky].items():
                for pTup in pTupS:
                    pD = pTup._asdict()
                    if "provReferenceId" in pD and pD["provReferenceId"] in refD:
                        pD["provData"] = refD[pD["provReferenceId"]]
                    rD.setdefault(ky, {}).setdefault(ky2, []).append(pD)
        return rD

    def __parseViewInfoList(self, pcKy, ky, infoDL, pcD):
        # ----
        #
        try:
            for infD in infoDL:
                referenceNumber = infD["ReferenceNumber"] if "ReferenceNumber" in infD else None
                reference = ",".join(infD["Reference"]) if "Reference" in infD else None
                #
                logger.debug("reference no %r reference %r", referenceNumber, reference)
                name = infD["Name"] if "Name" in infD else None
                if "Value" in infD:
                    if "StringWithMarkup" in infD["Value"]:
                        for sV in infD["Value"]["StringWithMarkup"]:
                            tupV = ExternalAnnotation(name=name, value=sV["String"], dataType="string", featureType=ky, provReferenceId=referenceNumber, provSourceName=reference)
                            pcD.setdefault(pcKy, {}).setdefault(ky, set()).update([tupV])
                    elif "Number" in infD["Value"]:
                        units = infD["Value"]["Unit"] if "Unit" in infD["Value"] else ""
                        for nV in infD["Value"]["Number"]:
                            tupV = ExternalAnnotation(name=name, value=nV, dataType="number", units=units, featureType=ky, provReferenceId=referenceNumber, provSourceName=reference)
                            # tV = str(str(sV) + " " + units).strip()
                            pcD.setdefault(pcKy, {}).setdefault(ky, set()).update([tupV])
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return pcD

    # ---
    def __parsePubChemRecord(self, recordD):
        """Parse selected content from a PubChem PUG API full RECORD response data.

        Args:
            recordD (dict): response from a PubChem PUG API record request with potentially multiple
                            compound records

        Returns:
            list : [{'cid': <pc_cid>, ....}, ...]

        """
        retL = []
        try:
            tDL = recordD["PC_Compounds"] if "PC_Compounds" in recordD else []
            for tD in tDL:
                rD = {}
                tId = self.__getKeyValue(tD, "id.id.cid")
                rD["cid"] = str(tId) if tId else None
                logger.debug("cId = %r", rD["cid"])
                elementCountD = defaultdict(int)
                elements = self.__getKeyValue(tD, "atoms.element")
                for el in elements:
                    elementCountD[str(el)] += 1
                rD["elementCounts"] = elementCountD
                rD["formalCharge"] = self.__getKeyValue(tD, "charge")
                pDL = self.__getKeyValue(tD, "props")
                for pD in pDL:
                    pName = self.__getKeyValue(pD, "urn.name")
                    pLabel = self.__getKeyValue(pD, "urn.label")
                    if pName == "XLogP3" and pLabel == "Log P":
                        rD["XlogP"] = self.__getKeyValue(pD, "value.fval")
                    elif pLabel == "Topological" and pName == "Polar Surface Area":
                        rD["polarSurfaceArea"] = self.__getKeyValue(pD, "value.fval")
                    elif pLabel == "IUPAC Name" and pName == "Preferred":
                        rD["iupacName"] = self.__getKeyValue(pD, "value.sval")
                    elif pLabel == "Molecular Formula":
                        rD["formula"] = self.__getKeyValue(pD, "value.sval")
                    elif pLabel == "InChIKey":
                        rD["InChiKey"] = self.__getKeyValue(pD, "value.sval")

                retL.append(rD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return retL

    def __parsePubChemXrefs(self, recordD):
        """Parse selected content from a PubChem PUG API xrefs response data.

        Args:
            recordD (dict): response from a PubChem PUG API record request with potentially multiple
                            records

        Returns:
            list : [{'cid': <pc_cid>, ....}, ...]

            {
             "InformationList": {
                "Information": [
                    {
                        "CID": 2244,
                        "RegistryID": [
                             ],
                        "RN": [
                        "11126-35-5",
                        ]
                    }
                }
            }
        """

        retL = []
        try:
            tDL = self.__getKeyValue(recordD, "InformationList.Information")
            for tD in tDL:
                rD = {}
                tId = self.__getKeyValue(tD, "CID")
                rD["cid"] = str(tId) if tId else None
                logger.debug("cId = %r", rD["cid"])
                #
                rnL = self.__getKeyValue(tD, "RN")
                rnD = {rn: True for rn in rnL} if rnL else {}
                rIdL = self.__getKeyValue(tD, "RegistryID")
                # This is a laundry list of IDs and it only possible to pull out  clearly prefixed idcodes.
                for rId in rIdL:
                    if rId in rnD:
                        rD.setdefault("CAS", []).append(rId)
                    if rId.startswith("CHEMBL"):
                        rD.setdefault("ChEMBL", []).append(rId)
                    if rId.startswith("CHEBI"):
                        rD.setdefault("ChEBI", []).append(rId)
                    if rId.startswith("HMDB"):
                        rD.setdefault("HMDB", []).append(rId)
                retL.append(rD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return retL

    def __parsePubChemProperties(self, recordD):
        """Parse selected content from a PubChem PUG API properties response data.

        Args:
            recordD (dict): response from a PubChem PUG API record request with potentially multiple
                            records

        Returns:
            list : [{'cid': <pc_cid>, ....}, ...]

                {
                "PropertyTable": {
                    "Properties": [
                        {
                            "CID": 2244,
                            "MolecularFormula": "C9H8O4",
                            "XLogP": 1.2,
                            "TPSA": 63.6,
                            "Volume3D": 136
                        }
                    ]
                }
                }
        """
        retL = []
        try:
            tDL = self.__getKeyValue(recordD, "PropertyTable.Properties")
            for tD in tDL:
                rD = {}
                tId = self.__getKeyValue(tD, "CID")
                rD["cid"] = str(tId) if tId else None
                logger.debug("cId = %r", rD["cid"])
                #
                rD["formula"] = self.__getKeyValue(tD, "MolecularFormula")
                rD["XLogP"] = self.__getKeyValue(tD, "XLogP")
                rD["polarSurfaceArea"] = self.__getKeyValue(tD, "TPSA")
                rD["volume3D"] = self.__getKeyValue(tD, "Volume3D")
                retL.append(rD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return retL

    def __parsePubChemSynonyms(self, recordD):
        """Parse selected content from a PubChem PUG API synonyms response data.

        Args:
            recordD (dict): response from a PubChem PUG API record request with potentially multiple
                            records

        Returns:
            list : [{'cid': <pc_cid>, ....}, ...]
        {
        "InformationList": {
            "Information": [
                {
                    "CID": 2244,
                    "Synonym": [
                    "aspirin",
                    "ACETYLSALICYLIC ACID",
                    "50-78-2",
                    "2-Acetoxybenzoic acid",
                    ]
                },

        """
        retL = []
        try:
            tDL = self.__getKeyValue(recordD, "InformationList.Information")
            for tD in tDL:
                rD = {}
                tId = self.__getKeyValue(tD, "CID")
                rD["cid"] = str(tId) if tId else None
                logger.debug("cId = %r", rD["cid"])
                rD["synonyms"] = self.__getKeyValue(tD, "Synonym")
                retL.append(rD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return retL

    def __parsePubChemClassifications(self, cId, recordD):
        """Parse selected content from a PubChem PUG API classification response data.

        Args:
            recordD (dict): response from a PubChem PUG API record request with potentially multiple
                            records

        Returns:
            list : [{'cid': <pc_cid>, ....}, ...]

            {
            "Hierarchies": {
                "Hierarchy": [
                    {
                    "SourceName": "WHO ATC",
                    "SourceID": "ATCTree",
                    "RootID": "root",
                    "HID": 79,
                    "Information": {
                    "Name": "ATC Code",
                    "Description": [
                        "In the World Health Organization (WHO) Anatomical Therapeutic Chemical (ATC) ..."
                    ],
                    "Comments": [
                        "Drugs are classified in groups at five different levels. The drugs are divided into fourteen main groups .."
                        "Medicinal products are classified according to the main therapeutic use of the main active ingredient, on the  ..."
                    ],
                    "URL": "https://www.whocc.no/atc_ddd_index/",
                    "HNID": 1949731,
                    "ChildID": [
                        "node_1",
                        "node_803",
                        "node_1105",
                        "node_1888",
                        "node_2369",
                        "node_2727",
                        "node_2844",
                        "node_3500",
                        "node_3931",
                        "node_4209",
                        "node_4938",
                        "node_5112",
                        "node_5565",
                        "node_5911"
                    ],
                    "HasCountsOfType": [
                        "CID",
                        "SID"
                    ],
                    "Counts": [
                        {
                            "Type": "CID",
                            "Count": 15837
                        },
                        {
                            "Type": "SID",
                            "Count": 178510
                        }
                    ]
                    },
                    "Node": [
                    {
                        "NodeID": "node_41",
                        "ParentID": [
                            "node_38"
                        ],
                        "Information": {
                            "Name": "A01AD05 - Acetylsalicylic acid",
                            "URL": "https://www.whocc.no/atc_ddd_index/?code=A01AD05",
                            "HNID": 1949772,
                            "Match": true,
                            "Counts": [
                                {
                                "Type": "CID",
                                "Count": 4
                                },
                                {
                                "Type": "SID",
                                "Count": 84
                                }
                            ]
                        }
                    },
                    {
                    "NodeID": "node_38",
                    "ParentID": [
                        "node_3"
                    ],
                    "Information": {
                        "Name": "A01AD - Other agents for local oral treatment",
                        "URL": "https://www.whocc.no/atc_ddd_index/?code=A01AD",
                        "HNID": 1949769,
                        "ChildID": [
                            "node_39",
                            "node_40",
                            "node_41",
                            "node_42",
                            "node_43",
                            "node_44",
                            "node_45"
                        ],
                        "Counts": [
                            {
                            "Type": "CID",
                            "Count": 20
                            },
                            {
                            "Type": "SID",
                            "Count": 303
                            }
                        ]
                    }
                },

        """
        retL = []
        rD = {}
        rD["cid"] = cId
        try:
            tDL = self.__getKeyValue(recordD, "Hierarchies.Hierarchy")
            for tD in tDL:
                sName = self.__getKeyValue(tD, "SourceName")
                sId = self.__getKeyValue(tD, "SourceID")
                logger.debug("Finding sourceName %rsourceId %r", sName, sId)
                #
                provName = self.__getKeyValue(tD, "Information.Name")
                provDescr = self.__getKeyValue(tD, "Information.Description")
                provComments = self.__getKeyValue(tD, "Information.Description")
                provUrl = self.__getKeyValue(tD, "Information.URL")
                provD = {"source": provName, "description": provDescr, "details": provComments, "url": provUrl}
                #
                logger.debug("++++>KEY (%r,%r)", sName, provName)
                if (sName, provName) in [
                    ("WHO ATC", "ATC Code"),
                    ("FDA Pharm Classes", "FDA Pharmacological Classification"),
                    ("IUPHAR/BPS Guide to PHARMACOLOGY", "Target Classification"),
                    # ("ChEMBL", "Target Tree"),
                    # ("ChemIDplus", "ChemIDplus Chemical Information Classification"),
                    ("ChEBI", "ChEBI Ontology"),
                    ("MeSH", "MeSH Tree"),
                ]:
                    nDL = self.__getKeyValue(tD, "Node")
                    cD = {}
                    for nD in nDL:
                        nodeId = self.__getKeyValue(nD, "NodeID")
                        parentL = self.__getKeyValue(nD, "ParentID")
                        childL = self.__getKeyValue(nD, "Information.ChildID")
                        childCount = len(childL) if childL else 0
                        nm = self.__getKeyValue(nD, "Information.Name")
                        descr = self.__getKeyValue(nD, "Information.Description")
                        url = self.__getKeyValue(nD, "Information.URL")
                        idCode = None
                        # --- Add a bit to separate unique identifiers -
                        try:
                            if provName in ["ATC Code"] and nm:
                                fL = nm.split("-")
                                idCode = fL[0].strip()
                                nm = fL[1].strip()
                            elif provName in ["MeSH Tree"] and url:
                                ff = url.split("/")
                                idCode = ff[-1].strip()
                            elif sName == "IUPHAR/BPS Guide to PHARMACOLOGY" and url:
                                ff = url.split("=")
                                idCode = ff[-1].strip()
                        except Exception as e:
                            logger.exception("%s failing with %s", provName, str(e))
                        # ---
                        cD[nodeId] = {"parentList": parentL, "childCount": childCount, "idCode": idCode, "name": nm, "description": descr, "url": url}
                    rD[sName] = {"data": cD, "provenance": provD}
                #
            retL.append(rD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        # --- Process the extracted classification data ----
        try:
            ontClassL = []
            for ret in retL:
                logger.debug("Classification keys %r", ret.keys())
                ontClassD = {}
                ontClassD["cid"] = cId
                for ont, oD in ret.items():
                    childD = {}
                    parentD = {}
                    #
                    if ont in ["cid"]:
                        continue
                    logger.debug("%s node count %d", ont, len(oD["data"]))
                    rootNodeL = []
                    for nm, nD in oD["data"].items():
                        if "parentList" in nD and len(nD["parentList"]) > 1:
                            logger.info("Multiple parents for %s nd %s", ont, nm)
                            continue
                        pN = nD["parentList"][0]
                        if pN not in oD["data"]:
                            rootNodeL.append(nm)
                        childD.setdefault(nm, []).append(pN)
                        parentD.setdefault(pN, []).append(nm)
                    #
                    logger.debug("%s root nodes (%d)", ont, len(rootNodeL))
                    #
                    for ch, prL in childD.items():
                        if len(prL) > 1:
                            logger.info("%s multiple parents for node %s %d", ont, ch, len(prL))
                    # -- enumerate leaf nodes - nodes without children
                    leafL = []
                    for nm in oD["data"]:
                        if nm not in parentD:
                            logger.debug("%s leaf node %r - %s", ont, nm, oD["data"][nm]["name"])
                            leafL.append(nm)
                    #
                    linD = {}
                    for nm in leafL:
                        pN = childD[nm][0] if nm in childD else None
                        linD.setdefault(nm, []).append(nm)
                        while pN and pN in oD["data"]:
                            linD[nm].append(pN)
                            pN = childD[pN][0] if pN in childD else None
                    #
                    for leafN, linL in linD.items():
                        logger.debug("%s: %r : %r %r", ont, leafN, linL, [oD["data"][nm]["name"] for nm in linL])
                    #
                    # Get the unique lineages -
                    #
                    uTupD = {}
                    for leafN, linL in linD.items():
                        linL.reverse()
                        logger.debug("%s linL (%d) %r", ont, len(linL), linL)
                        tupL = [
                            (
                                oD["data"][nm]["name"],
                                oD["data"][nm]["description"][0] if oD["data"][nm]["description"] else None,
                                oD["data"][nm]["idCode"] if oD["data"][nm]["idCode"] else None,
                            )
                            for nm in linL
                        ]
                        logger.debug("%s tupL %r", ont, tupL)
                        tt = tuple(tupL)
                        uTupD[tt] = len(tt)
                    #
                    # Expand this out to
                    for uTup, tLen in uTupD.items():
                        logger.debug("%s (%d) %r", ont, tLen, [u[0] for u in uTup])
                        linL = [{"name": u[0], "description": u[1], "idCode": u[2]} for u in uTup]
                        ontClassD.setdefault(ont, []).append(linL)
                    #
                ontClassL.append(ontClassD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ontClassL

    def __getKeyValues(self, dct, keyNames):
        """Return the tuple of values of corresponding to the input dictionary key names expressed in dot notation.

        Args:
            dct (dict): source dictionary object (nested)
            keyNames (list): list of dictionary keys in dot notatoin

        Returns:
            tuple: tuple of values corresponding to the input key names

        """
        rL = []
        try:
            for keyName in keyNames:
                rL.append(self.__getKeyValue(dct, keyName))
        except Exception as e:
            logger.exception("Failing for key names %r with %s", keyNames, str(e))

        return tuple(rL)

    def __getKeyValue(self, dct, keyName):
        """Return the value of the corresponding key expressed in dot notation in the input dictionary object (nested)."""
        try:
            kys = keyName.split(".")
            for key in kys:
                try:
                    dct = dct[key]
                except KeyError:
                    return None
            return dct
        except Exception as e:
            logger.exception("Failing for key %r with %s", keyName, str(e))

        return None
