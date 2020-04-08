##
# File:    PubChemUtils.py
# Date:    30-Mar-2020
#
# Updates:
#
##

import collections

import json
import logging

from rcsb.utils.io.UrlRequestUtil import UrlRequestUtil

try:
    from itertools import zip_longest
except ImportError:
    # Python 2
    from itertools import izip_longest as zip_longest

try:
    from urllib.parse import quote
except ImportError:
    from urllib2 import quote

logger = logging.getLogger(__name__)

ExternalAnnotationFields = ("featureType", "value", "description", "dataType", "units", "provReferenceId", "provSourceName")
ExternalAnnotation = collections.namedtuple("ExternalAnnotation", ExternalAnnotationFields, defaults=(None,) * len(ExternalAnnotationFields))

ChemicalIdentifierFields = ("idCode", "identifierSource", "identifierType", "identifier")
ChemicalIdentifier = collections.namedtuple("ChemicalIdentifier", ChemicalIdentifierFields, defaults=(None,) * len(ChemicalIdentifierFields))


class PubChemUtils(object):
    """ Manage search and fetch queries for PubChem compound, substance and assay data and related annotations.

    """

    def __init__(self, **kwargs):
        self.__urlPrimary = kwargs.get("urlPrimary", "https://pubchem.ncbi.nlm.nih.gov")
        #

    def fetchList(self, chemicalIdentifierList, searchType=None):
        """Fetch list of

           For tuple identifier types the first and last elements are used: (ky, ... , searchId).

        Args:
            chemicalIdentifierList (ChemicalIdentifier Tuple): list of chemical identifier tuples
                identifierType (string): search identifier type (cid|InChI|InChIKey|SMILES)
            searchType (string): search return type
        """
        retD = {}
        ok = True
        try:
            for chemicalIdentifier in chemicalIdentifierList:
                rOk, rD = self.fetch(chemicalIdentifier, searchType=searchType)
                if ok:
                    retD[chemicalIdentifier.idCode] = rD
                    ok = rOk and ok

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok, retD

    def fetch(self, chemicalIdentifier, searchType="lookup", returnType="record"):
        """Fetch data (returnType) for the input chemical identifier using the input searchType.

        Args:
            chemicalIdentifer (ChemicalIdentifier Tuple): chemical identifier tuple
                    identifierType (string): search identifier type (cid|InChI|InChIKey|SMILES)
            searchType (string): search type (lookup|fastsimilarity|fastidentity)
            returnType (string): object type to return (default: record)

        Return:
            (bool, dict): status, return object
        """
        ok = False
        if returnType == "view":
            ret, retCode = self.__doPugViewRequest(
                chemicalIdentifier.identifier, nameSpace=chemicalIdentifier.identifierType, domain="compound", requestType="GET", outputType="JSON"
            )
        else:
            ret, retCode = self.__doPugRequest(
                chemicalIdentifier.identifier,
                nameSpace=chemicalIdentifier.identifierType,
                domain="compound",
                searchType=searchType,
                returnType=returnType,
                requestType="POST",
                outputType="JSON",
            )
        ok = retCode in [200] and ret and len(ret) > 0
        return ok, ret

    def __doPugRequest(self, identifier, nameSpace="cid", domain="compound", searchType="lookup", returnType="record", requestType="GET", outputType="JSON"):
        """ Wrapper for PubChem PUG API requests

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
        response, retCode, ret = None, None, None
        try:
            baseUrl = self.__urlPrimary
            hL = []
            if outputType == "JSON":
                hL.append(("Accept", "application/json"))
            #
            pD = {}
            ureq = UrlRequestUtil()
            if domain == "compound":
                if nameSpace in ["cid", "name", "inchikey"] and returnType in ["record"] and searchType in ["lookup"] and requestType == "GET":
                    uId = quote(identifier.encode("utf8"))
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, uId, outputType])
                    ret, retCode = ureq.get(baseUrl, endPoint, pD, headers=hL)
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["record"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
                #
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["xrefs"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, "RegistryID,RN", outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["property"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, "MolecularFormula,XLogP,TPSA,Volume3D", outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["classification"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, outputType])
                    pD = {nameSpace: identifier}
                    ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
                elif nameSpace in ["cid", "name", "inchikey"] and returnType in ["synonyms"] and searchType in ["lookup"] and requestType == "POST":
                    endPoint = "/".join(["rest", "pug", domain, nameSpace, returnType, outputType])
                    pD = {nameSpace: identifier, "classification_type": "simple"}
                    ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
                #
                elif nameSpace in ["smiles", "inchi"] and returnType == "record" and searchType in ["fastidentity", "fastsimilarity_2d"]:
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

                    ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
            elif domain in ["substance", "assay"]:
                logger.error("Fetch not implemented for domain %s", domain)
            #
            if ret and outputType in ["JSON"]:
                response = json.loads(ret)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return response, retCode

    def __doPugViewRequest(self, identifier, nameSpace="cid", domain="compound", requestType="GET", outputType="JSON"):
        """ Wrapper for PubChem PUG_VIEW API requests

            Example:
                https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/1234/JSON

        """
        response, retCode, ret = None, None, None
        try:
            baseUrl = self.__urlPrimary
            hL = []
            if outputType == "JSON":
                hL.append(("Accept", "application/json"))
            #
            pD = {}
            ureq = UrlRequestUtil()
            if nameSpace in ["cid"] and requestType == "GET":
                uId = quote(identifier.encode("utf8"))
                endPoint = "/".join(["rest", "pug_view", "data", domain, uId, outputType])
                ret, retCode = ureq.get(baseUrl, endPoint, pD, headers=hL)
            elif nameSpace in ["cid"] and requestType == "POST":
                endPoint = "/".join(["rest", "pug_view", "data", domain, outputType])
                pD = {nameSpace: identifier}
                ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
            if ret and outputType in ["JSON"]:
                response = json.loads(ret)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return response, retCode

    def __makeSubLists(self, num, iterable):
        args = [iter(iterable)] * num
        return ([e for e in t if e is not None] for t in zip_longest(*args))

    def __makeSubListsWithPadding(self, num, iterable, padvalue=None):
        "__sublist(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
        return zip_longest(*[iter(iterable)] * num, fillvalue=padvalue)

    def traversePubChemCompoundView(self, vD):
        """
            "Record": {
                    "RecordType": "CID",
                    "RecordNumber": 123631,
                    "RecordTitle": "Gefitinib",
                    "Section": []

        """
        try:
            qD = vD["Record"] if "Record" in vD else {}
            logger.debug("keys %s", qD.keys())
            logger.debug("Record.RecordType %r", self.__getKeyValues(vD, ["Record.RecordType"]))
            logger.debug("Record.RecordNumber %r", self.__getKeyValues(vD, ["Record.RecordNumber"]))
            logger.debug("Record.RecordTitle %r", self.__getKeyValues(vD, ["Record.RecordTitle"]))

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

    def parsePubChemCompoundView(self, vD):
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
                                    }

        """
        pcD = {}
        sectionMapD = {
            "Computed Descriptors": "descriptors",
            "Other Identifiers": "identifiers",
            "Synonyns": "synonyms",
            "Computed Properties": "properties",
            "Experimental Properties": "properties",
        }
        mainSectionNameD = {
            "Names and Identifiers": {"Computed Descriptors": True, "Other Identifiers": True, "Synonyns": True,},
            "Chemical and Physical Properties": {"Computed Properties": True, "Experimental Properties": True},
        }
        #
        try:
            qD = vD["Record"] if "Record" in vD else {}
            logger.debug("keys %s", qD.keys())
            logger.debug("Record.RecordType %r", self.__getKeyValues(vD, ["Record.RecordType"]))
            logger.debug("Record.RecordNumber %r", self.__getKeyValues(vD, ["Record.RecordNumber"]))
            logger.debug("Record.RecordTitle %r", self.__getKeyValues(vD, ["Record.RecordTitle"]))

            for sectionD in qD["Section"]:
                logger.info("TOCHeading %r Description %r ...", sectionD["TOCHeading"], sectionD["Description"][:30])
                if "TOCHeading" in sectionD and sectionD["TOCHeading"] in mainSectionNameD:
                    logger.info("Reading section %r Description %r ...", sectionD["TOCHeading"], sectionD["Description"][:30])
                    if "Section" in sectionD:
                        for subSectionD in sectionD["Section"]:
                            logger.info("  >> subSection %r", subSectionD["TOCHeading"])

                            pcKy = sectionMapD[subSectionD["TOCHeading"]] if subSectionD["TOCHeading"] in sectionMapD else "misc"
                            # ---
                            if "Section" in subSectionD:
                                for subSubSectionD in subSectionD["Section"]:
                                    ky = subSubSectionD["TOCHeading"]
                                    logger.info("        >> subsubSection %r (%r)", ky, pcKy)
                                    try:
                                        for infD in subSubSectionD["Information"]:
                                            referenceNumber = infD["ReferenceNumber"] if "ReferenceNumber" in infD else None
                                            reference = ",".join(infD["Reference"]) if "Reference" in infD else None
                                            logger.info("reference no %r reference %r", referenceNumber, reference)
                                            if "Value" in infD:
                                                if "StringWithMarkup" in infD["Value"]:
                                                    for sV in infD["Value"]["StringWithMarkup"]:
                                                        tupV = ExternalAnnotation(
                                                            value=sV["String"], dataType="string", featureType=ky, provReferenceId=referenceNumber, provSourceName=reference
                                                        )
                                                        pcD.setdefault(pcKy, {}).setdefault(ky, set()).update([tupV])
                                                elif "Number" in infD["Value"]:
                                                    units = infD["Value"]["Unit"] if "Unit" in infD["Value"] else ""
                                                    for nV in infD["Value"]["Number"]:
                                                        tupV = ExternalAnnotation(
                                                            value=nV, dataType="number", units=units, featureType=ky, provReferenceId=referenceNumber, provSourceName=reference
                                                        )
                                                        # tV = str(str(sV) + " " + units).strip()
                                                        pcD.setdefault(pcKy, {}).setdefault(ky, set()).update([tupV])
                                    except Exception as e:
                                        logger.exception("Failing with %s", str(e))
                                        # pass
                        # ---
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return pcD

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
        """  Return the value of the corresponding key expressed in dot notation in the input dictionary object (nested).
        """
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
