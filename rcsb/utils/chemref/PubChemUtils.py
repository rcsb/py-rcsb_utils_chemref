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

FeatureLabel = collections.namedtuple("FeatureLabel", "type description id evidence reference original variation")


class PubChemUtils(object):
    """
    Manage search and fetch queries for PubChem compound, substance and assay data and related annotations.



    """

    def __init__(self, **kwargs):
        self.__urlPrimary = kwargs.get("urlPrimary", "https://pubchem.ncbi.nlm.nih.gov")
        #

    def fetchCompound(self, identifier, identifierType="cid", searchType=None):
        ok = False
        if searchType == "view":
            ret, retCode = self.__doPugViewRequest(identifier, nameSpace=identifierType, domain="compound", requestType="GET", outputType="JSON")
        else:
            ret, retCode = self.__doPugRequest(identifier, nameSpace=identifierType, domain="compound", searchType=searchType, requestType="POST", outputType="JSON")
        ok = retCode in [200] and ret and len(ret) > 0
        return ok, ret

    def __doPugRequest(self, identifier, nameSpace="cid", domain="compound", searchType=None, requestType="GET", outputType="JSON"):
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
            hL = [("Accept", "application/json")]
            #
            pD = {}
            ureq = UrlRequestUtil()
            if nameSpace in ["cid", "name", "inchikey"] and not searchType and requestType == "GET":
                uId = quote(identifier.encode("utf8"))
                endPoint = "/".join(["rest", "pug", domain, nameSpace, uId, outputType])
                ret, retCode = ureq.get(baseUrl, endPoint, pD, headers=hL)
            elif nameSpace in ["cid", "name", "inchikey"] and not searchType and requestType == "POST":
                endPoint = "/".join(["rest", "pug", domain, nameSpace, outputType])
                pD = {nameSpace: identifier}
                ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
            elif nameSpace in ["cid", "name", "inchikey"] and searchType in ["xrefs", "property"] and requestType == "POST":
                endPoint = "/".join(["rest", "pug", domain, nameSpace, searchType, "RegistryID,RN", outputType])
                pD = {nameSpace: identifier}
                ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
            elif nameSpace in ["smiles", "inchi",] and searchType in ["fastidentity", "fastsimilarity_2d"]:
                endPoint = "/".join(["rest", "pug", domain, searchType, nameSpace, outputType])
                pD = {nameSpace: identifier}
                ret, retCode = ureq.post(baseUrl, endPoint, pD, headers=hL)
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
            hL = [("Accept", "application/json")]
            #
            pD = {}
            ureq = UrlRequestUtil()
            if nameSpace in ["cid"] and requestType == "GET":
                uId = quote(identifier.encode("utf8"))
                endPoint = "/".join(["rest", "pug_view", "data", domain, uId, outputType])
                ret, retCode = ureq.get(baseUrl, endPoint, pD, headers=hL)
            elif nameSpace in ["cid",] and requestType == "POST":
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
                logger.debug("TOCHeading %r Description %r ...", tD["TOCHeading"], tD["Description"][:30])
                if "Section" in tD:
                    for ttD in tD["Section"]:
                        # logger.debug("ttD.keys() %r", ttD.keys())
                        if "Description" in ttD:
                            logger.debug("  >> TOCHeading %r Description %r ...", ttD["TOCHeading"], ttD["Description"][:30])
                        else:
                            logger.debug("  >> TOCHeading %r Description None ...", ttD["TOCHeading"])
                        if "Section" in ttD:
                            for tttD in ttD["Section"]:
                                if "Description" in tttD:
                                    logger.debug("        >> TOCHeading %r Description %r ...", tttD["TOCHeading"], tttD["Description"][:30])
                                else:
                                    logger.debug("        >> TOCHeading %r Description  None ...", tttD["TOCHeading"])
                                if "Section" in tttD:
                                    for ttttD in tttD["Section"]:
                                        if "Description" in ttttD:
                                            logger.debug("           >> TOCHeading %r Description %r ...", ttttD["TOCHeading"], ttttD["Description"][:30])
                                        else:
                                            logger.debug("           >> TOCHeading %r Description  None ...", ttttD["TOCHeading"])
                                        if "Section" in ttttD:
                                            for tttttD in ttttD["Section"]:
                                                if "Description" in tttttD:
                                                    logger.debug("              >> TOCHeading %r Description %r ...", tttttD["TOCHeading"], tttttD["Description"][:30])
                                                else:
                                                    logger.debug("              >> TOCHeading %r Description  None ...", tttttD["TOCHeading"])

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
        try:
            qD = vD["Record"] if "Record" in vD else {}
            logger.debug("keys %s", qD.keys())
            logger.debug("Record.RecordType %r", self.__getKeyValues(vD, ["Record.RecordType"]))
            logger.debug("Record.RecordNumber %r", self.__getKeyValues(vD, ["Record.RecordNumber"]))
            logger.debug("Record.RecordTitle %r", self.__getKeyValues(vD, ["Record.RecordTitle"]))

            for sectionD in qD["Section"]:
                logger.debug("TOCHeading %r Description %r ...", sectionD["TOCHeading"], sectionD["Description"][:30])
                #
                if sectionD["TOCHeading"] in ["Names and Identifiers"] and "Section" in sectionD:
                    for subSectionD in sectionD["Section"]:
                        logger.debug("  >> TOCHeading %r", subSectionD["TOCHeading"])
                        # ---
                        if subSectionD["TOCHeading"] in ["Computed Descriptors"] and "Section" in subSectionD:
                            for subSubSectionD in subSectionD["Section"]:
                                logger.debug("        >> TOCHeading %r", subSubSectionD["TOCHeading"])
                                if subSubSectionD["TOCHeading"] in ["IUPAC Name", "InChI", "InChI Key", "Canonical SMILES"]:
                                    ky = subSubSectionD["TOCHeading"]
                                    try:
                                        for infD in subSubSectionD["Information"]:
                                            pcD.setdefault(ky, []).append(infD["Value"]["StringWithMarkup"][0]["String"])
                                    except Exception:
                                        pass
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
