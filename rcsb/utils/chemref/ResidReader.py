##
# File:    ResidReader.py
# Author:  J. Westbrook
# Date:    17-Mar-2020
# Version: 0.001
#
# Update:
#
##
"""
Reader for the RESID repository data (RESIDUES.xml).
"""
import logging

logger = logging.getLogger(__name__)


class ResidReader(object):
    """Reader for the RESID repository data (RESIDUES.xml)."""

    def __init__(self):
        self.__ns = ""

    def read(self, xTree):
        rL = []
        xRoot = xTree.getroot()
        # <Database id="RESID" release="76.00" date="31-May-2018">
        version = xRoot.attrib.get("release", None)
        for el in xRoot:
            if el.tag == "Entry":
                rL.append(self.__processEntryElement(el))
        return version, rL

    def __processEntryElement(self, entryElement):
        """ """
        assert entryElement.tag == "{ns}Entry".format(ns=self.__ns)
        #
        doc = {
            "residCode": entryElement.findtext("{ns}Header/{ns}Code".format(ns=self.__ns)),
            "names": [name.text for name in entryElement.findall("{ns}Names/{ns}Name".format(ns=self.__ns))],
            "nameXrefs": [xref.text for xref in entryElement.findall("{ns}Names/{ns}Xref".format(ns=self.__ns))],
            "ontRefs": [xref.text for xref in entryElement.findall("{ns}SequenceCode/{ns}Xref".format(ns=self.__ns))],
            "genEnzymes": [xref.text for xref in entryElement.findall("{ns}GeneratingEnzyme/{ns}EnzymeName".format(ns=self.__ns))],
            "features": [feature.text for feature in entryElement.findall("{ns}Features/{ns}Feature".format(ns=self.__ns))],
        }
        return doc
