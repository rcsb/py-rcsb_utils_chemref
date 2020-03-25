##
# -*- coding: utf-8 -*-
#
# File:    PsiModProvider.py
# Author:  J. Westbrook
# Date:    19-Mar-2020
# Version: 0.001
#
# Update:
#
#
##
"""
Various utilities for extracting data from the PSI-MOD OBO files
and returning lineage details.
"""

import logging
import os
from collections import OrderedDict

import networkx
import obonet

from rcsb.utils.io.FileUtil import FileUtil

logger = logging.getLogger(__name__)


class PsiModProvider(object):
    """Various utilities for extracting data from the PSI-MOD OBO files
        and returning lineage details.
    """

    def __init__(self, **kwargs):
        urlTarget = kwargs.get("urlTarget", "http://data.bioontology.org/ontologies/PSIMOD/submissions/6/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb")
        dirPath = os.path.join(kwargs.get("cachePath", "."), "psimod")
        useCache = kwargs.get("useCache", True)
        self.__modGraph = self.__reload(urlTarget, dirPath, useCache=useCache)

    def testCache(self):
        if self.__modGraph:
            logger.info("Reading %d nodes and %d edges", len(self.__modGraph), self.__modGraph.number_of_edges())
            if networkx.is_directed_acyclic_graph(self.__modGraph) and len(self.__modGraph) > 1000:
                return True
        return False

    def getVersion(self):
        try:
            return self.__modGraph.graph.get("data-version", "1.013.0")
        except Exception:
            return None

    def exists(self, psimodId):
        try:
            return psimodId in self.__modGraph
        except Exception:
            return False

    def getNode(self, psimodId):
        try:
            return self.__modGraph[psimodId]
        except Exception:
            pass
        return None

    def getRootNodes(self):
        try:
            rootL = [n for n, d in self.__modGraph.out_degree() if d == 0]
            return rootL
        except Exception:
            pass
        return None

    def getName(self, psimodId):
        try:
            return self.__modGraph.nodes[psimodId]["name"]
        except Exception as e:
            logger.debug("Failing %r with %r", psimodId, str(e))
        return None

    def getAdjacentParents(self, psimodId):
        rL = []
        try:
            for child, parent, key in self.__modGraph.out_edges(psimodId, keys=True):
                logger.debug("%r %r - %r -> %r %r", child, self.getName(child), key, parent, self.getName(parent))
                rL.append((child, parent, key))
        except Exception:
            pass
        return rL

    def getPredecessors(self, psimodId):
        rL = []
        try:
            rL = [nd for nd in self.__modGraph.predecessors(psimodId)]
        except Exception:
            pass
        return rL

    def getSuccessors(self, psimodId):
        rL = []
        try:
            rL = [nd for nd in self.__modGraph.successors(psimodId)]
        except Exception:
            pass
        return rL

    def getDescendants(self, psimodId, includeSelf=True):
        linL = []
        try:
            if includeSelf:
                linL.append((psimodId, self.getName(psimodId)))
            for nd in networkx.descendants(self.__modGraph, psimodId):
                logger.debug("%r %r --> %r %r", psimodId, self.getName(psimodId), nd, self.getName(nd))
                linL.append((nd, self.getName(nd)))
        except Exception as e:
            logger.debug("Failing %s with %s", psimodId, str(e))
        return linL

    def getUniqueDescendants(self, psimodIdL, includeSelf=True):
        linL = []
        try:
            ndD = OrderedDict()
            for psimodId in psimodIdL:
                if includeSelf:
                    ndD[psimodId] = True
                for nd in networkx.descendants(self.__modGraph, psimodId):
                    ndD[nd] = True
            #
            for nd in ndD:
                linL.append((nd, self.getName(nd)))
        except Exception as e:
            logger.debug("Failing %s with %s", psimodId, str(e))
        return linL

    def getLineage(self, psiModId):
        linL = []
        lIdL = self.getUniqueDescendants([psiModId])
        for ii, lTup in enumerate(lIdL, 1):
            linL.append((lTup[0], lTup[1], ii))
        return linL

    def exportTreeNodeList(self, psimodIdL):
        """ For the input node list export full tree node list including parent nodes

        Args:
            psimodIdL (list): list of covered GO id lists
            includeSelf (bool, optional): include input and parent nodes. Defaults to True.

        Returns:
            list: [{'id': <> 'name': <name> 'parents': [<id>,<id>,...]}]
        """
        trL = []
        try:
            # Generate the full list of nodes and parents -
            ndS = set()
            for psimodId in psimodIdL:
                if not self.exists(psimodId):
                    logger.warning("%s not in current ontology", psimodId)
                    continue
                ndS.add(psimodId)
                for nd in networkx.descendants(self.__modGraph, psimodId):
                    ndS.add(nd)
            #
            for nd in ndS:
                # tupL = [(child, parent, key)]
                pIdL = [tup[1] for tup in self.getAdjacentParents(nd)]
                if not pIdL:
                    logger.info("Node %s (%s) has no parents", nd, self.getName(nd))
                    trL.append({"id": nd, "name": self.getName(nd)})
                else:
                    trL.append({"id": nd, "name": self.getName(nd), "parents": pIdL})

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return trL

    def __reload(self, urlTarget, dirPath, useCache=True):
        """ Reload input GO OBO ontology file and return a nx graph object.
'
        Returns:
            dictionary[psimodId] = {'name_list': ... , 'id_list': ... 'depth_list': ... }
        """
        goGraph = None
        #
        # mU = MarshalUtil()
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        oboFilePath = os.path.join(dirPath, fn)
        fU.mkdir(dirPath)
        #
        if not useCache:
            for fp in [oboFilePath]:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and fU.exists(oboFilePath):
            goGraph = obonet.read_obo(oboFilePath)
        else:
            logger.info("Fetching url %s to resource file %s", urlTarget, oboFilePath)
            ok = fU.get(urlTarget, oboFilePath)
            if ok:
                goGraph = obonet.read_obo(oboFilePath)
        if goGraph:
            logger.info("Reading %d nodes and %d edges", len(goGraph), goGraph.number_of_edges())
        else:
            logger.info("Go graph construction failing")
        #
        return goGraph
