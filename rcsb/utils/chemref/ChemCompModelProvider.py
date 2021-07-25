##
# File:    ChemCompModelProvider.py
# Author:  J. Westbrook
# Date:    1-Nov-2018
#
# Updates:
# 21-Jul-2021 jdw  Make this provider a subclass of StashableBase
##
"""
Utilities to read resource file containing compilation of CCDC models correspondences
for PDB chemical components.

"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class ChemCompModelProvider(StashableBase):
    """Utilities to read the resource file containing with a compilation of CCDC models correspondences
    for PDB chemical components.
    """

    def __init__(self, **kwargs):
        dirName = "chem_comp_model"
        cachePath = kwargs.get("cachePath", ".")
        super(ChemCompModelProvider, self).__init__(cachePath, [dirName])

        urlTarget = kwargs.get("urlTarget", "http://ftp.wwpdb.org/pub/pdb/data/component-models/complete/chem_comp_model.cif.gz")
        dirPath = os.path.join(cachePath, dirName)
        useCache = kwargs.get("useCache", True)
        mappingFileName = kwargs.get("mappingFileName", "ccdc_pdb_mapping.json")
        auditFileName = kwargs.get("mappingFileName", "ccdc_model_audit.json")
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__mappingD, self.__auditD = self.__reload(urlTarget, dirPath, mappingFileName, auditFileName, useCache=useCache)

    def testCache(self):
        if len(self.__mappingD) > 1000 and len(self.__auditD) > 1000:
            logger.info("chem_comp_model map length %d", len(self.__mappingD))
            return True
        return False

    def getMapping(self):
        return self.__mappingD

    def getAuditDetails(self):
        return self.__auditD

    def __reload(self, urlTarget, dirPath, mappingFileName, auditFileName, useCache=True):
        """Reload input mmCIF model mapping resource file and return a container list.

        Args:
            urlTarget (str): target url for resource file
            dirPath (str): path to the directory containing cache files
            mappingFileName (str): mapping file name
            auditFileName (str): audit details file name
            useCache (bool, optional): flag to use cached files. Defaults to True.

        Returns:
            (dict,dict): mapping dictionary (pdb_ccId -> CCDC details), audit details (ccId -> audit details)
        """
        mD = {}
        aD = {}
        #
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        filePath = os.path.join(dirPath, fn)
        mappingFilePath = os.path.join(dirPath, mappingFileName)
        auditFilePath = os.path.join(dirPath, auditFileName)
        self.__mU.mkdir(dirPath)
        #
        if not useCache:
            for fp in [filePath, mappingFilePath, auditFilePath]:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and fU.exists(mappingFilePath) and fU.exists(auditFilePath):
            mD = self.__mU.doImport(mappingFilePath, fmt="json")
            aD = self.__mU.doImport(auditFilePath, fmt="json")
        elif not useCache:
            ok = True
            if not (useCache and fU.exists(filePath)):
                logger.info("Fetching url %s for resource file %s", urlTarget, filePath)
                ok = fU.get(urlTarget, filePath)
            if ok:
                cL = self.__mU.doImport(filePath, fmt="mmcif")
                mD = self.__buildMapping(cL)
                ok = self.__mU.doExport(mappingFilePath, mD, fmt="json", indent=3)
                aD = self.__buildAuditDetails(cL)
                ok = self.__mU.doExport(auditFilePath, aD, fmt="json", indent=3)
        #
        return mD, aD

    def __buildMapping(self, cL):
        """Return a dictionary of identfier correspondences for CCDC/CSD models.

        data_M_010_00001
            #
            _pdbx_chem_comp_model.id        M_010_00001
            _pdbx_chem_comp_model.comp_id   010
            #
            _pdbx_chem_comp_model_reference.model_id   M_010_00001
            _pdbx_chem_comp_model_reference.db_name    CSD
            _pdbx_chem_comp_model_reference.db_code    YARXEW
            #
            loop_
            _pdbx_chem_comp_model_feature.model_id
            _pdbx_chem_comp_model_feature.feature_name
            _pdbx_chem_comp_model_feature.feature_value
            M_010_00001 experiment_temperature 123.0
            M_010_00001 publication_doi        10.1021/cg200971f
            M_010_00001 r_factor               3.91
            M_010_00001 all_atoms_have_sites   Y
        #
        """
        rD = {}
        for dataContainer in cL:
            logger.debug("Processing model %r", dataContainer.getName())
            cObj = dataContainer.getObj("pdbx_chem_comp_model")
            ccId = cObj.getValue("comp_id", 0)
            #
            if ccId not in rD:
                rD[ccId] = []

            tObj = dataContainer.getObj("pdbx_chem_comp_model_reference")
            for iRow in range(tObj.getRowCount()):
                dbName = tObj.getValue("db_name", iRow)
                dbCode = tObj.getValue("db_code", iRow)
                rD[ccId].append({"db_name": dbName, "db_code": dbCode})

        return rD

    def __buildAuditDetails(self, cL):
        """Return a dictionary of identfier correspondences for CCDC/CSD models.

                data_M_4YS_00001
                #
                _pdbx_chem_comp_model.id        M_4YS_00001
                _pdbx_chem_comp_model.comp_id   4YS
                #

                _pdbx_chem_comp_model_reference.model_id   M_4YS_00001
                _pdbx_chem_comp_model_reference.db_name    CSD
                _pdbx_chem_comp_model_reference.db_code    CIQKUJ
                #
                loop_
                _pdbx_chem_comp_model_feature.model_id
                _pdbx_chem_comp_model_feature.feature_name
                _pdbx_chem_comp_model_feature.feature_value
                M_4YS_00001 experiment_temperature '290 K'
                M_4YS_00001 publication_doi        10.1107/S1600536807057960
                M_4YS_00001 r_factor               4.230
                M_4YS_00001 csd_version            542
                #
                _pdbx_chem_comp_model_audit.model_id      M_4YS_00001
                _pdbx_chem_comp_model_audit.action_type   'Initial release'
                _pdbx_chem_comp_model_audit.date          2021-02-03

        #
        """
        rD = {}
        for dataContainer in cL:
            logger.debug("Processing model %r", dataContainer.getName())
            cObj = dataContainer.getObj("pdbx_chem_comp_model")
            ccId = cObj.getValue("comp_id", 0)
            modelId = cObj.getValue("id", 0)
            #
            tObj = dataContainer.getObj("pdbx_chem_comp_model_audit")
            aL = []
            for iRow in range(tObj.getRowCount()):
                auditAction = tObj.getValue("action_type", iRow)
                auditDate = tObj.getValue("date", iRow)
                aL.append({"audit_date": auditDate, "action_type": auditAction})

            tObj = dataContainer.getObj("pdbx_chem_comp_model_reference")
            dbName = dbCode = None
            for iRow in range(tObj.getRowCount()):
                dbName = tObj.getValue("db_name", iRow)
                dbCode = tObj.getValue("db_code", iRow)

            rD.setdefault(ccId, []).append({"model_id": modelId, "db_name": dbName, "db_code": dbCode, "audit_list": aL})

        return rD
