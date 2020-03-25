##
# File:    ChemCompModelProvider.py
# Author:  J. Westbrook
# Date:    1-Nov-2018
#
# Updates:
# http://ftp.wwpdb.org/pub/pdb/data/component-models/complete/chem_comp_model.cif.gz
##
"""
Utilities to read resource file containing compilation of CCDC models correspondences
for PDB chemical components.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChemCompModelProvider(object):
    """Utilities to read the resource file containing with a compilation of CCDC models correspondences
    for PDB chemical components.
    """

    def __init__(self, **kwargs):
        urlTarget = kwargs.get("urlTarget", "http://ftp.wwpdb.org/pub/pdb/data/component-models/complete/chem_comp_model.cif.gz")
        dirPath = os.path.join(kwargs.get("cachePath", "."), "chem_comp")
        useCache = kwargs.get("useCache", True)
        mappingFileName = kwargs.get("mappingFileName", "ccdc_pdb_mapping.json")
        #
        self.__mU = MarshalUtil(workPath=dirPath)
        self.__mappingD = self.__reload(urlTarget, dirPath, mappingFileName, useCache=useCache)

    def testCache(self):
        logger.info("Lengths map %d", len(self.__mappingD))
        if len(self.__mappingD) > 1000:
            return True
        return False

    def getMapping(self):
        return self.__mappingD

    def __reload(self, urlTarget, dirPath, mappingFileName, useCache=True):
        """Reload input mmCIF model mapping resource file and return a container list.

        Args:
            urlTarget (str): target url for resource file
            dirPath (str): path to the directory containing cache files
            mappingFileName (str): mapping file name
            useCache (bool, optional): flag to use cached files. Defaults to True.

        Returns:
            (dict): mapping dictionary (pdb_ccId -> CCDC details)
        """
        mD = {}
        #
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        filePath = os.path.join(dirPath, fn)
        mappingFilePath = os.path.join(dirPath, mappingFileName)
        self.__mU.mkdir(dirPath)
        #
        if not useCache:
            for fp in [filePath, mappingFilePath]:
                try:
                    os.remove(fp)
                except Exception:
                    pass
        #
        if useCache and fU.exists(mappingFilePath):
            mD = self.__mU.doImport(mappingFilePath, fmt="json")
        else:
            ok = True
            if not (useCache and fU.exists(filePath)):
                logger.info("Fetching url %s for resource file %s", urlTarget, filePath)
                ok = fU.get(urlTarget, filePath)
            if ok:
                cL = self.__mU.doImport(filePath, fmt="mmcif")
                mD = self.__buildMapping(cL)
                ok = self.__mU.doExport(mappingFilePath, mD, fmt="json", indent=3)
        #
        return mD

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
