##
#  File:           ChEMBLProvider.py
#  Date:           9-Nov-2020 jdw
#
#  Updated:
#  2-Dec-2020 jdw Add ChEMBL API access methods
##
"""
Accessors for ChEMBL small molecule data.

"""
# pylint: disable=line-too-long
import logging
import os.path
import time

from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.unichem import unichem_client

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logger = logging.getLogger(__name__)


class ChEMBLProvider:
    """Accessors for ChEMBL small molecule data."""

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirPath = os.path.join(cachePath, "CACHE", "ChEMBL")

        chemblDbUrl = kwargs.get("ChEMBLDbUrl", "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__retD = self.__reload(chemblDbUrl, self.__dirPath, useCache=useCache)
        #

    def testCache(self):
        if self.__retD and len(self.__retD) > 10000:
            logger.info("ChEMBL mappings (%d)", len(self.__retD))
            return True
        return False

    def getInChIKey(self, chemblId):
        rVal = None
        try:
            rVal = self.__retD[chemblId]["inchikey"]
        except Exception:
            pass
        return rVal

    def getSMILES(self, chemblId):
        rVal = None
        try:
            rVal = self.__retD[chemblId]["smiles"]
        except Exception:
            pass
        return rVal

    def __reload(self, chemblDbUrl, dirPath, useCache=True, fmt="json"):
        startTime = time.time()
        ok = False
        fU = FileUtil()
        fU.mkdir(dirPath)
        #
        baseVersion = 29
        # ChEMBL current version 27,...
        # template:  chembl_27_chemreps.txt.gz
        #
        inpFileName = "chembl_" + str(baseVersion) + ".fa.gz"
        #

        chemblDataPath = os.path.join(dirPath, "chembl_data." + fmt)

        #
        if useCache and self.__mU.exists(chemblDataPath):
            logger.info("useCache %r using %r", useCache, chemblDataPath)
            retD = self.__mU.doImport(chemblDataPath, fmt=fmt)
            ok = True
        else:
            #
            for vers in range(baseVersion, baseVersion + 10):
                inpFileName = "chembl_" + str(vers) + "_chemreps.txt.gz"
                inpPath = os.path.join(dirPath, inpFileName)
                url = os.path.join(chemblDbUrl, inpFileName)
                ok = fU.get(url, inpPath)
                logger.info("Fetching url %s path %s", url, inpPath)
                if ok:
                    break
            #
            logger.info("Completed fetches at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
            retD = self.__chemblParse(inpPath)
            self.__mU.doExport(chemblDataPath, retD, fmt=fmt, indent=3)
        return retD

    def __chemblParse(self, filePath):
        retD = {}
        rowL = self.__mU.doImport(filePath, fmt="tdd", rowFormat="list")
        logger.info("ChEMBL data length %d", len(rowL))
        for row in rowL[1:]:
            if len(row) < 4:
                logger.debug("Bad row %r", row)
                continue
            retD[row[0]] = {"smiles": row[1], "inchiKey": row[3]}
        return retD

    def getDrugData(self, moleculeChEMBLIdList):
        """Get drug data for the input ChEMBL molecule list.

        Args:
            moleculeChEMBLIdList (list): list of ChEMBL molecule identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {drug data}}

        """
        oD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(moleculeChEMBLIdList), chunkSize):
                drug = new_client.drug  # pylint: disable=no-member
                drug.set_format("json")
                mDL = drug.filter(molecule_chembl_id__in=moleculeChEMBLIdList[ii : ii + chunkSize])
                if mDL:
                    logger.info("mDL (%d)", len(mDL))
                    for mD in mDL:
                        oD.setdefault(mD["molecule_chembl_id"], []).append(mD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD

    def getMoleculeData(self, moleculeChEMBLIdList):
        """Get molecule data for the input ChEMBL molecule list.

        Args:
            moleculeChEMBLIdList (list): list of ChEMBL molecule identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {molecule data}}

        """
        oD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(moleculeChEMBLIdList), chunkSize):
                drug = new_client.molecule  # pylint: disable=no-member
                drug.set_format("json")
                mDL = drug.filter(molecule_chembl_id__in=moleculeChEMBLIdList[ii : ii + chunkSize])
                if mDL:
                    logger.info("mDL (%d)", len(mDL))
                    for mD in mDL:
                        oD.setdefault(mD["molecule_chembl_id"], []).append(mD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD

    def getMoleculeByInChIKey(self, inchiKeyList):
        """Get molecule data for the input InChI Key list.

        Args:
            inchiKeyList (list): list of InChI Key identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {molecule data}}

        """
        oD = {}
        chunkSize = 50
        try:
            for ii in range(0, len(inchiKeyList), chunkSize):
                drug = new_client.molecule  # pylint: disable=no-member
                drug.set_format("json")
                mDL = drug.get(inchiKeyList[ii : ii + chunkSize])
                if mDL:
                    logger.info("mDL (%d)", len(mDL))
                    for mD in mDL:
                        oD.setdefault(mD["molecule_chembl_id"], []).append(mD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD

    def getUniChemData(self, inchiKeyList):
        """Get UniChem data for the input InChiKey list.

        Args:
            InChIList (list): list of InChI key molecule identifiers

        Returns:
          (dict):  dictionary  {ChEMBId: {molecule data}}

        """
        mapD = {
            1: {"name": "chembl", "baseUrl": "https://www.ebi.ac.uk/chembl/", "entryUrl": "https://www.ebi.ac.uk/chembldb/compound/inspect/"},
            3: {"name": "pdb", "baseUrl": "http://www.ebi.ac.uk/pdbe/", "entryUrl": "http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/"},
            2: {"name": "drugbank", "baseUrl": "http://drugbank.ca/", "entryUrl": "http://www.drugbank.ca/drugs/"},
            5: {"name": "pubchem_dotf", "baseUrl": "http://pubchem.ncbi.nlm.nih.gov/sources/sources.cgi", "entryUrl": "http://pubchem.ncbi.nlm.nih.gov/substance/"},
            4: {"name": "gtopdb", "baseUrl": "http://www.guidetopharmacology.org", "entryUrl": "http://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId="},
            11: {"name": "ibm", "baseUrl": "http://www-935.ibm.com/services/us/gbs/bao/siip/nih/", "entryUrl": "http://www-935.ibm.com/services/us/gbs/bao/siip/nih/?sid="},
            6: {"name": "kegg_ligand", "baseUrl": "http://www.genome.jp/kegg/ligand.html", "entryUrl": "http://www.genome.jp/dbget-bin/www_bget?"},
            9: {"name": "zinc", "baseUrl": "http://zinc15.docking.org", "entryUrl": "http://zinc15.docking.org/substances/"},
            8: {"name": "nih_ncc", "baseUrl": "http://nihsmr.evotec.com/evotec/", "entryUrl": ""},
            10: {"name": "emolecules", "baseUrl": "https://www.emolecules.com/", "entryUrl": "https://www.emolecules.com/cgi-bin/more?vid="},
            12: {"name": "atlas", "baseUrl": "http://www.ebi.ac.uk/gxa/home", "entryUrl": "http://www.ebi.ac.uk/gxa/query?conditionQuery="},
            7: {"name": "chebi", "baseUrl": "http://www.ebi.ac.uk/chebi/downloadsForward.do", "entryUrl": "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI%3A"},
            14: {
                "name": "fdasrs",
                "baseUrl": "http://fdasis.nlm.nih.gov/srs/srs.jsp",
                "entryUrl": "http://fdasis.nlm.nih.gov/srs/ProxyServlet?mergeData=true&objectHandle=DBMaint&APPLICATION_NAME=fdasrs&actionHandle=default&nextPage=jsp/srs/ResultScreen.jsp&TXTSUPERLISTID=",
            },
            15: {"name": "surechembl", "baseUrl": "https://www.surechembl.org/search/", "entryUrl": "https://www.surechembl.org/chemical/"},
            21: {"name": "pubchem_tpharma", "baseUrl": "http://www.thomson-pharma.com/", "entryUrl": "http://pubchem.ncbi.nlm.nih.gov/substance/"},
            22: {"name": "pubchem", "baseUrl": "http://pubchem.ncbi.nlm.nih.gov", "entryUrl": "http://pubchem.ncbi.nlm.nih.gov/compound/"},
            27: {"name": "recon", "baseUrl": "https://vmh.uni.lu", "entryUrl": "https://vmh.uni.lu/"},
            28: {"name": "molport", "baseUrl": "https://www.molport.com/shop/index", "entryUrl": "https://www.molport.com/shop/molecule-link/"},
            31: {
                "name": "bindingdb",
                "baseUrl": "https://www.bindingdb.org/bind/index.jsp",
                "entryUrl": "http://www.bindingdb.org/bind/chemsearch/marvin/MolStructure.jsp?monomerid=",
            },
            41: {"name": "swisslipids", "baseUrl": "http://www.swisslipids.org/", "entryUrl": "http://www.swisslipids.org/"},
            29: {"name": "nikkaji", "baseUrl": "http://jglobal.jst.go.jp/en/", "entryUrl": "http://jglobal.jst.go.jp/en/redirect?Nikkaji_No="},
            32: {"name": "comptox", "baseUrl": "https://comptox.epa.gov/dashboard/", "entryUrl": "https://comptox.epa.gov/dashboard/"},
            33: {"name": "lipidmaps", "baseUrl": "http://www.lipidmaps.org", "entryUrl": "http://www.lipidmaps.org/data/LMSDRecord.php?LMID="},
            35: {"name": "carotenoiddb", "baseUrl": "http://carotenoiddb.jp/index.html", "entryUrl": "http://carotenoiddb.jp/Entries/"},
            36: {"name": "metabolights", "baseUrl": "http://www.ebi.ac.uk/metabolights/", "entryUrl": "http://www.ebi.ac.uk/metabolights/"},
            37: {"name": "brenda", "baseUrl": "https://www.brenda-enzymes.org/index.php", "entryUrl": "https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id="},
            17: {"name": "pharmgkb", "baseUrl": "https://www.pharmgkb.org", "entryUrl": "https://www.pharmgkb.org/drug/"},
            18: {"name": "hmdb", "baseUrl": "http://www.hmdb.ca", "entryUrl": "http://www.hmdb.ca/metabolites/"},
            24: {
                "name": "nmrshiftdb2",
                "baseUrl": "http://nmrshiftdb.nmr.uni-koeln.de/portal/media-type/html/user/anon/page/default.psml/js_pane/P-Home",
                "entryUrl": "http://nmrshiftdb.org/molecule/",
            },
            25: {"name": "lincs", "baseUrl": "http://www.lincsproject.org/", "entryUrl": "http://identifiers.org/lincs.smallmolecule/"},
            39: {"name": "chemicalbook", "baseUrl": "https://www.chemicalbook.com", "entryUrl": "https://www.chemicalbook.com/ChemicalProductProperty_EN_"},
            20: {"name": "selleck", "baseUrl": "http://www.selleckchem.com", "entryUrl": "http://www.selleckchem.com/products/"},
            23: {"name": "mcule", "baseUrl": "https://mcule.com", "entryUrl": "https://mcule.com/"},
            26: {"name": "actor", "baseUrl": "https://actor.epa.gov", "entryUrl": "http://actor.epa.gov/actor/chemical.xhtml?casrn="},
            34: {"name": "drugcentral", "baseUrl": "http://drugcentral.org", "entryUrl": "http://drugcentral.org/drugcard/"},
            38: {"name": "rhea", "baseUrl": "http://www.rhea-db.org", "entryUrl": "http://www.rhea-db.org/searchresults?q=CHEBI:"},
        }
        oD = {}
        try:
            for ky in inchiKeyList:
                unc = unichem_client  # pylint: disable=no-member
                # unc.set_format("json")
                uDL = unc.get(ky)
                if uDL:
                    qD = {}
                    for uD in uDL:
                        if "src_id" in uD and int(uD["src_id"]) in mapD:
                            qD[mapD[int(uD["src_id"])]["name"]] = uD["src_compound_id"]
                    if qD:
                        oD[ky] = qD

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return oD
