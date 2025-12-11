##
#  File:           RcsbLigandReferenceGenerator.py
#  Date:           2025-12-10 Chenghua Shao
#
#  Update:
##
"""
Generate ligand quality reference data that is used for ligand quality score computation.

"""
import csv
import logging
import json
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from rcsbapi.data import DataQuery as Query
from rcsbapi.data import ALL_STRUCTURES
from rcsbapi.config import config
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class RcsbLigandReferenceGenerator(StashableBase):
    """ This class generates ligand quality reference data by performing the following steps:
    Query ligand quality metrics from RCSB GraphQL API;
    Process the ligand quality data by filtering, aggregating, and formatting;
    Generate ligand reference data through PCA;
        Attributes:
            data: Stateful data updated with each step
        Methods:
            query(pdb_id): Fetch ligand quality metrics for given PDB IDs.
    """
    def __init__(self):
        self.data = None
        # Increase some settings to speed up querying
        config.DATA_API_MAX_CONCURRENT_REQUESTS = 15
        config.DATA_API_REQUESTS_PER_SECOND = 30
    
    def generate(self, pdb_ids: list[str] = []) -> "RcsbLigandReferenceGenerator":
        """
        Full pipeline to generate ligand quality reference data by running the steps of
        query -> filter -> reduce -> analyze.

        :param pdb_ids: List of specified PDB IDs, default to [] which leads to all PDB structures being queried.
        :return: The same object with updated self.data of ligand quality reference data.
        """
        self.query(pdb_ids).filter().reduce().analyze()
        return self

    def writeReference(self, output_file: str) -> bool:
        """
        Write the generated ligand quality reference data to a csv file.

        :param output_file: Path to the output csv file.
        :return: True upon successful write operation.
        """
        if not self.data:
            logger.error("No data to write. Please run generate() first.")
            return False
        if type(self.data) is not list:
            logger.error("Data format incorrect. Expected a list of dictionaries after generate()")
            return False
        if type(self.data[0]) is not dict:
            logger.error("Data format incorrect. Expected a list of dictionaries after generate()")
            return False
        fieldnames = self.data[0].keys()
        with open(output_file, mode="w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()  # write column headers
            writer.writerows(self.data)  # write all rows
        return True

    def query(self, pdb_ids: list[str] = []) -> "RcsbLigandReferenceGenerator":
        """
        Fetch ligand quality metrics for given PDB IDs, defaulting to all structures.

        :param pdb_ids: List of specified PDB IDs, default to [] which leads to all PDB structures being queried.
        :return: The same object with updated self.data of response from the GraphQL API containing ligand quality metrics.
        """
        if not pdb_ids:
            pdb_ids = ALL_STRUCTURES.get_all_ids("entries")
        logger.info("To query ligand quality metrics on %d PDB IDs", len(pdb_ids))
        query = Query(
            input_type="entries",
            input_ids=pdb_ids,
            return_data_list=[
                "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score.average_occupancy",
                "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score.mogul_bonds_RMSZ",
                "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score.mogul_angles_RMSZ",
                "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score.RSR",
                "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score.RSCC",
                "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_validation_score.completeness",
                "nonpolymer_entities.nonpolymer_entity_instances.rcsb_nonpolymer_instance_annotation.comp_id"
            ]
        )
        self.data = query.exec()
        logger.info("Queried ligand quality metrics for %d PDB IDs", len(self.data.get("data", {}).get("entries", [])))
        return self

    def filter(self) -> "RcsbLigandReferenceGenerator":
        """
        filter the queried ligand quality scores, and simplify the json data by using combined key of pdb-ligand.
        filter out entries without nonpolymer ligands;
        filter out entries without any ligand quality scores;
        filter out ligands whose quality scores do not meet thresholds;

        :return: The same object with updated self.data of filtered ligand quality scores with pdb-ligand as keys.
        """
        da = {}
        l_pdb_reviewed = []
        l_pdb_with_ligand = []
        l_pdb_with_ligand_filtered = []
        l_ligand_reviewed = []
        l_ligand_filtered = []
        for entry in self.data.get("data", {}).get("entries", []):
            pdb_id = entry.get("rcsb_id")
            l_pdb_reviewed.append(pdb_id)
            l_entity = entry.get("nonpolymer_entities", [])
            if not l_entity:
                continue
            for entity in l_entity:
                l_instance = entity.get("nonpolymer_entity_instances", [])
                if not l_instance:
                    continue
                l_pdb_with_ligand.append(pdb_id)
                for instance in l_instance:
                    l_annotation = instance.get("rcsb_nonpolymer_instance_annotation", [])
                    if not l_annotation:
                        continue
                    comp_id = l_annotation[0].get("comp_id", None)
                    if not comp_id:
                        continue
                    l_ligand_reviewed.append(comp_id)
                    l_score = instance.get("rcsb_nonpolymer_instance_validation_score", [])
                    if not l_score:
                        continue
                    score = self.filterScore(l_score[0])
                    if not score:
                        continue
                    l_pdb_with_ligand_filtered.append(pdb_id)
                    l_ligand_filtered.append(comp_id)
                    pdb_ligand = f"{pdb_id.upper()}-{comp_id.upper()}"
                    try:
                        da[pdb_ligand].append(score)
                    except KeyError:
                        da[pdb_ligand] = [score]
        logger.info("%s PDB IDs reviewed", len(set(l_pdb_reviewed)))
        logger.info("%s PDB IDs with nonpolymer ligands", len(set(l_pdb_with_ligand)))
        logger.info("%s PDB IDs with nonpolymer ligands after filtering", len(set(l_pdb_with_ligand_filtered)))
        logger.info("%s ligand instances reviewed, of them %s are unique ligand IDs", len(l_ligand_reviewed), len(set(l_ligand_reviewed)))
        logger.info("%s ligand instances kept after filtering, of them %s are unique ligand IDs", len(l_ligand_filtered), len(set(l_ligand_filtered)))
        self.data = da
        logger.info("%s unique pdb-ligand pairs kept after filtering", len(self.data))
        return self

    def filterScore(self, score: dict) -> dict:
        """
        Filter ligand quality scores to keep instances that meet all of the following criteria:
        0.9 <= average occupancy <= 1.0;
        0.9 <= completeness <= 1.0;
        0 < mogul_bonds_RMSZ <= 10.0 of reasonable range, values beyond 10 are exceptional cases;
        0 < mogul_angles_RMSZ <= 10.0 of reasonable range, values beyond 10 are exceptional cases;
        0 < RSR <= 1.0 of physcally meaningful range;
        0 < RSCC <= 1.0 of physically meaningful range;
        value 0 is theoretically possible for the quality metrics, but practically unlikely and sometimes
        indicate placeholder, hence excluded. 
        Instance with any value of None is excluded.
        Exceptional cases are not included in reference building.
        
        :param score: _description_
        :return: _description_
        """        
        for key, value in score.items():
            if key in ["average_occupancy", "completeness"]:  # ligand should be well occupied and complete
                if value is None:
                    return {}
                if value < 0.9 or value > 1:
                    return {}
            if key in ["mogul_bonds_RMSZ", "mogul_angles_RMSZ"]:  # ligand geometry should be reasonable
                if value is None:
                    return {}
                if value <=0 or value > 10:
                    return {}
            if key in ["RSR", "RSCC"]:  # ligand fit to density should be within the bounds
                if value is None:
                    return {}
                if value <=0 or value > 1:
                    return {}
        score_filtered = {}
        try:  # get values for the four primary scores that passed the above filters
            score_filtered["mogul_bonds_RMSZ"] = score["mogul_bonds_RMSZ"]
            score_filtered["mogul_angles_RMSZ"] = score["mogul_angles_RMSZ"]
            score_filtered["RSR"] = score["RSR"]
            score_filtered["RSCC"] = score["RSCC"]
        except KeyError:
            return {}
        return score_filtered

    def reduce(self) -> "RcsbLigandReferenceGenerator":
        """
        Reduce the data by combining multiple instances of the same ligand in the same PDB entry.
        This is to avoid the over-representation of ligands that appear hundreds of times in the same entry.
    
        :return: The same object with updated self.data of reduced ligand quality scores.
        """
        da = {}
        for pdb_ligand, l_score in self.data.items():
            n_instance = len(l_score)
            if n_instance == 1:
                da[pdb_ligand] = l_score[0]
            else:
                score_agg = {}
                for key in l_score[0].keys():
                    l_value = [instance[key] for instance in l_score]
                    array = np.array(l_value)
                    average = np.mean(array)
                    score_agg[key] = average
                da[pdb_ligand] = score_agg
        self.data = da
        logger.info("%s unique pdb-ligand pairs reduced by averaging all instances for each pair", len(self.data))
        return self

    def analyze(self):
        """
        Run PCA on the filtered and reduced ligand quality scores to generate reference data.
        The first principal component (PC1) is used as the ligand quality reference score.
    
        :return: A dictionary of ligand quality reference scores with pdb-ligand as keys.
        """
        # Prepare data for PCA, convert each score type into a separate list for matrix construction
        l_pdb_ligand = []
        l_mogul_bonds_RMSZ = []
        l_mogul_angles_RMSZ = []
        l_RSR = []
        l_RSCC = []
        for pdb_ligand, score in self.data.items():
            l_pdb_ligand.append(pdb_ligand)
            l_mogul_bonds_RMSZ.append(score["mogul_bonds_RMSZ"])
            l_mogul_angles_RMSZ.append(score["mogul_angles_RMSZ"])
            l_RSR.append(score["RSR"])
            l_RSCC.append(score["RSCC"])
        # Construct data matrices for geometry and fit scores
        matrix_geo = np.column_stack([l_mogul_bonds_RMSZ, l_mogul_angles_RMSZ])
        matrix_fit = np.column_stack([l_RSR, l_RSCC])
        # Run PCA separately on geometry and fit scores
        logging.info("To Run PCA on geometry scores")
        l_geo_pc1 = self.runPca(matrix_geo)
        logging.info("To Run PCA on fit scores")
        l_fit_pc1 = self.runPca(matrix_fit)
        # Combine the two PC1 scores into a final reference score
        da = []
        for i, pdb_ligand in enumerate(l_pdb_ligand):
            ref_score = {
                "pdb_ligand": pdb_ligand,
                "rsr": round(l_RSR[i], 4),
                "rscc": round(l_RSCC[i], 4),
                "mogul_bonds_rmsz": round(l_mogul_bonds_RMSZ[i], 4),
                "mogul_angles_rmsz": round(l_mogul_angles_RMSZ[i], 4),
                "fit_pc1": round(l_fit_pc1[i], 5),
                "geo_pc1": round(l_geo_pc1[i], 5)
            }
            da.append(ref_score)
        self.data = da
        logger.info("%s unique pdb-ligand pairs analzed by PCA", len(self.data))
        return self

    def runPca(self, X: np.ndarray) -> list:
        """
        Run PCA on numpy matrix and output the first principal components.
        Loadings and explained variance can be examined in log but not in outcome.
        Ensure lower scores correspond to better quality by adjusting the sign of PC1 if needed.

        :param X: numpy matrix of any dimension
        :return: first principal components as a list
        """        
        # Scale the data
        scaler = StandardScaler()
        X_std = scaler.fit_transform(X)
        # Perform PCA
        pca = PCA(n_components=1)
        principal_components = pca.fit_transform(X_std)
        logger.info("Explained variance ratios: %s", pca.explained_variance_ratio_)  # should be > 0.7(70%) for ligand geo and fit PC1
        loadings = pca.components_[0]
        logger.info("Resulting PCA loadings: %s", loadings)  # should be [sqrt(0.5), sqrt(0.5)] for scaled PCA or with minus sign on one component
        # Code below ensures loading on the first variable (RSR for fitting, and mogul_bonds_RMSZ for geometry) to be positive, 
        # This is because PCA direction is arbitrary, i.e. 0.707,-0.707 and -0.707,0.707 are equivalent, but we want to ensure 
        # lower scores correspond to better quality, which gives a correct direction for subsequent pecentile calculation.
        if loadings[0] < 0:
            loadings = -loadings
            logger.info("Adjusted PCA loadings: %s", loadings)
            principal_components = -principal_components
            logger.info("Change sign of the first principal components to ensure correct directionality.")
        return principal_components[:, 0].tolist()
