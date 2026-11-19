##
#  File:           testRcsbLigandReferenceGenerator.py
#  Date:           2025-12-10 Chenghua Shao
#
#  Update:
##
"""
Unit tests for RcsbLigandReferenceGenerator.py
"""
import json
import logging
import os
import platform
import resource
import time
import unittest
from rcsb.utils.chemref.RcsbLigandReferenceGenerator import RcsbLigandReferenceGenerator

HERE = os.path.dirname(os.path.abspath(__file__))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RcsbLigandReferenceGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        self.cRLRG = RcsbLigandReferenceGenerator()

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 10 ** 6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testQuery(self):
        """
        Test query on several PDB IDs.
        """
        pdb_ids = ["1C0T", "1DT4", "XXXX"]  # 1C0T with ligand, 1DT4 without ligand, XXXX invalid ID
        self.cRLRG.query(pdb_ids)
        self.assertTrue(self.cRLRG.data)
        response = self.cRLRG.data
        logger.info("query on %s with response %s", pdb_ids, response)
        self.assertIn("data", response)
        self.assertIn("entries", response["data"])
        l_pdb_id = []
        for entry in response["data"]["entries"]:
            pdb_id = entry["rcsb_id"]
            l_pdb_id.append(pdb_id)
            if pdb_id == "1DT4":
                self.assertIsNone(entry["nonpolymer_entities"])
                logger.info("No ligand entities for PDB ID %s as expected.", pdb_id)
            elif pdb_id == "1C0T":
                self.assertIsNotNone(entry["nonpolymer_entities"])
                logger.info("Ligand entities found for PDB ID %s as expected.", pdb_id)
        logger.info("PDB IDs with response: %s", l_pdb_id)
        self.assertNotIn("XXXX", l_pdb_id)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testQuery.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(response, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testQueryAll(self):
        """
        Test query on all PDB structures.
        """
        self.cRLRG.query()
        self.assertTrue(self.cRLRG.data)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testQueryAll.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(self.cRLRG.data, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testFilter(self):
        """
        Test filter function on several PDB IDs.
        """
        pdb_ids = ["1C0T", "1DT4", "6WJC", "4HHB"]
        self.cRLRG.query(pdb_ids)
        self.cRLRG.filter()
        data_filtered = self.cRLRG.data
        self.assertTrue(data_filtered)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testFilter.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(data_filtered, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testFilterAll(self):
        """
        Test filter function on all PDB structures.
        """
        self.cRLRG.query()
        self.cRLRG.filter()
        data_filtered = self.cRLRG.data
        self.assertTrue(data_filtered)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testFilterAll.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(data_filtered, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testReduce(self):
        """
        Test reduce function on several PDB IDs.
        """
        pdb_ids = ["1C0T", "1DT4", "6WJC", "4HHB"]
        self.cRLRG.query(pdb_ids)
        self.cRLRG.filter()
        self.cRLRG.reduce()
        data_reduced = self.cRLRG.data
        self.assertTrue(data_reduced)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testReduce.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(data_reduced, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testReduceAll(self):
        """
        Test reduce function on all PDB structures.
        """
        self.cRLRG.query()
        self.cRLRG.filter()
        self.cRLRG.reduce()
        data_reduced = self.cRLRG.data
        self.assertTrue(data_reduced)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testReduceAll.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(data_reduced, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testAnalyze(self):
        """
        Test analyze function on several PDB IDs.
        """
        pdb_ids = ["1C0T", "1DT4", "6WJC", "4HHB"]
        self.cRLRG.query(pdb_ids)
        self.cRLRG.filter()
        self.cRLRG.reduce()
        self.cRLRG.analyze()
        data_analyzed = self.cRLRG.data
        self.assertTrue(data_analyzed)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testAnalyze.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(data_analyzed, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testAnalyzeAll(self):
        """
        Test analyze function on all PDB structures.
        """
        self.cRLRG.query()
        self.cRLRG.filter()
        self.cRLRG.reduce()
        self.cRLRG.analyze()
        data_analyzed = self.cRLRG.data
        self.assertTrue(data_analyzed)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testAnalyzeAll.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(data_analyzed, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)

    def testGenerate(self):
        """
        Test generate function that runs the full pipeline on several PDB IDs.
        """
        pdb_ids = ["1C0T", "1DT4", "6WJC", "4HHB"]
        self.cRLRG.generate(pdb_ids)
        self.assertTrue(self.cRLRG.data)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testGenerate.json")
        with open(output_file, "w", encoding="utf-8") as file:
            json.dump(self.cRLRG.data, file, indent=2)
        logger.info("Wrote data to output file %s", output_file)
        # Write reference data to csv
        csv_output_file = os.path.join(HERE, "test-output", "ligand_score_reference_test.csv")
        self.assertTrue(self.cRLRG.writeReference(csv_output_file))
        logger.info("Wrote reference data to csv file %s", csv_output_file)

    def testGenerateAll(self):
        """
        Test generate function that runs the full pipeline on all PDB structures.
        """
        self.cRLRG.generate()
        self.assertTrue(self.cRLRG.data)
        # Write to output file
        output_file = os.path.join(HERE, "test-output", "RcsbLigandReferenceGenerator_testGenerate.json")
        with open(output_file, "w", encoding="utf-8") as f:
            json.dump(self.cRLRG.data, f, indent=2)
        logger.info("Wrote data to output file %s", output_file)
        # Write reference data to csv
        csv_output_file = os.path.join(HERE, "test-output", "ligand_score_reference.csv")
        self.assertTrue(self.cRLRG.writeReference(csv_output_file))
        logger.info("Wrote reference data to csv file %s", csv_output_file)


def genLigandRef():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testQuery"))
    # suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testQueryAll"))
    suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testFilter"))
    # suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testFilterAll"))
    suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testReduce"))
    # suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testReduceAll"))
    suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testAnalyze"))
    # suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testAnalyzeAll"))
    suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testGenerate"))
    # suiteSelect.addTest(RcsbLigandReferenceGeneratorTests("testGenerateAll"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = genLigandRef()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
