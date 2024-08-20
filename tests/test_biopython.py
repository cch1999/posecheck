import unittest
from posecheck.utils.biopython import load_biopython_structure, save_biopython_structure, ids_scriptly_increasing, reorder_ids
from Bio.PDB.Structure import Structure
import tempfile
import os

class TestBioPython(unittest.TestCase):
    def setUp(self):
        # Setup temporary PDB file
        self.temp_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        self.temp_pdb.write(b"ATOM      1  N   MET A   1      20.154  29.699  52.379  1.00 60.69           N\n")
        self.temp_pdb.close()

    def tearDown(self):
        # Cleanup
        os.unlink(self.temp_pdb.name)

    def test_load_biopython_structure(self):
        structure = load_biopython_structure(self.temp_pdb.name)
        self.assertIsInstance(structure, Structure)

    def test_save_biopython_structure(self):
        structure = load_biopython_structure(self.temp_pdb.name)
        save_biopython_structure(structure, self.temp_pdb.name)
        self.assertTrue(os.path.exists(self.temp_pdb.name))

    def test_ids_scriptly_increasing(self):
        structure = load_biopython_structure(self.temp_pdb.name)
        self.assertTrue(ids_scriptly_increasing(structure))

    def test_reorder_ids(self):
        structure = load_biopython_structure(self.temp_pdb.name)
        reordered_structure = reorder_ids(structure)
        self.assertTrue(ids_scriptly_increasing(reordered_structure))

if __name__ == '__main__':
    unittest.main()
