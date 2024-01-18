import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from posecheck.utils.strain import calculate_strain_energy
from posecheck.utils.constants import EXAMPLE_LIGAND_PATH, EXAMPLE_PDB_PATH
from posecheck.utils.loading import load_mols_from_sdf, load_protein_from_pdb




class StrainEnergyCalculationTest(unittest.TestCase):
    def test_get_strain_energy(self):
        # Example molecules
        #prot = load_protein_from_pdb(EXAMPLE_PDB_PATH)
        lig = load_mols_from_sdf(EXAMPLE_LIGAND_PATH)[0]

        # Calculate strain
        strain = calculate_strain_energy(lig)

        self.assertIsInstance(strain, float)


if __name__ == "__main__":
    unittest.main()
