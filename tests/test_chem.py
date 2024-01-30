import unittest

import datamol as dm
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import numpy as np

from posecheck.utils.chem import rmsd, has_radicals, remove_radicals
from posecheck.utils.constants import EXAMPLE_LIGAND_PATH, EXAMPLE_PDB_PATH


class TestChem(unittest.TestCase):
    def test_rmsd(self):
        mol1 = dm.read_sdf(EXAMPLE_LIGAND_PATH)[0]
        mol2 = dm.read_sdf(EXAMPLE_LIGAND_PATH)[0]
        result = rmsd(mol1, mol2)
        self.assertEqual(result, 0.0)

    def test_remove_radicals(self):
        """Test the remove_radicals function."""
        # Create a molecule with free radicals
        mol = Chem.MolFromSmiles("[O]")

        # Check if the molecule has radicals before removal
        self.assertTrue(has_radicals(mol))

        # Remove the radicals
        mol = remove_radicals(mol)

        # Check if the molecule has radicals after removal
        self.assertFalse(has_radicals(mol))

    def test_remove_radicals2(self):
        # Define a test molecule

        # smile for large molecule
        smiles = "CC(CC1=CC=C(C=C1)C(=O)O)C1=CC=C(C=C1)C(=O)O"

        mol = Chem.MolFromSmiles(smiles)
        # add generate conformer using UFF
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)

        # set radical
        for i, atom in enumerate(mol.GetAtoms()):
            if atom.GetSymbol() == "O":
                if atom.GetSymbol() == "C":
                    atom.SetNumRadicalElectrons(1)
                    break

        mol = remove_radicals(mol)


    def test_has_radicals(self):
        """Test the has_radicals function."""
        # Create a molecule without radicals
        mol = Chem.MolFromSmiles("CC")

        # Check if the molecule has radicals
        self.assertFalse(has_radicals(mol))

        # Create a molecule with radicals
        mol = Chem.MolFromSmiles("[O]")

        # Check if the molecule has radicals
        self.assertTrue(has_radicals(mol))


if __name__ == "__main__":
    unittest.main()
