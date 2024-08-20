import os
import subprocess
import unittest

from posecheck.utils.constants import EXAMPLE_PDB_PATH, REDUCE_PATH


class Testhydride(unittest.TestCase):
    def test_hydride_exists(self):
        command = 'hydride'
        exit_response = subprocess.run(
            command, capture_output=True, text=True, shell=True
        )
        self.assertEqual(exit_response.returncode, 1)

