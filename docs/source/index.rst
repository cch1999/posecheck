.. MolParse documentation master file, created by
   sphinx-quickstart on Wed Mar 13 13:02:45 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======================
PoseCheck Documentation
======================

**Warning: this documentation is WIP**

.. image:: _static/posecheck_logo.png
   :width: 100px
   :align: center


PoseCheck is a python package for checking the quality of 3D molecules poses from molecular generative models.

Installation
============

Install from PyPI

::

   git clone https://github.com/cch1999/posecheck.git
   cd posecheck

   pip install -e .
   pip install -r requirements.txt
   conda install -c mx reduce

Getting Started
===============

We provide a simple top level API to easily interact with the whole of the benchmark. Just define the `PoseCheck` object once at the top of your existing testing code and test molecules by loading them in iteratively. You can also use all the testing functions manually as well (see Docs for more info).

::
   from posecheck import PoseCheck

   # Initialize the PoseCheck object
   pc = PoseCheck()

   # Load a protein from a PDB file (will run reduce in the background)
   pc.load_protein_from_pdb("data/examples/1a2g.pdb")

   # Load ligands from an SDF file
   pc.load_ligands_from_sdf("data/examples/1a2g_ligand.sdf")
   # Alternatively, load RDKit molecules directly
   # pc.load_ligands_from_mols([rdmol])

   # Check for clashes
   clashes = pc.calculate_clashes()
   print(f"Number of clashes in example molecule: {clashes[0]}")

   # Check for strain
   strain = pc.calculate_strain_energy()
   print(f"Strain energy of example molecule: {strain[0]}")

   # Check for interactions
   interactions = pc.calculate_interactions()
   print(f"Interactions of example molecule: {interactions}")


Citation
===============

:: 
   @article{harris2023benchmarking,
  title={Benchmarking Generated Poses: How Rational is Structure-based Drug Design with Generative Models?},
  author={Harris, Charles and Didi, Kieran and Jamasb, Arian R and Joshi, Chaitanya K and Mathis, Simon V and Lio, Pietro and Blundell, Tom},
  journal={arXiv preprint arXiv:2308.07413},
  year={2023}
   }


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Home <self>
   Checks <checks>
   User API <api>
   Low-level API <low-level-api>