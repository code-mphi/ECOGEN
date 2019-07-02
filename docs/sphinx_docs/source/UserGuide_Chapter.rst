.. _Chap:UserGuide:

User's Guide
============

This user's guide is for the version 1.0 of ECOGEN package: ECOGEN_V1.0.zip. The package is available on the permanent link: https://github.com/code-mphi/ECOGEN.

The package includes several files or folders organized as described below:

- ECOGEN/src/ folder including 254 C++ source files in 23 subfolders.
- ECOGEN/libMeshes/ folder including examples of unstructured meshes in *.geo* format (gmsh files version 2).
- ECOGEN/libEOS/ folder including some possible parameters for Ideal Gas or Stiffened Gas Equation of State in XML files.
- ECOGEN/libTests folder including:
	- ECOGEN/libTests/referenceTestCases/ folder organized in a test cases library according the flow model (Euler Equations ECOGEN solver, Kapila's model for multiphase flow ECOGEN solver, Homogeneous Euler Equation ECOGEN solver, etc.)
	- 4 quick-manual XML files to create a new flow calculation with ECOGEN.
	- 1 reporting file for test cases.
- ECOGEN.xml main entry file to select running cases.
- Makefile: for compilation in Unix environment. This file may require some adaptation to the user's environment.
- LICENSE, COPYRIGHT and AUTHORS: Information files about authors and licensing.
- README_Developer: Information files for developers.
- ECOGEN_V1.0_userGuide.pdf: The present user's guide for ECOGEN. Include full description of input and output files to proceed to a run using ECOGEN.

ECOGEN should be compiled using a C++ compiler. It also required a functional system implementation of MPI library (not provided in this package). For prerequisities and installation, see the corresponding chapter :ref:`Chap:Start`

ECOGEN settings are managed via INPUT FILES only. The global INPUT FILES and OUTPUT FILES structure is depicted in :numref:`Fig:userGuide:overview`.

.. _Fig:userGuide:overview:

.. figure:: ./img/overview.png

	Structure of input and output files in ECOGEN

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   InputFiles_Chapter
   Materials
   OutputFiles_Chapter
   Tutorials_Chapter