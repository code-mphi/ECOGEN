Prerequisites
=============

ECOGEN must be compiled with C++. It also requires a functional system implementation of MPI library (not provided in this package). Depending on your operating system, you can follow the instructions below to set a full open-source installation.

Installing prerequisites on Ubuntu system
-----------------------------------------
ECOGEN requires two mandatory components to be installed on your Ubuntu system: A C++ compiler and an effective implementation of MPI.

Installing C++ compiler
~~~~~~~~~~~~~~~~~~~~~~~
Nothing is easier than installing C and C++ compiler on Ubuntu. In your terminal, just enter the following commands:

.. highlight:: console

::

	sudo apt-get update
	sudo apt-get upgrade
	sudo apt-get install gcc
	sudo apt-get install g++
	sudo apt-get install build-essential
	sudo apt-get install make

More information on the Ubuntu doc page https://doc.ubuntu-fr.org/gcc.

Installing openMPI
~~~~~~~~~~~~~~~~~~

Download the latest stable version of openMPI_ under compressed format. At the time this page was written, it corresponded to the compressed file: openmpi-4.0.1.tar.gz. Uncompresse and move it into the directory:

.. highlight:: console

::

	tar -xvf openmpi-4.0.1.tar.gz 
	cd openmpi-4.0.1/

Prepare the environment to use your favorite compiler:

.. highlight:: console

::

	export CC=gcc 
	export CXX=g++ 

Configure and proceed to the installation (you can choose a different directory). The "make" step should take some time (coffee time?):

.. highlight:: console

::

	./configure --prefix=/opt/openmpi 
	make 
	sudo make install 

Cleaning

.. highlight:: console

::

	cd .. 
	rm -rf openmpi-4.0.1/ 

Add openMPI library to the environment variable PATH (might be required to be root):

.. highlight:: console

:: 

	sudo echo 'export PATH=/opt/openmpi/bin:$PATH' >> /etc/bash.bashrc

Then source the file to take into consideration the modifications:

.. highlight:: console

::

	source /etc/bash.bashrc

If the installation succeeds you should be able to use the "mpicxx" command in your terminal. Then proceed to the download step below.

Download
========

The last ECOGEN version |version| can be downloaded from GitHub. The source files are available at the following address: https://github.com/code-mphi/ECOGEN/releases. 

The package includes:

* ECOGEN/src/ folder including C++ source files.
* ECOGEN/libMeshes/ folder including examples of unstructured meshes in *.geo* format (gmsh files version 2). See section :ref:`Sec:tuto:generatingMeshes` for details.
* ECOGEN/libEOS/ folder including some possible parameters for the equation-of-states in XML files. See section :ref:`Sec:IO:materials` for details.
* ECOGEN/libTests folder including:

	- ECOGEN/libTests/referenceTestCases/ folder organized as a test-case library according to the flow model (Euler-equation ECOGEN solver, pressure-velocity-equilibrium model (previously named Kapila's model) and velocity-equilibrium model for multiphase-flow ECOGEN solver, homogeneous-Euler-equation ECOGEN solver, etc.). A detailed list of available test cases is proposed in section :ref:`Chap:TestCases`.
	- 4 quick-manual XML files to create a new flow computation with ECOGEN.
* *ECOGEN.xml*: Main entry file to select running cases.
* *Makefile*: For compilation in Unix environment. This file may require some adaptation to the user's environment.
* *LICENSE*, *COPYRIGHT* and *AUTHORS*: Information files about authors and licensing.
* *README.md*: Information file.
* *ECOGEN_V1.0_documentation.pdf*: The full documentation for ECOGEN.

.. _Sec:installation:compileAndExecute:

Compilation/Execution on bash
=============================

Use the Makefile (can be adapted if necessary) to compile ECOGEN sources directly on bash (XX is the number of cores required for compilation):

.. highlight:: console

::

	make -j XX

Executing ECOGEN is really easy on bash (XX is the number of cores required for execution):

.. highlight:: console

::

	mpirun -np XX ECOGEN

Testing
=======

Once ECOGEN has been successfully compiled, the best way to test ECOGEN's installation is to run successively the two simple following commands:

.. highlight:: console

::

	./ECOGEN
	mpirun -np 2 ECOGEN

These will run the default test case included in the package two times:

* Once in sequential (single core). 
* Once in parallel using 2 cores.

These should print information in the terminal on the running default test case. If no error message appears, then your installation should be OK and you should be able to use ECOGEN for your own applications.

ECOGEN is including a given number of simple prebuild test cases. Each test can be used as a basis for a new one. Visit the tutorial section :ref:`Chap:Tutorials` for more information.

.. _openMPI: https://www.open-mpi.org/