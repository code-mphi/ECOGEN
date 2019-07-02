Prerequisities
==============

ECOGEN must be compiled with C++. It also requires a functional system implementation of MPI library (not provided in this package). Depending on your operating system, you can follow the instructions below to set a full open source installation:

Installing prerequisities on Ubuntu system
------------------------------------------
ECOGEN required two mandatory components to be installed on your Ubuntu system : a C++ compiler and an effective implementation of MPI.

Installing C++ compiler
~~~~~~~~~~~~~~~~~~~~~~~
Nothing is more easy than installing C and C++ compiler on Ubuntu. In your temrminal just enter the following commands:

.. highlight:: console

::

	sudo apt-get update
	sudo apt-get upgrade
	sudo apt-get install gcc
	sudo apt-get install g++
	sudo apt-get install build-essential

More information on the ubuntu doc page https://doc.ubuntu-fr.org/gcc

Installing openMPI
~~~~~~~~~~~~~~~~~~

Dowload the latest stable version of openMPI_ under compressed format. At the time this page is written, it corresponds to the compressed file : openmpi-4.0.1.tar.gz. Uncompresse and move into the directory:

.. highlight:: console

::

	tar -xvf openmpi-4.0.1.tar.gz 
	cd openmpi-4.0.1/

Prepare the environnement for using your favorite compiler:

.. highlight:: console

::

	export CC=gcc 
	export CXX=g++ 

Configure and proceed to the installation (you can choose a different directory):
.. highlight:: console

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

Modify the /etc/bash.bashrc by adding the line:

*export PATH=/opt/openmpi/bin:$PATH*

Then source the file to take into consideration the modifications:

.. highlight:: console

::

	source /etc/bash.bashrc

If the installation succeed you should use the mpicxx command in your terminal. Then proceed to the download step below.

Download
========

The last ECOGEN version can be downloaded from Git-hub. The source files are available at the following address: https://github.com/code-mphi/ECOGEN

Compilation/Execution on bash
=============================

Use the Makefile (can be adapted if necessary) to compile ECOGEN sources directly on bash (XX is the number of CPU required for compilation):

.. highlight:: console

::

	make -j XX

Executing ECOGEN is really easy on bash (XX is the number of CPU required for execution):

.. highlight:: console

::

	mpirun -np XX ECOGEN

ECOGEN is including a given number of simple prebuild test cases. Each test can be used as a basis for a new one. In the downloaded version, when executing ECOGEN, a simple default test case is running and you should see informations about the test output appering in the terminal output.


.. _openMPI: https://www.open-mpi.org/