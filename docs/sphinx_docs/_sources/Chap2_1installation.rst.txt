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

The simpliest way to install openMPI_ on Ubuntu is with the apt package manager:

.. highlight:: console

::

  sudo apt-get install libopenmpi-dev

Download
========

The last ECOGEN version |version| can be downloaded from GitHub. The source files are available at the following address: https://github.com/code-mphi/ECOGEN/releases.

The package includes:

* ECOGEN/src/ folder including C++ source files.
* ECOGEN/libMeshes/ folder including examples of unstructured meshes in *.geo* format (Gmsh files version 2). See section :ref:`Sec:tuto:generatingMeshes` for details.
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

Compilation on bash using the provided Makefile
===============================================

Compilation
-----------

Use the Makefile (can be adapted if necessary) to compile ECOGEN sources directly on bash (XX is the number of cores required for compilation):

.. highlight:: console

::

	make -j XX

It builds the ECOGEN application at root of the worktree. By default ECOGEN is built in `release` mode (`-O3 -DNDEBUG` compiler flags).

Predefined compilation modes
----------------------------

Other compilation modes are proposed in predefined Makefile targets:

   - `debug` mode for `-O0 -g` flags,
   - `coverage mode` for `-O3 --coverage` flags,
   - `profile mode` for `-O3 -pg` flags.

For example, to build ECOGEN with debug informations:

.. highlight:: console

::

   make -j XX debug


Adding compilation flags
------------------------

If needed, compilation flags and linker flags can be overloaded using, respectively, `CXXFLAGS` and `LDFLAGS` variables. For example, to build ECOGEN with `AddressSanitizer` (a memory error detector):

.. highlight:: console

::

    CXXFLAGS=-fsanitize=address LDFLAGS=-fsanitize=address make -j XX

Compiling ECOGEN library
------------------------

If one needs to use ECOGEN as a library, the provided Makefile also defines the `libECOGEN.so` target that produces a shared library of ECOGEN. To build it:

.. highlight:: console

::

    make -j XX libECOGEN.so

By default, the ECOGEN application built by calling `make` without providing an explicit target doesn't use the ecogen library. The ECOGEN application that calls the *libECOGEN.so* library can be built by calling the `exe_using_lib` target:

.. highlight:: console

::

    make -j XX exe_using_lib


Clean build files
-----------------

The `clean` target (`make clean` command) cleans temporary build files. The `distclean` target (`make distclean`) cleans all files created by the makefile (including *ECOGEN* and *libECOGEN.so* files).



Compilation using CMake tool
============================

Alternatively, more advanced users can choose to use `CMake <https://cmake.org/>`_ to build the project.

.. warning::

  To avoid any confusion, do not mix compilation using the provided Makefile and compilation using CMake in the same worktree.

Prerequisite
------------

Install CMake using your usual package manager. For example on Ubuntu:

.. highlight:: console

::

   apt install -y cmake

Compilation
-----------

1. Configuration step

   Generate the build files in the *build_dir* directory:

  .. highlight:: console

  ::

     cmake -S . -B build_dir

2. .. _cmake-compilation-step:

   Compilation step


   Build the project (XX is the number of cores required for compilation).

  .. highlight:: console

  ::

     cmake --build build_dir -j XX

The ECOGEN shared library (*libECOGEN.so*) and the ECOGEN executable (*ECOGEN*) are created inside the *build_dir* directory.

A static library can be built instead of the shared one by setting the *BUILD_SHARED_LIBS* CMake variable to *OFF* (`cmake -DBUILD_SHARED_LIBS=OFF ...`).

By default, the ECOGEN application doesn't use the ecogen library. The `ECOGEN_EXE_USE_LIB` cmake variable can be set to ON to build an ECOGEN application that calls the ecogen library:

.. highlight:: console

::

  cmake -S . -B build_dir -D ECOGEN_EXE_USE_LIB=ON
  cmake --build build_dir -j XX

Predefined compilation modes
----------------------------

CMake proposes predefined compilation modes through the setting of the `CMAKE_BUILD_TYPE` variable. Please refer to the `official documentation <https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html>`_ for more details. Typical values are `Debug` (`-g` compiler flag) at development step and `Release` (`-O3 -DNDEBUG` compiler flags) at production step:

.. highlight:: console

::

  cmake -S . -B build_dir -D CMAKE_BUILD_TYPE=Release

.. _enabling-hdf5:

Enabling HDF5 I/Os
==================

If needed, user can configure and build ECOGEN with `HDF5 I/Os <https://cmake.org/>`_ capabilities using the CMake build mode previously described.

Prerequisite
------------

Install **parallel HDF5**, either from source (see `HDF5 documentation <https://support.hdfgroup.org/documentation/hdf5/latest/md_release__docs_2_r_e_a_d_m_e___h_p_c.html#autotoc_md274>`_), or using a package manager.

For example, to install parallel HDF5 on *Ubuntu* with *OpenMPI* MPI library:

.. highlight:: console

::

  apt install -y libhdf5-openmpi-dev

CMake configuration
-------------------

If parallel HDF5 is installed on your system, CMake should detect it and enable HDF5 I/Os automatically at configuration step.
You can force CMake trigger an error if HDF5 is missing by enabling explicitly HDF5 support:

.. highlight:: console

::

  cmake -S . -B build_dir -DENABLE_HDF5=ON

You can then build your project as usually (see :ref:`cmake compilation step<cmake-compilation-step>`)

Executing ECOGEN
================

Executing ECOGEN is really easy on bash (XX is the number of cores required for execution):

.. highlight:: console

::

    mpirun -np XX <path-to-ecogen-binary>/ECOGEN

So, if ECOGEN has been built using provided makefile:

.. highlight:: console

::

    mpirun -np XX ./ECOGEN

If it has been built in `build_dir` directory using CMake:

.. highlight:: console

::

  mpirun -np XX ./build_dir/ECOGEN


ECOGEN may used colorized outputs using ANSI color escape sequences. To avoid garbage outputs in terminals
that doesn't support ANSI colors, you can set the NO_COLOR environment variable to 1 to disable output colorization:

.. highlight:: console

::

    NO_COLOR=1 mpirun -np XX <path-to-ecogen-binary>/ECOGEN

Testing
=======

Once ECOGEN has been successfully compiled, the best way to test ECOGEN's installation is to run successively the two simple following commands:

* if you built the project using the provided Makefile:

  .. highlight:: console

  ::

    ./ECOGEN
    mpirun -np 2 ECOGEN

* if you built it using CMake:

  .. highlight:: console

  ::

    ./build_dir/ECOGEN
    mpirun -np 2 ./build_dir/ECOGEN

These will run the default test case included in the package two times:

* Once in sequential (single core).
* Once in parallel using 2 cores.

These should print information in the terminal on the running default test case. If no error message appears, then your installation should be OK and you should be able to use ECOGEN for your own applications.

ECOGEN is including a given number of simple prebuild test cases. Each test can be used as a basis for a new one. Visit the tutorial section :ref:`Chap:Tutorials` for more information.

.. _openMPI: https://www.open-mpi.org/
