.. role:: xml(code)
  :language: xml

.. _Sec:tuto:runInAnotherFolder:

*************************************
Running a test case in another folder
*************************************

If one wants to execute a test located outside of the **ECOGEN** folder, it is now possible. To illustrate this tutorial, we will use a shock tube filled with air in 1D with the high pressure on the left side. This test will be in a folder called **test**.

Required files
==============

First, make sure that you still have the main **ECOGEN** folder, or this will not work. In the folder where you want to work, you need to have the same structure as the one for **ECOGEN**. So, for this example, you need to have:

- *ECOGEN.xml* with the test case associated.
- *libTests* containing the test cases.
- *libEOS* containing the fluids you want to use.

.. note:: 

  In the presented test case, the type of mesh used is Cartesian. However, if you want to use an unstructured mesh, you must have the *libMeshes* folder in the *test* folder.

For this example, the *ECOGEN.xml* file looks like:

.. code-block:: xml

  <?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
  <ecogen>
    <!-- <testCase>./libTests/euler/1D/shockTubes/HPLeft/</testCase> -->
  </ecogen>

As the fluid used here is air, the *libEOS* folder will contain *IG_air.xml*, and the *libTests* folder will contain the 4 input files that are detailed in :ref:`Chap:input`.

To sum up, you should have the following folder/file structure:

.. code-block:: console

  ├── ECOGEN
  │   ├── AUTHORS
  │   ├── COPYRIGHT
  │   ├── coverage_and_profile
  │   ├── docs
  │   ├── ECOGEN.xml
  │   ├── libEOS
  │   ├── libMeshes
  │   ├── libTests
  │   ├── LICENSE
  │   ├── Makefile
  │   ├── nonreg
  │   ├── README_Developer.md
  │   ├── README.md
  │   ├── scripts
  │   └── src
  └── test
      ├── ECOGEN.xml
      ├── libEOS
      └── libTests

Running a test in another folder
================================

As said before, the test used here is a shock tube filled with air, that is precisely detailed in :ref:`Sec:testcases:shocktube`.

.. _Fig:tutos:default:CI_2:

.. figure:: ./_static/tutos/default/CI.*
  :scale: 70%
  :align: center

  Initial condition for 1D, single-phase transport test case.

This test can be executed on a single core or on XX cores by the command:

.. code-block:: console

  mpirun -np XX ECOGEN path/folder/

So, if the folder named *test* is located in *Documents*, and one wants to run the test with 2 cores, the command is:

.. code-block:: console

  mpirun -np 2 ECOGEN ~/Documents/test/

A new folder *results* is created in the *test* folder, with another folder inside named *euler1DTShockTubeHPLeft* and where one can find output files of the test done.