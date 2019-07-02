.. role:: xml(code)
	:language: xml

Input Files
===========

The standard XML format is used for ECOGEN input files. This data format is using markups that gives to the input files some interesting flexibility. ECOGEN is using the `TinyXML-2`_  parser to read/write xml files.

More information about XML format can be found at: https://www.w3schools.com/xml/

ECOGEN package includes a sample of test cases. Each of them is independent, ready to use and can also be adapted and enriched to develop a new configuration. The test case to be run is chosen in the main input file: *ECOGEN.xml*. The minima structure of this file is:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<ecogen>
  	  <testCase>./libTests/referenceTestCases/euler/1D/transport/positiveVelocity/</testCase>
	</ecogen>

In this file, the :xml:`<testCase>` markup indicates the folder containing the test case to be run. It is then possible to run successively several cases by adding as many :xml:`<testCase>` markup as necessary.
Each folder indicated in a :xml:`<testCase>` markup must contain 4 input files:

- *mainV5.xml*
- *meshV5.xml*
- *modelV4.xml*
- *initialConditionsV4.xml*

Additional input files depending on the test case are necessary. They are placed in the following folders:

- **ECOGEN/libEOS/**: contains the files defining the material used in the test case.
- **ECOGEN/libMeshes/**: contains the files defining the mesh used in the test case.

In this section the structure of each input file is detailed. This information is also provided in a condensed form in the "handbook" files at the root folder of test cases library **ECOGEN/libTests/**. These files are named and constitute quick-reference manuals:

- *manualMainV5.xml*
- *manualMeshV5.xml*
- *manualModelV4.xml*
- *manualInitialConditionsV4.xml*


.. _`TinyXML-2`: http://www.grinninglizard.com/tinyxml2/

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Main
   Model
   Mesh
   InitialConditions


