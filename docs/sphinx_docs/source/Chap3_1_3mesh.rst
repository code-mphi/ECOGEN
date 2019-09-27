.. role:: xml(code)
	:language: xml

MeshV5.xml
==========

Input file *meshV5.xml* is necessary to specify the geometrical characteristics of the computational domain and the kind of mesh used. This file is **mandatory** and must be present in the folder of the current case. The minimalist content of this file is:

.. code-block:: xml

	<?xml version="1.0" encoding="UTF-8" standalone = "yes"?>
	<mesh>
	  <type structure="cartesian"/>
	  <cartesianMesh>
	    <dimensions x="1.e-1" y="5.e-2" z="1."/>
	    <numberCells x="50" y ="25" z="1"/>
	  </cartesianMesh>
	</mesh>
 
The  :xml:`<type>` markup is mandatory and specifies via the attribute :xml:`structure` the kind of mesh used in the current computation:

- *cartesian*: ECOGEN automatically generates its own Cartesian mesh. See section :ref:`Sec:input:cartesian` for details.
- *unstructured*: a specific external mesh generator (not provide with ECOGEN) must be used to generate the mesh.  See section :ref:`Sec:input:unstructured` for details.

.. _Sec:input:cartesian:

Cartesian mesh
--------------

.. code-block:: xml

  <cartesianMesh>
    <dimensions x="1.e-1" y="5.e-2" z="1."/>
    <numberCells x="50" y ="25" z="1"/>
  </cartesianMesh>

ECOGEN is able to automatically generate a cartesian mesh according the markup :xml:`<cartesianMesh>` which must contain the two following nodes:

- :xml:`<dimensions>`: Specifies the physical dimensions of the computational domain in each physical direction corresponding to the attributes: :xml:`x`, :xml:`y` and :xml:`z` (unit: m (SI)). Values must be real numbers.
- :xml:`<numberCells>`: Specify the number of cells in each direction corresponding to the :xml:`x`, :xml:`y` and :xml:`z` attributes. The values are integer numbers.

Optional Stretching
~~~~~~~~~~~~~~~~~~~

.. code-block:: xml

	<meshStretching>    <!-- Optionnal node -->
	  <XStretching>
	    <stretch startAt="0." endAt="0.5" factor="0.9" numberCells="20"/>
	    <stretch startAt="0.5" endAt="1." factor="1.1" numberCells="10"/>
	  </XStretching>
	</meshStretching>

Stretching can be set optionally adding the :xml:`<XStretching>` node to the :xml:`<cartesianMesh>` parent markup for X stretching. It should contain one or more :xml:`<stretch>` node(s) equipped with the following attributes:

- :xml:`startAt`: real number giving the beginning position of the stretched zone.
- :xml:`endAt`: real number giving the ending position of the stretched zone.
- :xml:`factor`: real number for the stretch factor (lower than 1 for shrinking, greater than 1 for stretching).
- :xml:`numberCells`: integer for the cell number attributed to the stretched zone.

**Remark:** 

1. Stretching can be set in each directions using :xml:`<YStretching>` and :xml:`<ZStretching>`. 
2. For each stretched direction, the sum of stretched zones should exactly recover the entire domain without overlaping, but the number of cell can differ than those precised in the :xml:`<numberCells>` initial markup.
3. A particular attention should be paid to the linking between zones that possibily present a bad quality.

Optional AMR
~~~~~~~~~~~~

.. code-block:: xml

	<AMR lvlMax="2" criteriaVar="0.2" varRho="true" varP="true" varU="false" varAlpha="false" xiSplit="0.11" xiJoin="0.11"/> <!-- Optionnal node -->

An efficient Adaptive Mesh refinement (AMR) technology  is embedded in ECOGEN. To do that *meshV5.xml* file must content the optional node :xml:`<AMR>` of the :xml:`<cartesianMesh>` markup to define the following attributes:

- :xml:`lvlMax`: integer to define the maximal number of refinements. 
- :xml:`criteriaVar`: real number that controls the detection of gradients for the location of the refinement.
- :xml:`varRho`, :xml:`varP`, :xml:`varU`, :xml:`varAlpha`: boolean (*true* or *false*) to select the flow quantity on which the gradient operator is applied to detect large gradient.
- :xml:`xiSplit`, :xml:`xiJoin`: normalized real numbers to control if a computational cell, selected by its high gradient value, must be refined or un-refined (values are in the range *0-1*.). 

**Remark:**

The global efficiency of the method is greatly depending on the chosen values for the :xml:`criteriaVar`, :xml:`xiSplit` and :xml:`xiJoin` attributes. These values depend on the physical problem and required a real *know-how*. More details about these criterion values can be found in :cite:`schmidmayer2019adaptive`.

.. _Sec:input:unstructured:

Unstructured mesh
-----------------

.. code-block:: xml

	<unstructuredMesh>
	  <file name="unstructured2D/testUS.msh"/>
	  <parallel GMSHPretraitement="true"/>  <!-- Optionnal node if multiCPU -->
	</unstructuredMesh>

When dealing with unstructured meshes, the :xml:`<unstructuredMesh>` markup **must be** present in the *meshV5.xml* input file and contains the following nodes:

- :xml:`<file>`: this **mandatory** node specifies the path of the mesh file via the attribute :xml:`name`. The file must be located in the folder **ECOGEN/libMeshes/**.
- :xml:`<modeParallele>` : This node is required only if the file mesh is a multi-CPU file. The attribute :xml:`GMSHPretraitement` can take the following values:
	- *true*: ECOGEN automatically splits the given mesh file in as many as necessary files according to the number of available CPUs.
	- *false*: do not redo the split of the given mesh (which has already been split in a precedent simulation).

**Remarks:**

1. The attribute :xml:`GMSHPretraitement` must be set as true if it is the first run with the given mesh file.
2. In the current version of ECOGEN, only mesh files generated with the opensource Gmsh_  software under file format *version 2* can be used.

Please refer to the section :ref:`Sec:tuto:generatingMeshes` for learning how to generate a mesh adapted to ECOGEN.

.. _Gmsh: http://gmsh.info/