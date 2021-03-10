.. role:: xml(code)
	:language: xml

.. _Sec:input:mesh:

MeshV5.xml
==========

Input file *meshV5.xml* is necessary to specify the geometrical characteristics of the computational domain and the type of mesh used. This file is **mandatory** and must be present in the folder of the current case. The minimalist content of this file is:

.. code-block:: xml

	<?xml version="1.0" encoding="UTF-8" standalone = "yes"?>
	<mesh>
	  <type structure="cartesian"/>
	  <cartesianMesh>
	    <dimensions x="1.e-1" y="5.e-2" z="1."/>
	    <numberCells x="50" y ="25" z="1"/>
	  </cartesianMesh>
	</mesh>
 
The  :xml:`<type>` markup is mandatory and specifies via the attribute :xml:`structure` the type of mesh used in the current computation:

- *cartesian*: ECOGEN automatically generates its own Cartesian mesh. See section :ref:`Sec:input:cartesian` for details.
- *unstructured*: A specific external mesh generator (not provide with ECOGEN) must be used to generate the mesh.  See section :ref:`Sec:input:unstructured` for details.

.. _Sec:input:cartesian:

Cartesian mesh
--------------

.. code-block:: xml

  <cartesianMesh>
    <dimensions x="1.e-1" y="5.e-2" z="1."/>
    <numberCells x="50" y ="25" z="1"/>
  </cartesianMesh>

ECOGEN is able to automatically generate a Cartesian mesh according to the markup :xml:`<cartesianMesh>` which must contain the two following nodes:

- :xml:`<dimensions>`: Specify the physical dimensions of the computational domain in each physical direction corresponding to the attributes :xml:`x`, :xml:`y` and :xml:`z` (unit: m (SI)). Values must be real numbers.
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

Stretching can be set optionally by adding the :xml:`<XStretching>` node to the :xml:`<cartesianMesh>` parent markup (in this example for stretching in the *x*-direction). It should contain one or more :xml:`<stretch>` node(s) equipped with the following attributes:

- :xml:`startAt`: Real number giving the starting position of the stretched region.
- :xml:`endAt`: Real number giving the ending position of the stretched region.
- :xml:`factor`: Real number for the stretch factor (lower than 1 for shrinking, greater than 1 for stretching).
- :xml:`numberCells`: Integer for the cell number attributed to the stretched region.

**Remarks:** 

1. Stretching can be set in each direction using :xml:`<XStretching>`, :xml:`<YStretching>` and :xml:`<ZStretching>`. 
2. For each stretched direction, the sum of stretched regions should exactly recover the entire domain without overlapping, but the number of cells can differ those specified in the :xml:`<numberCells>` initial markup.
3. A particular attention should be paid to the link between regions which can possibly present bad quality.

Optional AMR
~~~~~~~~~~~~

.. code-block:: xml

	<AMR lvlMax="2" criteriaVar="0.2" varRho="true" varP="true" varU="false" varAlpha="false" xiSplit="0.11" xiJoin="0.11"/> <!-- Optionnal node -->

An efficient Adaptive Mesh refinement (AMR) technology is embedded in ECOGEN :cite:`schmidmayer2019adaptive`. To use it, the *meshV5.xml* file must contain the optional node :xml:`<AMR>` of the :xml:`<cartesianMesh>` markup and define the following attributes:

- :xml:`lvlMax`: Integer to define the maximal number of refinements (levels).
- :xml:`criteriaVar`: Real number controlling the detection of gradients for the locations of refinement.
- :xml:`varRho`, :xml:`varP`, :xml:`varU`, :xml:`varAlpha`: Boolean (*true* or *false*) to select the flow quantities on which the gradient operator is applied to detect large gradient.
- :xml:`xiSplit`, :xml:`xiJoin`: Normalized real numbers to control if a computational cell, selected by its high gradient value, must be refined or un-refined (values are in the range *0-1*.).

**Remark:**

The global efficiency of the method is greatly depending on the chosen values for the :xml:`criteriaVar`, :xml:`xiSplit` and :xml:`xiJoin` attributes. These values depend on the physical problem and required a real *know-how*. More details about these criterion values can be found in :cite:`schmidmayer2019adaptive`.

.. _Sec:input:unstructured:

Unstructured mesh
-----------------

.. code-block:: xml

	<unstructuredMesh>
	  <file name="./libMeshes/unstructured2D/testUS.msh"/>
	  <parallel GMSHPretraitement="true"/>  <!-- Optionnal node if multi-core -->
	</unstructuredMesh>

When dealing with unstructured meshes, the :xml:`<unstructuredMesh>` markup **must be** present in the *meshV5.xml* input file and it contains the following nodes:

- :xml:`<file>`: This **mandatory** node specifies the path of the mesh file via the attribute :xml:`name`. The file must be located in the folder **ECOGEN/libMeshes/**.
- :xml:`<modeParallele>`: This node is required only if the mesh file is a multi-core file. The attribute :xml:`GMSHPretraitement` can take the following values:
		
	- *true*: ECOGEN automatically splits the given mesh file in as many as necessary files according to the number of available cores.
	- *false*: Do not redo the split of the given mesh (which has already been split in a previous simulation).

**Remarks:**

1. The attribute :xml:`GMSHPretraitement` must be set as *true* if this is the first run with the given mesh file.
2. In the current version |version| of ECOGEN, only mesh files generated with the opensource Gmsh_ software :cite:`geuzaine2009gmsh` under file format *version 2* can be used.

Please refer to the section :ref:`Sec:tuto:generatingMeshes` to learn how to generate a mesh adapted for ECOGEN.

.. _Gmsh: http://gmsh.info/