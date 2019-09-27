.. _Sec:tuto:generatingMeshes:

Generating Meshes
=================

ECOGEN is able to generate its own Cartesian meshes (see section :ref:`Sec:input:cartesian`), but for structured non Cartesian or unstructured meshes, as seen in section :ref:`Sec:input:unstructured`, an exeternal mesh software is required to generate mesh files.

Mesh files with Gmsh
--------------------

ECOGEN can use mesh files for mono- or multiprocessors computations generated with the opensource Gmsh_ software :cite:`geuzaine2009gmsh` with some specific precaution when editing the geometry file (.geo). Only the `MSH file format version 2`_ can be used in the current released version of ECOGEN. 

Download binaries of Gmsh in version 3.0.6 or lower : http://gmsh.info/bin/.

Here are the restrictions that should be used when generating a geometry with Gmsh_:

- Each part of the domain occupied by the fluid should correspond to a physical surface or a physical volume which is attribute to the value *10*.
- Each boundary condition must correspond to a physical line or a physical surface. The values are successively taken from *1* to maximum *9*. This is an important point that will be used to define boundary conditions physical treatment in the *initialConditionsV4.xml* input file described in section :ref:`Sec:input:InitialConditions`.

Below are presented examples:

Generating a nozzle structured mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider the geometry file of a simple nozzle depicted below:

.. _Fig:tutorials:nozzle_simple:

.. figure:: ./_static/tutos/gmsh/simpleNozzle.png

	Example of geometrical data file - nozzle2D_simple2.geo – for generating a mesh file at .msh format using Gmsh software and usable with ECOGEN.

This correponds to the geometry file available in `ECOGEN/libMeshes/nozzles/nozzle2D_simple2.geo`_.

The computational domain is a nozzle, the mesh is unstructured with quadrangles. In that case, one should take care that the fluid surface is defined by the value *10* and that 4 boundaries conditions are set with the following numerotation:

- Symetrical axis:  1 (condLimAxe =1)
- Wall: 		 	2 (condLimParoi =2)
- Inflow: 		 	3 (condLimEntree =3)
- Outflow: 	 		4 (condLimSortie =3)

The corresponding *initialConditionV4.xml* should then contains for example the following markup:

.. code-block:: xml

	<!-- LIST OF BOUNDARY CONDITIONS -->
	<boundaryConditions>
		<boundCond name="axe" type="wall" number="1" />
		<boundCond name="wall" type="wall" number="2" />
		<boundCond name="entrance" type="tank" number="3">
			<dataTank p0="2.e5" T0="187.5"/>
		</boundCond>
		<boundCond name="exit" type="outflow" number="4">
			<dataOutflow p0="1.e5"/>
		</boundCond>	
	</boundaryConditions>

That correponds to initialization of a nozzle connected to a tank on the left and to an outflow at imposed pressure to the right.

.. _Gmsh: http://gmsh.info/
.. _`MSH file format version 2`: http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
.. _`ECOGEN/libMeshes/nozzles/nozzle2D_simple2.geo`: https://github.com/code-mphi/ECOGEN/blob/master/libMeshes/nozzles/nozzle2D_simple2.geo