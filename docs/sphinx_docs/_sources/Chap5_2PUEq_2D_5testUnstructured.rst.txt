.. role:: xml(code)
  :language: xml

Test unstructured
==================

This is a case with an unstructured mesh generated with the open-source Gmsh_ software. Initial conditions are shown in :numref:`Fig:testCases:PUEq:testUnstructuredIni`. This test is referenced in *./libTests/referenceTestCases/PUEq/2D/testUnstructured/*. The corresponding uncommented line in *ECOGEN.xml* is: 

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/2D/testUnstructured/</testCase>


.. _Fig:testCases:PUEq:testUnstructuredIni:

.. figure:: ./_static/testCases/PUEq/testUnstructured/testUnstructuredIni.*
  :scale: 40%
  :align: center

  Initial conditions. Visualization using Paraview_ software.

The initial characteristics of the run are: 

+-----------------------------+----------------------+
| Characteristic              | Value                |
+=============================+======================+
| Initial mesh structure      | Unstructured         |
+-----------------------------+----------------------+
| Boundary conditions         | wall                 |
+-----------------------------+----------------------+
| Final solution time         | 2e-2 s               |
+-----------------------------+----------------------+
| Solution printing frequency | 4e-4 s               |
+-----------------------------+----------------------+

Results are shown in :numref:`Fig:testCases:PUEq:testUnstructuredAnim` and :numref:`Fig:testCases:PUEq:testUnstructuredAnim2`. There is a pressure gradient between the two regions, high pressure at t=0s is contained in red region. When the interface between the blue and the red regions disapears, a shock wave propagates, causing an increase in pressure in the blue region, and the red region relaxes.

.. _Fig:testCases:PUEq:testUnstructuredAnim:

.. figure:: ./_static/testCases/PUEq/testUnstructured/testUnstructuredAnim.*
  :scale: 40%
  :align: center

  Pressure evolution over time. Visualization using Paraview_ software.

.. figure:: ./_static/testCases/PUEq/testUnstructured/pressureToggle.*
  :scale: 40%
  :align: center

.. _Fig:testCases:PUEq:testUnstructuredAnim2:

.. figure:: ./_static/testCases/PUEq/testUnstructured/testUnstructuredAnim2.*
  :scale: 40%
  :align: center

  Mixture density over time. Visualization using Paraview_ software.

.. figure:: ./_static/testCases/PUEq/testUnstructured/densityToggle.*
  :scale: 40%
  :align: center


.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/
.. _Gmsh: http://gmsh.info/