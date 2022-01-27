.. role:: xml(code)
  :language: xml


Square to circle using symmetry
===============================

This case is a 2D Cartesian test case. A square domain composed of air, with inside another square domain composed of liquid (special for few surface-tension test cases). Initial conditions of this test are shown in :numref:`Fig:testCases:PUEq:squareToCircleSymmetryIni`. 
Input files for this test are available in *./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/squareToCircleSymmetry/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/squareToCircleSymmetry/</testCase>

.. _Fig:testCases:PUEq:squareToCircleSymmetryIni:

.. figure:: ./_static/testCases/PUEq/squareToCircleSymmetry/squareToCircleSymmetryIni.*
  :scale: 30%
  :align: center

  Initial conditions.

Computations are performed only on half of the domain because a symmetry is applied as indicated in the title of this test case.

The initial characteristics of the run are:

+------------------------------+---------------------------------+
| Characteristic               | Value                           |
+==============================+=================================+
| dimension                    | 75 cm x 75 cm                   |
+------------------------------+---------------------------------+
| Initial mesh size            | 20 x 10                         |
+------------------------------+---------------------------------+
| AMR max level                | 2                               |
+------------------------------+---------------------------------+
| boundary conditions          | outflow, symmetry               |
+------------------------------+---------------------------------+
| final solution time          | 0.5 s                           |
+------------------------------+---------------------------------+
| solution printing frequency  | 5 e-3 s                         |
+------------------------------+---------------------------------+
| precision                    | 2nd order (MC + Minmod + THINC) |
+------------------------------+---------------------------------+

Results are shown in :numref:`Fig:testCases:PUEq:squareToCircleSymmetryAnim`. When the interface between the blue and the white regions disapears, because of surface tension, the region that was initially square is deformed and tends to become circular with time. Note that it is necessary to have a stationary solution to obtain this circle and this requires a relatively long final solution time.

One can observe the distribution of computation on the CPUs (here 12) divided into different regions, according to AMR evolution. When the density variations are low or non-existent, the mesh is as coarse as possible (within the conditions given in meshV5.xml). On the contrary when the variations become significant, the mesh is refined. So, to distribute the work equally between the CPUs, highly refined regions are relatively small compared to coarse regions.

.. _Fig:testCases:PUEq:squareToCircleSymmetryAnim:

.. figure:: ./_static/testCases/PUEq/squareToCircleSymmetry/squareToCircleSymmetryAnim.*
  :scale: 45%
  :align: center

  CPU (upper half) and density (lower half) over time. Visualization using Paraview_ software. 

.. figure:: ./_static/testCases/PUEq/squareToCircleSymmetry/squareToCircleToggle.*
  :scale: 70%
  :align: center

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/