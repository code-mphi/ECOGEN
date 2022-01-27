.. role:: xml(code)
  :language: xml

Shock tubes
===========

Water--Air shock tube
---------------------
A shock tube between a high-pressure chamber filled with water and a low-pressure chamber filled with air is released. Input files for this test are available in *./libTests/referenceTestCases/PUEq/1D/shockTubes/interfaceWaterAir/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/1D/shockTubes/interfaceWaterAir/</testCase>

.. _Fig:testCases:PUEq:shockTubeWaterAirCI:

.. figure:: ./_static/testCases/PUEq/shockTubeWaterAir/schemaCI.png
  :scale: 70%
  :align: center

  Initial condition for 1D water--air shock tube.

The initial characteristics of the run are:

+------------------------------+---------------------------+
| Characteristic               | Value                     |
+==============================+===========================+
| dimension                    | 1 m                       |
+------------------------------+---------------------------+
| Initial mesh size            | 100 x 230                 |
+------------------------------+---------------------------+
| number of refinement level   | 4                         |
+------------------------------+---------------------------+
| diaphragm position           | 0.7 m                     |
+------------------------------+---------------------------+
| boundary conditions          | non-reflecting            |
+------------------------------+---------------------------+
| final solution time          | 0.240 ms                  |
+------------------------------+---------------------------+
| solution printing frequency  | 0.012 ms                  |
+------------------------------+---------------------------+
| precision                    | 2nd order (VanLeer/THINC) |
+------------------------------+---------------------------+

The AMR technique of :cite:`schmidmayer2019adaptive` is used with 4 refinement levels such that a maximum of 230 computational cells are used for this run.

.. _Fig:testCases:PUEq:shockTubeWaterAir:

.. figure:: ./_static/testCases/PUEq/shockTubeWaterAir/shockTubeWaterAir.*
  :scale: 50%
  :align: center

  Water--air shock tube. Visualization using Paraview_ software.

Epoxy--Spinel shock tube
------------------------

This test deals with shocks in mixture of materials. Epoxy and spinel are supposed mixed such that caracteristic times for wave propagation and drag effects are very small, allowing to consider the mixture as evolving in mechanical equilibrium. Input files for this test are available in *./libTests/referenceTestCases/PUEq/1D/shockTubes/epoxySpinel/*. The corresponding uncommented line in *ECOGEN.xml* is: 

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/1D/shockTubes/epoxySpinel/</testCase>

.. _Fig:testCases:PUEq:shockTubeEpoSpiCI:

.. figure:: ./_static/testCases/PUEq/shockTubeEpoSpi/schemaCI.png
  :scale: 70%
  :align: center

  Initial condition for mixture shock tube with epoxy and spinel.

The initial characteristics of the run are:

+------------------------------+---------------------+
| Characteristic               | Value               |
+==============================+=====================+
| dimension                    | 1 m                 |
+------------------------------+---------------------+
| Initial mesh size            | 200 x 237           |
+------------------------------+---------------------+
| number of refinement level   | 2                   |
+------------------------------+---------------------+
| diaphragm position           | 0.6 m               |
+------------------------------+---------------------+
| boundary conditions          | non-reflecting      |
+------------------------------+---------------------+
| final solution time          | 0.1 ms              |
+------------------------------+---------------------+
| solution printing frequency  | 0.025 ms            |
+------------------------------+---------------------+
| precision                    | 2nd order (VanLeer) |
+------------------------------+---------------------+

.. _Fig:testCases:PUEq:shockTubeEpoSpi:

.. figure:: ./_static/testCases/PUEq/shockTubeEpoSpi/shockTubeEpoSpi.*
  :scale: 50%
  :align: center

  Mixture shock tube with expoxy and spinel. Visualization using Paraview_ software.

Other shock-tube test case 
--------------------------

Another 1D shock-tube test is provided within the ECOGEN package and may be described in details later.

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/1D/shockTubes/interfaceWaterAirNASG/</testCase>

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/
