.. role:: xml(code)
  :language: xml

Square water explosion 
======================

This case is a 2D Cartesian test case.
The square domain is composed of low-pressure air and high-pressure water is present at its center.
Initial conditions are shown in :numref:`Fig:testCases:PUEq:squareWaterExplosion`.
This test is referenced in *./libTests/referenceTestCases/PUEq/2D/squareWaterExplosion/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/2D/squareWaterExplosion/</testCase>

.. _Fig:testCases:PUEq:squareWaterExplosion:

.. figure:: ./_static/testCases/PUEq/squareWaterExplosion/squareWaterExplosionIni.*
  :scale: 70%
  :align: center

  Initial condition before explosion.

The initial characteristics of the run are:

+------------------------------+---------------------+
| Characteristic               | Value               |
+==============================+=====================+
| dimension                    | 1 m                 |
+------------------------------+---------------------+
| Initial mesh size            | 80 x 80             |
+------------------------------+---------------------+
| boundary conditions          | non-reflecting      |
+------------------------------+---------------------+
| final solution time          | 10 ms               |
+------------------------------+---------------------+
| solution printing frequency  | 0.2 ms              |
+------------------------------+---------------------+
| precision                    | Minmod              |
+------------------------------+---------------------+

Results are shown in :numref:`Fig:testCases:PUEq:squareWaterExplosionAnim`. There is a pressure gradient between the two regions. When the interface between the blue (air) and red (water) regions disapears, shock waves propagate in the air and the red region relaxes. This causes an increase in pressure inside the blue region.

With the help of the density evolution, shown in :numref:`Fig:testCases:PUEq:squareWaterExplosionAnim` through a Schlieren visualization, it is easy to observe that there is an explosion between the water and the air regions.
The shock caused by the explosion propagates outward.
After that, one can observe from the center (where the square was) phenomena reminiscent of Richtmyer--Meshkov instabilities.


.. _Fig:testCases:PUEq:squareWaterExplosionAnim:

.. figure:: ./_static/testCases/PUEq/squareWaterExplosion/squareWaterExplosionAnim.*
  :scale: 40%
  :align: center

  Water explosion over time. Schlieren visualization using Paraview_ software.


.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/