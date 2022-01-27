.. role:: xml(code)
  :language: xml

Air shock on helium bubble
==========================

This is a 2D Cartesian test case representing the interaction between an air shock and a helium bubble.
Even before the shock, the bubble contains lighter fluid (helium) than the surrounding gas (air). 
Initial conditions are shown in :numref:`Fig:testCases:PUEq:HeliumAirIni`. This test is referenced in *./libTests/referenceTestCases/PUEq/2D/shockBubble/heliumAir/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/2D/shockBubble/heliumAir/</testCase>

.. _Fig:testCases:PUEq:HeliumAirIni:

.. figure:: ./_static/testCases/PUEq/shockBubble/heliumAir/heliumAirIni.*
  :scale: 70%
  :align: center

  Initial conditions.

The initial characteristics of the run are:

+------------------------------+-----------------------------+
| Characteristic               | Value                       |
+==============================+=============================+
| Initial mesh structure :     | Cartesian                   |
+------------------------------+-----------------------------+
| dimension                    | 0.3 m                       |
+------------------------------+-----------------------------+
| Initial mesh size            | 50 x 10                     |
+------------------------------+-----------------------------+
| boundary conditions          | non-reflecting and symmetry |
+------------------------------+-----------------------------+
| final solution time          | 700e-6 s                    |
+------------------------------+-----------------------------+
| solution printing frequency  | 700e-7 s                    |
+------------------------------+-----------------------------+
| precision                    | Van--Leer                   |
+------------------------------+-----------------------------+

Results are shown in :numref:`Fig:testCases:PUEq:HeliumAirAnim`.
At the beginning, the shock wave impacts the bubble and starts distorting it.
After the passage of the shock, the bubble is then surrounded by a post-shock state.
Kelvin-Helmholtz instabilities form at the bubble interface and then the bubble breaks up while moving in the same direction than the shock wave.
One can also note that the shock wave is deformed when passing through the helium bubble.
It somehow takes the shape of the bubble (similar to a parabola) and then it progressively comes back to a planar shape.

.. _Fig:testCases:PUEq:HeliumAirAnim:

.. figure:: ./_static/testCases/PUEq/shockBubble/heliumAir/heliumAirAnim.*
  :scale: 70%
  :align: center

  Evolution of a helium bubble impacted by a shock wave. Visualization using Paraview_ software.

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/