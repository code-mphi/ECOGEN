.. role:: xml(code)
  :language: xml

.. _Sec:testcases:shocktube:

******************
High-pressure left
******************

This is a classical single-phase shock tube filled with air. The test is available in the folder *./libTests/referenceTestCases/euler/1D/shockTubes/HPLeft/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/euler/1D/shockTubes/HPLeft/</testCase>

.. _Fig:testCases:Euler:shockTubeCI:

.. figure:: ./_static/testCases/Euler/shockTube/schema.jpg
  :scale: 70%
  :align: center

  Initial condition for single-phase shock tube.

The initial characteristics of the run are:

+-----------------------------+----------------------+
| Characteristic              | Value                |
+=============================+======================+
| dimension                   | 1 m                  |
+-----------------------------+----------------------+
| initial mesh size           | 100                  |
+-----------------------------+----------------------+
| AMR max level               | 3                    |
+-----------------------------+----------------------+
| diaphragm position          | 0.4 m                |
+-----------------------------+----------------------+
| boundary conditions         | non-reflecting       |
+-----------------------------+----------------------+
| final solution time         | 0.6 ms               |
+-----------------------------+----------------------+
| solution printing frequency | 0.06 ms              |
+-----------------------------+----------------------+
| precision                   | 2nd order (VanLeer)  |
+-----------------------------+----------------------+

Solution of this Riemann problem induces 3 waves: 

- expansion waves propagating in high-pressure chamber,
- a right-facing shock wave propagating in low-pressure chamber,
- a contact discontinuity.

These waves are clearly observable on the results:

.. _Fig:testCases:Euler:shockTube:

.. figure:: ./_static/testCases/Euler/shockTube/shockTube.*
  :scale: 70%
  :align: center

  Shock tube filled with air. Visualization using Paraview_ software.

This test is also equipped with 3 Eulerian sensors. For example, two sensors are positionned at :math:`x = 0.6 m` and :math:`x= 0.8 m`. They record the following pressures:

.. _Fig:testCases:Euler:sensors:

.. figure:: ./_static/testCases/Euler/shockTube/sensors.jpg
  :scale: 70%
  :align: center

  Pressure recorded by sensors at :math:`x = 0.6 m` (pink) and :math:`x = 0.8 m` (green). Visualization using gnuplot_ software.

The first sensor sees its pressure rising first because of the shock wave. Because this Riemann problem generates a supersonic flow after the shock wave, the tail of the expansion waves is seen by the sensor after 0.5 ms.

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/
