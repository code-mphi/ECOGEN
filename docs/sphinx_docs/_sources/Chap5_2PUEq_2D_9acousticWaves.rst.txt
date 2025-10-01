.. role:: xml(code)
  :language: xml

**************
Acoustic waves
**************

This test case is similar to the Euler one: :doc:`Single-phase acoustic waves <Chap5_1Euler_2D_5acousticWaves>`. However, a bubble is added in the middle of the domain. To recall, this is a plain pressure field where a source term is used to generate one-way acoustic waves. This work constitutes an adaptation of the work of Maeda & Colonius for single-phase flows :cite:`maeda2017source` to a multiphase system of equations. The waves here have a sinusoidal shape. 4 pulses are sent at a frequency of 300 kHz and have a pressure amplitude of 10 kPa.
This case is a 2D Cartesian test case. The domain is only composed of water. The source term is located at 10.e-3 m in *x* and *y* directions and pulses are directed toward the *x* direction.
This test is referenced in *./libTests/referenceTestCases/PUEq/2D/acousticWave/planeSinusoidalPulse/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/2D/acousticWave/planeSinusoidalPulse/</testCase>

The initial characteristics of the run are:

+-----------------------------+----------------------+
| Characteristic              | Value                |
+=============================+======================+
| Dimension                   | 40.e-3 m x 20.e-3 m  |
+-----------------------------+----------------------+
| Initial mesh size           | 200 x 100            |
+-----------------------------+----------------------+
| AMR max level               | 0                    |
+-----------------------------+----------------------+
| Boundary conditions         | non-reflecting       |
+-----------------------------+----------------------+
| Final solution time         | 25.e-6 s             |
+-----------------------------+----------------------+
| Solution printing frequency | 5.e-7 s              |
+-----------------------------+----------------------+
| Precision                   | 2nd order (MC)       |
+-----------------------------+----------------------+

Results are shown in :numref:`Fig:testCases:PUEq:AcousticWave:Anim` where one can observe the one-way pulses from a point source, therefore spreading cylindrically (2D here) but with a main direction (*x* direction). The pulses are reflected from the interactions with the bubble (in black).

.. _Fig:testCases:PUEq:AcousticWave:Anim:

.. figure:: ./_static/testCases/PUEq/acousticWave/sinusPulseAnimation.*
  :scale: 120%
  :align: center

  One-way sinusoidal acoustic waves in water interacting with a bubble. Visualization using Paraview_ software.


.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/
