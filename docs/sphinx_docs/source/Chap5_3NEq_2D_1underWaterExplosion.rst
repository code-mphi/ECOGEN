.. role:: xml(code)
  :language: xml

Under water explosion
======================

This case is a 2D Cartesian test.
The square domain is composed of air (blue), water (orange) and air cavity (red) as represented in :numref:`Fig:testCases:NEq:underWaterExplosion`.
This test is referenced in *./libTests/referenceTestCases/NEqDEM/2D/underWaterExplosion/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/NEqDEM/2D/underWaterExplosion/</testCase>

.. _Fig:testCases:NEq:underWaterExplosion:

.. figure:: ./_static/testCases/NEqDEM/2D/underWaterExplosion/initialConditions.*
  :scale: 40%
  :align: center

  Representation of density mixture initial conditions for 2D under water explosion.

.. math::
  :nowrap:

  \[
  \begin{array}{|c|c|c|c|}
      \hline
      & \text{Water} & \text{Air} & \text{Air cavity} \\ \hline
      \rho \, [\text{kg/m}^3] & 1000 & 1.225 & 1250 \\ \hline
      u \, [\text{m/s}] & 0 & 0 & 0 \\ \hline
      P \, [\text{Pa}] & 1.01325 \times 10^{5} & 1.01325 \times 10^{5} & 10^{9} \\ \hline
  \end{array}
  \]

  Initial conditions for 2D under water explosion.

The initial characteristics of the run are:

.. math::
  :nowrap:

  \begin{array}{|l|l|}
  \hline
  \textbf{Characteristic} & \textbf{Value} \\
  \hline
  \text{dimension } \mathrm{[m]} & 4 \times 3 \\
  \text{mesh size} & 100 \times 75 \\
  \text{boundary conditions} & \text{non-reflecting} \\
  \text{final solution time} & 2 \times 10^{-3} \ \mathrm{s}  \\
  \text{solution printing frequency} & 1 \times 10^{-4} \ \mathrm{s} \\
  \text{precision} & \text{1st order} \\
  \hline
  \end{array}

Results are shown in :numref:`Fig:testCases:NEq:gifSchlierenVolumeFraction`. Initially, the air cavity is at high pressure. When the interface between the air cavity and the water disapears, a shock wave propagates in the water and the air-cavity region relaxes. Part of this shock wave is transmitted in to the air region above.
With the help of the density evolution, shown through a Schlieren visualization, we can observe a kind of explosion between the water and the air regions.

.. _Fig:testCases:NEq:gifSchlierenVolumeFraction:

.. figure:: ./_static/testCases/NEqDEM/2D/underWaterExplosion/volumeFractionSchlieren.gif
  :align: center
  :width: 100%

  Schlieren visualization and volume fraction of under water explosion.


.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/
