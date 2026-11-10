.. role:: xml(code)
  :language: xml

Cavitation
===========

The phenomena of cavitation is explained in :ref:`pueq-cavitation`.

This section presents 3 different tests of cavitation, based on Baer-Nunziato model :cite:`baer`, without relaxations, with infinite velocity relaxation and with infinite velocity and pressure relaxations,
realised with the Discrete Equation Method and the HLLC Riemann solver.
The equation of state used is stiffened gas with :math:`\gamma = 1.4` and :math:`P_{\infty} = 0` for air and :math:`\gamma = 4.4` and :math:`P_{\infty} = 6 \times 10^{8}`. The input files for these tests can be found in *./libTests/referenceTestCases/NEqDEM/1D/cavitation/waterAir/* folder.
The line in *ECOGEN.xml* of the 2 tests is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/NEqDEM/1D/cavitation/waterAir</testCase>

**Remark:**
If one wants to simulate cavitation phenomena without infinite relaxation procedure and with the same initial conditions, he must comment the lines relative to relaxation of the *model.xml* file.

.. _Fig:testCases:NEqDEM:shockTube1D:interface:image:

.. figure:: ./_static/testCases/NEqDEM/1D/cavitation/initialConditions.png
  :scale: 100%
  :align: center

.. math::
  :nowrap:

  \[
  \begin{array}{|c|c|c|c|c|}
      \hline
      & \text{Left state} & \text{Left state} & \text{Right state} & \text{Right state} \\ \hline
      \text{Phase} & \text{Water} & \text{Air} & \text{Water} & \text{Air} \\ \hline
      \rho \, [\text{kg/m}^3] & 1000 & 10 & 1000 & 10 \\ \hline
      u \, [\text{m/s}] & -100 & -100 & 100 & 100 \\ \hline
      P \, [\text{Pa}] & 10^{5} & 10^{5} & 10^{5} & 10^{5} \\ \hline
  \end{array}
  \]

  Initial conditions for 1D cavitation tests cases.

The initial characteristics of the run are:

.. math::
  :nowrap:

  \begin{array}{|l|l|}
  \hline
  \textbf{Characteristic} & \textbf{Value} \\
  \hline
  \text{dimension} & 1 \text{ m} \\
  \text{mesh size} & 1000 \\
  \text{interface position} & 0.5 \text{ m} \\
  \text{boundary conditions} & \text{non-reflecting} \\
  \text{final solution time} & 200 \times 10^{-6} \ \mathrm{s}  \\
  \text{solution printing frequency} & 10 \times 10^{-6} \ \mathrm{s} \\
  \text{precision} & \text{1st order} \\
  \hline
  \end{array}

Legend of the tests is detailed below:

.. math::
  :nowrap:

  \begin{array}{|llll|}
  \hline
  \textbf{Legend} & \textbf{Model} & \textbf{Method} & \textbf{Infinite relaxations} \\
  \hline
  \textbf{DEM}    & 7 \text{ equations} & \text{Discrete Equation Method} & \text{Without relaxations} \\
  \textbf{DEM U}  & 7 \text{ equations} & \text{Discrete Equation Method} & \text{Velocity} \\
  \textbf{DEM UP} & 7 \text{ equations} & \text{Discrete Equation Method} & \text{Velocity and pressure} \\
  \hline
  \end{array}


.. _Fig:testCases:NEqDEM:shockTube1D:interface:gif:

.. container:: figures-row

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/volumeFractionDEM.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/volumeFractionDEMUP.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/pressureDEM.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/pressureDEMUP.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/densityDEM.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/densityDEMUP.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/velocityDEM.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/cavitation/velocityDEMUP.gif
    :align: center
    :width: 75%


.. _gnuplot: http://www.gnuplot.info/
