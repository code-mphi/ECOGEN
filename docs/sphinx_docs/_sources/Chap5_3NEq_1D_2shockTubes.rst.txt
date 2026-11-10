.. role:: xml(code)
  :language: xml

Shock tubes
===========

Water--Air interface
--------------------
This test case is important as it represents a limit case, realised with the Discrete Equation Method and the HLLC Riemann solver. Specifically, it is equivalent to solving the Euler equations with two pure fluids, resulting in an interface between the phases without mixture. Consequently, an exact solution exists, allowing comparison, as shown in :numref:`Fig:testCases:NEqDEM:shockTube1D:equalVolumeFraction:densityMixture`. The equation of state used is stiffened gas with :math:`\gamma = 1.4` and :math:`P_{\infty} = 0` for air and :math:`\gamma = 4.4` and :math:`P_{\infty} = 6 \times 10^{8}`. The input files for this test can be found in *./libTests/referenceTestCases/NEqDEM/1D/shockTubes/interfaceWaterAir/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/NEqDEM/1D/shockTubes/interfaceWaterAir/</testCase>

.. _Fig:testCases:NEqDEM:shockTube1D:interface:image:

.. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/interface/representation.png
  :scale: 100%
  :align: center

.. math::
  :nowrap:

  \[
  \begin{array}{|c|c|c|c|c|}
      \hline
      & \text{High pressure} & \text{High pressure} & \text{Low pressure} & \text{Low pressure} \\ \hline
      \text{Phase} & \text{Water} & \text{Air} & \text{Water} & \text{Air} \\ \hline
      \rho \, [\text{kg/m}^3] & 1000 & 50 & 1000 & 50 \\ \hline
      u \, [\text{m/s}] & 0 & 0 & 0 & 0 \\ \hline
      P \, [\text{Pa}] & 10^{9} & 10^{9} & 10^{5} & 10^{5} \\ \hline
  \end{array}
  \]

  Initial conditions for 1D shock tube of water-air interface.

The initial characteristics of the run are:

.. math::
  :nowrap:

  \begin{array}{|l|l|}
  \hline
  \textbf{Characteristic} & \textbf{Value} \\
  \hline
  \text{dimension} & 1 \text{ m} \\
  \text{mesh size} & 1000 \\
  \text{interface position} & 0.7 \text{ m} \\
  \text{boundary conditions} & \text{non-reflecting} \\
  \text{final solution time} & 241 \times 10^{-6} \ \mathrm{s}  \\
  \text{solution printing frequency} & 12 \times 10^{-6} \ \mathrm{s} \\
  \text{precision} & \text{1st order} \\
  \hline
  \end{array}

.. _Fig:testCases:NEqDEM:shockTube1D:equalVolumeFraction:densityMixture:

.. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/interface/densityMixture.png
  :scale: 33%
  :align: center

  Mixture density comparison between exact solution and Discrete Equation Method.

.. _Fig:testCases:NEqDEM:shockTube1D:interface:gif:

.. container:: figures-row

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/interface/volumeFraction.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/interface/velocity.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/interface/pressure.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/interface/density.gif
    :align: center
    :width: 75%


.. _gnuplot: http://www.gnuplot.info/

Water--Air mixture with equal volume fraction
---------------------------------------------

This test deals with shocks in mixture, realised with the Discrete Equation Method and the HLLC Riemann solver. Since we consider the volume fraction to be perfectly equal, this test is equivalent to solving the Euler equations twice with two pure fluids, allowing for a comparison with an exact solution. The outcome of this comparison is shown in :numref:`Fig:testCases:NEqDEM:shockTube1D:equalVolumeFraction:comparison`, along with another solver developed by Furfaro & Saurel :cite:`FURFARO2015159`. The equation of state used is stiffened gas with :math:`\gamma = 1.4` and :math:`P_{\infty} = 0` for air and :math:`\gamma = 4.4` and :math:`P_{\infty} = 6 \times 10^{8}` for water. Input files for this test are available in *./libTests/referenceTestCases/NEqDEM/1D/shockTubes/HPLeftWaterAirAlphaEqual_1e9/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/NEqDEM/1D/shockTubes/HPLeftWaterAirAlphaEqual_1e9/</testCase>

.. _Fig:testCases:NEqDEM:shockTube1D:equalVolumeFraction:image:

.. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaEqual/shockTubeWaterAirAlphaEqual.png
  :scale: 100%
  :align: center

.. math::
  :nowrap:

  \[
  \begin{array}{|c|c|c|c|c|}
      \hline
      & \text{High pressure} & \text{High pressure} & \text{Low pressure} & \text{Low pressure} \\ \hline
      \text{Phase} & \text{Water} & \text{Air} & \text{Water} & \text{Air} \\ \hline
      \rho \, [\text{kg/m}^3] & 1000 & 50 & 1000 & 50 \\ \hline
      u \, [\text{m/s}] & 0 & 0 & 0 & 0 \\ \hline
      P \, [\text{Pa}] & 10^{9} & 10^{9} & 10^{5} & 10^{5} \\ \hline
  \end{array}
  \]
  Initial conditions for 1D shock tube of water-air mixture with equal volume fraction.

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
  \text{final solution time} & 80 \times 10^{-6} \ \mathrm{s}  \\
  \text{solution printing frequency} & 4 \times 10^{-6} \ \mathrm{s} \\
  \text{precision} & \text{1st order} \\
  \hline
  \end{array}

.. _Fig:testCases:NEqDEM:shockTube1D:equalVolumeFraction:comparison:

.. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaEqual/densityComparison.png
  :align: center
  :width: 75%

  Phase density comparison between exact solution, Discrete Equation Method and Furfaro \& Saurel solver.

.. _Fig:testCases:NEqDEM:shockTube1D:equalVolumeFraction:gif:

.. container:: figures-row

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaEqual/volumeFraction.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaEqual/velocity.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaEqual/pressure.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaEqual/density.gif
    :align: center
    :width: 75%

Water--Air mixture with non-equal volume fraction
-------------------------------------------------

This test deals with shocks in mixture, realised with the Discrete Equation Method and the HLLC Riemann solver. The equation of state used is stiffened gas with :math:`\gamma = 1.4` and :math:`P_{\infty} = 0` for air and :math:`\gamma = 4.4` and :math:`P_{\infty} = 6 \times 10^{8}` for water. Input files for this test are available in *./libTests/referenceTestCases/NEqDEM/1D/shockTubes/HPLeftWaterAirAlphaDiff_2e8/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/NEqDEM/1D/shockTubes/HPLeftWaterAirAlphaDiff_2e8/</testCase>

.. _Fig:testCases:NEqDEM:shockTube1D:nonEqualVolumeFraction:image:

.. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaDiff/shockTubeWaterAirAlphaDiff.png
  :scale: 100%
  :align: center

.. math::
  :nowrap:

  \[
  \begin{array}{|c|c|c|c|c|}
      \hline
      & \text{High pressure} & \text{High pressure} & \text{Low pressure} & \text{Low pressure} \\ \hline
      \text{Phase} & \text{Water} & \text{Air} & \text{Water} & \text{Air} \\ \hline
      \rho \, [\text{kg/m}^3] & 1000 & 50 & 1000 & 50 \\ \hline
      u \, [\text{m/s}] & 0 & 0 & 0 & 0 \\ \hline
      P \, [\text{Pa}] & 2 \times 10^{8} & 2 \times 10^{8} & 10^{5} & 10^{5} \\ \hline
  \end{array}
  \]
  Initial conditions for 1D shock tube of water-air mixture with different volume fraction.

The initial characteristics of the run are:

.. math::
  :nowrap:

  \begin{array}{|l|l|}
  \hline
  \textbf{Characteristic} & \textbf{Value} \\
  \hline
  \text{dimension} & 1 \text{ m} \\
  \text{mesh size} & 400 \\
  \text{interface position} & 0.6 \text{ m} \\
  \text{boundary conditions} & \text{non-reflecting} \\
  \text{final solution time} & 20 \times 10^{-5} \ \mathrm{s}  \\
  \text{solution printing frequency} & 2 \times 10^{-5} \ \mathrm{s} \\
  \text{precision} & \text{1st order} \\
  \hline
  \end{array}

**Remark:** The default test case is realised with adaptive mesh refinement (AMR) with a maximum level of 3.

.. _Fig:testCases:NEqDEM:shockTube1D:nonEqualVolumeFraction:gif:

.. container:: figures-row

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaDiff/volumeFraction.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaDiff/velocity.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaDiff/pressure.gif
    :align: center
    :width: 75%

  .. figure:: ./_static/testCases/NEqDEM/1D/shockTubes/waterAirAlphaDiff/density.gif
    :align: center
    :width: 75%


.. _gnuplot: http://www.gnuplot.info/
