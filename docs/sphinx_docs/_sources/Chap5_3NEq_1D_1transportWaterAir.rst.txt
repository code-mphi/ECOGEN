.. role:: xml(code)
  :language: xml

Transport of water--air mixture
=================================

This test case, realised with the Discrete Equation Method and the HLLC Riemann solver, is important since it validates the capacity of the numerical method to treat advection of a mixture between fluids without creating spurious oscillations on pressure or velocity profiles. The equation of state used is stiffened gas with :math:`\gamma = 1.4` and :math:`P_{\infty} = 0` for air and :math:`\gamma = 4.4` and :math:`P_{\infty} = 6 \times 10^{8}` for water. Input files for this test are available in *./libTests/referenceTestCases/NEqDEM/1D/transport/waterAirAlphaDiff/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>libTests/referenceTestCases/NEqDEM/1D/transport/waterAirAlphaDiff/</testCase>

.. _Fig:testCases:NEqDEM:advection1D:image:

.. figure:: _static/testCases/NEqDEM/1D/transport/advection1D.png
  :scale: 100%
  :align: center


.. math::
  :nowrap:

  \[
  \begin{array}{|c|c|c|c|c|}
      \hline
      & \text{High pressure} & \text{High pressure} & \text{Low pressure} & \text{Low pressure} \\ \hline
      \text{Phase} & \text{Water} & \text{Air} & \text{Water} & \text{Air} \\ \hline
      \rho \, [\text{kg/m}^3] & 1000 & 100 & 100 & 1 \\ \hline
      u \, [\text{m/s}] & 20 & 20 & 20 & 20 \\ \hline
      P \, [\text{Pa}] & 10^{5} & 10^{5} & 10^{5} & 10^{5} \\ \hline
  \end{array}
  \]
  Initial conditions for 1D advection of water-air mixture.

The initial characteristics of the run are:

.. math::
  :nowrap:

  \begin{array}{|l|l|}
  \hline
  \textbf{Characteristic} & \textbf{Value} \\
  \hline
  \text{dimension} & 1 \text{ m} \\
  \text{mesh size} & 400 \\
  \text{interface position} & 0.4 \text{ m} \\
  \text{boundary conditions} & \text{non-reflecting} \\
  \text{final solution time} & 2 \times 10^{-2} \ \mathrm{s}  \\
  \text{solution printing frequency} & 1 \times 10^{-3} \ \mathrm{s} \\
  \text{precision} & \text{1st order} \\
  \hline
  \end{array}

.. _Fig:testCases:NEqDEM:advection1D:gif:

.. container:: figures-row

   .. figure:: ./_static/testCases/NEqDEM/1D/transport/volumeFraction.gif
      :align: center
      :width: 75%

   .. figure:: ./_static/testCases/NEqDEM/1D/transport/density.gif
      :align: center
      :width: 75%

   .. figure:: ./_static/testCases/NEqDEM/1D/transport/pressure.gif
      :align: center
      :width: 75%

   .. figure:: ./_static/testCases/NEqDEM/1D/transport/velocity.gif
      :align: center
      :width: 75%

.. _gnuplot: http://www.gnuplot.info/
