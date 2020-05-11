.. role:: xml(code)
  :language: xml

***********************
Single-phase test cases
***********************

Test cases presented in this section are dealing with single-phase compressible problems. In this part, ECOGEN solves the Euler equations :cite:`euler1757principes`:

.. math::
  :nowrap:
  
  \[
  \label{eqEuler}
  \begin{array}{l}
    \displaystyle \frac{\partial\rho}{\partial t}+ div(\rho {\mathbf{u}} )=0\\
    \displaystyle\frac{\partial\rho {\mathbf{u}}}{\partial t}+ div(\rho {\mathbf{u}} \otimes {\mathbf{u}} +p \mathbf{I})=\mathbf{0} \\ 
    \displaystyle \frac{\partial\rho E}{\partial t}+ div((\rho E+p){\mathbf{u}})=0
  \end{array} 
  \]

where :math:`\rho` represents the density, :math:`\mathbf{u}`` the velocity vector, :math:`p` the pressure and :math:`E = e +\frac{\mathbf{u}^2}{2}` the total energy, with :math:`e` the internal energy. 
The closure relation for this model is ensured by any convex equation of state (EOS) :math:`p = p(\rho,e)` (see section :ref:`Sec:IO:materials` for details about implemented EOS in ECOGEN).
Euler equations are solved thanks to an explicit finite-volume Godunov-like scheme :cite:`godunov79` that is coupled with HLLC approximate Riemann solver :cite:`toro2013riemann` for flux computation.

Advection test cases
====================
Basic advection test cases are proposed:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/euler/1D/transport/positiveVelocity/</testCase>
  <testCase>./libTests/referenceTestCases/euler/1D/transport/negativeVelocity/</testCase>
  <testCase>./libTests/referenceTestCases/euler/2D/transports/rectangleDiagonal/</testCase>

The first one is the default test case fully described in the tutorial section :ref:`Sec:tuto:default`. The second one is the reverse test advecting the contact discontinuity in the opposite direction. The third is described below.

rectangleDiagonal
-----------------
This test is a 2D Cartesian test case with advection of a rectangle of high-density air into low-density air environment. The computation uses the AMR technique of :cite:`schmidmayer2019adaptive`. This test is referenced in *./libTests/referenceTestCases/euler/2D/transports/rectangleDiagonal/*. A sketch of the initial conditions for this test is presented in figure :numref:`Fig:testCases:Euler:2DtransportCI`.

.. _Fig:testCases:Euler:2DtransportCI:

.. figure:: ./_static/testCases/Euler/transport2D/schema.jpg
  :scale: 70%
  :align: center

  Initial conditions for the 2D, single-phase transport test case.

The initial characteristics of the run are:

+-----------------------------+----------------------+
| Characteristic              | Value                |
+=============================+======================+
| dimensions                  | 1 m x 1 m            |
+-----------------------------+----------------------+
| initial mesh size           | 50 x 50              |
+-----------------------------+----------------------+
| AMR max level               | 2                    |
+-----------------------------+----------------------+
| rectangle position          | 0.2 m, 0.2 m         |
+-----------------------------+----------------------+
| boundary conditions         | non-reflecting       |
+-----------------------------+----------------------+
| final solution time         | 0.36 ms              |
+-----------------------------+----------------------+
| solution printing frequency | 0.036 ms             |
+-----------------------------+----------------------+
| precision                   | 2nd order (superbee) |
+-----------------------------+----------------------+

Because of the adaptative mesh refinement, the final number of computational cells is 5383. Results are shown on figure :numref:`Fig:testCases:Euler:2Dtransport`. The left picture shows the evolution of the mesh during the simulation (red color is for high gradients of density). Upper right picture represents density contours and lower right picture is to check for constant pressure preservation along the diagonal.

.. _Fig:testCases:Euler:2Dtransport:

.. figure:: ./_static/testCases/Euler/transport2D/transport2D.*
  :scale: 70%
  :align: center

  Advection of a high-density rectangle of air. Visualization using Paraview_ software.

Shock tubes
===========
Single-phase shock tubes are proposed in the following test cases:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/euler/1D/shockTubes/HPLeft/</testCase>
  <testCase>./libTests/referenceTestCases/euler/1D/shockTubes/HPRight/</testCase>
  <testCase>./libTests/referenceTestCases/euler/2D/HPCenter/</testCase>
  <testCase>./libTests/referenceTestCases/euler/2D/HPUnstructured/</testCase>

HPLeft
------
This is a classical single-phase shock tube filled with air. The test is available in the folder *./libTests/referenceTestCases/euler/1D/shockTubes/HPLeft/*.

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
| precision                   | 2nd order (vanleer)  |
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

The first sensor see its pressure rising first because of the shock wave. Because this Riemann problem generates a supersonic flow after the shock wave, the tail of the expansion waves is seen by the sensor after 0.5 ms.

Other test cases
================
Other tests are provided within the ECOGEN package. They will soon be described in details.

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/euler/2D/HPCenter/</testCase>
  <testCase>./libTests/referenceTestCases/euler/2D/HPUnstructured/</testCase>
  <testCase>./libTests/referenceTestCases/euler/2D/nozzles/tankWithShock/</testCase>
  <testCase>./libTests/referenceTestCases/euler/3D/LPCenter/</testCase>

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/