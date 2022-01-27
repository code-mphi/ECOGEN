.. role:: xml(code)
  :language: xml

******************
Rectangle diagonal
******************

This test is a 2D Cartesian test case with advection of a rectangle of high-density air into low-density air environment. The computation uses the AMR technique of :cite:`schmidmayer2019adaptive`. This test is referenced in *./libTests/referenceTestCases/euler/2D/transports/rectangleDiagonal/*. A sketch of the initial conditions for this test is presented in :numref:`Fig:testCases:Euler:2DtransportCI`.
The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/euler/2D/transports/rectangleDiagonal/</testCase>

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

Because of the adaptative mesh refinement, the final number of computational cells is 5383. Results are shown in :numref:`Fig:testCases:Euler:2Dtransport`. The left picture shows the evolution of the mesh during the simulation (red color is for high gradients of density). Upper right picture represents density contours and lower right picture is to check for constant pressure preservation along the diagonal.

.. _Fig:testCases:Euler:2Dtransport:

.. figure:: ./_static/testCases/Euler/transport2D/transport2D.*
  :scale: 70%
  :align: center

  Advection of a high-density rectangle of air. Visualization using Paraview_ software.


.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/