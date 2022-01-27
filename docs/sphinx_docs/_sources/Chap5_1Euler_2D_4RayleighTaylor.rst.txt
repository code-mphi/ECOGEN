.. role:: xml(code)
  :language: xml

****************
Rayleigh--Taylor
****************

This is an interface instability between two fluids of different densities which occurs when the lighter fluid is pushing the heavier fluid. Initial conditions are described in :numref:`Fig:testCases:Euler:RayleighTaylor`.
This case is a 2D Cartesian test case. The domain is only composed of gases called "light gas" and "heavy gas", with different densities. 
This test is referenced in *./libTests/referenceTestCases/euler/2D/RayleighTaylor/*. The corresponding uncommented line in *ECOGEN.xml* is: 

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/euler/2D/RayleighTaylor/</testCase>

The initial characteristics of the run are: 

+-----------------------------+----------------------+
| Characteristic              | Value                |
+=============================+======================+
| Dimension                   | 0.2 m x 1.2 m        |
+-----------------------------+----------------------+
| Initial mesh size           | 20 x 120             |
+-----------------------------+----------------------+
| AMR max level               | 2                    |
+-----------------------------+----------------------+
| Boundary conditions         | wall                 |
+-----------------------------+----------------------+
| Final solution time         | 5 s                  |
+-----------------------------+----------------------+
| Solution printing frequency | 0.50 s               |
+-----------------------------+----------------------+
| Precision                   | 2nd order (VanLeer)  |
+-----------------------------+----------------------+

.. _Fig:testCases:Euler:RayleighTaylor:

.. figure:: ./_static/testCases/Euler/RayleighTaylor/RTInitial.*
  :scale: 40%
  :align: center

  Initial gas conditions before instabilities.

Results are shown in :numref:`Fig:testCases:Euler:RayleighTaylor:Anim`. When the instability develops, irregularities (“dimples”) propagate downwards in Rayleigh--Taylor polyps that eventually mix. At the first stage, the perturbation amplitudes are small when compared to their wavelengths, then one observes the beginnings of the formation of the ubiquitous mushroom-shaped spikes (fluid structures of heavy fluid growing into light fluid) and bubbles (fluid structures of light fluid growing into heavy fluid). The growth of the mushroom structures continues in the second stage.

.. _Fig:testCases:Euler:RayleighTaylor:Anim:

.. figure:: ./_static/testCases/Euler/RayleighTaylor/RTAnimation.*
  :scale: 40%
  :align: center

  Rayleigh--Taylor instability. Visualization using Paraview_ software.


.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/