.. role:: xml(code)
  :language: xml

Droplet impact
===============

This case is a 2D Cartesian test case. A liquid droplet is moving upward in a box with solid walls and filled with air. Initial conditions of this test are shown in :numref:`Fig:testCases:PUEq:dropletImpactIni`. 
Input files for this test are available in *./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/dropletImpact/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/dropletImpact/</testCase>

.. _Fig:testCases:PUEq:dropletImpactIni:

.. figure:: ./_static/testCases/PUEq/dropletImpact/dropletImpactIni.png
  :scale: 30%
  :align: center

  Initial conditions. Visualization using Paraview_ software. 

The initial characteristics of the run are:

+------------------------------+-----------------------------+
| Characteristic               | Value                       |
+==============================+=============================+
| dimension                    | 10 mm x 5 mm                |
+------------------------------+-----------------------------+
| Initial mesh size            | 80 x 40                     |
+------------------------------+-----------------------------+
| AMR max level                | 2                           |
+------------------------------+-----------------------------+
| boundary conditions          | wall                        |
+------------------------------+-----------------------------+
| final solution time          | 8 e-4 s                     |
+------------------------------+-----------------------------+
| solution printing frequency  | 2 e-5 s                     |
+------------------------------+-----------------------------+
| precision                    | 2nd order (VanLeer + THINC) |
+------------------------------+-----------------------------+

Results are shown in :numref:`Fig:testCases:PUEq:dropletImpactAnim` and :numref:`Fig:testCases:PUEq:dropletImpactAnim2`. 
At the impact on the upper wall, the droplet is sliding and deforming, and thanks to surface tension it is then dividing in smaller droplets. Due to the gravity, some small droplets fall down while others are sliding down along the vertical wall. One can observe the rebounds of the droplets as well as the various droplet--droplet interactions leading to coalescence phenomena.

.. _Fig:testCases:PUEq:dropletImpactAnim:

.. figure:: ./_static/testCases/PUEq/dropletImpact/dropletImpactAnim.gif
  :scale: 100%
  :align: center

  Mixture density of a droplet impact over time. Visualization using Paraview_ software.

.. _Fig:testCases:PUEq:dropletImpactAnim2:

.. figure:: ./_static/testCases/PUEq/dropletImpact/dropletImpactAnim2.gif
  :scale: 100%
  :align: center

  AMR evolution during droplet impact. Visualization using Paraview_ software.

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/