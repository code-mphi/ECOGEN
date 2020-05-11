.. role:: xml(code)
  :language: xml

***************************************
Multiphase mechanical equilibrium flows
***************************************

Mechanical equilibrium flows are solved in ECOGEN using Kapila's model :cite:`kapila2001`. In the particular case of 2 phases involved and without any extra physics (surface tension, viscosity...), this model reads:

.. math::
  :nowrap:

  \begin{equation}
  \label{system_kapila}
  \left\{
  {\begin{array}{*{20}{l}}
    {\cfrac{{\partial {\alpha _1}}}{{\partial t}} + \mathbf{u} \cdot \nabla {\alpha _1}}&{ = K div( \mathbf{u} ),} \\ 
    {\cfrac{{\partial {\alpha _1}{\rho _1}}}{{\partial t}} +  div \left( {{\alpha _1}{\rho _1}\mathbf{u}} \right) } &{ = 0 ,} \\
    {\cfrac{{\partial {\alpha _2}{\rho _2}}}{{\partial t}} + div \left( {{\alpha _2}{\rho _2}\mathbf{u}} \right)}&{ = 0 ,} \\ 
    {\cfrac{{\partial \rho \mathbf{u}}}{{\partial t}} + div \left( {\rho \mathbf{u} \otimes \mathbf{u} + p \mathbf{I}} \right)}&{ = \mathbf{0} ,} \\ 
    {\cfrac{{\partial \rho E}}{{\partial t}} + div \left( {\left( {\rho E + p} \right) \mathbf{u}} \right)}&{ = 0 ,}
  \end{array}} \right.\
  \end{equation}

where subscripts :math:`1` and :math:`2` correspond to one of the two phases, respectively. :math:`\alpha_k` and :math:`\rho_k` are the volume fraction and density of phase :math:`k`. 

:math:`\rho = \sum\limits_{k} \alpha_k \rho_k`, :math:`\mathbf{u}`, :math:`p`, :math:`E = e + \cfrac{1}{2} \| \mathbf{u} \|^2` and :math:`e = \sum_k \alpha_k \rho_k e_k` are the mixture density, velocity, pressure, total energy and internal energy, respectively. 

The term :math:`K div (\mathbf{u})` accounts for the differences in the acoustic behavior of both phases or in other words, for the differences in expansion and compression of each phase in mixture regions. :math:`K` is given by:

.. math::
  :nowrap:

  \begin{equation*}
  K = \cfrac{\rho _2 s_2^2 - \rho _1 s_1^2}{\cfrac{\rho _2 s_2^2}{\alpha _2} + \cfrac{\rho _1 s_1^2}{\alpha _1}},
  \end{equation*}

:math:`s_k` being the speed of sound of phase :math:`k`.

This model is solved thanks to the numerical method presented in :cite:`relaxjcp`.

Advection test cases
====================

The code is provided with the following test cases for advections:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/kapila/1D/transports/interfaceWaterAir/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/2D/transportWaterSquareInAir/</testCase>

interfaceWaterAir
-----------------
This first test case is important since it validates the capacity of the numerical method to treat simple advection of an interface between pure fluid without creating spurious oscillations on pressure or velocity profiles. Input files for this test are available in *./libTests/referenceTestCases/kapila/1D/transports/interfaceWaterAir/*.

.. _Fig:testCases:Kapila:advection1DCI:

.. figure:: ./_static/testCases/Kapila/advection/advection1DCI.png
  :scale: 70%
  :align: center

  Initial condition for 1D advection of water--air interface.

The initial characteristics of the run are:

+-----------------------------+----------------+
| Characteristic              | Value          |
+=============================+================+
| dimension                   | 1 m            |
+-----------------------------+----------------+
| Mesh size                   | 800            |
+-----------------------------+----------------+
| interface position          | 0.3 m          |
+-----------------------------+----------------+
| boundary conditions         | non-reflecting |
+-----------------------------+----------------+
| final solution time         | 0.7 ms         |
+-----------------------------+----------------+
| solution printing frequency | 0.025 ms       |
+-----------------------------+----------------+
| precision                   | 1st order      |
+-----------------------------+----------------+

In the default test case, the computation is performed with a 1st order scheme. We compare this solution with those obtained using the 2nd order scheme with THINC limiter :cite:`shyue2014thinc`.

.. _Fig:testCases:Kapila:advection1D:

.. figure:: ./_static/testCases/Kapila/advection/advection1D.*
  :scale: 50%
  :align: center

  Advection of a water--air interface. Visualization using Paraview_ software.

Shock tubes
===========
The test cases relative to Kapila's model are those presented in :cite:`relaxjcp`. They are here reproduced using ECOGEN.

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/kapila/1D/shockTubes/interfaceWaterAir/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/1D/shockTubes/epoxySpinel/</testCase>

interfaceWaterAir shock tube
----------------------------
A shock tube between a high-pressure chamber filled with water and a low-pressure chamber filled with air is released. Input files for this test are available in *./libTests/referenceTestCases/kapila/1D/shockTubes/interfaceWaterAir/*.

.. _Fig:testCases:Kapila:shockTubeWaterAirCI:

.. figure:: ./_static/testCases/Kapila/shockTubeWaterAir/schemaCI.png
  :scale: 70%
  :align: center

  Initial condition for 1D water--air shock tube.

The initial characteristics of the run are:

+------------------------------+---------------------------+
| Characteristic               | Value                     |
+==============================+===========================+
| dimension                    | 1 m                       |
+------------------------------+---------------------------+
| Initial mesh size / max size | 100 / 230                 |
+------------------------------+---------------------------+
| number of refinement level   | 4                         |
+------------------------------+---------------------------+
| diaphragm position           | 0.7 m                     |
+------------------------------+---------------------------+
| boundary conditions          | non-reflecting            |
+------------------------------+---------------------------+
| final solution time          | 0.240 ms                  |
+------------------------------+---------------------------+
| solution printing frequency  | 0.012 ms                  |
+------------------------------+---------------------------+
| precision                    | 2nd order (Vanleer/THINC) |
+------------------------------+---------------------------+

AMR technique of :cite:`schmidmayer2019adaptive` is used with 4 refinement levels such that a maximum of 230 computational cells are used for this run.

.. _Fig:testCases:Kapila:shockTubeWaterAir:

.. figure:: ./_static/testCases/Kapila/shockTubeWaterAir/shockTubeWaterAir.*
  :scale: 50%
  :align: center

  Water--air shock tube. Visualization using Paraview_ software.

epoxySpinel
-----------
This test deals with shocks in mixture of materials. Epoxy and spinel are supposed mixed such that caracteristic times for wave propagation and drag effects are very small, allowing to consider the mixture as evolving in mechanical equilibrium. Input files for this test are available in *./libTests/referenceTestCases/kapila/1D/shockTubes/epoxySpinel/*.

.. _Fig:testCases:Kapila:shockTubeEpoSpiCI:

.. figure:: ./_static/testCases/Kapila/shockTubeEpoSpi/schemaCI.png
  :scale: 70%
  :align: center

  Initial condition for mixture shock tube with epoxy and spinel.

The initial characteristics of the run are:

+------------------------------+---------------------+
| Characteristic               | Value               |
+==============================+=====================+
| dimension                    | 1 m                 |
+------------------------------+---------------------+
| Initial mesh size / max size | 200 / 237           |
+------------------------------+---------------------+
| number of refinement level   | 2                   |
+------------------------------+---------------------+
| diaphragm position           | 0.6 m               |
+------------------------------+---------------------+
| boundary conditions          | non-reflecting      |
+------------------------------+---------------------+
| final solution time          | 0.1 ms              |
+------------------------------+---------------------+
| solution printing frequency  | 0.025 ms            |
+------------------------------+---------------------+
| precision                    | 2nd order (Vanleer) |
+------------------------------+---------------------+

.. _Fig:testCases:Kapila:shockTubeEpoSpi:

.. figure:: ./_static/testCases/Kapila/shockTubeEpoSpi/shockTubeEpoSpi.*
  :scale: 50%
  :align: center

  Mixture shock tube with expoxy and spinel. Visualization using Paraview_ software.

Other tests cases
=================
Other tests are provided with ECOGEN package. They will soon be described in details.

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/kapila/1D/cavitation/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/2D/transportWaterSquareInAir/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/2D/squareWaterExplosion/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/2D/shockBubble/heliumAir/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/2D/richtmyerMeshkov/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/2D/testUnstructured/</testCase>
  <testCase>./libTests/referenceTestCases/Kapila/AddPhysicalEffects/evap/evapShockTube/</testCase>
  <testCase>./libTests/referenceTestCases/Kapila/AddPhysicalEffects/evap/dodEvapShockTube/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/AddPhysicalEffects/surfaceTension/squareToCircle/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/AddPhysicalEffects/surfaceTension/squareToCircleSymmetry/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/AddPhysicalEffects/surfaceTension/waterCylinderInAir/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/AddPhysicalEffects/surfaceTension/waterDropletInAir/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/AddPhysicalEffects/surfaceTension/dropletImpact/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/AddPhysicalEffects/gravity/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/3D/unstructured/</testCase>
  <testCase>./libTests/referenceTestCases/kapila/3D/shockBubble/heliumAir/</testCase>


.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/