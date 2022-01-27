.. role:: xml(code)
  :language: xml

Cavitation 
==========

The phenomenon of cavitation consists in the disruption of continuity in the liquid where there is considerable local reduction of pressure. The formation of bubbles within liquids (cavitation) begins even in the presence of positive pressures that are equal to or close to the pressure of saturated vapor of the fluid at the given temperature.
The bubbles increase rapidly in size. Subsequently, when the bubbles enter a zone of higher pressure, they are reduced in size as a result of condensation of the vapors that they contain.
This process of condensation takes place fairly quickly, accompanied by local hydraulic shocks, the emission of sound, the destruction of material bonds and other undesirable phenomena.

In this test case, cavitation is only mimicked (no phase change added to the model). Indeed, the model takes into account the differences in the acoustic behavior of both phases or in other words, the differences in expansion and compression of each phase in mixture regions.
Initial conditions of this test are shown in :numref:`Fig:testCases:PUEq:cavitation`.
Input files for this test are available in *./libTests/referenceTestCases/PUEq/1D/cavitation/*. The corresponding uncommented line in *ECOGEN.xml* is:

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/PUEq/1D/cavitation/</testCase>


.. _Fig:testCases:PUEq:cavitation:

.. figure:: ./_static/testCases/PUEq/cavitation/initialCavitation.png
  :scale: 50%
  :align: center

  Initial condition before cavitation.

The initial characteristics of the run are:

+------------------------------+---------------------+
| Characteristic               | Value               |
+==============================+=====================+
| dimension                    | 1 m                 |
+------------------------------+---------------------+
| Initial mesh size            | 1000                |
+------------------------------+---------------------+
| boundary conditions          | non-reflecting      |
+------------------------------+---------------------+
| final solution time          | 1.85 ms             |
+------------------------------+---------------------+
| solution printing frequency  | 0.925 ms            |
+------------------------------+---------------------+
| precision                    | 2nd order (VanLeer) |
+------------------------------+---------------------+

Thanks to the evolution of the density (:numref:`Fig:testCases:PUEq:cavitationAnim`), one can observe that air in the mixture quickly expands.
One can also remark a peak in the center which is purely due to the initial discontinuous conditions. An initial smearing of the velocity discontinuity would prevent such an event.

.. _Fig:testCases:PUEq:cavitationAnim:

.. figure:: ./_static/testCases/PUEq/cavitation/cavitationAnim.*
  :scale: 100%
  :align: center

  Evolution of density. Visualization using Paraview_ software.

.. figure:: ./_static/testCases/PUEq/cavitation/cavitationColorbar.*
  :scale: 70%
  :align: center

.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/