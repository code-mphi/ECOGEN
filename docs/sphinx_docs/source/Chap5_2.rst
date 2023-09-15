.. role:: xml(code)
  :language: xml

*********************
Multiphase test cases
*********************

Test cases presented in this section are dealing with multiphase compressible problems.

.. toctree::

   ./Chap5_2PUEq.rst
   ./Chap5_2UEq.rst
   ./Chap5_2UEqTotE.rst
   ./Chap5_2PTUEq.rst
   ./Chap5_2eulerHomogeneous.rst


Other tests cases
=================
Other tests are provided with ECOGEN package and may be described in details later.

.. code-block:: xml

  <testCase>libTests/referenceTestCases/PUEq/2D/foil/noCavitation/</testCase>
  <testCase>libTests/referenceTestCases/PUEq/2D/foil/cavitation/</testCase>
  <testCase>libTests/referenceTestCases/PUEq/AddPhysicalEffects/phaseChange/evapExpansionTubeEquilibrium/</testCase>
  <testCase>libTests/referenceTestCases/PUEq/AddPhysicalEffects/phaseChange/condensation/</testCase>
  <testCase>libTests/referenceTestCases/PUEq/AddPhysicalEffects/phaseChange/evaporation/</testCase>
  <testCase>./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/waterCylinderInAir/</testCase>
  <testCase>./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/waterDropletInAir/</testCase>
  <testCase>./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/waterDropletInAir_axisym/</testCase>
  <testCase>./libTests/referenceTestCases/PUEq/AddPhysicalEffects/surfaceTension/dropletImpact_restart/</testCase>
  <testCase>./libTests/referenceTestCases/PUEq/AddPhysicalEffects/gravity/</testCase>



.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/

