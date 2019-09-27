.. _Chap:TestCases:

Test cases
==========

ECOGEN is provided including several test cases. Major part of these tests can be found in companion papers. The full list of available test cases in the package are listed in the conductor input file *ECOGEN.xml*. Description of some provided test cases is presented in this section. We will refer to the test case following the line format to uncomment in *ECOGEN.xml* input file.

.. admonition:: Example
  
  For example, the following line in *ECOGEN.xml* input file corresponds to the default test:

  .. code-block:: xml

    <testCase>./libTests/referenceTestCases/euler/1D/transport/positiveVelocity/</testCase>

  It is presented in details in the tutorial section :ref:`Sec:tuto:default` :

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Chap5_1Euler
   Chap5_2Kapila
   
.. Chap5_3thermalEq
.. Chap5_4eulerHomogene
