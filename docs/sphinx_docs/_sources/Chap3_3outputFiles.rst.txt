.. _Sec:IO:output:

Output Files
============

Writing the result files is done according to the user’s choice in *mainV5.xml* input file of the current test (see section :ref:`Sec:input:main`). For each run, the results are recorded in the specified folder **ECOGEN/results/XXX/** where XXX is the test case *name*.
One can select the following format: 

- *GNU*: Format in ASCII, results are given in column.
- *XML*: VTK file format (ASCII or binary).

The result files are placed within the **ECOGEN/results/XXX/datasets/** subfolder and their names follow the rules:

result(format)_CPU(proc)_TIME(instant).(ext)

where the variables are:

- (format): Data format (empty for ASCII, B64 for binary).
- (instant): Timestamp of the written result (depends on the selected frequency for writing).
- (proc): The core number where the results are from (in the case of a parallel simulation).
- (ext): Type of mesh.

Using GNU format
----------------
This format leads to result files with the *(ext)=out* extension. This format in column is very useful for a quick visualization of the results using the freeware `gnuplot`_ or any other tool. When this format is selected, ECOGEN automatically creates a script file *visualization.gnu* at the root of the result folder. This allows a very quick and efficient use for 1D runs.

XML VTK file format
-------------------
ECOGEN can provide output using XML files in VTK file format (ASCII or BINARY). Writing files in VTK file format leads to files with the extension:

- *(ext)=vtr*: Cartesian mesh.
- *(ext)=vtu*: Unstructured mesh.

ECOGEN also produces files named *collectionParaview.pvd* and *collectionVisIt.visit* (or *collectionParaviewB64.pvd* and *collectionParaviewB64.pvd* for BINARY) in order to load only one pack of files when the softwares `PARAVIEW`_ and `VisIt`_ are used, respectively.

Saving input files
------------------
In addition to the results, a copy of the input files of the current simulation is done in the subfolder **ECOGEN/results/XXX/savesInput/** to ensure a safe reproduction of the results.

Probes and cuts
---------------
All the output files linked to the probes and cuts of the simulation are placed within the subfolders **ECOGEN/results/XXX/probes/** and **ECOGEN/results/XXX/cuts/**, respectively.

Mesh information
----------------
When using AMR, mesh information is saved at a frequency specified in the *mainV5.xml* file to allow a restart of the simulation and the corresponding files are placed within the subfolders **ECOGEN/results/XXX/infoMesh/**.

Screen output
-------------
In real time, some data are available during the simulation directly at the terminal screen as soon as ECOGEN starts. Hereafter is an example of a screenshot:
 
.. figure:: ./_static/tutos/default/RunECOGEN_Logo.png
  :scale: 100%
  :align: center

  : Screenshot of the top of ECOGEN's default run console.

One reads:

- number of the result files,
- number of the last timestep,
- physical time in the simulation,
- value of the last timestep,
- CPU time since the beginning of the simulation,
- CPU time spent in AMR routines since the beginning of the simulation.

.. _gnuplot: http://www.gnuplot.info/
.. _PARAVIEW: https://www.paraview.org/
.. _VisIt: https://wci.llnl.gov/simulation/computer-codes/visit/