Output Files
============

Writing the result files is done according the user’s choice in *mainV5.xml* INPUT FILE of the current test (see section :ref:`Sec:input:main`). For each run, the results are recorded in the specified folder in **ECOGEN/results/XXX/** where XX is the test case *name*.
One can select the following format: 

- *GNU*: Format in ASCII, results are given in column.
- *XML*: VTK file format (ASCII or binary).

The name of the results files follows the rules:

result(format)_CPU(proc)_AMR(niveau)_TIME(instant). (ext)

that can select results files according:

- (format)	: data format (empty for ASCII, B64 for binary).
- (instant)	: time of writing results (depends on the selected frequency for writing.
- (niveau) 	: AMR level (in the case of an AMR simulation).
- (proc) 	: the number of the processor where the results are from (in the case of a parallel simulation).
- (ext) 	: kind of mesh.

Using GNU format
----------------
This format lead to results files with the (ext)=out extension. This format in colon is very useful for a quick visualization of the results when the freeware gnuplot  or any other tool. When this format is selected, ECOGEN automatically creates a script file visualisation.gnu at the root of the result folder. This allows a very quick and efficient use for 1D runs.

XML VTK file format
-------------------
ECOGEN can provide output using XML files in VTK file format (ASCII or BINARY). Writing files according VTK file format leads to files with the extension:

- (ext)=vtr	: cartesian mesh.
- (ext)=vtu	: unstructured mesh.

ECOGEN also produces a file named *collection.pvd* in order to load only one pack of files when the software PARAVIEW is used.

Saving input files
------------------
In addition to the results, a copy of the INPUT FILES of the current simulation is done in the subfolder **ECOGEN/results/dossierResExemple/savesEntrees/** to ensure a safe reproduction of the results. 

Screen output
-------------
In real time, some data are available during the simulation directly at the terminal screen as soon as ECOGEN starts. hereafter an example of a screenshot:
 
One reads:

- number of the results files.
- number of the last timestep.
- Real time in the simulation.
- Value of the last timestep.
- CPU time since the beginning of the simulation.
