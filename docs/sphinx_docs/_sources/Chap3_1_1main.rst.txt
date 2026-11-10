.. role:: xml(code)
	:language: xml

.. _Sec:input:main:

Main.xml
==========

*main.xml* is the main input file needed to define the test case. It **must be** in the current test case folder. Its minimal structure is:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<computationParam>
	  <run>example</run>
	  <outputMode format="VTK" binary="false"/>
	  <timeControlMode iterations="false">
	    <iterations number="3" iterFreq="1"/>
	    <physicalTime totalTime="8.e-3" timeFreq="8.e-4"/>
	  </timeControlMode>
	  <computationControl CFL="0.8"/>
	</computationParam>

It contains general parameters for the computation listed below. Some are mandatory, others are optional.

.. _Sec:input:main:runName:

Run name
--------

.. code-block:: xml

	<run>example</run>

The :xml:`<run>` markup is mandatory. It will be used to create the folder containing the test-case results in **ECOGEN/results/**.

Output format
-------------

.. code-block:: xml

	<outputMode format="VTK" binary="false" precision="10" reducedOutput="true"/>

The :xml:`<outputMode>` markup is mandatory. The user can choose the writing output format. Attributes are:

- :xml:`format`: Can take the value *GNU* (standard writing in column), *VTK* (`XML VTK format <https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html>`_) or *VTKHDF* (`HDF5 VTK format <https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html>`_, available only if ECOGEN has been built with :ref:`HDF5 support<enabling-hdf5>`)
- :xml:`binary`: Can take the value true or false. *Binary* (true) or *ASCII* (false) format can be chosen. *GNU* format only supports *false* value. *VTKHDF* format only supports *true* value.
- :xml:`precision`: Optional attribute. Floating-point precision in output files (number of digits). If not precised, uses the default stream precision (6). *VTKHDF* format only supports *8* (simple precision) or *16* (double precision) values and defaults to simple precision (8).
- :xml:`reducedOutput`: Optional attribute. Reduced number of output variables when possible (depends on the model). Only available for Euler and velocity/pressure-velocity equilibrium model. In the case of the Euler model, it toggles the output of the turbulent viscosity. Can take the value true or false. If not precised, default is false (complete output).
- :xml:`intSize`: Optional attribute. *VTKHDF* format only. Integer size to store mesh connectivity in output files (number of bits). *32* or *64* values can be chosen. If not precised, *32*-bit integers are saved.
- :xml:`transferMode`: Optional attribute. *VTKHDF* format only. Collective or independant transfer mode during data I/O operations. *H5FD_MPIO_COLLECTIVE* or *H5FD_MPIO_INDEPENDENT* values can be chosen. If not precised, *H5FD_MPIO_COLLECTIVE* mode is used.

Output result files will be placed in the folder **ECOGEN/results/** into a specific subfolder bearing the name of the run.

Time evolution control
----------------------

.. code-block:: xml

	  <timeControlMode iterations="false">
	    <iterations number="3" iterFreq="1"/>
	    <physicalTime totalTime="8.e-3" timeFreq="8.e-4"/>
	  </timeControlMode>

ECOGEN is a CFD tool based on an explicit integration scheme in time. The :xml:`<timeControlMode>` markup is mandatory and defines the temporal evolution of the current simulation. It contains the
:xml:`iterations` attribute that can take the two values:

-	*true*: The time control is done thanks to the total number of timesteps and the :xml:`<iterations>` node must be present.
-	*false*: The time control is done thanks to the final physical time and the :xml:`<physicalTime>` node must be present.

The :xml:`<iterations>` markup:

.. code-block:: xml

	<iterations number="3" iterFreq="1"/>

ECOGEN automatically computes the timestep value thanks to a numerical stability criterion (CFL  criterion). This markup is defined with following attributes:

-	:xml:`number`: Integer equal to the total number of temporal timesteps.
-	:xml:`iterFreq`: Integer equal to the writing frequency of the results (results are written every :xml:`iterFreq` timestep)

The :xml:`<physicalTime>` markup:

.. code-block:: xml

	<physicalTime totalTime="8.e-3" timeFreq="8.e-4"/>

If this markup is used, ECOGEN automatically determines the total amount of timestep to compute to reach the chosen physical time. Attributes are:

-	:xml:`totalTime`: Real number equal to the final physical time of the simulation (unit: s (SI)).
-	:xml:`timeFreq`: Real number equal to the writing frequency of the results (results are written every :xml:`timeFreq` seconds).

CFL criterion
-------------

.. code-block:: xml

	<computationControl CFL="0.8"/>

The :xml:`<computationControl>` markup is mandatory. It specifies the value of the attribute :xml:`CFL` which ensures the stability of the temporal integration scheme: The value (real number) must be less than 1.

Probes
------
It is possible to record over time flow variables at given locations in the computational domain. This is done by including to the *main.xml* input file the optional :xml:`<probe>` markup.

.. code-block:: xml

	<probe name="capteur1">
	  <vertex x="0.51" y="0.51" z="0.51"/>
	  <timeControl acqFreq="-1."/>       <!-- if negative or nul, recording at each time step -->
	</probe>

The two following nodes must be included in the :xml:`<probe>` markup:
-	:xml:`<vertex>`: Specify the location of the probe into the computational domain in each physical direction corresponding to the attributes :xml:`x`, :xml:`y` and :xml:`z` (unit: m (SI)). Values must be real numbers.
-	:xml:`<timeControl>`: Allow to specify the probe acquisition frequency (unit: s (SI)). If the value is set to zero or negative, flow values at the probe location are recorded at each time step.

Probe output-result files will be placed in the specific subfolder **ECOGEN/results/XXX/probes/** where *XXX* represents the name of the run specified in :xml:`<run>` markup.

**Remarks:**

1. Recording probes with a high frequency could have a significant impact on computation performances due to the computer memory time access. To prevent that, one should fix a reasonable acquisition frequency.
2. Several probes can be added simultaneously. For that, place as many as wanted :xml:`<probe>` markups in the *main.xml* input files.

Simulation resume option
-------------------------
ECOGEN can resume a run from previously written results. To do so, the :xml:`resumeFileNumber` attribute of the :xml:`<resumeSimulation>` markup must be specified. Furthermore, this markup can include an :xml:`AMRsaveFreq` attribute which gives the frequency at which the mesh information of an AMR simulation is saved. Therefore, an AMR simulation can only resume from a time coherent combination of result and mesh information files.

.. code-block:: xml

	<resumeSimulation resumeFileNumber="10" AMRsaveFreq="5"/>
