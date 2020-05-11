.. role:: xml(code)
	:language: xml

.. _Sec:input:main:

MainV5.xml
==========

*mainV5.xml* is the main input file needed to define the test case. It **must be** in the current test case folder. Its minimal structure is:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<computationParam>
	  <run>example</run>
	  <outputMode format="XML" binary="false" precision="10"/>
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

	<outputMode format="XML" binary="false" precision="10"/>

The :xml:`<outputMode>` markup is mandatory. The user can choose the writing output format. Attributes are:

- :xml:`format`: Can take the value *GNU* (standard writing in column) or *XML* (XML VTK format).
- :xml:`binary`: Can take the value true or false. *Binary* (true) or *ASCII* (false) format can be chosen.

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

Global accuracy order of the numerical scheme
---------------------------------------------
When it is possible, according to the mesh or to the flow model, ECOGEN can use a second-order scheme (based on MUSCL approach with a TVD slope limiter; see :cite:`schmidmayer2020ecogen` for details). In this case, the optional markup :xml:`<secondOrder>` can be inserted in the *mainV5.xml* input file as in the following example:

.. code-block:: xml

	<secondOrder>
	  <globalLimiter>minmod</globalLimiter>
	  <interfaceLimiter>superbee</interfaceLimiter>                            <!-- optionnal node -->
	  <globalVolumeFractionLimiter>minmod</globalVolumeFractionLimiter>        <!-- optionnal node -->
	  <interfaceVolumeFractionLimiter>thinc</interfaceVolumeFractionLimiter>   <!-- optionnal node -->
	</secondOrder>

The :xml:`<secondOrder>` markup must contain the node :xml:`<globalLimiter>`. The other nodes are optional. The slope limiters available in ECOGEN are the following: minmod :cite:`roe1986superbee`, vanleer :cite:`vanLeer1974`, vanalbada :cite:`vanAlbada1997`, mc :cite:`van1977towards`, superbee :cite:`roe1986superbee` and thinc :cite:`shyue2014thinc` (only for volume fraction at the interface location). Markup significations are:

- :xml:`<globalLimiter>`: Applied everywhere and on all variables unless it is overwritten by the following optional limiters.
- :xml:`<interfaceLimiter>`: Applied on all variables but only at the interface location. By default is equal to the global limiter.
- :xml:`<globalVolumeFractionLimiter>`: Applied everywhere but only on the volume-fraction and transport equations (THINC is only applied on the volume fraction) unless it is overwritten by the interface volume-fraction limiter. By default is equal to the global limiter.
- :xml:`<interfaceVolumeFractionLimiter>`: Applied only at the interface location and on the volume-fraction and transport equations (THINC is only applied on the volume fraction). By default is equal to the interface limiter.

Probes
------
It is possible to record over time flow variables at given locations in the computational domain. This is done by including to the *mainV5.xml* input file the optional :xml:`<probe>` markup.
 
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
2. Several probes can be added simultaneously. For that, place as many as wanted :xml:`<probe>` markups in the *mainV5.xml* input files.

Simulation restart option
-------------------------
ECOGEN can restart a run from previously written results. To do so, the :xml:`restartFileNumber` attribute of the :xml:`<restartSimulation>` markup must be specified. Furthermore, this markup can include an :xml:`AMRsaveFreq` attribute which gives the frequency at which the mesh information of an AMR simulation is saved. Therefore, an AMR simulation can only restart from a time coherent combination of result and mesh information files.

.. code-block:: xml

	<restartSimulation restartFileNumber="10" AMRsaveFreq="5"/>