.. role:: xml(code)
	:language: xml

.. _Sec:input:InitialConditions:

InitialConditionsV4.xml
=======================

 The *initialConditionsV4.xml* input file includes the initial conditions and the boundary conditions of the flow simulation. It is *mandatory* located in the folder of the current case. The typical structure of this file is:

 .. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<CI>
	  <!-- LIST OF GEOMETRICAL DOMAINS  -->
	  <physicalDomains> 
	    <domain name="base" state="leftChamber" type="entireDomain"/>
	  </physicalDomains>

	  <!-- LIST OF BOUNDARY CONDITIONS -->
	  <boundaryConditions>
	    <boundCond name="CLXm" type="abs" number="1"/>
	    <boundCond name="CLXp" type="abs" number="2"/>
	  </boundaryConditions>

	  <!--  LIST OF STATES  -->
	  <state name="leftChamber">
	    <material type="fluide" EOS="IG_air.xml">
	      <dataFluid alpha="0.5" density="1.0"/>  
	    </material>
	    <material type="fluide" EOS="SG_water.xml">
	      <dataFluid alpha="0.5" density="1000.0"/>
	    </material>
	    <mixture>
	      <dataMix pressure = "1.e5"/>
	      <velocity x="0." y="0." z="0."/>
	    </mixture>
	    <transport name="color" value="32."/>
	  </state>
	</CI>

Physical domains
----------------
The :xml:`<physicalDomains>` markup is mandatory. Some different initial conditions can be specified at different zones of the computational domain. This markup must contain as many nodes :xml:`<domain>` as necessary to correctly initialize the computational domain (overlaps are possible). A :xml:`<domain>` node contains the following attributes:
	
- :xml:`name`: a name for the domain. This name has no influence on the choices remaining in this file.
- :xml:`state`: This is the state of the fluid which will be specified with the :xml:`<state>` markup further in the file.
- :xml:`type`: to specify the kind of geometrical domain on which the state of the fluid must be attribute: :ref:`Sec:input:entireDomain`, :ref:`Sec:input:halfSpace`, :ref:`Sec:input:disc`, :ref:`Sec:input:rectangle`, :ref:`Sec:input:pavement`, :ref:`Sec:input:sphere`.

**Important remark:** 

The initial conditions are attributed on each domain by using a superposition principle. The order is important: in the case of overlapping, the last attributed data are considered in the flow computation. Hence, it is important to attribute at least the entire domain at the first-place thanks to the value entireDomain. 

According to the geometrical shape, additional information is required thanks to the use of the nodes among the list:

.. _Sec:input:entireDomain:

entireDomain
~~~~~~~~~~~~
Set the initial condition on the entire domain. No more information required.

.. code-block:: xml

	<domain name="base" state="leftChamber" type="entireDomain"/>

.. _Sec:input:halfSpace:

halfSpace
~~~~~~~~~
Set the initial condition on a half-domain. The node :xml:`<dataHalfSpace>` must be included with the following attributes:

- :xml:`axe`: can take the value *x*, *y* or *z*.
- :xml:`origin`: real number, indicates the location of the edge between the two subdomains on the specified axis.
- :xml:`direction`: can take the value positive or negative on the specified axis.

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="halfSpace">
	  <dataHalfSpace axe="x" origin="0.5" direction="positive"/>
	</domain>

.. _Sec:input:disc:

disc
~~~~
In 2D allows to define a disc on a plane, in 3D a cylinder with an infinite length is defined.  The node :xml:`<dataDisc>` must be added with the following attributes:

- Attributes :xml:`axe1`, :xml:`axe2`: The name of 2 axes to define the plane on which the disc is defined. Can take two different values among *x*, *y* or *z*.
- Attribute :xml:`radius`: Real number of the radius disc (unit: m (SI)).
- Node :xml:`<center>`: requires the attributes :xml:`x`, :xml:`y` et :xml:`z` giving the location of the center of the disc in the plan (axe1, axe2) in real numbers (unit: m (SI)).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="disc">
	  <dataDisc axe1="x" axe2="y" radius="0.5">
	    <center x="0." y="0." z="0."/>
	  </dataDisc>
	</domain>

.. _Sec:input:Rectangle:

Rectangle
~~~~~~~~~
In 2D allows to define a rectangle on a plane, in 3D a rectangular beam with an infinite length is defined. The node :xml:`< dataRectangle >` must be added with the following attributes:

- Attributes :xml:`axe1`, :xml:`axe2`: The name of 2 axes to define the plane on which the disc is defined.  Can take two different values among *x*, *y* or *z*.
- Attributes :xml:`lAxe1`, :xml:`lAxe2`: Length of both sides along (axe1,axe2).
- Node :xml:`<posInferiorVertex>`: equipped with the attributes :xml:`x`, :xml:`y` and :xml:`z`, real numbers giving the location of the inferior corner in the plane (axe1, axe2).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="rectangle">
	  <dataRectangle axe1="x" axe2="y" lAxe1="0.3" lAxe2="0.2">
	    <posInferiorVertex x="0.4" y="0.5" z="0."/>
	  </dataRectangle>
	</domain>

.. _Sec:input:Pavement:

Pavement
~~~~~~~~
Set the initial condition in a pavement. The additional node :xml:`<dataPavement>` must be added with the attributes:

- Attributes :xml:`lAxeX`, :xml:`lAxeY`, :xml:`lAxeZ`: Real numbers for length of each side of the pavement along axes (unit: m (SI)).
- Node :xml:`<posInferiorVertex>`: with the des attributes :xml:`x`, :xml:`y` and :xml:`z`, real numbers corresponding to the location of the inferior corner (unit: m (SI)).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="pavement">
	  <dataPavement lAxeX="1." lAxeY="1." lAxeZ="0.5">
	    <posInferiorVertex x="1." y="0.5" z="0.5"/>
	  </dataPavement>
	</domain>

.. _Sec:input:sphere:

sphere
~~~~~~
Set the initial condition in a sphere. The additional node :xml:`<dataSphere>` is required with the attributes or nodes:

- Attribute radius: real number giving the radius of the sphere (unit: m (SI)).
- Node :xml:`<center>`: with the attributes :xml:`x`, :xml:`y` et :xml:`z` real numbers giving the ocation on the ceter of the sphere (unit: m (SI)).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="sphere">
	  <dataSphere radius="0.5">
	    <center x="1." y="0.5" z="0.5"/>
	  </dataSphere>
	</domain>
 
Initializing using physical Identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Additional feature for geometrical domain: It is possible to use physicalIdentity number coming from mesh software to initialize a geometrical domain.

Example:

.. code-block:: xml

	<domain name="base" state="leftChamber" type="entireDomain" physicalEntity="10"/>

In this example the entire computation domain will be initialize accordingly to the physicalEntity 10.

Boundary conditions
-------------------
The :xml:`<boundaryConditions>` markup is mandatory. The boundary conditions are specified at the boundary of the computational domain. This markup must contain as many nodes :xml:`<boundCond>` as necessary to recover the entire boundary. Each :xml:`<boundCond>` node contains the following attributes:

- name: a name for the boundary condition. This name has no influence on the choices remaining in this file.
- type:  the kind of boundary condition, to choose among :ref:`Sec:input:abs`, :ref:`Sec:input:wall`, :ref:`Sec:input:tank`, :ref:`Sec:input:outflow` 
- numero: integer number that correspond to the number of the boundary.

According :xml:`<type>`, additional information is required thanks the use of the nodes among the list:  

.. _Sec:input:abs:

abs
~~~
The numerical treatment corresponds to an outgoing flow with no wave reflection. No more information required.

.. code-block:: xml

	<boundCond name="exit" type="abs" number="1" />

.. _Sec:input:wall:

wall
~~~~
The numerical treatment corresponds to a wall boundary condition. No more information required.

.. code-block:: xml

	<boundCond name="wall" type="wall" number="3" />

.. _Sec:input:tank:

tank
~~~~
The numerical treatment corresponds to the link between the boundary with an infinite tank. An infinite tank is characterized by a null velocity while pressure and temperature are constant). :xml:`<tank>` requires the :xml:`<dataTank>` node with the following attributes:

- :xml:`p0`: Stagnation pressure, real number (unit: Pa(SI)).
- :xml:`T0`: Stagnation pressure, real number (unit: K (SI)).
- An additional :xml:`<fluidsProp>` node is necessary to define the presence of each phase in the tank. It must contain as many nodes :xml:`<dataFluid>` as the number of phases in the flow simulation and contains the attributes:
	- :xml:`EOS`: the name of the file corresponding to the choice of the EOS for the phase in the tank. This file must correspond to the one specified in modelV4.xml input file for every fluid.
	- :xml:`alpha`: The volume fraction of the fluid in the tank, real number in the range ]0.,1.[

.. code-block:: xml

	<boundCond name="entrance" type="tank" number="3">
	  <dataTank p0="4.e6" T0="93.3"/>
	  <fluidsProp>
	    <dataFluid EOS="IG_oxyVap.xml" alpha="0.0001"/>
	    <dataFluid EOS="SG_oxyLiq.xml" alpha="0.9999"/>
	  </fluidsProp>
	</boundCond>

.. _Sec:input:outflow:

outflow
~~~~~~~
In the case of a subsonic flow, the pressure is set equal to the ambient pressure at the boundary. The additional :xml:`<dataOutflow>` node is required with the attributes:

- p0: outside pressure, real number (unit: Pa(SI)).

.. code-block:: xml

	<boundCond name="exit" type="outflow" number="5">
	  <dataOutflow p0="1.e5">
	    <transport name="color" value="1.e-6"/>
	  </dataOutflow>
	</boundCond>

**Important remark**

The choice of the boundary condition number is made according to the kind of mesh given in meshV5.xml input file according to the following rules:

- cartesian: The boundaries are ordered and labeled from 1 to 6 (in 3D) according to:
	1. boundary condition at the minimal x location
	2. boundary condition at the maximal x location
	3. boundary condition at the minimal y location
	4. boundary condition at the maximal y location
	5. boundary condition at the minimal z location
	6. boundary condition at the maximal z location
- unStructured: When an unstructured mesh is used, the number of the boundary condition must correspond to the number specified in the mesh file .geo (see example in section :ref:`Sec:tuto:generatingMeshes`).
 
**Remark**

The boundary conditions are dependent on the flow model specified in modelV4.xml input file. Some boundary conditions may be not available for the flow model considered.


Mechanical and thermodynamical states of the fluid
--------------------------------------------------
For each physical domain in the :xml:`<physicalDomains>` markup, a fluid state must correspond. It implies an additional :xml:`<state>` markup for each state of fluid. This :xml:`<state>` markup contains:

- as many :xml:`<material>` nodes as the number of phases involved in the simulation. 
- A :xml:`<mixture>` node is required if a multiphase model is used.

Each :xml:`<material>` node corresponds to a phase and contains the following attributes or nodes:

- Attribute :xml:`type`: Only the value *fluide* is available in the current ECOGEN version.
- Attribute :xml:`EOS`: the name of the file corresponding to the fluid Equation of State parameters. This file must correspond to the one specified in *modelV4.xml* input file for each phase (see section :ref:`Sec:input:FlowModel`).
- Node :xml:`<dataFluid>`: contains data related to the considered state of the fluid in the current phase. 

This last node :xml:`<dataFluid>` as well as the :xml:`<mixture>` node are dependent on the flow model according to:

.. _Sec:input:euler:

euler
~~~~~
Single phase flow. In this case, the :xml:`<mixture>` node is absent and the :xml:`<dataFluid>` node contains the following attributes or nodes:

- Attribute :xml:`temperature`: Initial temperature of the fluid, real number (unit: K (SI)).
- Attribute :xml:`pressure`: Initial pressure of the fluid, real number (unit: Pa(SI)).
- Node :xml:`<velocity>`: with :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector, real numbers (unit: m/s (SI)).


.. code-block:: xml

	<material type="fluide" EOS="IG_air.xml">
	  <dataFluid density="10.0" pressure="1.e5">
	    <velocity x="1000." y="1000." z="0."/>
	  </dataFluid>
	</material>

.. _Sec:input:Kapila:

Kapila
~~~~~~
Multiphase flow at pressure and velocity equilibrium (same velocity and same pressure for every phase). Each :xml:`<dataFluid>` node corresponds to a phase with the following attributes:

- :xml:`alpha`: Volume fraction of the phase, real number in the range ]0.,1.[.
- :xml:`density`: Initial specific mass of the fluid, real number (unit: kg/m3 (SI)) or :xml:`temperature`: Initial temperature (unit: K).

Moreover, in this case, the :xml:`<mixture>` node contains:

- the :xml:`<dataMix>` node with :xml:`pressure` attribute for initial pressure of the fluid, real number (unit: Pa(SI)).
- the :xml:`<velocity>` node with :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector, real number (unit: m/s (SI)).

.. code-block:: xml

	<material type="fluide" EOS="IG_air.xml">
	  <dataFluid alpha="0.5" density="1.0"/>  
	</material>
	<material type="fluide" EOS="SG_water.xml">
	  <dataFluid alpha="0.5" density="1000.0"/>   
	</material>
	<mixture>
	  <dataMix pressure = "1.e5"/>
	  <velocity x="0." y="0." z="0."/>
	</mixture>

.. _Sec:input:ThermalEq:

ThermalEq
~~~~~~~~~
Multiphase flow at pressure, velocity and thermal equilibrium (same velocity, same pressure and same temperature for every phase). In this case, every :xml:`<dataFluid>` node corresponds to a phase with only one attribute :xml:`alpha` setting the volume fraction real number in the range ]0.,1.[.

The :xml:`<mixture>` node contains the following attributes and nodes:

- the :xml:`<dataMix>` node with temperature and pressure attributes for initial temperature of the fluid, real number (unit: K (SI)) and initial pressure, real number (unit: Pa).
- Attribute :xml:`pressure`: Initial pressure of the mixture, real number (unit: Pa(SI)).
- Node :xml:`<velocity>`: with :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector of the mixture, real numbers (unit : m/s (SI)).

.. code-block:: xml

	<material type="fluide" EOS="IG_waterVap.xml">
	  <dataFluid alpha="0.2"/>
	</material>
	<material type="fluide" EOS="SG_waterLiq.xml">
	  <dataFluid alpha="0.8"/>
	</material>
	<mixture>
	  <dataMix pressure = "1.e5" temperature ="300."/>
	  <velocity x="0." y="0." z="0."/>
	</mixture>

.. _Sec:input:EulerHomogeneous:

EulerHomogeneous
~~~~~~~~~~~~~~~~
Multiphase flow at mechanical and thermodynamical equilibrium. In this case, every :xml:`<dataFluid>` node corresponds to a phase with only one attribute :xml:`alpha` setting the volume fraction real number in the range ]0.,1.[.

Moreover, in this case, the :xml:`<mixture>` contains the following attributes and nodes:

- the :xml:`<dataMix>` node with :xml:`pressure` attribute for initial pressure of the mixture, real number (unit: Pa(SI)).
- Node :xml:`<velocity>`: with :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector of the mixture, real numbers (unit: m/s (SI)).

.. code-block:: xml

	<material type="fluide" EOS="SG_waterLiq.xml">
	  <dataFluid alpha="0.99"/>  
	</material>
	<material type="fluide" EOS="IG_waterVap.xml">
	  <dataFluid alpha="0.01"/>   
	</material>
	<mixture>
	  <dataMix pressure = "1.e6"/>
	  <velocity x="0." y="0." z="0."/>
	</mixture>

**Remark**

Be careful to set the volume fraction in the range ]0,.1.[ as well as the sum over the phases equal to 1.
