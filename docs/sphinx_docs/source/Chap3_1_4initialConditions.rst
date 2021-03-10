.. role:: xml(code)
	:language: xml

.. _Sec:input:InitialConditions:

InitialConditionsV4.xml
=======================

 The *initialConditionsV4.xml* input file includes the initial conditions and the boundary conditions of the flow simulation. It is **mandatory** located in the folder of the current case. The typical structure of this file is:

 .. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<CI>
	  <!-- LIST OF GEOMETRICAL DOMAINS  -->
	  <physicalDomains> 
	    <domain name="base" state="leftChamber" type="entireDomain"/>
	  </physicalDomains>

	  <!-- LIST OF BOUNDARY CONDITIONS -->
	  <boundaryConditions>
	    <boundCond name="CLXm" type="nonReflecting" number="1"/>
	    <boundCond name="CLXp" type="nonReflecting" number="2"/>
	  </boundaryConditions>

	  <!--  LIST OF STATES  -->
	  <state name="leftChamber">
	    <material type="fluid" EOS="IG_air.xml">
	      <dataFluid alpha="0.5" density="1.0"/>  
	    </material>
	    <material type="fluid" EOS="SG_water.xml">
	      <dataFluid alpha="0.5" density="1000.0"/>
	    </material>
	    <mixture>
	      <dataMix pressure = "1.e5"/>
	      <velocity x="0." y="0." z="0."/>
	    </mixture>
	    <transport name="color" value="32."/>
	  </state>
	</CI>

.. _Sec:input:physicalDomains:

Physical domains
----------------
The :xml:`<physicalDomains>` markup is mandatory. Different initial conditions can be specified for different regions of the computational domain. This markup must contain as many nodes :xml:`<domain>` as necessary to correctly initialize the computational domain (overlaps are possible). A :xml:`<domain>` node contains the following attributes:
	
- :xml:`name`: A name for the domain. This name has no influence on the choices remaining in this file.
- :xml:`state`: This is the state of the fluid which will be specified with the :xml:`<state>` markup further in the file.
- :xml:`type`: Specify the type of geometric domain to which the state of the fluid will be assigned; :ref:`Sec:input:EntireDomain`, :ref:`Sec:input:HalfSpace`, :ref:`Sec:input:Disc`, :ref:`Sec:input:Rectangle`, :ref:`Sec:input:Ellipse`, :ref:`Sec:input:Cuboid`, :ref:`Sec:input:Sphere`, :ref:`Sec:input:Ellipsoid`, :ref:`Sec:input:Cylinder`.

**Important remark:** 

The initial conditions are attributed on each domain by using an overlapping principle. The order is therefore important: in the case of overlapping, the last attributed data are considered in the flow computation. Hence, it is important to attribute at least the entire domain at the first place thanks to the value *entireDomain*.

Depending on the geometrical shape, additional information is required through the use of the following nodes.

.. _Sec:input:EntireDomain:

EntireDomain
~~~~~~~~~~~~
Set the initial condition on the entire domain. No more information required.

.. code-block:: xml

	<domain name="base" state="leftChamber" type="entireDomain"/>

.. _Sec:input:HalfSpace:

HalfSpace
~~~~~~~~~
Set the initial condition on a half domain. The node :xml:`<dataHalfSpace>` must be included with the following attributes:

- :xml:`axis`: Can take the value *x*, *y* or *z*.
- :xml:`origin`: Real number indicating the location of the edge between the two subdomains on the specified axis.
- :xml:`direction`: Can take a *positive* or *negative* value on the specified axis.

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="halfSpace">
	  <dataHalfSpace axis="x" origin="0.5" direction="positive"/>
	</domain>

.. _Sec:input:Disc:

Disc
~~~~
In 2D, a disc is defined on a plane, while in 3D, a cylinder with an infinite length is defined. The node :xml:`<dataDisc>` must be added with the following attributes:

- :xml:`axis1` and :xml:`axis2`: The name of the 2 axes to define the plane on which the disc is defined. Can take two different values among *x*, *y* and *z*.
- :xml:`radius`: Real number indicating the disc radius (unit: m (SI)).
- Node :xml:`<center>`: Require the attributes :xml:`x`, :xml:`y` and :xml:`z`, as real numbers (unit: m (SI)), giving the location of the center of the disc in the plane (axis1, axis2).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="disc">
	  <dataDisc axis1="x" axis2="y" radius="0.5">
	    <center x="0." y="0." z="0."/>
	  </dataDisc>
	</domain>

.. _Sec:input:Rectangle:

Rectangle
~~~~~~~~~
In 2D, a rectangle is defined on a plane, while in 3D, a rectangular beam with an infinite length is defined. The node :xml:`< dataRectangle >` must be added with the following attributes:

- :xml:`axis1` and :xml:`axis2`: The name of the 2 axes to define the plane on which the rectangle is defined. Can take two different values among *x*, *y* and *z*.
- :xml:`lAxis1` and :xml:`lAxis2`: Real number indicating the length of both sides along (axis1, axis2).
- Node :xml:`<posInferiorVertex>`: Require the attributes :xml:`x`, :xml:`y` and :xml:`z`, as real numbers (unit: m (SI)), giving the location of the inferior corner in the plane (axis1, axis2).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="rectangle">
	  <dataRectangle axis1="x" axis2="y" lAxis1="0.3" lAxis2="0.2">
	    <posInferiorVertex x="0.4" y="0.5" z="0."/>
	  </dataRectangle>
	</domain>

.. _Sec:input:Ellipse:

Ellipse
~~~~~~~
In 2D, an ellipse is defined on a plane, while in 3D, an ellipsoid with an infinite length is defined. The node :xml:`<dataEllipse>` must be added with the following attributes:

- :xml:`axis1` and :xml:`axis2`: The name of the 2 axes to define the plane on which the ellipse is defined. Can take two different values among *x*, *y* and *z*.
- :xml:`radius1` and :xml:`radius2`: Real numbers indicating the ellipse radii (unit: m (SI)) along the corresponding axes.
- Node :xml:`<center>`: Require the attributes :xml:`x`, :xml:`y` and :xml:`z`, as real numbers (unit: m (SI)), giving the location of the center of the ellipse in the plane (axis1, axis2).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="ellipse">
	  <dataEllipse axis1="x" axis2="y" radius1="1." radius2="1.5">
	    <center x="0." y="0." z="0."/>
	  </dataEllipse>
	</domain>

.. _Sec:input:Cuboid:

Cuboid
~~~~~~
Set the initial condition of a cuboid. The additional node :xml:`<dataCuboid>` must be added with the attributes:

- :xml:`lAxisX`, :xml:`lAxisY` and :xml:`lAxisZ`: Real numbers for length of each side of the cuboid along axes (unit: m (SI)).
- Node :xml:`<posInferiorVertex>`: With the attributes :xml:`x`, :xml:`y` and :xml:`z`, real numbers corresponding to the location of the inferior corner (unit: m (SI)).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="cuboid">
	  <dataCuboid lAxisX="1." lAxisY="1." lAxisZ="0.5">
	    <posInferiorVertex x="1." y="0.5" z="0.5"/>
	  </dataCuboid>
	</domain>

.. _Sec:input:Sphere:

Sphere
~~~~~~
Set the initial condition of a sphere. The additional node :xml:`<dataSphere>` is required with the attributes:

- :xml:`radius`: Real number giving the radius of the sphere (unit: m (SI)).
- Node :xml:`<center>`: With the attributes :xml:`x`, :xml:`y` and :xml:`z`, real numbers giving the location on the center of the sphere (unit: m (SI)).

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="sphere">
	  <dataSphere radius="0.5">
	    <center x="1." y="0.5" z="0.5"/>
	  </dataSphere>
	</domain>

.. _Sec:input:Ellipsoid:

Ellipsoid
~~~~~~~~~
Set the initial condition of an ellipsoid. The additional node :xml:`<dataEllipsoid>` is required with the attributes:

- :xml:`axis1`, :xml:`axis2` and :xml:`axis3`: The name and therefore order of the 3 axes on which the ellipsoid is defined. Can take values among *x*, *y* and *z*.
- :xml:`radius1`, :xml:`radius2` and :xml:`radius3`: Real numbers indicating the ellipsoid radii (unit: m (SI)) along the corresponding axes.
- Node :xml:`<center>`: Require the attributes :xml:`x`, :xml:`y` and :xml:`z`, as real numbers (unit: m (SI)), giving the location of the center of the ellipsoid.

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="ellipsoid">
	  <dataEllipsoid axis1="x" axis2="y" axis3="z" radius1="1." radius2="1.5" radius3="1.5">
	    <center x="0." y="0." z="0."/>
	  </dataEllipsoid>
	</domain>

.. _Sec:input:Cylinder:

Cylinder
~~~~~~~~
Set the initial condition of a cylinder. The additional node :xml:`<dataCylinder>` is required with the attributes:

- :xml:`axis1` and :xml:`axis2`: The name of the 2 axes to define the plane on which the disc surface of the cylinder is defined. Can take two different values among *x*, *y* and *z*.
- :xml:`radius`: Real number indicating the disc-surface radius (unit: m (SI)).
- :xml:`length`: Real number indicating the length of the cylinder (unit: m (SI)).
- Node :xml:`<center>`: Require the attributes :xml:`x`, :xml:`y` and :xml:`z`, as real numbers (unit: m (SI)), giving the location of the center of the cylinder.

.. code-block:: xml

	<domain name="HP"  state="rightChamber" type="cylinder">
	  <dataCylinder axis1="x" axis2="y" radius="0.5" length="1.">
	    <center x="0." y="0." z="0."/>
	  </dataCylinder>
	</domain>

Initializing using physical identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
An additional possible feature for the geometric domain is to use :xml:`physicalIdentity` number coming from mesh software to initialize a geometrical domain.

Example:

.. code-block:: xml

	<domain name="base" state="leftChamber" type="entireDomain" physicalEntity="10"/>

In this example, the entire computation domain will be initialized accordingly to the :xml:`physicalIdentity` 10 from the mesh file.

.. _Sec:input:boundaryConditions:

Boundary conditions
-------------------
The :xml:`<boundaryConditions>` markup is mandatory. The boundary conditions are specified at the boundary of the computational domain. This markup must contain as many nodes :xml:`<boundCond>` as necessary to cover the entire boundary. Each :xml:`<boundCond>` node contains the following attributes:

- :xml:`name`: A name for the boundary condition. This name has no influence on the choices remaining in this file.
- :xml:`type`:  The type of boundary condition, to choose among :ref:`Sec:input:NonReflecting`, :ref:`Sec:input:Symmetry`, :ref:`Sec:input:Wall`, :ref:`Sec:input:Injection`, :ref:`Sec:input:Outflow` and :ref:`Sec:input:Tank`.
- :xml:`number`: Integer corresponding to the identifier of the boundary.

Depending on the :xml:`<type>`, additional information is required through the use of the following nodes.

.. _Sec:input:NonReflecting:

Non-reflecting
~~~~~~~~~~~~~~
The numerical treatment corresponds to an in- or out-going flow without any wave reflection. No more information required.

.. code-block:: xml

	<boundCond name="exit" type="nonReflecting" number="1" />

.. _Sec:input:Symmetry:

Symmetry
~~~~~~~~
The numerical treatment corresponds to a symmetry condition. No more information required.

.. code-block:: xml

	<boundCond name="symmetry" type="symmetry" number="2" />

.. _Sec:input:Wall:

Wall
~~~~
The numerical treatment corresponds to a wall boundary condition. No more information required.

.. code-block:: xml

	<boundCond name="wall" type="wall" number="3" />

.. _Sec:input:Injection:

Injection
~~~~~~~~~
The numerical treatment corresponds to the link between the boundary with a inflow. The inflow is characterized by an incoming mass-flow rate at a given thermodynamical state. :xml:`injection` requires the :xml:`<dataInjection>` node with the following attributes:

- :xml:`m0`: Incoming mass-flow rate, real number (unit: kg/s.m-2 (SI)).
- Node :xml:`<dataFluid>` for each phase: It must contain as many nodes :xml:`<dataFluid>` as the number of phases in the flow simulation and each contains the attributes:

	- :xml:`EOS`: The name of the file corresponding to the choice of the EOS for the phase in the tank. This file must correspond to the one specified in *modelV4.xml* input file for every fluid.
	- :xml:`density`: The density of the fluid incoming, real number (unit: kg/m3 (SI)).
	- :xml:`pressure`: The pressure of the fluid incoming, real number (unit: Pa (SI)).
	- :xml:`alpha`: The volume fraction of the fluid incoming, real number in the range ]0.,1.[.

.. code-block:: xml

	<boundCond name="entrance" type="injection" number="1">
	  <dataInjection m0="2.413092"/>
	  <dataFluid EOS="IG_air.xml" density="4.04e-3" pressure="2.57404e2" alpha="0.1"/>
	  <dataFluid EOS="SG_water.xml" density="1000." pressure="2.57404e2" alpha="0.9"/>
	</boundCond>

.. _Sec:input:Outflow:

Outflow
~~~~~~~
In the case of a subsonic flow, the pressure is set equal to the ambient (distant) pressure at the boundary. The additional :xml:`<dataOutflow>` node is required with the attributes:

- :xml:`p0`: Outside pressure, real number (unit: Pa (SI)).
- Node :xml:`<transport>`: This node is also required for each transport equation used.

.. code-block:: xml

	<boundCond name="exit" type="outflow" number="5">
	  <dataOutflow p0="1.e5">
	    <transport name="color" value="1.e-6"/>
	  </dataOutflow>
	</boundCond>

.. _Sec:input:Tank:

Tank
~~~~
The numerical treatment corresponds to the link between the boundary with an infinite tank. An infinite tank is characterized by a null velocity while pressure and temperature are constant. :xml:`tank` requires the :xml:`<dataTank>` node with the following attributes:

- :xml:`p0`: Stagnation pressure, real number (unit: Pa (SI)).
- :xml:`T0`: Stagnation temperature, real number (unit: K (SI)).
- Node :xml:`<fluidsProp>`: Necessary to define the presence of each phase in the tank. It must contain as many nodes :xml:`<dataFluid>` as the number of phases in the flow simulation and each contains the attributes:

	- :xml:`EOS`: The name of the file corresponding to the choice of the EOS for the phase in the tank. This file must correspond to the one specified in *modelV4.xml* input file for every fluid.
	- :xml:`alpha`: The volume fraction of the fluid in the tank, real number in the range ]0.,1.[.

.. code-block:: xml

	<boundCond name="entrance" type="tank" number="3">
	  <dataTank p0="4.e6" T0="93.3"/>
	  <fluidsProp>
	    <dataFluid EOS="IG_oxyVap.xml" alpha="0.0001"/>
	    <dataFluid EOS="SG_oxyLiq.xml" alpha="0.9999"/>
	  </fluidsProp>
	</boundCond>

**Important remark**

The choice of the boundary-condition number is made according to the type of mesh given in *meshV5.xml* input file and it follows the rules:

- Cartesian: The boundaries are ordered and labeled from 1 to 6 (in 3D) according to:

	1. boundary condition at the minimal x location,
	2. boundary condition at the maximal x location,
	3. boundary condition at the minimal y location,
	4. boundary condition at the maximal y location,
	5. boundary condition at the minimal z location,
	6. boundary condition at the maximal z location.

- UnStructured: When an unstructured mesh is used, the number of the boundary condition must correspond to the number specified in the mesh file .geo (see example in section :ref:`Sec:tuto:generatingMeshes`).
 
**Remark**

The boundary conditions are dependent on the flow model specified in *modelV4.xml* input file. Some boundary conditions may be not available for the flow model considered.


Mechanical and thermodynamical states of the fluid
--------------------------------------------------
For each physical domain in the :xml:`<physicalDomains>` markup, a fluid state must correspond. It implies an additional :xml:`<state>` markup for each state of fluid. This :xml:`<state>` markup contains:

- As many :xml:`<material>` nodes as the number of phases involved in the simulation. 
- A :xml:`<mixture>` node is required if a multiphase model is used.

Each :xml:`<material>` node corresponds to a phase and contains the following attributes or nodes:

- Attribute :xml:`type`: Only the value *fluid* is available in the current ECOGEN version.
- Attribute :xml:`EOS`: The name of the file corresponding to the fluid equation-of-state parameters. This file must correspond to the one specified in *modelV4.xml* input file for each phase (see section :ref:`Sec:input:FlowModel`).
- Node :xml:`<dataFluid>`: Contain data related to the considered state of the fluid in the current phase. 

This last node :xml:`<dataFluid>` as well as the :xml:`<mixture>` node are dependent on the flow model according to:

.. _Sec:input:Euler:

Euler
~~~~~
Single phase flow. In this case, the :xml:`<mixture>` node is absent and the :xml:`<dataFluid>` node contains the following attributes or nodes:

- Attribute :xml:`temperature`: Initial temperature of the fluid, real number (unit: K (SI)).
- Attribute :xml:`pressure`: Initial pressure of the fluid, real number (unit: Pa(SI)).
- Node :xml:`<velocity>`: With :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector, real numbers (unit: m/s (SI)).

.. code-block:: xml

	<material type="fluid" EOS="IG_air.xml">
	  <dataFluid density="10.0" pressure="1.e5">
	    <velocity x="1000." y="1000." z="0."/>
	  </dataFluid>
	</material>

.. _Sec:input:UEq:

UEq (previously named MultiP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Multiphase flow at velocity equilibrium (same velocity for every phase within each cell). Each :xml:`<dataFluid>` node corresponds to a phase with the following attributes:

- :xml:`alpha`: Volume fraction of the phase, real number in the range ]0.,1.[. The range can be increased to [0.;1.] if the option *alphaNull* is turned on (*true*) in the *modelV4.xml* input file.
- :xml:`density` or :xml:`temperature`: Initial density or temperature of the fluid, real number (unit: kg/m3 or K (SI)), respectively.
- :xml:`pressure`: Initial pressure of the fluid, real number (unit: Pa (SI)).

Moreover, in this case, the :xml:`<mixture>` node contains:

- Node :xml:`<velocity>`: With :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector, real numbers (unit: m/s (SI)).

.. code-block:: xml

	<material type="fluid" EOS="IG_air.xml">
	  <dataFluid alpha="1." density="50." pressure="1.e5"/>
	</material>
	<material type="fluid" EOS="SG_water.xml">
	  <dataFluid alpha="0." density="1000.0" pressure="1.e5"/>
	</material>
	<mixture>
	  <velocity x="0." y="0." z="0."/>
	</mixture>

.. _Sec:input:PUEq:

PUEq (previously named Kapila)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Multiphase flow at pressure and velocity equilibrium (same velocity and pressure for every phase within each cell). Each :xml:`<dataFluid>` node corresponds to a phase with the following attributes:

- :xml:`alpha`: Volume fraction of the phase, real number in the range ]0.,1.[. The range can be increased to [0.;1.] if the option *alphaNull* is turned on (*true*) in the *modelV4.xml* input file.
- :xml:`density` or :xml:`temperature`: Initial density or temperature of the fluid, real number (unit: kg/m3 or K (SI)), respectively.

Moreover, in this case, the :xml:`<mixture>` node contains:

- Node :xml:`<dataMix>`: With :xml:`pressure` attribute for initial pressure of the fluid, real number (unit: Pa (SI)).
- Node :xml:`<velocity>`: With :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector, real numbers (unit: m/s (SI)).

.. code-block:: xml

	<material type="fluid" EOS="IG_air.xml">
	  <dataFluid alpha="0.5" density="1.0"/>  
	</material>
	<material type="fluid" EOS="SG_water.xml">
	  <dataFluid alpha="0.5" density="1000.0"/>   
	</material>
	<mixture>
	  <dataMix pressure = "1.e5"/>
	  <velocity x="0." y="0." z="0."/>
	</mixture>

.. _Sec:input:PTUEq:

PTUEq (previously named ThermalEq)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Multiphase flow at pressure, velocity and thermal equilibrium (same velocity, pressure and temperature for every phase within each cell). In this case, every :xml:`<dataFluid>` node corresponds to a phase with only one attribute :xml:`alpha` setting the volume-fraction real number in the range ]0.,1.[.

The :xml:`<mixture>` node contains the following attributes and nodes:

- Node :xml:`<dataMix>`: With :xml:`temperature` and :xml:`pressure` attributes for initial temperature and pressure of the fluid, real numbers (unit: K and Pa (SI)), respectively.
- Node :xml:`<velocity>`: With :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector of the mixture, real numbers (unit: m/s (SI)).

.. code-block:: xml

	<material type="fluid" EOS="IG_waterVap.xml">
	  <dataFluid alpha="0.2"/>
	</material>
	<material type="fluid" EOS="SG_waterLiq.xml">
	  <dataFluid alpha="0.8"/>
	</material>
	<mixture>
	  <dataMix pressure = "1.e5" temperature ="300."/>
	  <velocity x="0." y="0." z="0."/>
	</mixture>

.. _Sec:input:EulerHomogeneous:

EulerHomogeneous
~~~~~~~~~~~~~~~~
Multiphase flow at mechanical and thermodynamical equilibrium. In this case, every :xml:`<dataFluid>` node corresponds to a phase with only one attribute :xml:`alpha` setting the volume-fraction real number in the range ]0.,1.[.

Moreover, in this case, the :xml:`<mixture>` contains the following attributes and nodes:

- Node :xml:`<dataMix>`: With :xml:`pressure` attribute for initial pressure of the mixture, real number (unit: Pa (SI)).
- Node :xml:`<velocity>`: With :xml:`x`, :xml:`y` and :xml:`z` attributes setting the initial values for the components of the velocity vector of the mixture, real numbers (unit: m/s (SI)).

.. code-block:: xml

	<material type="fluid" EOS="SG_waterLiq.xml">
	  <dataFluid alpha="0.99"/>  
	</material>
	<material type="fluid" EOS="IG_waterVap.xml">
	  <dataFluid alpha="0.01"/>   
	</material>
	<mixture>
	  <dataMix pressure = "1.e6"/>
	  <velocity x="0." y="0." z="0."/>
	</mixture>

**Remark**

Be careful to set the volume fraction in the range ]0,.1.[ (unless you use the option *alphaNull=true* for pressure-velocity- and velocity-equilibrium models) as well as the its sum over all the phases equal to 1.
