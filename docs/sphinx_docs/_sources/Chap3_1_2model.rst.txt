.. role:: xml(code)
	:language: xml

.. _Sec:input:model:

Model.xml
===========

Fluid mechanics models used in the computation are specified in the *model.xml* file. It is **mandatory** located in the folder of the current case. A typical form of this file is:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<model>
	  <flowModel name="PressureVelocityEq" numberPhases="2" alphaNull="false"/>
	  <EOS name="IG_air.xml"/>
	  <EOS name="SG_water.xml"/>
	</model>

.. _Sec:input:FlowModel:

Flow model
----------

.. code-block:: xml

	<flowModel name="PressureVelocityEq" numberPhases="2" numberTransports="1" alphaNull="false"/>

The :xml:`<flowModel>` markup is **mandatory** to specify the mathematical model to solve during the computation. This markup may contain the following attributes:

- :xml:`name`: Name of the mathematical flow model. This attribute can take the values: *Euler* :cite:`euler1757principes`, *PressureVelocityEq* (previously named *Kapila*) :cite:`relaxjcp`:cite:`kapila2001`, *VelocityEq* (previously named *multiP*) :cite:`schmidmayer2023UEq`, *VelocityEqTotEnergy*, *NonEqDEM* :cite:`absau2003`, *NonEq* :cite:`FURFARO2015159`, *TemperaturePressureVelocityEq* (previously named *ThermalEq*) :cite:`martelotboiling`, *EulerHomogeneous*, *ExtendedNavierStokesKorteweg1D* :cite:`DidierlaurentInReview`.
- :xml:`numberPhases`: Integer number corresponding to the number of phases present in the simulation. The total amount of equations is related to this number. This attribute is not necessary when the values of :xml:`name` are *Euler* or *EulerHomogeneous*.
- :xml:`numberTransports`: This attribute is optional and is set to 0 as default. It can be used to add specific variables advected in the flow (color functions).
- :xml:`alphaNull`: For *PressureVelocityEq*, *VelocityEq*, *VelocityEqTotEnergy* and *NonEqDEM* models, the volume fraction can either be null (*true*) or not (*false*) and this choice is determined by this parameter. Default value is *false*.

**Remark:**

If *EulerHomogeneous* model is chosen, two other attributes are to be used: :xml:`liquid` and :xml:`vapor` to specify the number corresponding to the liquid phase and the vapor phase. They are phase 0 for the first and 1 for the second.

.. code-block:: xml

	<flowModel name="EulerHomogeneous" liquid="0" vapor="1"/>
	<EOS name="SG_waterLiq.xml"/>
	<EOS name="IG_waterVap.xml"/>

If *ExtendedNavierStokesKorteweg1D* model is chosen, other attributes are needed, in accordance with the work of Didierlaurent et al. :cite:`DidierlaurentInReview`. From this generic model, one can choose different equations of state to obtain particular models lying in different fields: Defocusing Nonlinear Schrödinger equation :cite:`Dhaouadi2019`, dynamics of a bistable tape spring :cite:`Bourgeois2020`, or the dynamics of phase transition by considering the van der Waals EOS :cite:`DidierlaurentInReview`.

.. code-block:: xml

	<flowModel name="ExtendedNavierStokesKorteweg1D" lambda="84000" beta="2.15e-14" tau="1." eps="1."/>
	<EOS name="VDW_dodecane.xml"/>

If *NonEqDEM* model is selected, two additional attributes may be specified. :xml:`RiemannSolver` to determine which Riemann solver will be used during simulation, either *HLLC* or *Exact*, and :xml:`contactAreas` to specify which contact area algorithm will be used, either *Symmetric* or *NonSymmetric*.

.. code-block:: xml

	<flowModel name="NonEqDEM" numberPhases="2" alphaNull="false" RiemannSolver="HLLC" contactAreas="Symmetric"/>
	<flowModel name="NonEqDEM" numberPhases="2" alphaNull="false" RiemannSolver="Exact" contactAreas="NonSymmetric"/>

If *NonEq* model is selected, an additional attribute must be set: :xml:`solver` to specify the solver of the Riemann problem. The HLLAC (HLLC-type) solver presented in :cite:`LOCHON2016733`, the HLLC4waves solver presented in :cite:`FURFARO2015159`, the HLL- or Rusanov-type solvers presented in :cite:`sab99` (with different closure laws) are implemented in ECOGEN to solve the total-non-equilibrium model :cite:`sab99`. They are referenced in Section :ref:`Sec_model_NonEq`.

.. code-block:: xml

	<flowModel name="NonEq" numberPhases="2" alphaNull="false" solver="HLLAC"/>
	<flowModel name="NonEq" numberPhases="2" alphaNull="false" solver="HLLC4waves"/>
	<flowModel name="NonEq" numberPhases="2" alphaNull="false" solver="HLL"/>

If *PressureVelocityEq* model is selected, an additional attribute might be set: :xml:`threshold`. If the volume fraction of a phase is superior to 1 - threshold, the pressure-relaxation procedure won't be applied, otherwise, it will be applied. If this attribute is not specified, the default value is 1.e-5.

.. code-block:: xml

	<flowModel name="PressureVelocityEq" numberPhases="2" numberTransports="0" alphaNull="false" threshold="1.e-6"/>

Equations of state (EOS)
------------------------

.. code-block:: xml

	<EOS name="IG_air.xml"/>

The *model.xml* input file **must contain** as many :xml:`<EOS>` markups as number of phases specified in the :ref:`Sec:input:FlowModel` markup. Each phase is described thanks to relations and parameters (see Section :ref:`Sec:IO:materials` for more details). The values of these parameters are specified in a separate file: the name attribute contains the name of this file which must be placed in the **ECOGEN/libEOS/** folder. Some fluid files are already present in the ECOGEN package.

.. _Sec:input:Transport:

Global accuracy order of the numerical scheme
---------------------------------------------

By default, the numerical scheme of the fluid-flow solver is first order in time and space.
ECOGEN can use a second-order scheme in time and space using the MUSCL method with TVD limiters (see :cite:`toro2013riemann` for example for an overview of the method).
Depending on whether the mesh is **Cartesian** (with/without AMR) or **unstructured**, the MUSCL method used and its parameters are different.

Cartesian mesh (with/without AMR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For a cartesian mesh (with/without AMR), second-order scheme is based on a TVD slope limiter; see :cite:`schmidmayer2020ecogen` for details.
In this case, the optional markup :xml:`<secondOrder>` can be inserted in the *model.xml* input file as in the following example:

.. code-block:: xml

	<secondOrder>
	  <globalLimiter>minmod</globalLimiter>
	  <interfaceLimiter>superbee</interfaceLimiter>                            <!-- optional node -->
	  <globalVolumeFractionLimiter>minmod</globalVolumeFractionLimiter>        <!-- optional node -->
	  <interfaceVolumeFractionLimiter>thinc</interfaceVolumeFractionLimiter>   <!-- optional node -->
	</secondOrder>

The :xml:`<secondOrder>` markup must contain the node :xml:`<globalLimiter>`. The other nodes are optional. The slope limiters available in ECOGEN are the following: minmod :cite:`roe1986superbee`, vanleer :cite:`vanLeer1974`, vanalbada :cite:`vanAlbada1997`, mc :cite:`van1977towards`, superbee :cite:`roe1986superbee`, koren :cite:`vreugdenhil1993koren` and thinc :cite:`shyue2014thinc` (only for volume fraction at the interface location). Markup significations are:

- :xml:`<globalLimiter>`: Applied everywhere and on all variables unless it is overwritten by the following optional limiters.
- :xml:`<interfaceLimiter>`: Applied on all variables but only at the interface location. By default is equal to the global limiter.
- :xml:`<globalVolumeFractionLimiter>`: Applied everywhere but only on the volume-fraction and transport equations (THINC is only applied on the volume fraction) unless it is overwritten by the interface volume-fraction limiter. By default is equal to the global limiter.
- :xml:`<interfaceVolumeFractionLimiter>`: Applied only at the interface location and on the volume-fraction and transport equations (THINC is only applied on the volume fraction). By default is equal to the interface limiter.

Unstructured mesh
~~~~~~~~~~~~~~~~~
For an unstructured mesh, second-order scheme is based on a TVD gradient limitation; see :cite:`chiapolino2017sharpening` for details.
In this case, the chosen limiter is a global limiter applied everywhere and on all variables.
Currently the gradient limiters are the following: minmod :cite:`roe1986superbee` and superbee :cite:`roe1986superbee`.
The optional markup :xml:`<secondOrder>` can be inserted in the *model.xml* input file as in the following example:

.. code-block:: xml

  <secondOrder>
    <globalLimiter>minmod</globalLimiter>
  </secondOrder>

.. important::

  Since the MUSCL method relies on gradients computation, the gradient method chosen should be compatible with unstructured meshes (see the markup :xml:`<gradient>`).

Gradient method
---------------

In ECOGEN, to compute gradients, it is possible to use:

- finite-difference-like gradient on **cartesian** mesh (without AMR) (:xml:`<method>finite-difference</method>`)
- Green-Gauss gradient on **cartesian** and **unstructured** mesh (:xml:`<method>green-gauss</method>`).

By default (without the XML markup :xml:`<gradient>`), the gradients will be computed using the Green-Gauss scheme.
To define explicitly the gradient method, one can use:

.. code-block:: xml

  <gradient>
    <method>finite-difference</method> <!-- For finite difference like gradient scheme -->
  </gradient>

.. note::

	The Green-Gauss method for the gradient computation must be selected when:

	- AMR is used,
	- periodic boundary conditions are used,
	- unstructured meshes, in case second-order scheme is set and/or additional physics are used (see Section :ref:`Sec:input:additionalPhysics`).

Advected additional variables
-----------------------------

.. code-block:: xml

	<transport name="color"/>

The *model.xml* input file **must contain** as many :xml:`<transport>` markups as number of transports specified in the :ref:`Sec:input:FlowModel` markup. Each transported variable is described by its name. The default number of advected variable is 0.

Relaxation procedures
---------------------

An additional markup :xml:`<relaxation>` may be used to impose some specific equilibrium between the phases depending on the flow model used. The attribute :xml:`type` specifies the type of equilibrium:

- *U*: A velocity equilibrium is imposed at every location of the flow. The attribute :xml:`speed` can be added to specify the speed at which the relaxation operates (only *infinite* :cite:`sab99` for now). Default is *infinite*.

.. code-block:: xml

	<relaxation type="U" speed="infinite"/>

- *P*: A pressure equilibrium is imposed at every location of the flow. The attribute :xml:`speed` can be added to specify the speed at which the relaxation operates (*infinite* :cite:`relaxjcp` or *finite* :cite:`schmidmayer2023UEq, biasiori2023phenomenological`). Default is *infinite*. In this case, an additional attribute might be set: :xml:`threshold` to determine when the pressure-relaxation procedure is bypassed. If a volume fraction exceeds 1 - threshold, the procedure is not applied; otherwise, it is executed. If a finite relaxation is chosen, the :xml:`rate` and :xml:`solver` have to be specified. The rate is a real number while the solver can either be *Euler* :cite:`schmidmayer2023UEq` or *LSODA* :cite:`hindmarsh1983odepack, petzold1983LSODA`.

.. code-block:: xml

	<relaxation type="P" speed="infinite" threshold="1.e-6"/>

.. code-block:: xml

	<relaxation type="P" speed="finite" rate="1.e-5" solver="Euler"/>

- *PT*: Both pressure and thermal equilibrium are imposed at every location of the flow. It does not require additional attributes.

.. code-block:: xml

	<relaxation type="PT"/>

- *PTMu*: A thermodynamical equilibrium is imposed at every location of the flow. It must be associated with the node :xml:`<dataPTMu>` whose attributes are :xml:`liquid` and :xml:`vapor` to specify the name of the EOS of the liquid and the vapor phase. Hereafter the complete node when *PTMu* is used:

.. code-block:: xml

	<relaxation type="PTMu">
	  <dataPTMu liquid="SG_waterLiq.xml" vapor="IG_waterVap.xml"/>
	</relaxation>

Source terms
------------

The additional :xml:`<sourceTerms>` markup can be used to numerically integrate some source terms in the equations. The attribute :xml:`type` selects the source term:

- *heating*: Related to a thermal energy heating/cooling. This attribute requires the :xml:`<dataHeating>` node with the attribute :xml:`volumeHeatPower`; a real number corresponding to the power by volume unit added to the flow (unit: W/m3 (SI)).

.. code-block:: xml

	<sourceTerms type="heating" order="EULER">
	  <dataHeating volumeHeatPower="1.e6"/>
	</sourceTerms>

- *gravity*: If the gravity is considered. The node :xml:`<dataGravity>` must be present with the attributes :xml:`x`, :xml:`y` and :xml:`z` giving the coordinates for the gravity acceleration vector in real numbers (unit: m/s2 (SI)).

.. code-block:: xml

	<sourceTerms type="gravity" order="EULER">
	  <gravity x="0." y="-9.81" z="0."/>
	</sourceTerms>

- *MRF*: For a simulation in a moving reference frame. Allow to compute solution in a rotating frame :cite:`caze2024modeling`. The node :xml:`<omega>` requires the attributes :xml:`x`, :xml:`y` and :xml:`z` giving the coordinates for the rotating vector in real numbers (unit: rad/s (SI)). The node :xml:`<timeToOmega>` is optional and allow to specify a progressing acceleration (linear) to the final rotating velocity (requires the attribute :xml:`tf` for acceleration time).

.. code-block:: xml

	<sourceTerms type="MRF" order="RK4">
	  <omega x="0." y="0." z="1."/>
	  <timeToOmega tf="1.e-3"/>  <!-- Optional: If activated, the angular velocity increase linearly to omega in during tf -->
	</sourceTerms>

- *Acoustic waves*: Acoustic waves can be modelled with sinusoidal or Gaussian profils. This follows the work of :cite:`maeda2017source` for single-phase flows and has been extended for multiphase system of equations. The node :xml:`typeAcousticWave` requires you to mention either a sinusoidal or a Gaussian type of pulses. For both, nodes :xml:`sourceLocation` and :xml:`pulseDirection` are required and have the attributes :xml:`x`, :xml:`y` and :xml:`z` giving the coordinates of the source and the direction of the pulses in real numbers (unit: m (SI) and none), respectively. For :xml:`pulseDirection`, the sum of the attributes must be equal to 1. Additionally, the :xml:`pressureAmplitude` (unit: Pa (SI)) and :xml:`delayedTime` (unit: s (SI)) of the pulses are required for both. The sinusoidal pulses require the :xml:`numberOfPulse` (can be a real number) and :xml:`frequency` (unit: Hz (SI)), while the Gaussian pulses require the :xml:`pulseWidth` (unit: m (SI)). XML blocks can be written for sinusoidal or Gaussian pulses as:

.. code-block:: xml

	<sourceTerms type="acousticWave" order="EULER">
	  <typeAcousticWave type="planeSinusoidalPulse">
	    <dataPlaneSinusoidalPulse pressureAmplitude="10000." delayedTime="5.e-7" numberOfPulse ="4." frequency="3.e5">
	      <sourceLocation x="10.e-3" y="10.e-3" z="0.5"/>
	      <pulseDirection x="1." y="0." z="0."/>
	    </dataPlaneSinusoidalPulse>
	  </typeAcousticWave>
	</sourceTerms>

.. code-block:: xml

	<sourceTerms type="acousticWave" order="EULER">
	  <typeAcousticWave type="planeGaussianPulse">
	    <dataPlaneGaussianPulse pressureAmplitude="10." delayedTime="20.e-6" pulseWidth="5.e-6">
	      <sourceLocation x="10.03e-3" y="0." z="0."/> <!-- Must be placed at: cell_center + 1/4 * cell_size_in_wanted_direction -->
	      <pulseDirection x="1." y="0." z="0."/>
	    </dataPlaneGaussianPulse>
	  </typeAcousticWave>
	</sourceTerms>

Source-term integration can be done up to 4th order. Available integration schemes are:

- *EULER* (1st order)
- Runge--Kutta 2 *RK2* (2nd order)
- Runge--Kutta 4 *RK4* (4th order)

Scheme selection is done through :xml:`order` attribute.

Symmetry terms
--------------

Both cylindrical (2D) and spherical (1D) symmetries are implemented. The additional :xml:`<symmetryTerm>` markup can be used. It requires the attribute :xml:`type` that can take the value *cylindrical* or *spherical*. It also requires an additional node to specify the radial axis:

- Cylindrical:

.. code-block:: xml

	<symmetryTerm type="cylindrical">
	  <dataSymCyl radialAxis="X"/>
	</symmetryTerm>

- Spherical:

.. code-block:: xml

	<symmetryTerm type="spherical">
	  <dataSymSpher radialAxis="X"/>
	</symmetryTerm>

.. _Sec:input:additionalPhysics:

Additional physics
------------------

Depending on the model chosen in section :ref:`Sec:input:FlowModel`, additional physical effects can be added. This is the case for surface tension, viscosity, turbulence and conductive heat transfers. These additional physical effects are obtained thanks to the additional markup :xml:`additionalPhysics` with the attribute :xml:`type` that can take different value according to the chosen effect.

Surface tension
~~~~~~~~~~~~~~~
This physical effect is obtained by using the type *surface tension* and its treatment is detailed in :cite:`schmidmayer2017capillary`. Then it requires the node :xml:`dataSurfaceTension` with following attributes:

- :xml:`transport`: This is the name of advected variable used as color function for surface-tension terms. This advected variable has been precised in the section :ref:`Sec:input:Transport`. The name should be the same.
- :xml:`sigma`: A real number for the surface-tension coefficient (unit: N/m (SI)).

.. code-block:: xml

	<additionalPhysics type="surfaceTension" >
	  <dataSurfaceTension transport="color" sigma="72.e-3"/>
	</additionalPhysics>

Viscosity
~~~~~~~~~
This physical effect is obtained by using the type *viscosity* and it needs the attribute :xml:`mu` of the :xml:`<physicalParameters>` to be set in the EOS files (:ref:`Sec:IO:materials`).

.. code-block:: xml

	<additionalPhysics type="viscosity"/>

Turbulence: static Smagorinsky
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This physical effect is obtained by using the type *turbulence*. One can choose the turbulence model by changing the :xml:`<model>` node. The current implemented LES model is the static Smagorinsky :cite:`smagorinsky1963`. It needs the attributes :xml:`C_s`, :xml:`C_I` and :xml:`Pr_T` of the :xml:`<smagorinsky>` node to be set in the model file (:ref:`Sec:input:Model`). The value of :xml:`C_s` depends on the problem, but here are some common values: :xml:`C_s=0.16` :cite:`pinomartin2000`, :xml:`C_s=0.0085` :cite:`kotov2016`. Additionally, it needs the attributes :xml:`mu` and :xml:`lambda` of the :xml:`<physicalParameters>` to be set in the EOS files (:ref:`Sec:IO:materials`) to include the molecular viscosity and conductivity effects.
Finally, one can set the attribute :xml:`ducros` to :xml:`true/false` in order to activate/deactivate the Ducros et al. flow sensor :cite:`ducros1999`. This feature is currently available for the Euler flow model.

.. code-block:: xml

	<additionalPhysics type="turbulence">
		<turbulence model = "smagorinsky">
			<smagorinsky C_s="0.12" C_I="0.09" Pr_T="0.71" ducros="false"/>
		</turbulence>
	</additionalPhysics>

.. warning::
	When enabling the turbulence model, one must not enable viscosity and conductivity, as they are already taken into account. Setting the attributes :xml:`mu` and :xml:`lambda` of the :xml:`<physicalParameters>` is enough.

Strang splitting (only for Sources and Additional physics)
----------------------------------------------------------

For a system of balance laws of this form:

.. math::
  :nowrap:

  \begin{align}
    &\cfrac{\partial \mathbf{Q}}{\partial t} + \textrm{div}( \mathbf{F}(\mathbf{Q}) ) + \mathbf{G}(\mathbf{Q}) \textrm{div}(\mathbf{u}) = \mathbf{S}(\mathbf{Q}),
  \end{align}

it is common to split the hyperbolic part and the source part:

.. math::
  :nowrap:

  \begin{array}{*{20}{l}}
    \mathscr{H}(t) : \quad \cfrac{\partial \mathbf{Q}}{\partial t} + \textrm{div}( \mathbf{F}(\mathbf{Q}) ) + \mathbf{G}(\mathbf{Q}) \textrm{div}(\mathbf{u}) = 0,\\
    \mathscr{S}(t) : \quad \cfrac{\partial \mathbf{Q}}{\partial t} = \mathbf{S}(\mathbf{Q}),
  \end{array}

where :math:`\mathscr{H}(t)` is the hyperbolic operator and :math:`\mathscr{S}(t)` the source operator. Each operator is discretized by an appropriate numerical method of any order of accuracy. Such splitting procedure leads to a splitting error independant of the numerical discretization of the operators. To obtain the solution at the next step, one have to combine the numerical solutions of each operator in some ordered sequence. Depending on the sequence, the splitting error (and thus the accuracy) change. For a first order accurate in time, the Godunov splitting can be used:

.. math::
  :nowrap:

  \begin{align}
    &\mathbf{Q}^{n+1} = \mathscr{S}(\Delta t) \, \mathscr{H}(\Delta t) \, \mathbf{Q}^{n}
  \end{align}

For a second order accurate in time, the Strang splitting can be considered :cite:`Strang1968` :cite:`Leveque2002`:

.. math::
  :nowrap:

  \begin{align}
    &\mathbf{Q}^{n+1} = \mathscr{S} \left( \cfrac{\Delta t}{2} \right) \, \mathscr{H}(\Delta t) \, \mathscr{S} \left( \cfrac{\Delta t}{2} \right) \, \mathbf{Q}^{n}
  \end{align}

If the nondiscretized operators :math:`\mathscr{H}` and :math:`\mathscr{S}` commute :math:`(\mathscr{H} \mathscr{S} - \mathscr{S} \mathscr{F} = 0)`, then the Godunov splitting is also second order accurate in time. By default, the Godunov splitting is used in ECOGEN. To use the Strang splitting, put the following line in model.xml:

.. code-block:: xml

	<splitting type="strang"/>

Others
~~~~~~

In dev...

Low-Mach preconditioning
------------------------

In case of low-speed flows, a low-Mach preconditioning can be applied to the selected flow model.
The method is detailed in :cite:`lemartelot2013lowmach`.
Currently this option is compatible with the flow models: Euler, UEq and PUEq.

.. code-block:: xml

  <lowMach state="true"/>

Note that the minimum reference Mach number can be chosen if required (default is 0.01).
When the local Mach number of the flow is lower than this threshold, it is replaced by the threshold value in order to keep a reasonable computation time.

.. code-block:: xml

  <lowMach state="true" machRefMin="5.e-2"/>

For an example of use, see the test case presented in the Section :ref:`Sec:tests:euler:2d:nozzleLowMach`.

1D geometry with smooth cross section variation
-----------------------------------------------
In case of a 2D geometry with smooth cross section variation, it is possible to use a 1D geometry and reproduce 2D effects with this attribute as detailed in the first appendix of :cite:`lemartelot2013lowmach`.
This option is only compatible with unstructured meshes defined with flow direction along the X-axis.
Boundary condition contributions in other directions than the X-axis should use *nullFlux* boundary condition (see :ref:`Sec:input:InitialConditions`).
Note that this option is available for the flow models: Euler, UEq, PUEq and EulerHomogeneous.
For an example of use, see the test case presented in the Section :ref:`Sec:tests:euler:2d:nozzleLowMach`.

.. code-block:: xml

  <geometry smoothCrossSection1d="true"/>
