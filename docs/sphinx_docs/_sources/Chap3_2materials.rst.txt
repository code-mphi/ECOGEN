.. _Sec:IO:materials:

Materials
=========

For each phase, an equation of state is required. The data for a given phase are gathered in a file. Every file must be in the folder **ECOGEN/libEOS/** and must follow rules depending on the EOS. Six equations of state are implemented in ECOGEN:

- :ref:`Sec:input:IdealGas`: For gaseous phase only.
- :ref:`Sec:input:StiffenedGas`: For condensed matter (liquid, solid) in a pressure range where the compressible assumption is reasonable.
- :ref:`Sec:input:NobleAbelStiffenedGas`: For condensed matter (liquid, solid) subject to phase change.

Only for the Extended-Navier-Stokes-Kortewed 1D model:

- :ref:`Sec:input:NonlinearShrodinger`: For quantum fluids described by the Defocusing non-linear Schrödinger equation.
- :ref:`Sec:input:Polynomial` (non-convex): For a bistable bar.
- :ref:`Sec:input:VanDerWaals` (non-convex): For liquid and gaseous phase.

.. _Sec:input:IdealGas:

Ideal Gas
---------
It is a simple relation linking the pressure :math:`p`, the temperature :math:`T` and the density :math:`\rho` or the specific volume :math:`v=1/\rho`:

.. math::

	pv=rT

:math:`r` is the universal gas constant divided by the molar weight of the gas (unit: :math:`J/(K.kg)` (SI)). The Gibbs relation as well as Maxwell relations yield the definition of the internal energy:

.. math::

	e(p,v) = \frac{p}{\gamma-1} v + e_{ref}

One can easily obtain the entropy:

.. math::

	s(p,T) = c_v \ln{\left(\frac{T^{\gamma}}{p^{\gamma-1}}\right)} + s_{ref}

The following parameters are assumed constant and :math:`\gamma`, :math:`c_v`, :math:`e_{ref}` and :math:`s_{ref}` must be specified in the material file, as in the following example:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<parametersEOS>
		<EOS fichier="IG_air.xml" type="IG"/>
		<parameters
			gamma="1.4"
			cv="800" cv_save="717.46"
			energyRef="0."
			entropyRef="0.">
		</parameters>
		<physicalParameters
			mu="1.85e-5"
			lambda="0.0262">
		</physicalParameters>
	</parametersEOS>

.. _Sec:input:StiffenedGas:

Stiffened Gas
-------------
This equation of state was first presented in 1971 by Harlow & Amdsen :cite:`harlow1968stiffenedGas`:

.. math::

	e = \left(\frac{p + \gamma p_{\infty}}{\gamma-1}\right) v + e_{ref}

This EOS takes into account the molecular attraction within condensed matter. This quite simple EOS is largely used in numerical simulation based on diffuse interface methods because it is able to reproduce shock relation as well as saturation curve for a liquid--vapor mixture when the phase transition is modeled.
Hereafter, some useful relations for the specific volume, the enthalpy and the entropy:

.. math::

	v(p,T) &= \frac{ (\gamma-1) c_v T}{p + p_{\infty}} \\
	h(T) &= \gamma c_v T + e_{ref} \\
	s(p,T) &= c_v \ln{\left(\frac{T^{\gamma}}{\left(p+p_{\infty}\right)^{\gamma-1}} \right)} + s_{ref}

A useful description of these relations can be found in :cite:`SGEOS`.

The ideal-gas law is recovered if :math:`p_{\infty}=0`. The following parameters are assumed constant :math:`\gamma`, :math:`p_{\infty}`, :math:`c_v`, :math:`e_{ref}` and :math:`s_{ref}`. These parameters must be specified in the material file, as in the following example:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<parametersEOS>
		<EOS fichier="SG_water.xml" type="SG"/>
		<parameters
			gamma="4.4"
			pInf="6.e8"
			cv="1000.0"
			energyRef="0."
			entropyRef="0.">
		</parameters>
		<physicalParameters
			mu="1.e-3"
			lambda="0.6">
		</physicalParameters>
	</parametersEOS>

.. _Sec:input:NobleAbelStiffenedGas:

Noble-Abel Stiffened Gas
------------------------

This equation of state is an improved version of the stiffened-gas equation of state in which a parameter of covolume `b` is introduced. This parameter allows to take into account the repulsive short distance effects linked to molecular motion. The Noble-Abel Stiffened-Gas (NASG) EOS was developed in 2016 by Le Métayer & Saurel :cite:`le2016noble` and it reads:

.. math::

	e(p,v) = \frac{p+\gamma p_{\infty}}{\gamma-1}(v-b) + e_{ref}

As the stiffened-gas EOS, the NASG EOS is able to take into account thermal agitation and attractive effects of condensed matter. This EOS is particularly interesting when phase transition is modeled because it is able to reproduce more accurately saturation curves of a liquid--vapor mixture than the traditional stiffened-gas EOS.
Hereafter, some useful relations for the specific volume, the enthalpy and the entropy:

.. math::

	v(p,T) &= \frac{ (\gamma-1) c_v T}{p + p_{\infty}} + b \\
	h(p,T) &= \gamma c_v T + + b p +e_{ref} \\
	s(p,T) &= c_v \ln{\left(\frac{T^{\gamma}}{\left(p+p_{\infty}\right)^{\gamma-1}} \right)} + s_{ref}

A complete description of the obtention of these relations is given in :cite:`le2016noble`.

The stiffened-gas EOS is recovered if :math:`b=0` and the ideal-gas law is recovered if :math:`b = p_{\infty} = 0`. The parameters :math:`\gamma`, :math:`p_{\infty}`, :math:`b`, :math:`c_v`, :math:`e_{ref}` and :math:`s_{ref}` are assumed constant. These parameters must be specified in the material file, as in the following example:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<parametersEOS>
		<EOS fichier="NASG_dodLiq.xml" type="NASG"/>
		<parameters
			gamma="1.09"
			pInf="1.159e8"
			b="7.51e-4"
			cv="2393"
			energyRef="-794696."
			entropyRef="0.">
		</parameters>
		<physicalParameters
			mu="1.e-3"
			lambda="0.6">
		</physicalParameters>
	</parametersEOS>

.. _Sec:input:NonlinearShrodinger:

Defocusing non-linear Schrödinger
---------------------------------

The hydrodynamic form of the Defocusing Nonlinear Schrödinger equation :cite:`Dhaouadi2019` can be obtain by choosing:

.. math::

	e_h(\rho) = \cfrac{\rho^2}{2}, \quad \mathcal{K} = \cfrac{1}{4 \, \rho}.

No parameters are specified for this EOS:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<parametersEOS>
		<EOS type="NLS"/>
		<parameters no_params="1."/>
	</parametersEOS>


.. _Sec:input:Polynomial:

Polynomial
----------

A bistable Ericksen model is used to describe the dynamic of a bistable tape spring :cite:`Bourgeois2020`. The potential energy is:

.. math::

	e_h(\rho) = a \, \epsilon^4 + b \, \epsilon^3 + c \, \epsilon^2, \quad \epsilon(\rho) := \cfrac{1}{\rho} - 1, \quad \mathcal{K} = \cfrac{1}{\rho^3}.

In the following example, the parameters :math:`a,b,c` are chosen from :cite:`Bourgeois2020` in order to reproduce a double well potential:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
	<parametersEOS>
		<EOS type="Polynomial"/>
		<parameters
			a="0.14583333333333334"
			b="-1.1666666666666667"
			c="2.4583333333333335">
		</parameters>
	</parametersEOS>

.. _Sec:input:VanDerWaals:

Van der Waals
-------------

To describe both the liquid and the vapor through a single EOS, we can choose the van der Waals EOS :cite:`DidierlaurentInReview`:

.. math::

	 &p(\rho,T) = \cfrac{\rho \, r \, T}{1 - \rho \, b} - a \, \rho^2,\\
	 &e_h(\rho,T) = C_v \, T - a \, \rho,\\[6pt]
	 &s(\rho,T) = C_v \, \log{T} + r \, \log(1/\rho - b).

The parameters :math:`a,b,C_v` are defined at the critical point of the fluid. The constants :math:`T_c,P_c,M` are given. For this EOS, one can choose :math:`\mathcal{K}` as a function of :math:`\rho` or as a constant. Here is an example:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
		<parametersEOS>
		<EOS type="VDW"/>
		<parameters
			Tc="658.2"
			Pc="18.2e5"
			M="170.34"
			c_v="3200."
			kappa="2.7778e-6">
		</parameters>
	</parametersEOS>

Here, a value for :math:`\mathcal{K}` has been specified. If you want :math:`\mathcal{K}` as a function of :math:`\rho`, set the value to zero:

.. code-block:: xml

	<?xml version = "1.0" encoding = "UTF-8" standalone = "yes"?>
		<parametersEOS>
		<EOS type="VDW"/>
		<parameters
			Tc="658.2"
			Pc="18.2e5"
			M="170.34"
			c_v="3200."
			kappa="0.">
		</parameters>
	</parametersEOS>

By default, the choice :math:`\mathcal{K} = \cfrac{1}{4 \, \rho}` is made.
