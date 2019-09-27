.. _Sec:IO:materials:

Materials
=========

For each phase, an Equation of State is required. The data for a given phase are gathered in a file. Every file must be in the folder **ECOGEN/libEOS/** and must follow rules depending on the EOS. Two equations of states are implemented in ECOGEN:

- *Ideal gas*: for gaseous phase only.
- *Stiffened gas*: for condensed matter (liquid, solid) in a pressure range where the compressible assumption is reasonable. 

Ideal Gas
---------
It is a simple relation linking the pressure (:math:`p`), the temperature (:math:`T`) and the density (:math:`\rho`) or the specific volume (:math:`v=1/ρ`): 

.. math::
	
	pv=rT 

:math:`r` is the universal gas constant divided by the molar weight of the gas (unit: :math:`J/(K.kg)` (SI)). The Gibbs relation as well as Maxwell relations yield the definition of the internal energy: 

.. math:: 

	e(p,v) = \frac{p}{\gamma-1} v + e_{ref}

One can easily obtain the entropy: 

.. math::

	s(p,T) = c_v \ln{\left(\frac{T^{\gamma}}{p^{\gamma-1}}\right)} + s_{ref}

The following parameters are assumed constant and :math:`\gamma`, :math:`c_v`, :math:`e_{ref}` and :math:`s_{ref}` must be specified in the material file, as in the following example :

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

Stiffened Gas
-------------
This Equation of state was first presented in 1971 by Harlow & Amdsen :cite:`harlow1968stiffenedGas`:

.. math::
	
	e = \left(\frac{p + \gamma p_{\infty}}{\gamma-1}\right) v + e_{ref}

This EOS takes into account for the molecular attraction within condensed matter. This quite simple EOS is largely used in numerical simulation based on diffuse interface methods because it is able to reproduce shock relation as well as saturation curve for a liquid-vapor mixture when the phase transition is modeled.
Hereafter some useful relations for the specific volume, the enthalpy and the entropy:

.. math::

	v(p,T) &= \frac{ (\gamma-1) c_v T}{p + p_{\infty}} \\
	v(p,T) &= \frac{ (\gamma-1) c_v T}{p + p_{\infty}} \\
	h(T) &= \gamma c_v T + e_{ref} \\
	s(p,T) &= c_v \ln{\left(\frac{T^{\gamma}}{\left(p+p_{\infty}\right)^{\gamma-1}} \right)} + s_{ref}

A usefull description of these relations can be found in :cite:`SGEOS`.

The ideal gas law is recovered if :math:`p_{\infty}=0`. The following parameters are assumed constant :math:`\gamma`, :math:`p_{\infty}`, :math:`c_v`, :math:`e_{ref}` and :math:`s_{ref}`. These parameters must be specified in the material file, as in the following example:

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

