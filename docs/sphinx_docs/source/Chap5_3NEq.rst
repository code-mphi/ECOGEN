.. role:: xml(code)
  :language: xml

.. _Sec_model_NonEq:

**************************
Non-equilibrium model
**************************

Total-non-equilibrium flows are solved in ECOGEN using the total-non-equilibrium model :cite:`sab99`. For *N* phases involved and without any extra physics (surface tension, viscosity, heat flux ...), this model reads:

.. math::
  :nowrap:

  \begin{equation}
  \label{system_NEq}
  \left\{
  {\begin{array}{*{20}{l}}
    \cfrac{\partial \alpha_{k}}{\partial t} + \mathbf{u}_{I} \cdot \nabla \alpha_{k} = R_{\alpha_{k}} , \\
		\cfrac{\partial \alpha_{k} \rho_{k}}{\partial t} + \nabla \cdot \left( \alpha_{k} \rho_{k} \mathbf{u}_{k} \right) = R_{\mathrm{masse}_{k}} , \\
		\cfrac{\partial \alpha_{k} \rho_{k} \mathbf{u}_{k}}{\partial t} + \nabla \cdot \left( \alpha_{k} \rho_{k} \mathbf{u}_{k} \otimes \mathbf{u}_{k} \right) = p_{I} \nabla \alpha_{k} + \boldsymbol{R}_{\mathrm{mvt}_{k}},  \\
		\cfrac{\partial \alpha_{k} \rho_{k} E_{k}}{\partial t} + \nabla \cdot \left( \alpha_{k} \rho_{k} E_{k} \mathbf{u}_{k} \right) = p_{I} \mathbf{u}_{I} \cdot \nabla \alpha_{k} + R_{{E}_{k}},
  \end{array}} \right.\
  \end{equation}

where :math:`\alpha_k`, :math:`\rho_k`, :math:`p_k`, :math:`\mathbf{u_k}` and :math:`E_k` are the volume fraction, density, pressure, velocity and total energy of each phase, respectively, and for which :math:`k` indicates the phase index.
The mixture density, pressure and velocity are

.. math::
  :nowrap:

  \begin{equation}
    \rho = \sum_{k=1}^N \alpha _k \rho_k  \quad \text{,} \quad p = \sum_{k=1}^N \alpha _k p_k \quad \text{and} \quad \mathbf{u} = \sum_{k=1}^N Y_k u_k,
  \end{equation}

with :math:`Y_k` the mass fractions define as

.. math::
  :nowrap:

  \begin{equation}
    Y_k = \frac{\alpha_k \rho_k}{\rho} .
  \end{equation}

The total energy per phase is

.. math::
  :nowrap:

  \begin{equation}
    E_k = e_k + \frac{1}{2} \| \mathbf{u_k} \|^2,
  \end{equation}

:math:`e` is the mixture specific internal energy

.. math::
  :nowrap:

  \begin{equation}
    e = \sum_{k=1}^N Y_k e_k \left( \rho_k , p_k \right) .
  \end{equation}

:math:`e_k \left( \rho_k , p_k \right)` is defined via the equation of state (EOS), such as of perfect gas or stiffened gas.
:math:`p_{I}` and :math:`\mathbf{u}_{I}` denote the pressure and velocity at the interfaces. Their definition depends on the Riemann solver employed to compute the flux at each cell interface. Finally, the coefficients :math:`R` correspond to relaxation source terms.

This model can be solved thanks to the Discrete Equation Method (DEM) presented in :cite:`absau2003`, or alternatively by the HLLAC solver presented in :cite:`LOCHON2016733`, the HLLC4waves solver presented in :cite:`FURFARO2015159`, the HLL- or Rusanov-type solvers presented in :cite:`sab99` with different closure laws.

**Remarks:**
  - The DEM can be applied with two different Riemann-problem solvers: an HLLC and an exact solver. By default, the contact surface definition adopted is that of Paula et al. :cite:`paula2023robust`, which is symmetric. Alternatively, the formulation proposed by Chinnayya :cite:`chinnayya2002construction`, which is non-symmetric, extended to *N* phases, can be chosen.
  - The HLLAC solvers must be used with the gas phase placed first, followed by the solid phase.
  - The HLL/Rusanov solver can be used with symmetric (default) or Baer & Nunziato (commented out) closure law. To ensure numerical stability, the CFL number must be set to 0.5 or lower.

.. toctree::

   ./Chap5_3NEq_1D.rst
   ./Chap5_3NEq_2D.rst
