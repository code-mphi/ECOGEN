.. role:: xml(code)
  :language: xml

****************************************
Extended-Navier-Stokes-Korteweg-1D model
****************************************

The 1D Extended-Navier-Stokes-Korteweg (ENSK) model is solved in ECOGEN. This model allows to consider compressible flows with capillarity and liquid-vapor phase transition :cite:`DidierlaurentInReview`. The model writes:

.. math::
  :nowrap:

  \begin{equation}
  \label{system_ENSK}
  \left\{
  {\begin{array}{*{20}{l}}
    \cfrac{\partial}{\partial t} \rho + \cfrac{\partial}{\partial x} ( \rho u ) = 0, \\
    \cfrac{\partial}{\partial t} (\rho \, u) + \cfrac{\partial }{\partial x} ( \rho \, u^2 - \sigma_e ) = 0,\\
    \cfrac{\partial }{\partial t} \left(\cfrac{\rho}{\eta} \right) + \cfrac{\partial }{\partial x} \left( \cfrac{\rho}{\eta} \, u \right) = - \cfrac{\rho}{\eta^2} \, \psi,\\
    \cfrac{\partial }{\partial t} (\rho \, \psi) + \cfrac{\partial }{\partial x} \left( \rho \, \psi \, u - \varepsilon^2 \, \mathcal{K} \, \phi \right) = \cfrac{\lambda}{\beta} \left( 1 - \cfrac{\eta}{\rho} \right) - \cfrac{\tau \, \mu}{\beta \, \rho^2} \psi,\\
    \cfrac{\partial }{\partial t} \phi + \cfrac{\partial }{\partial x} ( \phi \, u - \psi ) = 0,\\
    \cfrac{\partial }{\partial t} (\rho \, E_e) + \cfrac{\partial }{\partial x} ( \rho \, u \, E_e - \sigma_e \, u - \varepsilon^2 \,\mathcal{K} \, \phi \, \psi) = 0,
  \end{array}} \right.\
  \end{equation}


where :math:`\rho`, :math:`u` and :math:`E_e` are the density, velocity and extended total energy, respectively. The extended stress :math:`\sigma_e` and the extended total energy are:

.. math::
  :nowrap:

  \begin{align}
    &\sigma_e = - \left( p + \varepsilon^2 \, (\rho \, \mathcal{K}^{\prime} + \mathcal{K} ) \cfrac{\phi^2}{2} + \lambda \, \eta \, \left( 1 - \cfrac{\eta}{\rho} \right) \right),\\
    &E_e = \cfrac{\beta}{2} \, \psi^2 + \cfrac{1}{2} u^2 + e_h + \cfrac{\varepsilon^2}{2} \cfrac{\mathcal{K}}{\rho} \, \phi^2 + \cfrac{\lambda}{2} \left( 1 - \cfrac{\eta}{\rho} \right)^2,
  \end{align}

where :math:`p` is the thermodynamic pressure and :math:`e_h` the internal energy. The parameters :math:`\varepsilon` and :math:`\tau` control the magnitude of the capillarity and of the viscosity, respectively. The extended variables :math:`\eta`, :math:`\phi` and :math:`\psi` are additional variables close to :math:`\rho`, :math:`\cfrac{\partial \rho }{\partial x}`, :math:`\cfrac{d \rho }{d t}`, provided that the extended parameters :math:`(\lambda,\beta)` are well-chosen :cite:`Didierlaurent2025` :cite:`DidierlaurentInReview`. For :math:`\lambda` sufficiently large and :math:`\beta` sufficiently small, the ENSK model provides a good hyperbolic approximation of the Navier-Stokes-Korteweg model:

.. math::
  :nowrap:

  \begin{equation}
  \label{system_NSK}
  \left\{
  {\begin{array}{*{20}{l}}
    \cfrac{\partial }{\partial t} \rho + \cfrac{\partial }{\partial x} (\rho u) = 0, \\
    \cfrac{\partial }{\partial t} (\rho \, u) + \cfrac{\partial }{\partial x} ( \rho \, u^2 - \sigma ) = 0,\\
    \cfrac{\partial }{\partial t} (\rho \, E) + \cfrac{\partial }{\partial x} \left( \rho \, u \, E - \sigma \, u - \varepsilon^2 \, \mathcal{K} \, \cfrac{\partial \rho }{\partial x} \, \cfrac{d \rho }{d t} \right) = 0,
  \end{array}} \right.\
  \end{equation}

where

.. math::
  :nowrap:

  \begin{align}
    &\sigma = - \left( p - \varepsilon^2 \,\left[ \cfrac{1}{2} (\rho \, \mathcal{K}^{\prime} - \mathcal{K} ) \left(\cfrac{\partial \rho }{\partial x}\right)^2 + \rho \, \mathcal{K} \, \cfrac{\partial^2 \rho }{\partial x^2} \right] - \tau \, \mu \, \cfrac{\partial u }{\partial x} \right),\\
    &E = \cfrac{1}{2} u^2 + e_h + \cfrac{\varepsilon^2}{2} \cfrac{\mathcal{K}}{\rho} \, \left(\cfrac{\partial \rho }{\partial x}\right)^2.
  \end{align}

Equations of state
==================

From this generic model, one can choose different equations of state to obtain particular models lying in different fields. For example, the hydrodynamic form of the Defocusing Nonlinear Schrödinger (NLS) equation :cite:`Dhaouadi2019` can be obtain by choosing:

.. math::
  :nowrap:

  \begin{align}
    &\textrm{Defocusing NLS EOS:} \quad e_h = \cfrac{\rho^2}{2}, \quad \mathcal{K} = \cfrac{1}{4 \, \rho}, \quad \varepsilon \quad \textrm{arbitrary}, \quad \tau = 0.\\
  \end{align}

It is also possible to describe the dynamics of a bistable tape spring :cite:`Bourgeois2020` by choosing:

.. math::
  :nowrap:

  \begin{align}
    &\textrm{Polynomial EOS:} \quad e_h = a \, \epsilon^4 + b \, \epsilon^3 + c \, \epsilon^2, \quad \epsilon := \cfrac{1}{\rho} - 1 \quad \mathcal{K} = \cfrac{1}{\rho^3}, \quad \varepsilon \quad \textrm{arbitrary}, \quad \tau = 0.\\
  \end{align}

Finally, one can also describe the dynamics of phase transition by considering the van der Waals EOS :cite:`DidierlaurentInReview`:

.. math::
  :nowrap:

  \begin{align}
    &\textrm{Van der Waals EOS:} \quad e_h = C_v \, T - a \, \rho, \quad \left\{
  {\begin{array}{*{20}{l}}
    \mathcal{K} = \cfrac{1}{4 \, \rho}, \quad (\varepsilon,\tau) \quad \textrm{arbitrary},\\
    \mathcal{K} = C, \quad (\varepsilon,\tau) \quad \textrm{arbitrary},\\
  \end{array}} \right.\ \quad \mu = 2 \, \rho \, \sqrt{ \rho \, \mathcal{K}}.
  \end{align}

In this last EOS, :math:`\mathcal{K}` can be considered either as a function of the density or as a constant :math:`C`.


Tests cases
===========
Tests are provided with ECOGEN package and may be described in details later.

.. code-block:: xml

  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/shockTubes/nonlinearSchrodingerEos/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/shockTubes/nonlinearSchrodingerEosAmr/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/shockTubes/polynomialEos/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/shockTubes/polynomialEosAmr/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/shockTubes/vanDerWaalsTestEos/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/shockTubes/vanDerWaalsDodecaneEos/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/cavitations/vanDerWaalsTestEos/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/impacts/vanDerWaalsTestEos/</testCase>
  <testCase>./libTests/referenceTestCases/ExtendedNavierStokesKorteweg1D/travellingWave/vanDerWaalsTestEos/</testCase>
