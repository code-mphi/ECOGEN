.. role:: xml(code)
  :language: xml

**************************
Velocity-equilibrium model
**************************

Velocity-equilibrium and mechanical-equilibrium flows are solved in ECOGEN using the velocity-equilibrium model :cite:`schmidmayer2021UEq`. For *N* phases involved and without any extra physics (surface tension, viscosity...), this model reads:

.. math::
  :nowrap:

  \begin{equation}
  \label{system_UEq}
  \left\{
  {\begin{array}{*{20}{l}}
    \cfrac{\partial \alpha_k}{\partial t} + \mathbf{u} \cdot \nabla \alpha_k & = \delta p_k , \\ 
    \cfrac{\partial \alpha_k \rho_k}{\partial t} + \nabla \cdot \left( \alpha_k \rho_k \mathbf{u} \right) & = 0 , \\
    \cfrac{\partial \rho \mathbf{u}}{\partial t} + \nabla \cdot \left( \rho \mathbf{u} \otimes \mathbf{u} + p \mathbf{I} \right) & = \mathbf{0} , \\ 
    \cfrac{\partial \alpha_k \rho_k e_k}{\partial t} + \nabla \cdot \left( \alpha_k \rho_k e_k \mathbf{u} \right) + \alpha_k p_k \nabla \cdot \mathbf{u} & = - p_I \delta p_k ,
  \end{array}} \right.\
  \end{equation}

where :math:`\alpha_k`, :math:`\rho_k`, :math:`p_k` and :math:`e_k` are the volume fraction, density, pressure and internal energy of each phase, respectively, and for which :math:`k` indicates the phase index.
The mixture density and pressure are

.. math::
  :nowrap:

  \begin{equation}
    \rho = \sum_{k=1}^N \alpha _k \rho_k  \quad \text{and} \quad p = \sum_{k=1}^N \alpha _k p_k ,
  \end{equation}

while the mixture total energy is 

.. math::
  :nowrap:

  \begin{equation}
    E = e + \frac{1}{2} \| \mathbf{u} \|^2,
  \end{equation}

where :math:`e` is the mixture specific internal energy

.. math::
  :nowrap:

  \begin{equation}
    e = \sum_{k=1}^N Y_k e_k \left( \rho_k , p_k \right) .
  \end{equation}

:math:`e_k \left( \rho_k , p_k \right)` is defined via an equation of state (EOS) and :math:`Y_k` are the mass fractions

.. math::
  :nowrap:

  \begin{equation}
    Y_k = \frac{\alpha_k \rho_k}{\rho} .
  \end{equation}

The relaxation of pressures between the phases is

.. math::
  :nowrap:

  \begin{equation}
    \delta p_k = \sum_{j \neq k}^N \mu_{k,j} \left( p_k - p_j \right) ,
  \end{equation}

where :math:`j` are phases different from :math:`k` and :math:`\mu_{k,j}` are the pressure-relaxation coefficients related to the :math:`k`--:math:`j` interactions. Herein, the pressure-relaxation coefficient :math:`\mu` is considered the same for each phase combination.
The interfacial pressure is defined as

.. math::
  :nowrap:

  \begin{equation}
    p_I = \cfrac{\sum_k^N \left( p_k \sum_{j \neq k}^N z_j \right)}{\sum_k^N z_k} ,
  \end{equation}

where :math:`z_k = \rho_k c_k` and :math:`c_k` are the acoustic impedance and speed of sound of the phase :math:`k`, respectively.

Since pressures are in disequilibrium here, the total energy equation of the mixture is replaced by the internal-energy equation for each phase. Nevertheless, conservation of the mixture total energy can be written in its usual form

.. math::
  :nowrap:

  \begin{equation}
    \frac{\partial \rho E}{\partial t} + \nabla \cdot \left[ \left( \rho E + p \right) \mathbf{u} \right] = 0 .
  \end{equation}

We note that this equation is redundant when the internal energy equations are also computed. However, in practice, we include it in our computations to ensure that the total energy is numerically conserved, and thus preserve a correct treatment of shock waves.

This model is solved thanks to the numerical method presented in :cite:`schmidmayer2021UEq` where infinite as well as finite pressure-relaxation rates are possible.

Tests cases
===========
Tests are provided with ECOGEN package and may be described in details later.

.. code-block:: xml

  <testCase>libTests/referenceTestCases/UEq/1D/transports/interfaceWaterAir/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/shockTubes/interfaceAirHelium/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/shockTubes/interfaceWaterAir/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/shockTubes/interfaceWaterAirNASG/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/shockTubes/epoxySpinel/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/shockTubes/mixtures/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/mixture/waterAir/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/shockOnInterface/sharpInterfaceWaterAir/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/shockOnInterface/diffusedInterfaceWaterAir/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/cavitation/</testCase>
  <testCase>libTests/referenceTestCases/UEq/1D/sphericalCollapse/Pratio1427/</testCase>
  <testCase>libTests/referenceTestCases/UEq/2D/nonSphericalCollapseNearWall/</testCase>
  <testCase>libTests/referenceTestCases/UEq/2D/squareToCircleSymmetry/</testCase>