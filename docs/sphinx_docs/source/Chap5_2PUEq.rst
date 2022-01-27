.. role:: xml(code)
  :language: xml

****************************************
Pressure- and velocity-equilibrium model
****************************************

Mechanical-equilibrium flows are solved in ECOGEN using the pressure-velocity-equilibrium model (previously named Kapila's model) :cite:`kapila2001`. In the particular case of 2 phases involved and without any extra physics (surface tension, viscosity...), this model reads:

.. math::
  :nowrap:

  \begin{equation}
  \label{system_PUEq}
  \left\{
  {\begin{array}{*{20}{l}}
    {\cfrac{{\partial {\alpha _1}}}{{\partial t}} + \mathbf{u} \cdot \nabla {\alpha _1}}&{ = K div( \mathbf{u} ),} \\ 
    {\cfrac{{\partial {\alpha _1}{\rho _1}}}{{\partial t}} +  div \left( {{\alpha _1}{\rho _1}\mathbf{u}} \right) } &{ = 0 ,} \\
    {\cfrac{{\partial {\alpha _2}{\rho _2}}}{{\partial t}} + div \left( {{\alpha _2}{\rho _2}\mathbf{u}} \right)}&{ = 0 ,} \\ 
    {\cfrac{{\partial \rho \mathbf{u}}}{{\partial t}} + div \left( {\rho \mathbf{u} \otimes \mathbf{u} + p \mathbf{I}} \right)}&{ = \mathbf{0} ,} \\ 
    {\cfrac{{\partial \rho E}}{{\partial t}} + div \left( {\left( {\rho E + p} \right) \mathbf{u}} \right)}&{ = 0 ,}
  \end{array}} \right.\
  \end{equation}

where subscripts :math:`1` and :math:`2` correspond to one of the two phases, respectively. :math:`\alpha_k` and :math:`\rho_k` are the volume fraction and density of phase :math:`k`. 

:math:`\rho = \sum\limits_{k} \alpha_k \rho_k`, :math:`\mathbf{u}`, :math:`p`, :math:`E = e + \cfrac{1}{2} \| \mathbf{u} \|^2` and :math:`e = \sum_k \alpha_k \rho_k e_k` are the mixture density, velocity, pressure, total energy and internal energy, respectively. 

The term :math:`K div (\mathbf{u})` accounts for the differences in the acoustic behavior of both phases or in other words, for the differences in expansion and compression of each phase in mixture regions. :math:`K` is given by:

.. math::
  :nowrap:

  \begin{equation*}
  K = \cfrac{\rho _2 s_2^2 - \rho _1 s_1^2}{\cfrac{\rho _2 s_2^2}{\alpha _2} + \cfrac{\rho _1 s_1^2}{\alpha _1}},
  \end{equation*}

:math:`s_k` being the speed of sound of phase :math:`k`.

This model is solved thanks to the numerical method presented in :cite:`relaxjcp`.

**Remark:**
This model can also be solved thanks to the numerical method presented in :cite:`schmidmayer2021UEq` (velocity-equilibrium model) with an infinite pressure relaxation.

.. toctree::

   ./Chap5_2PUEq_1D.rst
   ./Chap5_2PUEq_2D.rst
