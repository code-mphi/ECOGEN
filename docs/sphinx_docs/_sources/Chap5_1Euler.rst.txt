.. role:: xml(code)
  :language: xml

***********************
Single-phase test cases
***********************

Test cases presented in this section are dealing with single-phase compressible problems. In this part, ECOGEN solves the Euler equations :cite:`euler1757principes`:

.. math::
  :nowrap:
  
  \[
  \label{eqEuler}
  \begin{array}{l}
    \displaystyle \frac{\partial\rho}{\partial t}+ div(\rho {\mathbf{u}} )=0\\
    \displaystyle\frac{\partial\rho {\mathbf{u}}}{\partial t}+ div(\rho {\mathbf{u}} \otimes {\mathbf{u}} +p \mathbf{I})=\mathbf{0} \\ 
    \displaystyle \frac{\partial\rho E}{\partial t}+ div((\rho E+p){\mathbf{u}})=0
  \end{array} 
  \]

where :math:`\rho` represents the density, :math:`\mathbf{u}` the velocity vector, :math:`p` the pressure and :math:`E = e +\frac{\mathbf{u}^2}{2}` the total energy, with :math:`e` the internal energy. 
The closure relation for this model is ensured by any convex equation of state (EOS) :math:`p = p(\rho,e)` (see section :ref:`Sec:IO:materials` for details about implemented EOS in ECOGEN).
Euler equations are solved thanks to an explicit finite-volume Godunov-like scheme :cite:`godunov79` that is coupled with HLLC approximate Riemann solver :cite:`toro2013riemann` for flux computation.


.. toctree::
   :maxdepth: 10
   :caption: Contents:

   ./Chap5_1Euler_1D.rst
   ./Chap5_1Euler_2D.rst
   ./Chap5_1Euler_3D.rst
.. _Paraview: https://www.paraview.org/
.. _gnuplot: http://www.gnuplot.info/