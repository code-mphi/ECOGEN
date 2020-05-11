Introduction
============

ECOGEN is a CFD

- multi-model (single phase, multiphase with or without equilibriums),
- multi-physics (thermal transfers, viscosity, surface tension, mass transfers),
- multi-mesh (Cartesian, unstructured, AMR),
- multi-core,

plateform written in C++ object-oriented-programming language. It is dedicated to numerical simulation of compressible multiphase flows. It has the vocation to share academics research in the multiphase flow field in direction to other academics but also for industrials, students, etc.

ECOGEN stands for:

- **E**\ volutive: Made for easier future developments.
- **C**\ ompressible: Dedicated to compressible flows.
- **O**\ pen-source: Distributed under `GNU GPLv3 License`_.
- **G**\ enuine: Uses the "Diffuse Interface Method" (DIM).
- **E**\ asy: Simple to install and use (C++ compiler and MPI).
- **N**\ -phase: Liquids, vapors, inert gases and/or reactives.

.. _`GNU GPLv3 License`: http://www.gnu.org/licenses

What kind of physical problems?
-------------------------------

ECOGEN is designed for the following:

- Solving interface problems between pure or multicomponent fluids and mixtures of multiphase flows.
- Treating surface tension, heat and mass transfers for cavitation, evaporating and condensing flows.
- Computing wave propagation in strongly unsteady situations using a specific Adaptive Mesh Refinement.
- Computing on unstructured grids to simulate complex geometries.
- Parallel computing using open MPI libraries.

What about the engine?
----------------------

ECOGEN is a receptacle of a story of diffuse-interface-method (DIM) theory that started in the late 90s. DIM summarizes more than 20 years of researches on multiphase flow modelling with the goal to develop mathematical models as well as their associated numerical methods.

What is a diffuse interface? In DIM theory, the interfaces between pure phases are captured as diffuse numerical regions, meaning that one goes continuously from one phase to another.

.. _Fig:introduction:diffInterface:

.. figure:: ./_static/intro/diffInterface.png

  1D extraction of the diffuse interface region of a water droplet.

This way is possible thanks to a thermodynamics consistency. Then, the flow solution does no longer requires interface tracking algorithms: It becomes easy to simulate complex topological-shape evolutions between miscible or non-miscible fluids. Moreover, pressure waves (shock waves, acoustic waves) can propagate and interact properly in the whole flow.

The base research on DIM is now matured enough to propose ECOGEN, a numerical tool that can be largely cast and use to solve industrial as well as research multiphase flow problems.
