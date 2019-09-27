Introduction
============

ECOGEN is a CFD plateform written in C++ object oriented programming langage. It is dedicated to numerical simulation of compressible multiphase flows. It has the vocation to share academics researches in the multiphase flow field in direction to ohter academics but also for industrials, students, etc.

- multi-models (single phase, multiphase with or without equilibrium)
- multi-physics (thermal transfers, viscosity, surface tension, mass transfers)
- multi-meshes (Cartesian, unstructured, AMR)
- multi-CPU

ECOGEN stands for:

- **E**\ volutive: makes easier future developpements
- **C**\ ompressible: dedicated to compressible flows
- **O**\ penSource: Distributed under GPLv3 Licence
- **G**\ enuine: Uses the "Diffuse Interface Method" (DIM)
- **E**\asy: simple to install and use (C++ compiler and MPI)
- **N**-phase: liquids, vapors, inert gases and/or reactives 

What kind of physical problems ?
--------------------------------

ECOGEN is designed for following applications:

- Solving interface problems between pure or multicomponents fluids and mixtures of multiphase flows.
- Treating surface tension, heat and mass transfers for evaporating and condensing flow, cavitation.
- Computing wave propagation in strongly unsteady situations using a specific Adaptative Mesh Reffinement .
- Computing on unstructured grids to simulate complex geometries.
- Parallel computing using open MPI libraries.

What about the engine ?
-----------------------

ECOGEN is a receptacle of a story of diffuse interface method (DIM) theory that started in the late 90s. DIM summarizes more than 20 years of researches on multiphase flow modelling with the objective to develop mathematical models as well as their associated numerical methods.

What is a diffuse interface? In DIM theory, the interfaces between pure phases are captured as diffuse numerical zones meaning that one goes continuously from one phase to another.

.. _Fig:introduction:diffInterface:

.. figure:: ./_static/intro/diffInterface.png

  1D extraction of the diffuse interface zone of a water droplet.

This way is possible thanks to a thermodynamical consistency. Then, the flow solution does no longer requires interface tracking algorithms: It became easy to simulate complex topological shape evolutions between miscible or non miscible fluid. Moreover, pressure waves (shock waves, acoustic waves) can propagate and interact properly in the whole flow.

The basic research on DIM is now matured enough to propose ECOGEN, a numerical tool that can be largely cast and use to solve industrial as well as research multiphase flow problems.

