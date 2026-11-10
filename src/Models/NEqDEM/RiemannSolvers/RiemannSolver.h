//
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-.
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| |
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | |
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  |
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)|
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_)
//      (__)              (_)      (__)     (__)     (__)
//      Official webSite: https://code-mphi.github.io/ECOGEN/
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names
//  are listed in the copyright file included with this source
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published
//  by the Free Software Foundation, either version 3 of the License,
//  or (at your option) any later version.
//
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

#include "../../Phase.h"
#include "../../Euler/FluxEuler.h"

//! \class     RiemannSolver
//! \brief     RiemannSolver abstract class for the non-equilibrium system of equations

class RiemannSolver
{
  public:
    //! \brief  constructor of the Riemann Solver
    RiemannSolver();
    //! \brief  virtual destructor of the Riemann Solver
    virtual ~RiemannSolver();
    //! \brief  call the solver to solve the Riemann problem, either HLLC or Exact
    virtual void solver(const Phase& /*phaseLeft*/,
                        const Phase& /*phaseRight*/,
                        const double& /*dxLeft*/,
                        const double& /*dxRight*/,
                        double& /*dtMax*/,
                        double& /*pStar*/,
                        FluxEuler& /*m_fluxEuler*/)
    {
      Errors::errorMessage("solver not available for required RiemannSolver");
    };
    //! \brief  call the solver to solve the Riemann problem at a wall BC, HLLC (Exact not implemented)
    virtual void solverWall(const Phase& /*phaseLeft*/, const double& /*dxLeft*/, double& /*dtMax*/, double& /*pStar*/, FluxEuler& /*m_fluxEuler*/)
    {
      Errors::errorMessage("solverWall not available for required RiemannSolver");
    };
};

#endif // RIEMANNSOLVER_H
