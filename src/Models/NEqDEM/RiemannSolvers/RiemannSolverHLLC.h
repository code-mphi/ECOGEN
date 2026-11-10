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

#ifndef RIEMANNSOLVERHLLC_H
#define RIEMANNSOLVERHLLC_H

#include "RiemannSolver.h"

//! \class     RiemannSolverHLLC
//! \brief     RiemannSolverHLLC class for the non-equilibrium system of equations

class RiemannSolverHLLC : public RiemannSolver
{
  public:
    //! \brief  constructor of the Riemann Solver HLLC
    RiemannSolverHLLC();
    //! \brief  destructor of the Riemann Solver HLLC
    ~RiemannSolverHLLC() override;
    //! \brief  call the HLLC Riemann solver to solve Riemann Problem at the interface of two cells
    void solver(const Phase& phaseLeft,
                const Phase& phaseRight,
                const double& dxLeft,
                const double& dxRight,
                double& dtMax,
                double& pStar,
                FluxEuler& m_fluxEuler) override;
    //! \brief  call the HLLC Riemann solver to solve Riemann Problem at a wall BC
    void solverWall(const Phase& phaseLeft, const double& dxLeft, double& dtMax, double& pStar, FluxEuler& m_fluxEuler) override;
};

#endif // RIEMANNSOLVERHLLC_H
