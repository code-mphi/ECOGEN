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

#ifndef RIEMANNSOLVEREXACT_H
#define RIEMANNSOLVEREXACT_H

#include "RiemannSolver.h"

//! \class     RiemannSolverExact
//! \brief     RiemannSolverExact class for the non-equilibrium system of equations

class RiemannSolverExact : public RiemannSolver
{
  public:
    //! \brief  constructor of the exact Riemann Solver
    RiemannSolverExact();
    //! \brief  destructor of the exact Riemann Solver
    ~RiemannSolverExact() override;
    //! \details Solve a discontinuity problem between two cells by using an exact Riemann problem solver
    void solver(const Phase& phaseLeft,
                const Phase& phaseRight,
                const double& dxLeft,
                const double& dxRight,
                double& dtMax,
                double& pStar,
                FluxEuler& m_fluxEuler) override;
    //! \details Compute wave speed inside Newton algorithm in order to find star pression
    void computeWaveSpeed(
      double p0, double rho0, double soundSpeed0, double gamma, double pInf, double targetPressure, double& waveSpeed, double& waveSpeedDerivative);
    //! \details Compute a coefficient ratio according to the pressure ratio and gas constant
    void computeCoefficientRatio(double pressureRatio, double gamma, double& coefficientRatio);
};

#endif // RIEMANNSOLVEREXACT_H
