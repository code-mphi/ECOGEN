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

#ifndef BOUNDCONDPERIODICO2CARTESIAN_H
#define BOUNDCONDPERIODICO2CARTESIAN_H

#include "BoundCondPeriodic.h"

class BoundCondPeriodicO2Cartesian : public BoundCondPeriodic
{
  public:
    BoundCondPeriodicO2Cartesian(const BoundCondPeriodicO2Cartesian& Source, const int& lvl = 0); //Copy ctor (useful for AMR)
    BoundCondPeriodicO2Cartesian(int numPhysique);
    ~BoundCondPeriodicO2Cartesian() override;

    void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces) override;
    void allocateSlopes(int& /*allocateSlopeLocal*/) override;
    void computeSlopes(Prim /*type*/ = vecPhases) override;
    /*!< Specific Riemann problem for 2nd order */
    void solveRiemann(double& dtMax,
                      Limiter& globalLimiter,
                      Limiter& interfaceLimiter,
                      Limiter& globalVolumeFractionLimiter,
                      Limiter& interfaceVolumeFractionLimiter,
                      Prim type) override;

    int whoAmI() const override { return TypeBC::PERIODIC; };

    //Accessors
    Phase* getSlopesPhase(const int& phaseNumber) const override;
    Mixture* getSlopesMixture() const override;
    Transport* getSlopesTransport(const int& numberTransport) const override;

    //For AMR method
    void creerCellInterfaceChild() override; /*!< Create a child cell interface (not initialized) */

  protected:
    Phase** m_vecPhasesSlopes;        /*!< vector des slopes des phases */
    Mixture* m_mixtureSlopes;         /*!< vector des slopes de mixture */
    Transport* m_vecTransportsSlopes; /*!< vector des slopes des transports */
};

#endif // BOUNDCONDPERIODICO2CARTESIAN_H
