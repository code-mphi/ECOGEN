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

#ifndef BOUNDCONDMOVINGWALLO2CARTESIAN_H
#define BOUNDCONDMOVINGWALLO2CARTESIAN_H

#include "BoundCondMovingWall.h"

class BoundCondMovingWallO2Cartesian : public BoundCondMovingWall
{
  public:
    BoundCondMovingWallO2Cartesian(const BoundCondMovingWallO2Cartesian& Source, const int& lvl = 0); //Copy ctor (useful for AMR)
    BoundCondMovingWallO2Cartesian(int numPhysique, tinyxml2::XMLElement* element, std::string fileName);
    BoundCondMovingWallO2Cartesian(int numPhysique);
    ~BoundCondMovingWallO2Cartesian() override;

    void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces) override;
    void allocateSlopes(int& /*allocateSlopeLocal*/) override;
    void computeSlopes(Prim type = vecPhases) override;
    void solveRiemann(double& dtMax,
                      Limiter& globalLimiter,
                      Limiter& interfaceLimiter,
                      Limiter& globalVolumeFractionLimiter,
                      Limiter& interfaceVolumeFractionLimiter,
                      Prim type = vecPhases) override;

    int whoAmI() const override { return TypeBC::MOVINGWALL; };

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
  private:
};

#endif // BOUNDCONDMOVINGWALLO2CARTESIAN_H
