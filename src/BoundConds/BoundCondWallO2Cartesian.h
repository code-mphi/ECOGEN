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

#ifndef BOUNDCONDWALLO2CARTESIAN_H
#define BOUNDCONDWALLO2CARTESIAN_H

#include "BoundCondWall.h"


class BoundCondWallO2Cartesian : public BoundCondWall
{
  public:
    BoundCondWallO2Cartesian(const BoundCondWallO2Cartesian& Source, const int& lvl = 0); //Copy ctor (useful for AMR)
    BoundCondWallO2Cartesian(int numPhysique, tinyxml2::XMLElement* element, std::string fileName);
    BoundCondWallO2Cartesian(int numPhysique);
    virtual ~BoundCondWallO2Cartesian();

    virtual void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces);
    virtual void allocateSlopes(int& /*allocateSlopeLocal*/);
    virtual void computeSlopes(Prim type = vecPhases);
    virtual void solveRiemann(double& dtMax, Limiter& globalLimiter, Limiter& interfaceLimiter, Limiter& globalVolumeFractionLimiter, Limiter& interfaceVolumeFractionLimiter, Prim type = vecPhases);

    virtual int whoAmI() const { return WALL; };

    //Accessors
    virtual Phase* getSlopesPhase(const int& phaseNumber) const;
    virtual Mixture* getSlopesMixture() const;
    virtual Transport* getSlopesTransport(const int& numberTransport) const;

    //For AMR method
    virtual void creerCellInterfaceChild();  /*!< Create a child cell interface (not initialized) */

  protected:
    Phase** m_vecPhasesSlopes;         /*!< vector des slopes des phases */
    Mixture* m_mixtureSlopes;          /*!< vector des slopes de mixture */
    Transport* m_vecTransportsSlopes;  /*!< vector des slopes des transports */
  private:
};

#endif // BOUNDCONDWALLO2CARTESIAN_H
