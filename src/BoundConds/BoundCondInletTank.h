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

#ifndef BOUNDCONDINLETTANK_H
#define BOUNDCONDINLETTANK_H

#include "BoundCond.h"

class BoundCondInletTank : public BoundCond
{
  public:
    BoundCondInletTank(int numPhysique, tinyxml2::XMLElement* element, const int& numbPhases, const int& numbTransports, std::vector<std::string> nameTransports, Eos** eos, std::string fileName = "Unknown file");
    BoundCondInletTank(const BoundCondInletTank &Source, const int& lvl = 0); //Copy ctor (useful for AMR)
    virtual ~BoundCondInletTank();

    virtual void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces);
    virtual void solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax);
    virtual void solveRiemannTransportBoundary(Cell& cellLeft) const;

    virtual int whoAmI() const { return INLETTANK; };
    virtual void printInfo();

    //For AMR method
    virtual void creerCellInterfaceChild();  /*!< Create a child cell interface (non initialize) */

  protected:
  private:
    double* m_ak0;
    double* m_Yk0;
    double* m_rhok0;
    double m_p0;
    double m_T0;
    double* m_valueTransport;
};

#endif // BOUNDCONDINLETTANK_H
