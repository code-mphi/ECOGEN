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

#ifndef BOUNDCONDINJ_H
#define BOUNDCONDINJ_H

#include "BoundCond.h"


class BoundCondInj : public BoundCond
{
  public:
    BoundCondInj(int numPhysique, tinyxml2::XMLElement* element, int& numberPhases, int& numberTransports, std::vector<std::string> nameTransports, Eos** eos, std::string fileName = "Unknown file");
    BoundCondInj(const BoundCondInj &Source, const int& lvl = 0); //Copy ctor (useful for AMR)
    virtual ~BoundCondInj();

    virtual void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces);
    virtual void solveRiemannBoundary(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax);
    virtual void solveRiemannTransportBoundary(Cell& cellLeft, const int& numberTransports) const;

    virtual int whoAmI() const { return INJ; };
    virtual void printInfo();

    //For AMR method
    virtual void creerCellInterfaceChild();  //!< Create a child cell interface (not initialized) */

  protected:
  private:
    int m_numberPhases;
    double m_m0; //!< Specific massflow (kg.s-1.m-2)
    double m_T0;
    double* m_ak0;
    double* m_rhok0;
    double* m_pk0;
    int m_numberTransports;
    double* m_valueTransport;
};

#endif // BOUNDCONDINJ_H
