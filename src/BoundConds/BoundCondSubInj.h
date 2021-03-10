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

#ifndef BOUNDCONDSUBINJ_H
#define BOUNDCONDSUBINJ_H

#include "BoundCond.h"

// This boundary condition is specific for Euler model

class BoundCondSubInj : public BoundCond
{
public:
  BoundCondSubInj(int numPhysique, tinyxml2::XMLElement* element, int& numberPhases, std::string fileName = "Unknown file");
  BoundCondSubInj(const BoundCondSubInj& Source, const int& lvl = 0); // Copy ctor (useful for AMR)
  ~BoundCondSubInj();

  virtual void createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces);
  virtual void solveRiemannBoundary(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax);

  virtual int whoAmI() const { return SUBINJ; };
  virtual void printInfo();

  // For AMR method
  virtual void creerCellInterfaceChild();  /*!< Create a child cell interface (uninitialized) */

private:
  double m_m0; //!< Inflow specific massflow
  double m_T0; //!< Inflow temperature
};

#endif