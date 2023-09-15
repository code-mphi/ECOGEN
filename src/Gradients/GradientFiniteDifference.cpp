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

#include "GradientFiniteDifference.h"

//****************************************************************************

GradientFiniteDifference::GradientFiniteDifference()
{
  m_type = TypeGrad::FD;
}

//****************************************************************************

GradientFiniteDifference::~GradientFiniteDifference() {}

//****************************************************************************

void GradientFiniteDifference::computeGradients(Cell* cell, std::vector<Coord>& grads, const std::vector<Variable>& nameVariables, const std::vector<int>& numPhases)
{
  int typeCellInterface(0);
  double cg(0.), cd(0.), gradCellInterface(0.);
  double distance(0.), distanceX(0.), distanceY(0.), distanceZ(0.);
  double sumDistanceX(0.), sumDistanceY(0.), sumDistanceZ(0.);
  for (unsigned int g = 0; g < grads.size(); g++) { grads[g] = 0.; }

  for (int b = 0; b < cell->getCellInterfacesSize(); b++) {
    if (!cell->getCellInterface(b)->getSplit()) {
      typeCellInterface = cell->getCellInterface(b)->whoAmI();
      if (typeCellInterface == 0) //Cell interface of type CellInterface/O2 (inner)
      {
        // Sum for each cell interface with ponderation using distance in each direction
        // then the cell gradient is normalized by sum of distances.
        distanceX = std::fabs(cell->getCellInterface(b)->getCellLeft()->distanceX(cell->getCellInterface(b)->getCellRight()));
        distanceY = std::fabs(cell->getCellInterface(b)->getCellLeft()->distanceY(cell->getCellInterface(b)->getCellRight()));
        distanceZ = std::fabs(cell->getCellInterface(b)->getCellLeft()->distanceZ(cell->getCellInterface(b)->getCellRight()));
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = cell->getCellInterface(b)->getCellLeft()->distance(cell->getCellInterface(b)->getCellRight());

        // Extracting left and right variables values for each cell interface
        // and calculus of the gradients normal to the face
        for (unsigned int g = 0; g < grads.size(); g++) {
          cg = cell->getCellInterface(b)->getCellLeft()->selectScalar(nameVariables[g], numPhases[g]);
          cd = cell->getCellInterface(b)->getCellRight()->selectScalar(nameVariables[g], numPhases[g]);
          gradCellInterface = (cd - cg) / distance;

          // Projection in the absolute system of coordinate
          grads[g].setXYZ(grads[g].getX() + cell->getCellInterface(b)->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                          grads[g].getY() + cell->getCellInterface(b)->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                          grads[g].getZ() + cell->getCellInterface(b)->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
        }
      }
      else {
        distanceX = std::fabs(cell->distanceX(cell->getCellInterface(b))) * 2.;
        distanceY = std::fabs(cell->distanceY(cell->getCellInterface(b))) * 2.;
        distanceZ = std::fabs(cell->distanceZ(cell->getCellInterface(b))) * 2.;
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = cell->distance(cell->getCellInterface(b));
          
        for (unsigned int g = 0; g < grads.size(); g++) {
          if (nameVariables[g] == Variable::velocityU || nameVariables[g] == Variable::velocityV || nameVariables[g] == Variable::velocityW) {
            // Extracting left variables values
            // and calculus of the gradient normal to the face
            cg = cell->getCellInterface(b)->getCellLeft()->selectScalar(nameVariables[g], numPhases[g]);
            gradCellInterface = - cg / distance;

            if (typeCellInterface == SYMMETRY) {
              // Multiplication of the gradient by the normal direction to guarantee symmetry
              if (nameVariables[g] == Variable::velocityU) { gradCellInterface = gradCellInterface * std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getX()); }
              if (nameVariables[g] == Variable::velocityV) { gradCellInterface = gradCellInterface * std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getY()); }
              if (nameVariables[g] == Variable::velocityW) { gradCellInterface = gradCellInterface * std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getZ()); }

              // Projection in the absolute system of coordinate
              grads[g].setXYZ(grads[g].getX() + cell->getCellInterface(b)->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                              grads[g].getY() + cell->getCellInterface(b)->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                              grads[g].getZ() + cell->getCellInterface(b)->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
            }
            else if (typeCellInterface == WALL) {
              // Projection in the absolute system of coordinate
              grads[g].setXYZ(grads[g].getX() + cell->getCellInterface(b)->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                              grads[g].getY() + cell->getCellInterface(b)->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                              grads[g].getZ() + cell->getCellInterface(b)->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
            }
          }
        }
      }
    }
  }

  // Verifications in multiD
  if (sumDistanceX <= 1.e-12) { sumDistanceX = 1.; }
  if (sumDistanceY <= 1.e-12) { sumDistanceY = 1.; }
  if (sumDistanceZ <= 1.e-12) { sumDistanceZ = 1.; }

  // Final normalized gradient on the cell
  for (unsigned int g = 0; g < grads.size(); g++) {
    grads[g].setXYZ(grads[g].getX() / sumDistanceX, grads[g].getY() / sumDistanceY, grads[g].getZ() / sumDistanceZ);
  }
}