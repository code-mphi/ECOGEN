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

#include "GradientGreenGauss.h"

//****************************************************************************

GradientGreenGauss::GradientGreenGauss() 
{
  m_type = GG;
}

//****************************************************************************

GradientGreenGauss::~GradientGreenGauss() {}

//****************************************************************************

void GradientGreenGauss::addGradInterface(Coord &grad, CellInterface &cellInterface, Cell &cell, double const& faceValue)
{
  // Since face normal is constructed through mesh vertex there is no indication of outward direction.
  // Outward normal direction is defined according to cell center. 

  Coord vFaceToElt(0.);
  vFaceToElt.setFromSubtractedVectors(cellInterface.getFace()->getPos(), cell.getElement()->getPosition());

  double scalarProduct(Coord::scalarProduct(cellInterface.getFace()->getNormal(), vFaceToElt));

  if (scalarProduct < 0) {
    grad.setXYZ(grad.getX() + faceValue * cellInterface.getFace()->getSurface() * cellInterface.getFace()->getNormal().getX(), 
                grad.getY() + faceValue * cellInterface.getFace()->getSurface() * cellInterface.getFace()->getNormal().getY(),
                grad.getZ() + faceValue * cellInterface.getFace()->getSurface() * cellInterface.getFace()->getNormal().getZ());
  }
  else {
    grad.setXYZ(grad.getX() - faceValue * cellInterface.getFace()->getSurface() * cellInterface.getFace()->getNormal().getX(), 
                grad.getY() - faceValue * cellInterface.getFace()->getSurface() * cellInterface.getFace()->getNormal().getY(),
                grad.getZ() - faceValue * cellInterface.getFace()->getSurface() * cellInterface.getFace()->getNormal().getZ());
  }
}

//****************************************************************************

Coord GradientGreenGauss::computeGradient(Cell* cell, Variable nameVariable, int numPhase)
{
// Compute gradient of cell using Green-Gauss method
  int typeCellInterface(0);
  double wl(0.), wr(0.), wf(0.);
  double dl(0.), dr(0.);
  Coord grad(0.);

  for (int b = 0; b < cell->getCellInterfacesSize(); b++) {
    if (!cell->getCellInterface(b)->getSplit()) 
    {
      typeCellInterface = cell->getCellInterface(b)->whoAmI();
      if (typeCellInterface == 0) // Inner interface (O1 or O2)
      {
        // Distances
        dl = cell->getCellInterface(b)->distance(cell->getCellInterface(b)->getCellGauche());
        dr = cell->getCellInterface(b)->distance(cell->getCellInterface(b)->getCellDroite());

        // Extracting left and right variables values
        wl = cell->getCellInterface(b)->getCellGauche()->selectScalar(nameVariable, numPhase);
        wr = cell->getCellInterface(b)->getCellDroite()->selectScalar(nameVariable, numPhase);

        // Interpolation on face using barycenter
        wf = (wl * dr + wr * dl) / (dl + dr);

        // Build gradient of the face
        this->addGradInterface(grad, *cell->getCellInterface(b), *cell, wf);
      }
      else // Boundary
      {
        // Extracting left variable value
        wl = cell->getCellInterface(b)->getCellGauche()->selectScalar(nameVariable, numPhase);

        // For gradients other than velocity the face value is equals to the cell value 
        // Hence the gradient between the face and the center of the cell is null leading to adiatic condition for temperature for example.
        // For velocity it is assumed that the flow is viscous leading to null value on a wall boundary for all components.

        if (nameVariable == velocityU || nameVariable == velocityV || nameVariable == velocityW) 
        {
          switch (typeCellInterface) 
          {
            case SYMMETRY:
              // Multiplication of the gradient by the normal direction to guarantee symmetry
              if (nameVariable == velocityU) { wf = wl * (1. - std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getX())); }
              if (nameVariable == velocityV) { wf = wl * (1. - std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getY())); }
              if (nameVariable == velocityW) { wf = wl * (1. - std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getZ())); }
              break;

            case WALL:
                if (!cell->getCellInterface(b)->isMRFWall()) wf = 0.;
                else {
                  Coord velocityWall(0.);
                  velocityWall = cell->getCellInterface(b)->getWallRotationalVelocityMRF().cross(cell->getCellInterface(b)->getFace()->getPos());
                  if (nameVariable == velocityU) { wf = velocityWall.getX(); }
                  if (nameVariable == velocityV) { wf = velocityWall.getY(); }
                  if (nameVariable == velocityW) { wf = velocityWall.getZ(); }
                }
                break;

            default:
              wf = wl; // Null gradient on this interface
          }
        }
        else {
          wf = wl; // Null gradient on this interface
        }

        // Build gradient of the face
        this->addGradInterface(grad, *cell->getCellInterface(b), *cell, wf);
      }
    }
  }
  
  // Cell gradient average
  grad /= cell->getElement()->getVolume();

  return grad;
}

//****************************************************************************

void GradientGreenGauss::computeGradient(Cell* cell, std::vector<Coord>& grads, std::vector<Variable>& nameVariables, std::vector<int>& numPhases)
{
  // Compute gradient of cell using Green-Gauss method
  int typeCellInterface(0);
  double wl(0.), wr(0.), wf(0.);
  double dl(0.), dr(0.);
  for (unsigned int g = 0; g < grads.size(); g++) { grads[g] = 0.; }

  for (int b = 0; b < cell->getCellInterfacesSize(); b++) {
    if (!cell->getCellInterface(b)->getSplit()) 
    {
      typeCellInterface = cell->getCellInterface(b)->whoAmI();
      if (typeCellInterface == 0) // Inner interface (O1 or O2)
      {
        // Distances
        dl = cell->getCellInterface(b)->distance(cell->getCellInterface(b)->getCellGauche());
        dr = cell->getCellInterface(b)->distance(cell->getCellInterface(b)->getCellDroite());

        for (unsigned int g = 0; g < grads.size(); g++) 
        {
          // Extracting left and right variables values
          wl = cell->getCellInterface(b)->getCellGauche()->selectScalar(nameVariables[g], numPhases[g]);
          wr = cell->getCellInterface(b)->getCellDroite()->selectScalar(nameVariables[g], numPhases[g]);

          // Interpolation on face using barycenter
          wf = (wl * dr + wr * dl) / (dl + dr);

          // Build gradient of the face
          this->addGradInterface(grads[g], *cell->getCellInterface(b), *cell, wf);
        }
      }
      else // Boundary
      {
        for (unsigned int g = 0; g < grads.size(); g++) 
        { 
          // Extracting left variable value
          wl = cell->getCellInterface(b)->getCellGauche()->selectScalar(nameVariables[g], numPhases[g]);

          // For gradients other than velocity the face value is equals to the cell value 
          // Hence the gradient between the face and the center of the cell is null leading to adiatic condition for temperature for example.
          // For velocity it is assumed that the flow is viscous leading to null value on a wall boundary for all components.

          if (nameVariables[g] == velocityU || nameVariables[g] == velocityV || nameVariables[g] == velocityW) 
          {
            switch (typeCellInterface) 
            {
              case SYMMETRY:
                // Multiplication of the gradient by the normal direction to guarantee symmetry
                if (nameVariables[g] == velocityU) { wf = wl * (1. - std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getX())); }
                if (nameVariables[g] == velocityV) { wf = wl * (1. - std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getY())); }
                if (nameVariables[g] == velocityW) { wf = wl * (1. - std::fabs(cell->getCellInterface(b)->getFace()->getNormal().getZ())); }
                break;

              case WALL:
                if (!cell->getCellInterface(b)->isMRFWall()) wf = 0.;
                else {
                  Coord velocityWall(0.);
                  velocityWall = cell->getCellInterface(b)->getWallRotationalVelocityMRF().cross(cell->getCellInterface(b)->getFace()->getPos());
                  if (nameVariables[g] == velocityU) { wf = velocityWall.getX(); }
                  if (nameVariables[g] == velocityV) { wf = velocityWall.getY(); }
                  if (nameVariables[g] == velocityW) { wf = velocityWall.getZ(); }
                }
                break;

              default:
                wf = wl; // Null gradient on this interface
            }
          }
          else {
            wf = wl; // Null gradient on this interface
          }

          // Build gradient of the face
          this->addGradInterface(grads[g], *cell->getCellInterface(b), *cell, wf);
        }
      }
    }
  }
  
  // Cell gradient average
  for (unsigned int g = 0; g < grads.size(); g++) 
  {
    grads[g] /= cell->getElement()->getVolume();
  }
}