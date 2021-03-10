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

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APEConductivity.h"

//***********************************************************************

APEConductivity::APEConductivity(int& numberQPA, Eos** eos)
{
  m_lambda = eos[0]->getLambda();
  m_numQPA = numberQPA++;
}

//***********************************************************************

APEConductivity::~APEConductivity() {}

//***********************************************************************

void APEConductivity::addQuantityAddPhys(Cell* cell)
{
  cell->getVecQuantitiesAddPhys().push_back(new QAPConductivity(this, cell->getNumberPhases()));
}

//***********************************************************************

void APEConductivity::solveFluxAddPhys(CellInterface* cellInterface, const int& /*numberPhases*/)
{
  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Copy and projection on orientation axis attached to the edge of gradients of left and right cells
  m_gradTLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(0);
  m_gradTRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGrad(0);

  m_gradTLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradTRight.localProjection(m_normal, m_tangent, m_binormal);

  this->solveFluxConductivityInner(m_gradTLeft, m_gradTRight);

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APEConductivity::solveFluxAddPhysBoundary(CellInterface* cellInterface, const int& /*numberPhases*/)
{
  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
  m_gradTLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(0);
  m_gradTLeft.localProjection(m_normal, m_tangent, m_binormal);

  int typeCellInterface(cellInterface->whoAmI());
  if (typeCellInterface == NONREFLECTING) {
    this->solveFluxConductivityNonReflecting(m_gradTLeft);
  }
  else if (typeCellInterface == WALL) {
    int heatWalltype(cellInterface->whoAmIHeat());
    if (heatWalltype == IMPOSEDTEMP) {
      this->solveFluxConductivityWallImposedTemp(cellInterface);
    }
    else if (heatWalltype == IMPOSEDFLUX) {
      this->solveFluxConductivityWallImposedFlux(cellInterface);
    }
    else {
      this->solveFluxConductivityOther(); // Adiabatic wall, i.e. null flux
    }
  }
  else { 
    this->solveFluxConductivityOther();
  }
  // etc... Boundaries not fully taken into account yet for conductivity, pay attention

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APEConductivity::solveFluxConductivityInner(const Coord& gradTLeft, const Coord& gradTRight) const
{
  //Data of the cell interface
  double dTdx;
  dTdx = (gradTLeft.getX() + gradTRight.getX()) / 2.;

  // Compute heat fluxes
  static_cast<FluxEuler*> (fluxBuff)->m_masse = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_energ = - m_lambda * dTdx;
}

//***********************************************************************

void APEConductivity::solveFluxConductivityNonReflecting(const Coord& gradTLeft) const
{
  this->solveFluxConductivityInner(gradTLeft, gradTLeft);
}

//***********************************************************************

void APEConductivity::solveFluxConductivityWallImposedTemp(CellInterface *cellInterface)
{
  // Retrieve local temperatures
  double temperatureLeft(cellInterface->getCellGauche()->getPhase(0)->getTemperature());
  double tempWall(cellInterface->getBoundaryHeatQuantity());

  // Distances cell center/cell interface
  double distLeft(cellInterface->getCellGauche()->distance(cellInterface));

  // Compute the temperature gradient locally
  double dTdx((tempWall - temperatureLeft) / distLeft);

  // Compute heat fluxes
  static_cast<FluxEuler*> (fluxBuff)->m_masse = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_energ = - m_lambda * dTdx;
}

//***********************************************************************

void APEConductivity::solveFluxConductivityWallImposedFlux(CellInterface* cellInterface)
{
  // Retrieve local imposed flux
  double flux(cellInterface->getBoundaryHeatQuantity());

  // Compute heat fluxes
  static_cast<FluxEuler*> (fluxBuff)->m_masse = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_energ = flux;
}

//***********************************************************************

void APEConductivity::solveFluxConductivityOther() const
{
  // To avoid bug when not managed or for adiabatic wall
  static_cast<FluxEuler*> (fluxBuff)->m_masse = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxEuler*> (fluxBuff)->m_energ = 0.;
}

//***********************************************************************

void APEConductivity::communicationsAddPhys(const int& /*numberPhases*/, const int& dim, const int& lvl)
{
  parallel.communicationsVector(QPA, dim, lvl, m_numQPA, 0); //m_gradT
}

//***********************************************************************