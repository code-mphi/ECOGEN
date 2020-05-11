//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
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

//! \file      APKConductivity.cpp
//! \author    K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APKConductivity.h"

//***********************************************************************

APKConductivity::APKConductivity() {}

//***********************************************************************

APKConductivity::APKConductivity(int& numberQPA, Eos** eos, int &numberPhases, std::string nameFile)
{
  m_lambdak = new double[numberPhases];
  for (int k = 0; k < numberPhases; k++) {
    m_lambdak[k] = eos[k]->getLambda();
  }
  m_numQPA = numberQPA++;
}

//***********************************************************************

APKConductivity::~APKConductivity(){ delete[] m_lambdak; }

//***********************************************************************

void APKConductivity::addQuantityAddPhys(Cell *cell)
{
  cell->getVecQuantitiesAddPhys().push_back(new QAPConductivity(this,cell->getNumberPhases()));
}

//***********************************************************************

void APKConductivity::solveFluxAddPhys(CellInterface *cellInterface, const int &numberPhases)
{
  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Reset of fluxBuffKapila
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxKapila*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = 0.;

  for (int numPhase = 0; numPhase < numberPhases; numPhase++) {
    // Copy and projection on orientation axis attached to the edge of gradients of left and right cells
    m_gradTkLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(numPhase);
    m_gradTkRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGrad(numPhase);
    m_gradTkLeft.localProjection(m_normal, m_tangent, m_binormal);
    m_gradTkRight.localProjection(m_normal, m_tangent, m_binormal);

    // Extraction of alphak
    double alphakLeft = cellInterface->getCellGauche()->getPhase(numPhase)->getAlpha();
    double alphakRight = cellInterface->getCellDroite()->getPhase(numPhase)->getAlpha();

    this->solveFluxConductivityInner(m_gradTkLeft, m_gradTkRight, alphakLeft, alphakRight, numPhase);
  }

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKConductivity::solveFluxAddPhysBoundary(CellInterface *cellInterface, const int &numberPhases)
{
  //KS//DEV// On ne fait rien aux limites avec la conductivite pour le moment, a gerer un jour

  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Reset of fluxBuffKapila (allow to then do the sum of conductivity effects for the different phases combinations)
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxKapila*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = 0.;

  for (int numPhase = 0; numPhase < numberPhases; numPhase++) {
    // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
    m_gradTkLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(numPhase);
    m_gradTkLeft.localProjection(m_normal, m_tangent, m_binormal);

    // Extraction of alphak
    double alphakLeft = cellInterface->getCellGauche()->getPhase(numPhase)->getAlpha();

    int typeCellInterface = cellInterface->whoAmI();
    if (typeCellInterface == 1) { this->solveFluxConductivityNonReflecting(m_gradTkLeft, alphakLeft, numPhase); }
    else if (typeCellInterface == 2 || typeCellInterface == 6) { this->solveFluxConductivityWall(m_gradTkLeft, alphakLeft, numPhase); }
    else if (typeCellInterface == 3) { this->solveFluxConductivityOutflow(m_gradTkLeft, alphakLeft, numPhase); }
    else if (typeCellInterface == 4) { this->solveFluxConductivityInflow(m_gradTkLeft, alphakLeft, numPhase); }
    else { this->solveFluxConductivityOther(m_gradTkLeft, alphakLeft, numPhase); }
    // etc... Boundaries not taken into account yet for conductivity, pay attention
  }

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKConductivity::solveFluxConductivityInner(Coord &gradTkLeft, Coord &gradTkRight, double &alphakL, double &alphakR, int &numPhase) const
{
  //Data of the cell interface
  double dTkdx, alphak;
  dTkdx = (gradTkLeft.getX() + gradTkRight.getX()) / 2.;
  alphak = (alphakL + alphakR) / 2.;

  //Writing of conductive terms on each equation of fluxBuffKapila
  static_cast<FluxKapila*> (fluxBuff)->m_energ[numPhase] = -alphak*m_lambdak[numPhase] * dTkdx;
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture   += -alphak*m_lambdak[numPhase] * dTkdx;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityNonReflecting(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  this->solveFluxConductivityInner(gradTkLeft, gradTkLeft, alphakL, alphakL, numPhase);
}

//***********************************************************************

void APKConductivity::solveFluxConductivityWall(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example
  //KS//DEV// A faire pour couche limite thermique !!! ...

  // To avoid bug when not manage
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityOutflow(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityInflow(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityOther(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example
  std::cout << "Conductive boundary not manage" << std::endl;

  // To avoid bug when not manage
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::communicationsAddPhys(int numberPhases, const int &dim, const int &lvl)
{
  for (int k = 0; k < numberPhases; k++) {
		parallel.communicationsVector(QPA, dim, lvl, m_numQPA, k);
	}
}

//***********************************************************************
