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
//! \version   1.0
//! \date      December 20 2017

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APKConductivity.h"

using namespace std;

//***********************************************************************

APKConductivity::APKConductivity() {}

//***********************************************************************

APKConductivity::APKConductivity(int& numberQPA, Eos** eos, int &numberPhases, string nameFile)
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

void APKConductivity::solveFluxAddPhys(CellInterface *cellBound, const int &numberPhases)
{
  m_normal = cellBound->getFace()->getNormal();
  m_tangent = cellBound->getFace()->getTangent();
  m_binormal = cellBound->getFace()->getBinormal();

  // Reset of fluxBufferKapila
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = 0.;
  }
  fluxBufferKapila->m_qdm = 0.;
  fluxBufferKapila->m_energMixture = 0.;

  for (int numPhase = 0; numPhase < numberPhases; numPhase++) {
    // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
    m_gradTkLeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(numPhase);
    m_gradTkRight = cellBound->getCellDroite()->getQPA(m_numQPA)->getGrad(numPhase);
    m_gradTkLeft.localProjection(m_normal, m_tangent, m_binormal);
    m_gradTkRight.localProjection(m_normal, m_tangent, m_binormal);

    // Extraction of alphak
    double alphakLeft = cellBound->getCellGauche()->getPhase(numPhase)->getAlpha();
    double alphakRight = cellBound->getCellDroite()->getPhase(numPhase)->getAlpha();

    this->solveFluxConductivityInner(m_gradTkLeft, m_gradTkRight, alphakLeft, alphakRight, numPhase);
  }

  // Flux projection on the absolute orientation axes
  cellBound->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKConductivity::solveFluxAddPhysBoundary(CellInterface *cellBound, const int &numberPhases)
{
  m_normal = cellBound->getFace()->getNormal();
  m_tangent = cellBound->getFace()->getTangent();
  m_binormal = cellBound->getFace()->getBinormal();

  // Reset of fluxBufferKapila (allow to then do the sum of conductivity effects for the different phases combinations)
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = 0.;
  }
  fluxBufferKapila->m_qdm = 0.;
  fluxBufferKapila->m_energMixture = 0.;

  for (int numPhase = 0; numPhase < numberPhases; numPhase++) {
    // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
    m_gradTkLeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(numPhase);
    m_gradTkLeft.localProjection(m_normal, m_tangent, m_binormal);

    // Extraction of alphak
    double alphakLeft = cellBound->getCellGauche()->getPhase(numPhase)->getAlpha();

    int typeBord = cellBound->whoAmI();
    if (typeBord == 1) { this->solveFluxConductivityAbs(m_gradTkLeft, alphakLeft, numPhase); }
    else if (typeBord == 2 || typeBord == 6) { this->solveFluxConductivityWall(m_gradTkLeft, alphakLeft, numPhase); }
    else if (typeBord == 3) { this->solveFluxConductivityOutflow(m_gradTkLeft, alphakLeft, numPhase); }
    else if (typeBord == 4) { this->solveFluxConductivityInflow(m_gradTkLeft, alphakLeft, numPhase); }
    else { this->solveFluxConductivityOther(m_gradTkLeft, alphakLeft, numPhase); }
    // etc... Boundaries not taken into account yet for conductivity, pay attention
  }

  // Flux projection on the absolute orientation axes
  cellBound->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKConductivity::solveFluxConductivityInner(Coord &gradTkLeft, Coord &gradTkRight, double &alphakL, double &alphakR, int &numPhase) const
{
  //Data of the cell boundary
  double dTkdx, alphak;
  dTkdx = (gradTkLeft.getX() + gradTkRight.getX()) / 2.;
  alphak = (alphakL + alphakR) / 2.;

  //Writing of conductive terms on each equation of fluxTempXXX
  fluxBufferKapila->m_energ[numPhase] += -alphak*m_lambdak[numPhase]* dTkdx;
  fluxBufferKapila->m_energMixture += -alphak*m_lambdak[numPhase] * dTkdx;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityAbs(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  this->solveFluxConductivityInner(gradTkLeft, gradTkLeft, alphakL, alphakL, numPhase);
}

//***********************************************************************

void APKConductivity::solveFluxConductivityWall(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + 0.);
  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + 0.);
  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + 0.);
  fluxBufferKapila->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityOutflow(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + 0.);
  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + 0.);
  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + 0.);
  fluxBufferKapila->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityInflow(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + 0.);
  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + 0.);
  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + 0.);
  fluxBufferKapila->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::solveFluxConductivityOther(Coord &gradTkLeft, double &alphakL, int &numPhase) const
{
  //Not manage at the moment, just an example
  cout << "Conductive boundary not manage" << endl;

  // To avoid bug when not manage
  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + 0.);
  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + 0.);
  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + 0.);
  fluxBufferKapila->m_energMixture += 0.;
}

//***********************************************************************

void APKConductivity::communicationsAddPhys(Cell **cells, const int &dim)
{
  for (int k = 0; k < cells[0]->getNumberPhases(); k++) {
    parallel.communicationsVector(cells, "QPA", dim, m_numQPA, k);
  }
}

//***********************************************************************

void APKConductivity::communicationsAddPhysAMR(Cell **cells, const int &dim, const int &lvl)
{
	for (int k = 0; k < cells[0]->getNumberPhases(); k++) {
		parallel.communicationsVectorAMR(cells, "QPA", dim, lvl, m_numQPA, k);
	}
}

//***********************************************************************