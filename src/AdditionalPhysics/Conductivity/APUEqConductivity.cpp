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

#include "APUEqConductivity.h"

//***********************************************************************

APUEqConductivity::APUEqConductivity() {}

//***********************************************************************

APUEqConductivity::APUEqConductivity(int& numberQPA, Eos** eos, const int& numbPhases)
{
  m_lambdak = new double[numbPhases];
  for (int k = 0; k < numbPhases; k++) {
    m_lambdak[k] = eos[k]->getLambda();
  }
  m_numQPA = numberQPA++;
}

//***********************************************************************

APUEqConductivity::~APUEqConductivity(){ delete[] m_lambdak; }

//***********************************************************************

void APUEqConductivity::addQuantityAddPhys(Cell* cell)
{
  cell->getVecQuantitiesAddPhys().push_back(new QAPConductivity(this));
}

//***********************************************************************

void APUEqConductivity::solveFluxAddPhys(CellInterface* cellInterface)
{
  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Reset of fluxBuffUEq
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum = 0.;
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;

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

void APUEqConductivity::solveFluxAddPhysBoundary(CellInterface* cellInterface)
{
  //KS//DEV// On ne fait rien aux limites avec la conductivite pour le moment, a gerer un jour

  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Reset of fluxBuffUEq (allow to then do the sum of conductivity effects for the different phases combinations)
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxUEq*> (fluxBuff)->m_momentum = 0.;
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;

  for (int numPhase = 0; numPhase < numberPhases; numPhase++) {
    // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
    m_gradTkLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(numPhase);
    m_gradTkLeft.localProjection(m_normal, m_tangent, m_binormal);

    // Extraction of alphak
    double alphakLeft = cellInterface->getCellGauche()->getPhase(numPhase)->getAlpha();

    int typeCellInterface = cellInterface->whoAmI();
    if (typeCellInterface == NONREFLECTING) { this->solveFluxConductivityNonReflecting(m_gradTkLeft, alphakLeft, numPhase); }
    else if (typeCellInterface == WALL || typeCellInterface == SYMMETRY) { this->solveFluxConductivityWall(); }
    else if (typeCellInterface == OUTFLOW) { this->solveFluxConductivityOutflow(); }
    else if (typeCellInterface == INJ) { this->solveFluxConductivityInflow(); }
    else { this->solveFluxConductivityOther(); }
    // etc... Boundaries not taken into account yet for conductivity, pay attention
  }

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APUEqConductivity::solveFluxConductivityInner(const Coord& gradTkLeft, const Coord& gradTkRight, const double& alphakL, const double& alphakR, const int& numPhase) const
{
  //Data of the cell interface
  double dTkdx, alphak;
  dTkdx = (gradTkLeft.getX() + gradTkRight.getX()) / 2.;
  alphak = (alphakL + alphakR) / 2.;

  //Writing of conductive terms on each equation of fluxBuffUEq
  static_cast<FluxUEq*> (fluxBuff)->m_energ[numPhase] = -alphak*m_lambdak[numPhase] * dTkdx;
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture   += -alphak*m_lambdak[numPhase] * dTkdx;
}

//***********************************************************************

void APUEqConductivity::solveFluxConductivityNonReflecting(const Coord& gradTkLeft, const double& alphakL, const int& numPhase) const
{
  this->solveFluxConductivityInner(gradTkLeft, gradTkLeft, alphakL, alphakL, numPhase);
}

//***********************************************************************

void APUEqConductivity::solveFluxConductivityWall() const
{
  //Not manage at the moment, just an example
  //KS//DEV// A faire pour couche limite thermique !!! ...

  // To avoid bug when not manage
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getX() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getY() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getZ() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APUEqConductivity::solveFluxConductivityOutflow() const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getX() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getY() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getZ() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APUEqConductivity::solveFluxConductivityInflow() const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getX() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getY() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getZ() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APUEqConductivity::solveFluxConductivityOther() const
{
  //Not manage at the moment, just an example
  std::cout << "Conductive boundary not manage" << std::endl;

  // To avoid bug when not manage
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getX() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getY() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(static_cast<FluxUEq*> (fluxBuff)->m_momentum.getZ() + 0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APUEqConductivity::communicationsAddPhys(const int& dim, const int& lvl)
{
  for (int k = 0; k < numberPhases; k++) {
		parallel.communicationsVector(QPA, dim, lvl, m_numQPA, k);
	}
}

//***********************************************************************
