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

//! \file      APKViscosity.cpp
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APKViscosity.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

APKViscosity::APKViscosity(){}

//***********************************************************************

APKViscosity::APKViscosity(int& numberQPA, Eos** eos, int &numberPhases, string fileName)
{
  m_muk = new double[numberPhases];
  for (int k = 0; k < numberPhases; k++) {
    m_muk[k] = eos[k]->getMu();
  }
  m_numQPA = numberQPA++;
}

//***********************************************************************

APKViscosity::~APKViscosity(){ delete[] m_muk; }


//***********************************************************************

void APKViscosity::addQuantityAddPhys(Cell *cell)
{
  cell->getVecQuantitiesAddPhys().push_back(new QAPViscosity(this));
}

//***********************************************************************

void APKViscosity::solveFluxAddPhys(CellInterface *cellBound, const int &numberPhases)
{
  // Copy velocities and gradients of left and right cells
  m_velocityLeft = cellBound->getCellGauche()->getMixture()->getVelocity();
  m_velocityRight = cellBound->getCellDroite()->getMixture()->getVelocity();

  m_gradULeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(1);
  m_gradURight = cellBound->getCellDroite()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(2);
  m_gradVRight = cellBound->getCellDroite()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(3);
  m_gradWRight = cellBound->getCellDroite()->getQPA(m_numQPA)->getGrad(3);

  // Compute the mixture mu on left and right
  double muMixLeft(0.), muMixRight(0.);
  for (int k = 0; k < numberPhases; k++) {
    muMixLeft += cellBound->getCellGauche()->getPhase(k)->getAlpha()*m_muk[k];
    muMixRight += cellBound->getCellDroite()->getPhase(k)->getAlpha()*m_muk[k];
  }

  m_normal = cellBound->getFace()->getNormal();
  m_tangent = cellBound->getFace()->getTangent();
  m_binormal = cellBound->getFace()->getBinormal();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_velocityRight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradULeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradURight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVRight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWRight.localProjection(m_normal, m_tangent, m_binormal);

  // Distances cells/cellBoundary for weighting on the flux
  double distLeft = cellBound->getCellGauche()->distance(cellBound);
  double distRight = cellBound->getCellDroite()->distance(cellBound);

  this->solveFluxViscosityInner(m_velocityLeft, m_velocityRight, m_gradULeft, m_gradURight,
    m_gradVLeft, m_gradVRight, m_gradWLeft, m_gradWRight, muMixLeft, muMixRight, distLeft, distRight, numberPhases);

  // Flux projection on the absolute orientation axes
  cellBound->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKViscosity::solveFluxAddPhysBoundary(CellInterface *cellBound, const int &numberPhases)
{
  // Copy velocities and gradients of left and right cells
  m_velocityLeft = cellBound->getCellGauche()->getMixture()->getVelocity();
  m_gradULeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft = cellBound->getCellGauche()->getQPA(m_numQPA)->getGrad(3);

  // Compute the mixture mu on left and right
  double muMixLeft(0.);
  for (int k = 0; k < numberPhases; k++) {
    muMixLeft += cellBound->getCellGauche()->getPhase(k)->getAlpha()*m_muk[k];
  }

  m_normal = cellBound->getFace()->getNormal();
  m_tangent = cellBound->getFace()->getTangent();
  m_binormal = cellBound->getFace()->getBinormal();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradULeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWLeft.localProjection(m_normal, m_tangent, m_binormal);

  // Distances cells/cellBoundary for weighting on the flux
  double distLeft = cellBound->getCellGauche()->distance(cellBound);

  int typeCellBound = cellBound->whoAmI();
  if (typeCellBound == 1 || typeCellBound == 3 || typeCellBound == 4 || typeCellBound == 5 || typeCellBound == 6) {
    // Cell boundary of type Abs, Inlet, Outlet, Res or Symmetry
    this->solveFluxViscosityAbs(m_velocityLeft, m_gradULeft, m_gradVLeft, m_gradWLeft, muMixLeft, distLeft, numberPhases);
  }
  else if (typeCellBound == 2) {
    // Cell boundary of type Wall
    this->solveFluxViscosityWall(m_velocityLeft, muMixLeft, distLeft, numberPhases);
  }
  else { this->solveFluxViscosityOther(m_velocityLeft, m_gradULeft, m_gradVLeft, m_gradWLeft, muMixLeft, distLeft, numberPhases); }

  // Flux projection on the absolute orientation axes
  cellBound->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKViscosity::solveFluxViscosityInner(Coord &velocityLeft, Coord &velocityRight, Coord &gradULeft, Coord &gradURight,
  Coord &gradVLeft, Coord &gradVRight, Coord &gradWLeft, Coord &gradWRight, double &muMixLeft, double &muMixRight, double &distLeft, double &distRight, int numberPhases) const
{
	//Extraction of data
	double uL, vL, wL, uR, vR, wR;
	double du1L, du2L, du3L, du1R, du2R, du3R;
  double dv1L, dv2L, dv1R, dv2R;
  double dw1L, dw3L, dw1R, dw3R;
	uL = velocityLeft.getX();
	vL = velocityLeft.getY();
	wL = velocityLeft.getZ();
	uR = velocityRight.getX();
	vR = velocityRight.getY();
	wR = velocityRight.getZ();

	du1L = gradULeft.getX();
  du2L = gradULeft.getY();
  du3L = gradULeft.getZ();
  du1R = gradURight.getX();
  du2R = gradURight.getY();
  du3R = gradURight.getZ();

  dv1L = gradVLeft.getX();
  dv2L = gradVLeft.getY();
  dv1R = gradVRight.getX();
  dv2R = gradVRight.getY();

  dw1L = gradWLeft.getX();
  dw3L = gradWLeft.getZ();
  dw1R = gradWRight.getX();
  dw3R = gradWRight.getZ();

	//Data of the cell boundary
  double u, v, w;
  double du1, du2, du3;
  double dv1, dv2;
  double dw1, dw3;
  double muMel;
	u = (uL*distRight + uR*distLeft) / (distLeft + distRight);
	v = (vL*distRight + vR*distLeft) / (distLeft + distRight);
	w = (wL*distRight + wR*distLeft) / (distLeft + distRight);

	du1 = (du1L*distRight + du1R*distLeft) / (distLeft + distRight);
  du2 = (du2L*distRight + du2R*distLeft) / (distLeft + distRight);
  du3 = (du3L*distRight + du3R*distLeft) / (distLeft + distRight);

  dv1 = (dv1L*distRight + dv1R*distLeft) / (distLeft + distRight);
  dv2 = (dv2L*distRight + dv2R*distLeft) / (distLeft + distRight);

  dw1 = (dw1L*distRight + dw1R*distLeft) / (distLeft + distRight);
  dw3 = (dw3L*distRight + dw3R*distLeft) / (distLeft + distRight);

  muMel = (muMixLeft*distRight + muMixRight*distLeft) / (distLeft + distRight);

	//Writing of viscous terms on each equation of fluxTempXXX
	for (int k = 0; k<numberPhases; k++)
	{
	  fluxBufferKapila->m_alpha[k] = 0.;
	  fluxBufferKapila->m_masse[k] = 0.;
	  fluxBufferKapila->m_energ[k] = 0.;
	}
  fluxBufferKapila->m_qdm.setX(-muMel / 3. * (4.*du1 - 2.*(dv2 + dw3)));
  fluxBufferKapila->m_qdm.setY(-muMel * (dv1 + du2));
  fluxBufferKapila->m_qdm.setZ(-muMel * (dw1 + du3));
  fluxBufferKapila->m_energMixture = -muMel * (4. / 3.*du1*u + (dv1 + du2)*v + (dw1 + du3)*w - 2. / 3.*(dv2 + dw3)*u);
}

//***********************************************************************

void APKViscosity::solveFluxViscosityAbs(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, double &distLeft, int numberPhases) const
{
  this->solveFluxViscosityInner(velocityLeft, velocityLeft, gradULeft, gradULeft,
    gradVLeft, gradVLeft, gradWLeft, gradWLeft, muMixLeft, muMixLeft, distLeft, distLeft, numberPhases);
}

//***********************************************************************

void APKViscosity::solveFluxViscosityWall(Coord &velocityLeft, double &muMixLeft, double &distLeft, int numberPhases) const
{
  double du1 = -velocityLeft.getX() / distLeft;
  double dv1 = -velocityLeft.getY() / distLeft;
  double dw1 = -velocityLeft.getZ() / distLeft;

  for (int k = 0; k<numberPhases; k++)
  {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = 0.;
  }
  fluxBufferKapila->m_qdm.setX(-muMixLeft / 3. * 4. * du1 );
  fluxBufferKapila->m_qdm.setY(-muMixLeft * dv1);
  fluxBufferKapila->m_qdm.setZ(-muMixLeft * dw1);
  fluxBufferKapila->m_energMixture = 0.;
}

//***********************************************************************

void APKViscosity::solveFluxViscosityOther(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, double &distLeft, int numberPhases) const
{
  //Not manage at the moment, just an example
  cout << "Viscous boundary not manage" << endl;

  // To avoid bug when not manage
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = 0.;
  }
  fluxBufferKapila->m_qdm = 0.;
  fluxBufferKapila->m_energMixture = 0.;
}

//***********************************************************************

void APKViscosity::addNonCons(Cell *cell, const int &numberPhases)
{
  double du1 = cell->getQPA(m_numQPA)->getGrad(1).getX();
  double dv2 = cell->getQPA(m_numQPA)->getGrad(2).getY();
  double dw3 = cell->getQPA(m_numQPA)->getGrad(3).getZ();
  double termeNonCons = 2.*(du1*du1 + dv2*dv2 + dw3*dw3) - 2. / 3.*(du1 + dv2 + dw3)*(du1 + dv2 + dw3);

  for (int k = 0; k<numberPhases; k++) {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = cell->getPhase(k)->getAlpha()*m_muk[k] * termeNonCons;
  }
  fluxBufferKapila->m_qdm = 0.;
  fluxBufferKapila->m_energMixture = 0.;

  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKViscosity::communicationsAddPhys(Cell **cells, const int &dim)
{
  parallel.communicationsVector(cells, "QPA", dim, m_numQPA, 1); //m_gradU
  parallel.communicationsVector(cells, "QPA", dim, m_numQPA, 2); //m_gradV
  parallel.communicationsVector(cells, "QPA", dim, m_numQPA, 3); //m_gradW
}

//***********************************************************************

void APKViscosity::communicationsAddPhysAMR(Cell **cells, const int &dim, const int &lvl)
{
	parallel.communicationsVectorAMR(cells, "QPA", dim, lvl, m_numQPA, 1); //m_gradU
	parallel.communicationsVectorAMR(cells, "QPA", dim, lvl, m_numQPA, 2); //m_gradV
	parallel.communicationsVectorAMR(cells, "QPA", dim, lvl, m_numQPA, 3); //m_gradW
}

//***********************************************************************
