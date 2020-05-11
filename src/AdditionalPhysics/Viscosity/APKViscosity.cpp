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
//! \author    K. Schmidmayer, J. Caze
//! \version   1.1
//! \date      November 19 2019

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APKViscosity.h"

using namespace tinyxml2;

//***********************************************************************

APKViscosity::APKViscosity(){}

//***********************************************************************

APKViscosity::APKViscosity(int& numberQPA, Eos** eos, int &numberPhases, std::string fileName)
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

void APKViscosity::solveFluxAddPhys(CellInterface *cellInterface, const int &numberPhases)
{
  // Copy velocities and gradients of left and right cells
  m_velocityLeft = cellInterface->getCellGauche()->getMixture()->getVelocity();
  m_velocityRight = cellInterface->getCellDroite()->getMixture()->getVelocity();

  m_gradULeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(1);
  m_gradURight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(2);
  m_gradVRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(3);
  m_gradWRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGrad(3);

  // Compute the mixture mu on left and right
  double muMixLeft(0.), muMixRight(0.);
  for (int k = 0; k < numberPhases; k++) {
    muMixLeft += cellInterface->getCellGauche()->getPhase(k)->getAlpha()*m_muk[k];
    muMixRight += cellInterface->getCellDroite()->getPhase(k)->getAlpha()*m_muk[k];
  }

  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_velocityRight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradULeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradURight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVRight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWRight.localProjection(m_normal, m_tangent, m_binormal);

  this->solveFluxViscosityInner(m_velocityLeft, m_velocityRight, m_gradULeft, m_gradURight,
    m_gradVLeft, m_gradVRight, m_gradWLeft, m_gradWRight, muMixLeft, muMixRight, numberPhases);

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKViscosity::solveFluxAddPhysBoundary(CellInterface *cellInterface, const int &numberPhases)
{
  ////KS//DEV// BC Injection, Tank, Outflow to do

  // Copy velocities and gradients of left and right cells
  m_velocityLeft = cellInterface->getCellGauche()->getMixture()->getVelocity();
  m_gradULeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(1);
  m_gradVLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(2);
  m_gradWLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGrad(3);

  // Compute the mixture mu on left and right
  double muMixLeft(0.);
  for (int k = 0; k < numberPhases; k++) {
    muMixLeft += cellInterface->getCellGauche()->getPhase(k)->getAlpha()*m_muk[k];
  }

  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradULeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWLeft.localProjection(m_normal, m_tangent, m_binormal);

  // Distances cells/cell interfaces for weighting on the flux
  double distLeft = cellInterface->getCellGauche()->distance(cellInterface);

  int typeCellInterface = cellInterface->whoAmI();
  if (typeCellInterface == 1 || typeCellInterface == 3 || typeCellInterface == 4 || typeCellInterface == 5 || typeCellInterface == 6) {
    // Cell interface of type NonReflecting, Outflow, Injection, Tank or Symmetry
    this->solveFluxViscosityNonReflecting(m_velocityLeft, m_gradULeft, m_gradVLeft, m_gradWLeft, muMixLeft, numberPhases);
  }
  else if (typeCellInterface == 2) {
    // Cell interface of type Wall
    this->solveFluxViscosityWall(m_velocityLeft, muMixLeft, distLeft, numberPhases);
  }
  else { this->solveFluxViscosityOther(m_velocityLeft, m_gradULeft, m_gradVLeft, m_gradWLeft, muMixLeft, numberPhases); }

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKViscosity::solveFluxViscosityInner(Coord &velocityLeft, Coord &velocityRight, Coord &gradULeft, Coord &gradURight,
  Coord &gradVLeft, Coord &gradVRight, Coord &gradWLeft, Coord &gradWRight, double &muMixLeft, double &muMixRight, int numberPhases) const
{
	//Extraction of data
	double uL, vL, wL, uR, vR, wR;
	double dudxL, dudyL, dudzL, dudxR, dudyR, dudzR;
  double dvdxL, dvdyL, dvdxR, dvdyR;
  double dwdxL, dwdzL, dwdxR, dwdzR;
	uL = velocityLeft.getX();
	vL = velocityLeft.getY();
	wL = velocityLeft.getZ();
	uR = velocityRight.getX();
	vR = velocityRight.getY();
	wR = velocityRight.getZ();

	dudxL = gradULeft.getX();
  dudyL = gradULeft.getY();
  dudzL = gradULeft.getZ();
  dudxR = gradURight.getX();
  dudyR = gradURight.getY();
  dudzR = gradURight.getZ();

  dvdxL = gradVLeft.getX();
  dvdyL = gradVLeft.getY();
  dvdxR = gradVRight.getX();
  dvdyR = gradVRight.getY();

  dwdxL = gradWLeft.getX();
  dwdzL = gradWLeft.getZ();
  dwdxR = gradWRight.getX();
  dwdzR = gradWRight.getZ();

	//Data of the cell interface
  double u, v, w;
  double dudx, dudy, dudz;
  double dvdx, dvdy;
  double dwdx, dwdz;
  double muMel;
	u = (uL + uR) / 2.;
	v = (vL + vR) / 2.;
	w = (wL + wR) / 2.;

	dudx = (dudxL + dudxR) / 2.;
  dudy = (dudyL + dudyR) / 2.;
  dudz = (dudzL + dudzR) / 2.;

  dvdx = (dvdxL + dvdxR) / 2.;
  dvdy = (dvdyL + dvdyR) / 2.;

  dwdx = (dwdxL + dwdxR) / 2.;
  dwdz = (dwdzL + dwdzR) / 2.;

  muMel = (muMixLeft + muMixRight) / 2.;

	//Writing of viscous terms on each equation of fluxBuffKapila
	for (int k = 0; k<numberPhases; k++)
	{
	  static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
	  static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
	  static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
	}
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(-muMel / 3. * (4.*dudx - 2.*(dvdy + dwdz)));
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(-muMel * (dvdx + dudy));
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(-muMel * (dwdx + dudz));
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = -muMel * (4./3.*dudx*u + (dvdx + dudy)*v + (dwdx + dudz)*w - 2./3.*(dvdy + dwdz)*u);
}

//***********************************************************************

void APKViscosity::solveFluxViscosityNonReflecting(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, int numberPhases) const
{
  this->solveFluxViscosityInner(velocityLeft, velocityLeft, gradULeft, gradULeft,
    gradVLeft, gradVLeft, gradWLeft, gradWLeft, muMixLeft, muMixLeft, numberPhases);
}

//***********************************************************************

void APKViscosity::solveFluxViscosityWall(Coord &velocityLeft, double &muMixLeft, double &distLeft, int numberPhases) const
{
  double dudx = -velocityLeft.getX() / distLeft;
  double dvdx = -velocityLeft.getY() / distLeft;
  double dwdx = -velocityLeft.getZ() / distLeft;

  for (int k = 0; k<numberPhases; k++)
  {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(-muMixLeft / 3. * 4. * dudx );
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(-muMixLeft * dvdx);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(-muMixLeft * dwdx);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = 0.;
}

//***********************************************************************

void APKViscosity::solveFluxViscosityOther(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, int numberPhases) const
{
  //Not manage at the moment, just an example
  std::cout << "Viscous boundary not manage" << std::endl;

  // To avoid bug when not manage
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxKapila*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = 0.;
}

//***********************************************************************

void APKViscosity::addNonCons(Cell *cell, const int &numberPhases)
{
  double dudx = cell->getQPA(m_numQPA)->getGrad(1).getX();
  double dudy = cell->getQPA(m_numQPA)->getGrad(1).getY();
  double dudz = cell->getQPA(m_numQPA)->getGrad(1).getZ();
  double dvdx = cell->getQPA(m_numQPA)->getGrad(2).getX();
  double dvdy = cell->getQPA(m_numQPA)->getGrad(2).getY();
  double dvdz = cell->getQPA(m_numQPA)->getGrad(2).getZ();
  double dwdx = cell->getQPA(m_numQPA)->getGrad(3).getX();
  double dwdy = cell->getQPA(m_numQPA)->getGrad(3).getY();
  double dwdz = cell->getQPA(m_numQPA)->getGrad(3).getZ();
  double termeNonCons = - 2./3.*(dudx + dvdy + dwdz)*(dudx + dvdy + dwdz)
                        + 2.*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
                        + (dudy+dvdx)*(dudy+dvdx) + (dudz+dwdx)*(dudz+dwdx) + (dvdz+dwdy)*(dvdz+dwdy);

  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = cell->getPhase(k)->getAlpha()*m_muk[k] * termeNonCons;
  }
  static_cast<FluxKapila*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = 0.;

  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKViscosity::addSymmetricTermsRadialAxisOnX(Cell *cell, const int &numberPhases)
{
  //Extraction of data
  double r = cell->getPosition().getX();
  double dudx = cell->getQPA(m_numQPA)->getGrad(1).getX();
  double dudy = cell->getQPA(m_numQPA)->getGrad(1).getY();
  double dvdx = cell->getQPA(m_numQPA)->getGrad(2).getX();
  double dvdy = cell->getQPA(m_numQPA)->getGrad(2).getY();

  //Compute the mixture mu
  double muMix(0.);
  for (int k = 0; k < numberPhases; k++) {
    muMix += cell->getPhase(k)->getAlpha()*m_muk[k];
  }

  //Writing of symmetrical viscous terms on each equation of fluxBuffKapila
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(muMix * 2. * dudx / r);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(muMix * (dvdx + dudy) / r);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = muMix * (1./3.*(4.*dudx - 2.*dvdy)*cell->getMixture()->getU() + (dudy+dvdx)*cell->getMixture()->getV()) / r;

  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKViscosity::addSymmetricTermsRadialAxisOnY(Cell *cell, const int &numberPhases)
{
  //Extraction of data
  double r = cell->getPosition().getY();
  double dudx = cell->getQPA(m_numQPA)->getGrad(1).getX();
  double dudy = cell->getQPA(m_numQPA)->getGrad(1).getY();
  double dvdx = cell->getQPA(m_numQPA)->getGrad(2).getX();
  double dvdy = cell->getQPA(m_numQPA)->getGrad(2).getY();

  //Compute the mixture mu
  double muMix(0.);
  for (int k = 0; k < numberPhases; k++) {
    muMix += cell->getPhase(k)->getAlpha()*m_muk[k];
  }

  //Writing of symmetrical viscous terms on each equation of fluxBuffKapila
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(muMix * (dvdx + dudy) / r);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(muMix * 2. * dvdy / r);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = muMix * (1./3.*(4.*dvdy - 2.*dudx)*cell->getMixture()->getV() + (dudy+dvdx)*cell->getMixture()->getU()) / r;
  
  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKViscosity::communicationsAddPhys(int numberPhases, const int &dim, const int &lvl)
{
	parallel.communicationsVector(QPA, dim, lvl, m_numQPA, 1); //m_gradU
	parallel.communicationsVector(QPA, dim, lvl, m_numQPA, 2); //m_gradV
	parallel.communicationsVector(QPA, dim, lvl, m_numQPA, 3); //m_gradW
}

//***********************************************************************
