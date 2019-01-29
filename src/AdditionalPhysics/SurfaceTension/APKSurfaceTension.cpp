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

//! \file      APKSurfaceTension.cpp
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APKSurfaceTension.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

APKSurfaceTension::APKSurfaceTension(){}

//***********************************************************************
/*!
*  APKSurfaceTension constructor from a read in XML format
*  ex : <dataSurfaceTension transport="color1" sigma="800."/>
*/
APKSurfaceTension::APKSurfaceTension(XMLElement *element, int& numberQPA, vector<string> nameTransports, vector<string> namePhases, string nameFichier)
{
  //Collecting attributes
  //---------------------
  XMLElement *sousElement(element->FirstChildElement("dataSurfaceTension"));
  if (sousElement == NULL) throw ErrorXMLElement("dataSurfaceTension", nameFichier, __FILE__, __LINE__);
  XMLError error;
  //Name of the associated transport equation
  m_nameTransportAssociated = sousElement->Attribute("transport");
  if (m_nameTransportAssociated == "") throw ErrorXMLAttribut("transport", nameFichier, __FILE__, __LINE__);
  //We directly associate with the corresponding number in m_vecTransport of cells
  unsigned int t(0);
  for (t = 0; t < nameTransports.size(); t++) {
    if (m_nameTransportAssociated == nameTransports[t]) { break; }
  }
  if (t != nameTransports.size()) { m_numTransportAssociated= t; }
  else { Errors::errorMessage("Transport equation not found for surface tension"); }

  //Value of surface tension coefficient
  error = sousElement->QueryDoubleAttribute("sigma", &m_sigma);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("sigma", nameFichier, __FILE__, __LINE__);

  //Reinitialization of transport with a given fluid
  //------------------------------------------------
  m_reinitializationActivated = false;
  XMLElement *sousElement2(element->FirstChildElement("reinitializationTransport"));
  if (sousElement2 != NULL) {
    m_reinitializationActivated = true;
    //Name du fluide associee
    m_namePhaseAssociated = sousElement2->Attribute("phase");
    if (m_namePhaseAssociated == "") throw ErrorXMLAttribut("phase", nameFichier, __FILE__, __LINE__);
    //On associe tout de suite avec le number correspondant dans m_vecPhases des cells
    unsigned int p(0);
    for (p = 0; p < namePhases.size(); p++) {
      if (m_namePhaseAssociated == namePhases[p]) { break; }
    }
    if (p != namePhases.size()) { m_numPhaseAssociated = p; }
    else { Errors::errorMessage("Equation of phase not found for surface tension (reinitialization)"); }
  }

  //Matching with the Quantities of additional physics
  m_numQPAGradC = numberQPA++;
}

//***********************************************************************

APKSurfaceTension::~APKSurfaceTension(){}

//***********************************************************************

void APKSurfaceTension::addQuantityAddPhys(Cell *cell)
{
  cell->getVecQuantitiesAddPhys().push_back(new QAPSurfaceTension(this));
}

//***********************************************************************

double APKSurfaceTension::computeEnergyAddPhys(QuantitiesAddPhys* QPA)
{
  double energyCap = m_sigma*QPA->getGrad().norm(); //prepare volume capillary energy (sigma*gradC)
  return energyCap;
}

//***********************************************************************

void APKSurfaceTension::solveFluxAddPhys(CellInterface *cellBound, const int& numberPhases)
{
  m_normal = cellBound->getFace()->getNormal();
  m_tangent = cellBound->getFace()->getTangent();
  m_binormal = cellBound->getFace()->getBinormal();

  // Copy and projection on orientation axes attached to the edge of velocities of left and right cells
  m_velocityLeft = cellBound->getCellGauche()->getMixture()->getVelocity();
  m_velocityRight = cellBound->getCellDroite()->getMixture()->getVelocity();
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_velocityRight.localProjection(m_normal, m_tangent, m_binormal);

  // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
  m_gradCLeft = cellBound->getCellGauche()->getQPA(m_numQPAGradC)->getGrad();
  m_gradCRight = cellBound->getCellDroite()->getQPA(m_numQPAGradC)->getGrad();
  m_gradCLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradCRight.localProjection(m_normal, m_tangent, m_binormal);

  // Reset of fluxBufferKapila
  fluxBufferKapila->setToZero(numberPhases);

  this->solveFluxSurfaceTensionInner(m_velocityLeft, m_velocityRight, m_gradCLeft, m_gradCRight);

  // Flux projection on the absolute orientation axes
  cellBound->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKSurfaceTension::solveFluxAddPhysBoundary(CellInterface *cellBound, const int &numberPhases)
{
  m_normal = cellBound->getFace()->getNormal();
  m_tangent = cellBound->getFace()->getTangent();
  m_binormal = cellBound->getFace()->getBinormal();

  // Copy and projection on orientation axes attached to the edge of velocities of left and right cells
  m_velocityLeft = cellBound->getCellGauche()->getMixture()->getVelocity();
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  
  // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
  m_gradCLeft = cellBound->getCellGauche()->getQPA(m_numQPAGradC)->getGrad();
  m_gradCLeft.localProjection(m_normal, m_tangent, m_binormal);

  // Reset of fluxBufferKapila (allow to then do the sum of surface-tension effects for the different phases combinations)
  fluxBufferKapila->setToZero(numberPhases);

  int typeBord = cellBound->whoAmI();
  if (typeBord == 1) { this->solveFluxSurfaceTensionAbs(m_velocityLeft, m_gradCLeft); }
  else if (typeBord == 2 || typeBord == 6) { this->solveFluxSurfaceTensionWall(m_gradCLeft); }
  else if (typeBord == 3) { this->solveFluxSurfaceTensionOutflow(m_velocityLeft, m_gradCLeft); }
  else if (typeBord == 4) { this->solveFluxSurfaceTensionInflow(m_velocityLeft, m_gradCLeft); }
  else { this->solveFluxSurfaceTensionOther(m_velocityLeft, m_gradCLeft); }
  // etc... Boundaries not taken into account yet for surface tension, pay attention

  // Flux projection on the absolute orientation axes
  cellBound->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionInner(Coord &velocityLeft, Coord &velocityRight, Coord &gradCLeft, Coord &gradCRight) const
{
	//Extraction of data
	double uL, vL, wL, uR, vR, wR;
	double w1L, w2L, w3L, normWL, w1R, w2R, w3R, normWR;
	uL = velocityLeft.getX();
	vL = velocityLeft.getY();
	wL = velocityLeft.getZ();
	uR = velocityRight.getX();
	vR = velocityRight.getY();
	wR = velocityRight.getZ();
	w1L = gradCLeft.getX();
	w2L = gradCLeft.getY();
	w3L = gradCLeft.getZ();
	normWL = gradCLeft.norm();
	w1R = gradCRight.getX();
	w2R = gradCRight.getY();
	w3R = gradCRight.getZ();
	normWR = gradCRight.norm();

	//Data of the cell boundary
	double u, v, w, w1, w2, w3, normW;
	u = (uL + uR) / 2.;
	v = (vL + vR) / 2.;
	w = (wL + wR) / 2.;
	w1 = (w1L + w1R) / 2.;
	w2 = (w2L + w2R) / 2.;
	w3 = (w3L + w3R) / 2.;
	normW = (normWL + normWR) / 2.;

	//Writing of surface-tension terms on each equation of fluxBufferKapila
	if (normW > 1.e-6) {
	  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + m_sigma*(w1*w1 / normW - normW));
	  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + m_sigma *w1*w2 / normW);
	  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + m_sigma *w1*w3 / normW);
    fluxBufferKapila->m_energMixture += m_sigma / normW*(w1*w1*u + w1*w2*v + w1*w3*w);
	}
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionAbs(Coord &velocityLeft, Coord &gradCLeft) const
{
  this->solveFluxSurfaceTensionInner(velocityLeft, velocityLeft, gradCLeft, gradCLeft);
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionWall(Coord &gradCLeft) const
{
  //Considered as a symmetry and not as a wall with a triple vertex yet

  //Data of the cell boundary
  double normW;
  normW = gradCLeft.norm();

  //Writing of surface-tension terms on each equation of fluxBufferKapila
  if (normW > 1.e-6) {
    fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() - m_sigma*normW);
  }
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionOutflow(Coord &velocityLeft, Coord &gradCLeft) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + 0.);
  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + 0.);
  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + 0.);
  fluxBufferKapila->m_energMixture += 0.;
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionInflow(Coord &velocityLeft, Coord &gradCLeft) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + 0.);
  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + 0.);
  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + 0.);
  fluxBufferKapila->m_energMixture += 0.;
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionOther(Coord &velocityLeft, Coord &gradCLeft) const
{
  //Not manage at the moment, just an example
  cout << "Surface-tension boundary not manage" << endl;

  // To avoid bug when not manage
  fluxBufferKapila->m_qdm.setX(fluxBufferKapila->m_qdm.getX() + 0.);
  fluxBufferKapila->m_qdm.setY(fluxBufferKapila->m_qdm.getY() + 0.);
  fluxBufferKapila->m_qdm.setZ(fluxBufferKapila->m_qdm.getZ() + 0.);
  fluxBufferKapila->m_energMixture += 0.;
}

//***********************************************************************

void APKSurfaceTension::addSymmetricTermsRadialAxeOnX(Cell *cell, const int &numberPhases)
{
  //Extraction of data
  double r, w1, w2, normW;
  r = cell->getPosition().getY();
  w1 = cell->getQPA(m_numQPAGradC)->getGrad().getX();
  w2 = cell->getQPA(m_numQPAGradC)->getGrad().getY();
  normW = cell->getQPA(m_numQPAGradC)->getGrad().norm();

  //Writing of symmetrical surface-tension terms on each equation of fluxBufferKapila
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = 0.;
  }
  if (normW > 1.e-6) {
    fluxBufferKapila->m_qdm.setX(-m_sigma * w1*w1 / normW / r);
    fluxBufferKapila->m_qdm.setY(-m_sigma * w1*w2 / normW / r);
    fluxBufferKapila->m_energMixture = -m_sigma / normW * (w1*w1*cell->getMixture()->getU() + w1*w2*cell->getMixture()->getV()) / r;
  }
  else {
    fluxBufferKapila->m_qdm = 0.;
    fluxBufferKapila->m_energMixture = 0.;
  }

  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKSurfaceTension::addSymmetricTermsRadialAxeOnY(Cell *cell, const int &numberPhases)
{
  //Extraction of data
  double r, w1, w2, normW;
  r = cell->getPosition().getY();
  w1 = cell->getQPA(m_numQPAGradC)->getGrad().getX();
  w2 = cell->getQPA(m_numQPAGradC)->getGrad().getY();
  normW = cell->getQPA(m_numQPAGradC)->getGrad().norm();

  //Writing of symmetrical surface-tension terms on each equation of fluxBufferKapila
  for (int k = 0; k<numberPhases; k++) {
    fluxBufferKapila->m_alpha[k] = 0.;
    fluxBufferKapila->m_masse[k] = 0.;
    fluxBufferKapila->m_energ[k] = 0.;
  }
  if (normW > 1.e-6) {
    fluxBufferKapila->m_qdm.setX(-m_sigma * w1*w2 / normW / r);
    fluxBufferKapila->m_qdm.setY(-m_sigma * w2*w2 / normW / r);
    fluxBufferKapila->m_energMixture = -m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r;

    //if ((abs(-m_sigma * w1*w2 / normW / r) + abs(-m_sigma * w2*w2 / normW / r)) > 2*(abs(cell->getCons()->getQdm().getX()) + abs(cell->getCons()->getQdm().getY()))) {
    //  fluxBufferKapila->m_qdm.setX((-m_sigma * w1*w2 / normW / r) * 2*(abs(cell->getCons()->getQdm().getX()) + abs(cell->getCons()->getQdm().getY())) / (abs(-m_sigma * w1*w2 / normW / r) + abs(-m_sigma * w2*w2 / normW / r)));
    //  fluxBufferKapila->m_qdm.setY((-m_sigma * w2*w2 / normW / r) * 2*(abs(cell->getCons()->getQdm().getX()) + abs(cell->getCons()->getQdm().getY())) / (abs(-m_sigma * w1*w2 / normW / r) + abs(-m_sigma * w2*w2 / normW / r)));
    //}
    //else {
    //  fluxBufferKapila->m_qdm.setX(-m_sigma * w1*w2 / normW / r);
    //  fluxBufferKapila->m_qdm.setY(-m_sigma * w2*w2 / normW / r);
    //}
    //if (abs(-m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r) > 2*abs(cell->getCons()->getEnergyMix())) {
    //  fluxBufferKapila->m_energMixture = (-m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r) * 2*abs(cell->getCons()->getEnergyMix()) / abs(-m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r);
    //}
    //else {
    //  fluxBufferKapila->m_energMixture = -m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r;
    //}
  }
  else {
    fluxBufferKapila->m_qdm = 0.;
    fluxBufferKapila->m_energMixture = 0.;
  }
  
  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKSurfaceTension::reinitializeColorFunction(vector<Cell *> *cellsLvl, int &lvl)
{
  if (m_reinitializationActivated) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) { cellsLvl[lvl][i]->reinitializeColorFunction(m_numTransportAssociated, m_numPhaseAssociated); }
    }
  }
}

//***********************************************************************

void APKSurfaceTension::communicationsAddPhys(Cell **cells, const int &dim)
{
  parallel.communicationsVector(cells, "QPA", dim, m_numQPAGradC);
}

//***********************************************************************

void APKSurfaceTension::communicationsAddPhysAMR(Cell **cells, const int &dim, const int &lvl)
{
	parallel.communicationsVectorAMR(cells, "QPA", dim, lvl, m_numQPAGradC);
}

//***********************************************************************

int APKSurfaceTension::getNumTransportAssociated() const
{
  return m_numTransportAssociated;
}

//***********************************************************************