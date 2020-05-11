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
//! \author    K. Schmidmayer, J. Caze
//! \version   1.1
//! \date      November 19 2019

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APKSurfaceTension.h"

using namespace tinyxml2;

//***********************************************************************

APKSurfaceTension::APKSurfaceTension(){}

//***********************************************************************
/*!
*  APKSurfaceTension constructor from a read in XML format
*  ex : <dataSurfaceTension transport="color1" sigma="800."/>
*/
APKSurfaceTension::APKSurfaceTension(XMLElement *element, int& numberQPA, std::vector<std::string> nameTransports, std::vector<std::string> namePhases, std::string nameFichier)
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
    //Name of associated fluid
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

void APKSurfaceTension::solveFluxAddPhys(CellInterface *cellInterface, const int& numberPhases)
{
  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Copy and projection on orientation axes attached to the edge of velocities of left and right cells
  m_velocityLeft = cellInterface->getCellGauche()->getMixture()->getVelocity();
  m_velocityRight = cellInterface->getCellDroite()->getMixture()->getVelocity();
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_velocityRight.localProjection(m_normal, m_tangent, m_binormal);

  // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
  m_gradCLeft = cellInterface->getCellGauche()->getQPA(m_numQPAGradC)->getGrad();
  m_gradCRight = cellInterface->getCellDroite()->getQPA(m_numQPAGradC)->getGrad();
  m_gradCLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradCRight.localProjection(m_normal, m_tangent, m_binormal);

  // Reset of fluxBuffKapila
  static_cast<FluxKapila*> (fluxBuff)->setToZero(numberPhases);

  this->solveFluxSurfaceTensionInner(m_velocityLeft, m_velocityRight, m_gradCLeft, m_gradCRight);

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APKSurfaceTension::solveFluxAddPhysBoundary(CellInterface *cellInterface, const int &numberPhases)
{
  //KS//DEV// Nothing special is done at the boundaries with surface tension right now (considered as symmetry)

  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Copy and projection on orientation axes attached to the edge of velocities of left and right cells
  m_velocityLeft = cellInterface->getCellGauche()->getMixture()->getVelocity();
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  
  // Copy and projection on orientation axes attached to the edge of gradients of left and right cells
  m_gradCLeft = cellInterface->getCellGauche()->getQPA(m_numQPAGradC)->getGrad();
  m_gradCLeft.localProjection(m_normal, m_tangent, m_binormal);

  // Reset of fluxBuffKapila (allow to then do the sum of surface-tension effects for the different phases combinations)
  static_cast<FluxKapila*> (fluxBuff)->setToZero(numberPhases);

  int typeCellInterface = cellInterface->whoAmI();
  if (typeCellInterface == 1) { this->solveFluxSurfaceTensionNonReflecting(m_velocityLeft, m_gradCLeft); } //Non-reflecting
  else if (typeCellInterface == 2 || typeCellInterface == 6) { this->solveFluxSurfaceTensionWall(m_gradCLeft); } //Wall or Symmetry
  else if (typeCellInterface == 3) { this->solveFluxSurfaceTensionOutflow(m_velocityLeft, m_gradCLeft); } //Outflow
  else if (typeCellInterface == 4) { this->solveFluxSurfaceTensionInflow(m_velocityLeft, m_gradCLeft); } //Injection
  else { this->solveFluxSurfaceTensionOther(m_velocityLeft, m_gradCLeft); } //Tank or else
  // etc... Boundaries not taken into account yet for surface tension, pay attention

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
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

	//Data of the cell interface
	double u, v, w, w1, w2, w3, normW;
	u = (uL + uR) / 2.;
	v = (vL + vR) / 2.;
	w = (wL + wR) / 2.;
	w1 = (w1L + w1R) / 2.;
	w2 = (w2L + w2R) / 2.;
	w3 = (w3L + w3R) / 2.;
	normW = (normWL + normWR) / 2.;

	//Writing of surface-tension terms on each equation of fluxBuffKapila
	if (normW > 1.e-6) {
	  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + m_sigma*(w1*w1 / normW - normW));
	  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + m_sigma *w1*w2 / normW);
	  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + m_sigma *w1*w3 / normW);
    static_cast<FluxKapila*> (fluxBuff)->m_energMixture += m_sigma / normW*(w1*w1*u + w1*w2*v + w1*w3*w);
	}
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionNonReflecting(Coord &velocityLeft, Coord &gradCLeft) const
{
  this->solveFluxSurfaceTensionInner(velocityLeft, velocityLeft, gradCLeft, gradCLeft);
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionWall(Coord &gradCLeft) const
{
  //Considered as a symmetry and not as a wall with a triple vertex yet

  //Data of the cell interface
  double normW;
  normW = gradCLeft.norm();

  //Writing of surface-tension terms on each equation of fluxBuffKapila
  if (normW > 1.e-6) {
    static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() - m_sigma*normW);
  }
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionOutflow(Coord &velocityLeft, Coord &gradCLeft) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionInflow(Coord &velocityLeft, Coord &gradCLeft) const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APKSurfaceTension::solveFluxSurfaceTensionOther(Coord &velocityLeft, Coord &gradCLeft) const
{
  //Not manage at the moment, just an example
  std::cout << "Surface-tension boundary not manage" << std::endl;

  // To avoid bug when not manage
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getX() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getY() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setZ(static_cast<FluxKapila*> (fluxBuff)->m_qdm.getZ() + 0.);
  static_cast<FluxKapila*> (fluxBuff)->m_energMixture += 0.;
}

//***********************************************************************

void APKSurfaceTension::addSymmetricTermsRadialAxisOnX(Cell *cell, const int &numberPhases)
{
  //Extraction of data
  double r, w1, w2, normW;
  r = cell->getPosition().getX();
  w1 = cell->getQPA(m_numQPAGradC)->getGrad().getX();
  w2 = cell->getQPA(m_numQPAGradC)->getGrad().getY();
  normW = cell->getQPA(m_numQPAGradC)->getGrad().norm();

  //Writing of symmetrical surface-tension terms on each equation of fluxBuffKapila
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  if (normW > 1.e-6) {
    static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(-m_sigma * w1*w1 / normW / r);
    static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(-m_sigma * w1*w2 / normW / r);
    static_cast<FluxKapila*> (fluxBuff)->m_energMixture = -m_sigma / normW * (w1*w1*cell->getMixture()->getU() + w1*w2*cell->getMixture()->getV()) / r;
  }
  else {
    static_cast<FluxKapila*> (fluxBuff)->m_qdm = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energMixture = 0.;
  }

  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKSurfaceTension::addSymmetricTermsRadialAxisOnY(Cell *cell, const int &numberPhases)
{
  //Extraction of data
  double r, w1, w2, normW;
  r = cell->getPosition().getY();
  w1 = cell->getQPA(m_numQPAGradC)->getGrad().getX();
  w2 = cell->getQPA(m_numQPAGradC)->getGrad().getY();
  normW = cell->getQPA(m_numQPAGradC)->getGrad().norm();

  //Writing of symmetrical surface-tension terms on each equation of fluxBuffKapila
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxKapila*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_masse[k] = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energ[k] = 0.;
  }
  if (normW > 1.e-6) {
    static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(-m_sigma * w1*w2 / normW / r);
    static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(-m_sigma * w2*w2 / normW / r);
    static_cast<FluxKapila*> (fluxBuff)->m_energMixture = -m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r;

    //if ((std::fabs(-m_sigma * w1*w2 / normW / r) + std::fabs(-m_sigma * w2*w2 / normW / r)) > 2*(std::fabs(cell->getCons()->getQdm().getX()) + std::fabs(cell->getCons()->getQdm().getY()))) {
    //  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX((-m_sigma * w1*w2 / normW / r) * 2*(std::fabs(cell->getCons()->getQdm().getX()) + std::fabs(cell->getCons()->getQdm().getY())) / (std::fabs(-m_sigma * w1*w2 / normW / r) + std::fabs(-m_sigma * w2*w2 / normW / r)));
    //  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY((-m_sigma * w2*w2 / normW / r) * 2*(std::fabs(cell->getCons()->getQdm().getX()) + std::fabs(cell->getCons()->getQdm().getY())) / (std::fabs(-m_sigma * w1*w2 / normW / r) + std::fabs(-m_sigma * w2*w2 / normW / r)));
    //}
    //else {
    //  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setX(-m_sigma * w1*w2 / normW / r);
    //  static_cast<FluxKapila*> (fluxBuff)->m_qdm.setY(-m_sigma * w2*w2 / normW / r);
    //}
    //if (std::fabs(-m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r) > 2*std::fabs(cell->getCons()->getEnergyMix())) {
    //  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = (-m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r) * 2*std::fabs(cell->getCons()->getEnergyMix()) / std::fabs(-m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r);
    //}
    //else {
    //  static_cast<FluxKapila*> (fluxBuff)->m_energMixture = -m_sigma / normW * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / r;
    //}
  }
  else {
    static_cast<FluxKapila*> (fluxBuff)->m_qdm = 0.;
    static_cast<FluxKapila*> (fluxBuff)->m_energMixture = 0.;
  }
  
  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APKSurfaceTension::reinitializeColorFunction(std::vector<Cell *> *cellsLvl, int &lvl)
{
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
    if (!cellsLvl[lvl][i]->getSplit()) { cellsLvl[lvl][i]->reinitializeColorFunction(m_numTransportAssociated, m_numPhaseAssociated); }
  }
}

//***********************************************************************

void APKSurfaceTension::communicationsAddPhys(int numberPhases, const int &dim, const int &lvl)
{
	parallel.communicationsVector(QPA, dim, lvl, m_numQPAGradC);
}

//***********************************************************************
