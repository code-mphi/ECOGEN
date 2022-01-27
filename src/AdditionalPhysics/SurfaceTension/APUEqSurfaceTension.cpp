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

#include "APUEqSurfaceTension.h"

using namespace tinyxml2;

//***********************************************************************
/*!
*  APUEqSurfaceTension constructor from a read in XML format
*  ex : <dataSurfaceTension transport="color1" sigma="800."/>
*/
APUEqSurfaceTension::APUEqSurfaceTension(XMLElement* element, int& numberQPA, std::vector<std::string> nameTransports, std::vector<std::string> namePhases, std::string nameFichier)
{
  //Collecting attributes
  //---------------------
  XMLElement* sousElement(element->FirstChildElement("dataSurfaceTension"));
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
  XMLElement* sousElement2(element->FirstChildElement("reinitializationTransport"));
  if (sousElement2 != NULL) {
    m_reinitializationActivated = true;
    //Name of the associated fluid
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

APUEqSurfaceTension::~APUEqSurfaceTension(){}

//***********************************************************************

void APUEqSurfaceTension::addQuantityAddPhys(Cell* cell)
{
  cell->getVecQuantitiesAddPhys().push_back(new QAPSurfaceTension(this));
}

//***********************************************************************

double APUEqSurfaceTension::computeEnergyAddPhys(QuantitiesAddPhys* QPA)
{
  double energyCap = m_sigma*QPA->getGrad().norm(); //prepare volume capillary energy (sigma*gradC)
  return energyCap;
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxAddPhys(CellInterface* cellInterface)
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

  // Reset of fluxBuffUEq
  static_cast<FluxUEq*> (fluxBuff)->setToZero();

  this->solveFluxSurfaceTensionInner(m_velocityLeft, m_velocityRight, m_gradCLeft, m_gradCRight);

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxAddPhysBoundary(CellInterface* cellInterface)
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

  // Reset of fluxBuffUEq (allow to then do the sum of surface-tension effects for the different phases combinations)
  static_cast<FluxUEq*> (fluxBuff)->setToZero();

  int typeCellInterface = cellInterface->whoAmI();
  if (typeCellInterface == NONREFLECTING) { this->solveFluxSurfaceTensionNonReflecting(m_velocityLeft, m_gradCLeft); }
  else if (typeCellInterface == WALL || typeCellInterface == SYMMETRY) { this->solveFluxSurfaceTensionWall(m_gradCLeft); }
  else if (typeCellInterface == OUTFLOW) { this->solveFluxSurfaceTensionOutflow(); }
  else if (typeCellInterface == INJ) { this->solveFluxSurfaceTensionInflow(); }
  else { this->solveFluxSurfaceTensionOther(); } //Else
  // etc... Boundaries not taken into account yet for surface tension, pay attention

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxSurfaceTensionInner(const Coord& velocityLeft, const Coord& velocityRight, const Coord& gradCLeft, const Coord& gradCRight) const
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

	//Writing of surface-tension terms on each equation of fluxBuffUEq
	if (normW > 1.e-6) {
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(- m_sigma*(w2*w2 + w3*w3) / normW);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(m_sigma *w1*w2 / normW);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(m_sigma *w1*w3 / normW);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = m_sigma * (w1*w1*u + w1*w2*v + w1*w3*w) / normW;
	}
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxSurfaceTensionNonReflecting(const Coord& velocityLeft, const Coord& gradCLeft) const
{
  this->solveFluxSurfaceTensionInner(velocityLeft, velocityLeft, gradCLeft, gradCLeft);
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxSurfaceTensionWall(const Coord& gradCLeft) const
{
  //Considered as a symmetry and not as a wall with a triple vertex yet

  //Data of the cell interface
  double normW, w2, w3;
  normW = gradCLeft.norm();
  w2 = gradCLeft.getY();
  w3 = gradCLeft.getZ();

  //Writing of surface-tension terms on each equation of fluxBuffUEq
  if (normW > 1.e-6) {
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(- m_sigma*(w2*w2 + w3*w3) / normW);
  }
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxSurfaceTensionOutflow() const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxSurfaceTensionInflow() const
{
  //Not manage at the moment, just an example

  // To avoid bug when not manage
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;
}

//***********************************************************************

void APUEqSurfaceTension::solveFluxSurfaceTensionOther() const
{
  //Not manage at the moment, just an example
  std::cout << "Surface-tension boundary not manage" << std::endl;

  // To avoid bug when not manage
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_momentum.setZ(0.);
  static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;
}

//***********************************************************************

void APUEqSurfaceTension::addSymmetricTermsRadialAxisOnX(Cell* cell)
{
  //Extraction of data
  double r, w1, w2, normW;
  r = cell->getPosition().getX();
  w1 = cell->getQPA(m_numQPAGradC)->getGrad().getX();
  w2 = cell->getQPA(m_numQPAGradC)->getGrad().getY();
  normW = cell->getQPA(m_numQPAGradC)->getGrad().norm();

  //Writing of symmetrical surface-tension terms on each equation of fluxBuffUEq
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  if (normW > 1.e-6) {
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(-m_sigma * w1*w1 / normW / r);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(-m_sigma * w1*w2 / normW / r);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = -m_sigma * (w1*w1*cell->getMixture()->getU() + w1*w2*cell->getMixture()->getV()) / normW / r;
  }
  else {
    static_cast<FluxUEq*> (fluxBuff)->m_momentum = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;
  }

  cell->getCons()->addFlux(1.);
}

//***********************************************************************

void APUEqSurfaceTension::addSymmetricTermsRadialAxisOnY(Cell* cell)
{
  //Extraction of data
  double r, w1, w2, normW;
  r = cell->getPosition().getY();
  w1 = cell->getQPA(m_numQPAGradC)->getGrad().getX();
  w2 = cell->getQPA(m_numQPAGradC)->getGrad().getY();
  normW = cell->getQPA(m_numQPAGradC)->getGrad().norm();

  //Writing of symmetrical surface-tension terms on each equation of fluxBuffUEq
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxUEq*> (fluxBuff)->m_alpha[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_mass[k] = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energ[k] = 0.;
  }
  if (normW > 1.e-6) {
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setX(-m_sigma * w1*w2 / normW / r);
    static_cast<FluxUEq*> (fluxBuff)->m_momentum.setY(-m_sigma * w2*w2 / normW / r);
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = -m_sigma * (w1*w2*cell->getMixture()->getU() + w2*w2*cell->getMixture()->getV()) / normW / r;
  }
  else {
    static_cast<FluxUEq*> (fluxBuff)->m_momentum = 0.;
    static_cast<FluxUEq*> (fluxBuff)->m_energMixture = 0.;
  }
  
  cell->getCons()->addFlux(1.);
}

//***********************************************************************

void APUEqSurfaceTension::reinitializeColorFunction(std::vector<Cell*>* cellsLvl, const int& lvl)
{
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
    if (!cellsLvl[lvl][i]->getSplit()) { cellsLvl[lvl][i]->reinitializeColorFunction(m_numTransportAssociated, m_numPhaseAssociated); }
  }
}

//***********************************************************************

void APUEqSurfaceTension::communicationsAddPhys(const int& dim, const int& lvl)
{
	parallel.communicationsVector(QPA, dim, lvl, m_numQPAGradC);
}

//***********************************************************************
