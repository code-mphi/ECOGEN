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

//! \file      SymmetryCylindrical.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017

#include "SymCylindrical.h"

using namespace tinyxml2;

//***********************************************************************

SymCylindrical::SymCylindrical() {}

//***********************************************************************
/*!
*  Cylindrical symmetry constructor from a read in XML format
*  ex : <dataSymCyl RadialAxis="X"/>
*/
SymCylindrical::SymCylindrical(XMLElement *element, std::string nameFile)
{
  XMLElement *sousElement(element->FirstChildElement("dataSymCyl"));
  if (sousElement == NULL) throw ErrorXMLElement("dataSymCyl", nameFile, __FILE__, __LINE__);
  //Attributes collecting
  //---------------------
  //Applicated axis
  std::string axis(sousElement->Attribute("radialAxis"));
  Tools::uppercase(axis);
  if (axis == "X") { m_radialAxis = X; }
  else if (axis == "Y") { m_radialAxis = Y; }
  else if (axis == "Z") { m_radialAxis = Z; }
  else { throw ErrorXMLAttribut("radialAxis", nameFile, __FILE__, __LINE__); }
}

//***********************************************************************

SymCylindrical::~SymCylindrical() {}

//***********************************************************************

void SymCylindrical::addSymmetricTerms(Cell *cell, const int &numberPhases, Prim type)
{
  double r(0.), v(0.);
  if (numberPhases > 1) { // Multiphase models
	  switch (m_radialAxis) {
		case X: r = cell->getPosition().getX(); v = cell->getMixture(type)->getU(); break;
		case Y: r = cell->getPosition().getY(); v = cell->getMixture(type)->getV(); break;
		case Z: r = cell->getPosition().getZ(); v = cell->getMixture(type)->getW(); break;
		default: Errors::errorMessage("Name of the axis is unknown in SymCylindrical::addSymmetricTerms");
	  }
	  cell->getCons()->addSymmetricTerms(cell->getPhases(type), cell->getMixture(type), numberPhases, r, v);
  }
  else { // Euler monophasic model //JC//Q// Why MixEuler is not linked to PhaseEuler ? Like getMixture() redirects to getPhase(0)  
	  switch (m_radialAxis) {
		case X: r = cell->getPosition().getX(); v = cell->getPhase(0)->getU(); break;
		case Y: r = cell->getPosition().getY(); v = cell->getPhase(0)->getV(); break;
		case Z: r = cell->getPosition().getZ(); v = cell->getPhase(0)->getW(); break;
		default: Errors::errorMessage("Name of the axis is unknown in SymCylindrical::addSymmetricTerms");
	  }
	  cell->getCons()->addSymmetricTerms(cell->getPhases(type), cell->getMixture(type), numberPhases, r, v);
  }
}

//***********************************************************************

void SymCylindrical::addSymmetricTermsAddPhys(Cell *cell, const int &numberPhases, AddPhys &addPhys)
{
  switch (m_radialAxis) {
  case X: addPhys.addSymmetricTermsRadialAxisOnX(cell, numberPhases); break;
  case Y: addPhys.addSymmetricTermsRadialAxisOnY(cell, numberPhases); break;
  default: Errors::errorMessage("Name of the axis is unknown in SymCylindrical::addSymmetricTermsAddPhys");
  }
}

//***********************************************************************
