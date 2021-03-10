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

#include "BoundCondSubInj.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondSubInj::BoundCondSubInj(int numPhysique, XMLElement* element, int& numberPhases, std::string fileName) : BoundCond(numPhysique)
{
  XMLElement* subElement(element->FirstChildElement("dataInjection"));
  if (subElement == NULL) throw ErrorXMLElement("dataInjection", fileName, __FILE__, __LINE__);

  // Reading specific massflow (kg.s-1.m-2)
  // --------------------------------------
  XMLError error;
  error = subElement->QueryDoubleAttribute("m0", &m_m0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("m0", fileName, __FILE__, __LINE__);
  m_m0 = - m_m0; // Sign change due to right side Riemann solver convention

  // Reading temperature of fluid
  // ----------------------------
  error = subElement->QueryDoubleAttribute("T0", &m_T0);
  if (error != tinyxml2::XML_NO_ERROR) throw ErrorXMLAttribut("T0", fileName, __FILE__, __LINE__);

  // Boundary condition specific to Euler model
  if (numberPhases > 1) { 
    Errors::errorMessage("Subsonic inflow is only available for Euler model");;
    throw ErrorXMLLimite(fileName, __FILE__, __LINE__); 
  }
}

//****************************************************************************

BoundCondSubInj::BoundCondSubInj(const BoundCondSubInj& Source, const int& lvl) : BoundCond(Source, lvl)
{
  m_m0 = Source.m_m0;
  m_T0 = Source.m_T0;
}

//****************************************************************************

BoundCondSubInj::~BoundCondSubInj()
{
}

//****************************************************************************

void BoundCondSubInj::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondSubInj(*(this)));
}

//****************************************************************************

void BoundCondSubInj::solveRiemannBoundary(Cell& cellLeft, const int& numberPhases, const double& dxLeft, double& dtMax)
{
  m_mod->solveRiemannSubInj(cellLeft, numberPhases, dxLeft, dtMax, m_m0, m_T0, m_massflow, m_powerFlux);
}

//****************************************************************************

void BoundCondSubInj::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_m0 << std::endl;
  std::cout << m_T0 << std::endl;
}

//****************************************************************************
//******************************AMR Method***********************************
//****************************************************************************

void BoundCondSubInj::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondSubInj(*this, m_lvl + 1));
}

//****************************************************************************