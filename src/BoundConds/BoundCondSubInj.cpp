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

BoundCondSubInj::BoundCondSubInj(int numPhysique, XMLElement* element, const int& numbPhases, std::string fileName) : BoundCond(numPhysique)
{
  m_Tk0 = new double[numbPhases];
  m_ak0 = new double[numbPhases];
  
  XMLElement* subElement(element->FirstChildElement("dataInjection"));
  if (subElement == NULL) throw ErrorXMLElement("dataInjection", fileName, __FILE__, __LINE__);

  // Reading specific massflow (kg.s-1.m-2)
  // --------------------------------------
  XMLError error;
  error = subElement->QueryDoubleAttribute("m0", &m_m0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("m0", fileName, __FILE__, __LINE__);
  m_m0 = - m_m0; // Sign change due to right side Riemann solver convention

  // Reading volume fraction and temperature of the phases
  // -----------------------------------------------------
  XMLElement* fluid(element->FirstChildElement("dataFluid"));
  for (int k = 0; k < numbPhases; k++) {
    // Attributes reading
    error = fluid->QueryDoubleAttribute("temperature", &m_Tk0[k]);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("temperature", fileName, __FILE__, __LINE__);
    if (numbPhases > 1) {
      error = fluid->QueryDoubleAttribute("alpha", &m_ak0[k]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("alpha", fileName, __FILE__, __LINE__);
    }
    else { m_ak0[0] = 1.; }
    fluid = fluid->NextSiblingElement("dataFluid");
  }
  bool notIsotherm(false);
  for (int k = 0; k < numbPhases; k++) {
    if (std::fabs(m_Tk0[0] - m_Tk0[k]) > 1.e-3) { notIsotherm = true; break; }
  }
  if (notIsotherm) throw ErrorXMLAttribut("temperature", fileName, __FILE__, __LINE__);
}

//****************************************************************************

BoundCondSubInj::BoundCondSubInj(const BoundCondSubInj& Source, const int& lvl) : BoundCond(Source, lvl)
{
  m_Tk0 = new double[numberPhases];
  m_ak0 = new double[numberPhases];
  m_m0 = Source.m_m0;

  for (int k = 0; k < numberPhases; k++)
  {
    m_Tk0[k] = Source.m_Tk0[k];
    m_ak0[k] = Source.m_ak0[k];
  }
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

void BoundCondSubInj::solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax)
{
  model->solveRiemannSubInj(cellLeft, dxLeft, dtMax, m_m0, m_Tk0, m_ak0, m_boundData);
}

//****************************************************************************

void BoundCondSubInj::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_m0 << std::endl;
  std::cout << m_Tk0[0] << std::endl;
}

//****************************************************************************
//******************************AMR Method***********************************
//****************************************************************************

void BoundCondSubInj::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondSubInj(*this, m_lvl + 1));
}

//****************************************************************************