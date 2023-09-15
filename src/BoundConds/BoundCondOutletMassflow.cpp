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

#include "BoundCondOutletMassflow.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondOutletMassflow::BoundCondOutletMassflow(int numPhysique, XMLElement* element, std::string fileName) :
  BoundCond(numPhysique)
{
  XMLElement* sousElement(element->FirstChildElement("dataOutletMassflow"));
  if (sousElement == NULL) throw ErrorXMLElement("dataOutletMassflow", fileName, __FILE__, __LINE__);
  
  // Reading specific massflow (kg.s-1.m-2)
  // --------------------------------------
  XMLError error;
  error = sousElement->QueryDoubleAttribute("m0", &m_m0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("m0", fileName, __FILE__, __LINE__);
}

//****************************************************************************

BoundCondOutletMassflow::BoundCondOutletMassflow(const BoundCondOutletMassflow& Source, const int& lvl) : BoundCond(Source, lvl)
{
  m_m0 = Source.m_m0;
}

//****************************************************************************

BoundCondOutletMassflow::~BoundCondOutletMassflow()
{
}

//****************************************************************************

void BoundCondOutletMassflow::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondOutletMassflow(*(this)));
}

//****************************************************************************

void BoundCondOutletMassflow::solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax)
{
  model->solveRiemannOutletMassflow(cellLeft, dxLeft, dtMax, m_m0, m_boundData);
}

//****************************************************************************

void BoundCondOutletMassflow::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_m0 << std::endl;
}

//***************************************************************************
//******************************AMR Method***********************************
//***************************************************************************

void BoundCondOutletMassflow::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondOutletMassflow(*this, m_lvl + 1));
}

//****************************************************************************