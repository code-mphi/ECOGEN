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

#include "BoundCondOutflow.h"

using namespace tinyxml2;

//****************************************************************************

BoundCondOutflow::BoundCondOutflow(int numPhysique, XMLElement* element, const int& numbTransports, std::vector<std::string> nameTransports, std::string fileName) :
  BoundCond(numPhysique)
{
  //Reading imposed outflow pressure
  XMLElement* sousElement(element->FirstChildElement("dataOutflow"));
  if (sousElement == NULL) throw ErrorXMLElement("dataOutflow", fileName, __FILE__, __LINE__);
  //Reading attributes
  //------------------
  XMLError error;
  error = sousElement->QueryDoubleAttribute("p0", &m_p0);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("p0", fileName, __FILE__, __LINE__);

  //Reading transports
  int couleurTrouvee(0);
  m_valueTransport = new double[numbTransports];
  XMLElement* elementTransport(sousElement->FirstChildElement("transport"));
  std::string nameTransport;
  while (elementTransport != NULL)
  {
    nameTransport = elementTransport->Attribute("name");
    if (nameTransport == "") throw ErrorXMLAttribut("name", fileName, __FILE__, __LINE__);
    int e(0);
    for (e = 0; e < numbTransports; e++) {
      if (nameTransport == nameTransports[e]) { break; }
    }
    if (e != numbTransports) {
      error = elementTransport->QueryDoubleAttribute("value", &m_valueTransport[e]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("value", fileName, __FILE__, __LINE__);
      couleurTrouvee++;
    }
    //Next transport
    elementTransport = elementTransport->NextSiblingElement("transport");
  }
  if (numbTransports > couleurTrouvee) throw ErrorXMLAttribut("Not enough transport equations in BC inj", fileName, __FILE__, __LINE__);
}

//****************************************************************************

BoundCondOutflow::BoundCondOutflow(const BoundCondOutflow& Source, const int& lvl) : BoundCond(Source, lvl)
{
  m_p0 = Source.m_p0;
  m_valueTransport = new double[numberTransports];
  for (int k = 0; k < numberTransports; k++) {
    m_valueTransport[k] = Source.m_valueTransport[k];
  }
}

//****************************************************************************

BoundCondOutflow::~BoundCondOutflow()
{
  delete[] m_valueTransport;
}

//****************************************************************************

void BoundCondOutflow::createBoundary(TypeMeshContainer<CellInterface*>& cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondOutflow(*(this)));
}

//****************************************************************************

void BoundCondOutflow::solveRiemannBoundary(Cell& cellLeft, const double& dxLeft, double& dtMax)
{
  model->solveRiemannOutflow(cellLeft, dxLeft, dtMax, m_p0, m_boundData);
}

//****************************************************************************

void BoundCondOutflow::solveRiemannTransportBoundary(Cell& cellLeft) const
{
	model->solveRiemannTransportOutflow(cellLeft, m_valueTransport);
}

//****************************************************************************

void BoundCondOutflow::printInfo()
{
  std::cout << m_numPhysique << std::endl;
  std::cout << m_p0 << std::endl;
}

//***************************************************************************
//******************************AMR Method***********************************
//***************************************************************************

void BoundCondOutflow::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondOutflow(*this, m_lvl + 1));
}

//****************************************************************************