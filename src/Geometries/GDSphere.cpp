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

#include "GDSphere.h"

using namespace tinyxml2;

//***************************************************************

GDSphere::GDSphere(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, XMLElement* element, const int& physicalEntity, std::string fileName) :
GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{
  XMLElement* sousElement(element->FirstChildElement("dataSphere"));
  if (sousElement == NULL) throw ErrorXMLElement("dataSphere", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //sphere radius
  error = sousElement->QueryDoubleAttribute("radius", &m_radius);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("radius", fileName, __FILE__, __LINE__);
  //center position
  double x(0.), y(0.), z(0.);
  XMLElement* coin(sousElement->FirstChildElement("center"));
  if (coin == NULL) throw ErrorXMLElement("center", fileName, __FILE__, __LINE__);
  error = coin->QueryDoubleAttribute("x", &x);
  error = coin->QueryDoubleAttribute("y", &y);
  error = coin->QueryDoubleAttribute("z", &z);
  m_centerPos.setXYZ(x, y, z);
}

//***************************************************************

GDSphere::~GDSphere(){}

//***************************************************************

bool GDSphere::belong(Coord& posElement, const int& /*lvl*/) const
{
  double sum;
  sum = std::pow(posElement.getX() - m_centerPos.getX(), 2.)
        + std::pow(posElement.getY() - m_centerPos.getY(), 2.)
        + std::pow(posElement.getZ() - m_centerPos.getZ(), 2.);
  if (sum <= m_radius*m_radius) return true;
  return false;
}