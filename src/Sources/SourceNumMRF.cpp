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

#include "SourceNumMRF.h"

using namespace tinyxml2;

//***********************************************************************

SourceNumMRF::SourceNumMRF(XMLElement* element, int order, int physicalEntity, std::string fileName) : SourceNum(order, physicalEntity)
{
  XMLElement* sousElement(element->FirstChildElement("omega"));
  if (sousElement == NULL) throw ErrorXMLElement("omega", fileName, __FILE__, __LINE__);
  //Collecting attributes
  //---------------------
  XMLError error;
  double x, y, z;
  //Angular velocity components (rad/s)
  error = sousElement->QueryDoubleAttribute("x", &x);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
  error = sousElement->QueryDoubleAttribute("y", &y);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
  error = sousElement->QueryDoubleAttribute("z", &z);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);
  m_omega.setXYZ(x, y, z);
  m_incr = 1.;
  m_tf = 0.;

  //Optional time to increase linearly omega
  sousElement = element->FirstChildElement("timeToOmega");
  if (sousElement != NULL) {
    error = sousElement->QueryDoubleAttribute("tf", &m_tf);
    m_incr = 0.;
  }
}

//***********************************************************************

SourceNumMRF::~SourceNumMRF()
{}

//***********************************************************************

void SourceNumMRF::prepSourceTerms(Cell* cell, const int& i)
{
  sourceCons[i]->prepSourceTermsMRF(cell, m_incr*m_omega); 
}

//***********************************************************************

void SourceNumMRF::sourceEvolution(const double& time)
{
  if (m_incr < 1.) {
    m_incr = std::min(time/m_tf,1.);
  }
}

//***********************************************************************

Coord SourceNumMRF::computeAbsVelocity(const Coord& relVelocity, const Coord& position)
{
  Coord absVelocity(m_incr*m_omega);
  absVelocity = relVelocity + absVelocity.cross(position);
  return absVelocity;
}

//***********************************************************************