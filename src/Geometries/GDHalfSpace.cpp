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

//! \file      GDHalfSpace.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 19 2017

#include "GDHalfSpace.h"

using namespace std;
using namespace tinyxml2;

//***************************************************************

GDHalfSpace::GDHalfSpace(string name, vector<Phase*> vecPhases, Mixture *mixture, vector<Transport> vecTransports, XMLElement *element, const int &physicalEntity, string fileName) :
  GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{
  XMLElement *sousElement(element->FirstChildElement("dataHalfSpace"));
  if (sousElement == NULL) throw ErrorXMLElement("dataHalfSpace", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //Origin
  error = sousElement->QueryDoubleAttribute("origin", &m_position);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("origin", fileName, __FILE__, __LINE__);
  //Axe
  string axe(sousElement->Attribute("axe"));
  Tools::uppercase(axe);
  if      (axe == "X"){ m_axe = X; }
  else if (axe == "Y"){ m_axe = Y; }
  else if (axe == "Z"){ m_axe = Z; }
  else { throw ErrorXMLAttribut("axe", fileName, __FILE__, __LINE__); }
  //Direction
  string direction(sousElement->Attribute("direction"));
  Tools::uppercase(direction);
  if      (direction == "POSITIVE"){ m_direction = 1; }
  else if (direction == "NEGATIVE"){ m_direction = -1; }
  else { throw ErrorXMLAttribut("direction", fileName, __FILE__, __LINE__); }
}

//***************************************************************

GDHalfSpace::~GDHalfSpace(){}

//***************************************************************

bool GDHalfSpace::belong(Coord &posElement, const int &lvl) const
{
  bool result(false);
  switch (m_axe)
  {
  case X:
    if (m_direction >= 0){if (posElement.getX() >= m_position) result = true;}
    else{ if (posElement.getX() <= m_position) result = true; }
    break;
  case Y:
    if (m_direction >= 0){ if (posElement.getY() >= m_position) result = true; }
    else{ if (posElement.getY() <= m_position) result = true; }
    break;
  case Z:
    if (m_direction >= 0){ if (posElement.getZ() >= m_position) result = true; }
    else{ if (posElement.getZ() <= m_position) result = true; }
    break;
  default:
    result = false; //false if not belongs to
    break;
  }
  return result;
}

//***************************************************************