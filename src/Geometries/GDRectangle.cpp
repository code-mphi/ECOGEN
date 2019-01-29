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

//! \file      GDRectangle.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 19 2017

#include <vector>
#include "GDRectangle.h"
#include <iostream>

using namespace std;
using namespace tinyxml2;

//***************************************************************

GDRectangle::GDRectangle(string name, vector<Phase*> vecPhases, Mixture *mixture, vector<Transport> vecTransports, XMLElement *element, const int &physicalEntity, string fileName) :
GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{
  XMLElement *sousElement(element->FirstChildElement("dataRectangle"));
  if (sousElement == NULL) throw ErrorXMLElement("dataRectangle", fileName, __FILE__, __LINE__);
  //Attributes reading
  //--------------------------
  XMLError error;
  //width along axe 1
  error = sousElement->QueryDoubleAttribute("lAxe1", &m_lAxe1);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lAxe1", fileName, __FILE__, __LINE__);
  //width along axe2
  error = sousElement->QueryDoubleAttribute("lAxe2", &m_lAxe2);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lAxe2", fileName, __FILE__, __LINE__);
  //Axe1
  string axe(sousElement->Attribute("axe1"));
  Tools::uppercase(axe);
  if (axe == "X"){ m_axe1 = X; }
  else if (axe == "Y"){ m_axe1 = Y; }
  else if (axe == "Z"){ m_axe1 = Z; }
  else { throw ErrorXMLAttribut("axe1", fileName, __FILE__, __LINE__); }
  //Axe2
  axe = sousElement->Attribute("axe2");
  Tools::uppercase(axe);
  if (axe == "X"){ m_axe2 = X; }
  else if (axe == "Y"){ m_axe2 = Y; }
  else if (axe == "Z"){ m_axe2 = Z; }
  else { throw ErrorXMLAttribut("axe2", fileName, __FILE__, __LINE__); }
  //Inferior vertex position
  double x(0.), y(0.), z(0.);
  XMLElement *coin(sousElement->FirstChildElement("posInferiorVertex"));
  if (coin == NULL) throw ErrorXMLElement("posInferiorVertex", fileName, __FILE__, __LINE__);
  error = coin->QueryDoubleAttribute("x", &x);
  error = coin->QueryDoubleAttribute("y", &y);
  error = coin->QueryDoubleAttribute("z", &z);
  m_posLeftBottom.setXYZ(x, y, z);
}

//***************************************************************

GDRectangle::~GDRectangle(){}

//***************************************************************

bool GDRectangle::belong(Coord &posElement, const int &lvl) const
{
  double somme(0.);
  vector<Axe> axes;
  axes.push_back(m_axe1);
  axes.push_back(m_axe2);
  vector<double> lengths;
  lengths.push_back(m_lAxe1);
  lengths.push_back(m_lAxe2);

  for (unsigned int i = 0; i < axes.size(); i++)
  {
    switch (axes[i])
    {
    case X:
      if (posElement.getX() - m_posLeftBottom.getX()<0) return false;
      if (posElement.getX() - m_posLeftBottom.getX()>lengths[i]) return false;
      break;
    case Y:
      if (posElement.getY() - m_posLeftBottom.getY()<0) return false;
      if (posElement.getY() - m_posLeftBottom.getY()>lengths[i]) return false;
      break;
    case Z:
      if (posElement.getZ() - m_posLeftBottom.getZ()<0) return false;
      if (posElement.getZ() - m_posLeftBottom.getZ()>lengths[i]) return false;
      break;
    }
  }
  return true;
}