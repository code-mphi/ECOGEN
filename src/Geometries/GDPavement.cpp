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

//! \file      GDPavement.cpp
//! \author    F. Petitpas
//! \version   1.0
//! \date      December 19 2017

#include <vector>
#include "GDPavement.h"
#include <iostream>

using namespace std;
using namespace tinyxml2;

//***************************************************************

GDPavement::GDPavement(string name, vector<Phase*> vecPhases, Mixture *mixture, vector<Transport> vecTransports, XMLElement *element, const int &physicalEntity, string fileName) :
GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity)
{
  XMLElement *sousElement(element->FirstChildElement("dataPavement"));
  if (sousElement == NULL) throw ErrorXMLElement("dataPavement", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  //width along X axis
  error = sousElement->QueryDoubleAttribute("lAxeX", &m_lX);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lAxeX", fileName, __FILE__, __LINE__);
  //width along Y axis
  error = sousElement->QueryDoubleAttribute("lAxeY", &m_lY);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lAxeX", fileName, __FILE__, __LINE__);
  //width along Z axis
  error = sousElement->QueryDoubleAttribute("lAxeZ", &m_lZ);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lAxeZ", fileName, __FILE__, __LINE__);
  //inferior vertex of pavement
  double x(0.), y(0.), z(0.);
  XMLElement *coin(sousElement->FirstChildElement("posInferiorVertex"));
  if (coin == NULL) throw ErrorXMLElement("posInferiorVertex", fileName, __FILE__, __LINE__);
  error = coin->QueryDoubleAttribute("x", &x);
  error = coin->QueryDoubleAttribute("y", &y);
  error = coin->QueryDoubleAttribute("z", &z);
  m_posXmYmZm.setXYZ(x, y, z);
}

//***************************************************************

GDPavement::~GDPavement(){}

//***************************************************************

bool GDPavement::belong(Coord &posElement, const int &lvl) const
{
  if ((posElement.getX() - m_posXmYmZm.getX())<0) return false;
  if ((posElement.getX() - m_posXmYmZm.getX())>m_lX) return false;
  if ((posElement.getY() - m_posXmYmZm.getY())<0) return false;
  if ((posElement.getY() - m_posXmYmZm.getY())>m_lY) return false;
  if ((posElement.getZ() - m_posXmYmZm.getZ())<0) return false;
  if ((posElement.getZ() - m_posXmYmZm.getZ())>m_lZ) return false;
  return true;
}