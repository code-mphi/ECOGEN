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

//! \file      SourceMRF.cpp
//! \author    F. Petitpas, J. Caze
//! \version   1.0
//! \date      October 29 2019

#include "SourceMRF.h"

using namespace tinyxml2;

//***********************************************************************

SourceMRF::SourceMRF(){}

//***********************************************************************

SourceMRF::SourceMRF(XMLElement *element, int order, std::string fileName) : Source(order)
{
  XMLElement *sousElement(element->FirstChildElement("omega"));
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

SourceMRF::~SourceMRF()
{}

//***********************************************************************

void SourceMRF::prepSourceTerms(Cell *cell, const int &numberPhases, const double &dt, const int i)
{
  sourceCons[i]->prepSourceTermsMRF(cell, dt, numberPhases, m_incr*m_omega); 
}

//***********************************************************************

void SourceMRF::sourceEvolution(const double &time)
{
  if (m_incr < 1.) {
    m_incr = std::min(time/m_tf,1.);
  }
}

//***********************************************************************

Coord SourceMRF::computeAbsVelocity(const Coord relVelocity, const Coord position)
{
  Coord absVelocity(m_incr*m_omega);
  absVelocity = relVelocity + absVelocity.cross(position);
  return absVelocity;
}

//***********************************************************************