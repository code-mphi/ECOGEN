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

#include <cmath>
#include "MixPUEq.h"

using namespace tinyxml2;

//***************************************************************************

MixPUEq::MixPUEq() : MixUEq(){}

//***************************************************************************

MixPUEq::MixPUEq(XMLElement* state, std::string fileName) : MixUEq()
{
  XMLElement* sousElement(state->FirstChildElement("mixture"));
  if (sousElement == NULL) throw ErrorXMLElement("mixture", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;
  XMLElement* dataMix(sousElement->FirstChildElement("dataMix"));
  if (dataMix == NULL) throw ErrorXMLElement("dataMix", fileName, __FILE__, __LINE__);
  //pressure
  error = dataMix->QueryDoubleAttribute("pressure", &m_pressure);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("pressure", fileName, __FILE__, __LINE__);
  //velocity
  XMLElement* velocity(sousElement->FirstChildElement("velocity"));
  if (velocity == NULL) throw ErrorXMLElement("velocity", fileName, __FILE__, __LINE__);
  double velocityX(0.), velocityY(0.), velocityZ(0.);
  error = velocity->QueryDoubleAttribute("x", &velocityX);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
  error = velocity->QueryDoubleAttribute("y", &velocityY);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
  error = velocity->QueryDoubleAttribute("z", &velocityZ);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);
  m_velocity.setXYZ(velocityX, velocityY, velocityZ);
}

//***************************************************************************

MixPUEq::~MixPUEq(){}

//***************************************************************************

void MixPUEq::allocateAndCopyMixture(Mixture** mixture)
{
  *mixture = new MixPUEq(*this);
}

//***************************************************************************