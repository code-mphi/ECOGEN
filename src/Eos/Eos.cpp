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

#include "Eos.h"

using namespace tinyxml2;

double epsilonAlphaNull;

//***********************************************************************

Eos::Eos(int& number) :
m_number(number), m_mu(-1.), m_lambda(-1.)
{
  number++;
}

//***********************************************************************

Eos::~Eos(){}

//***********************************************************************

void Eos::readPhysicalParameter(XMLNode *element)
{
  XMLError error;

  XMLElement* sousElement(element->FirstChildElement("physicalParameters"));
  if (sousElement != NULL) {
    //Recuperation des donnees
    error = sousElement->QueryDoubleAttribute("mu", &m_mu);
    if (error != XML_NO_ERROR) m_mu = -1.;
    error = sousElement->QueryDoubleAttribute("lambda", &m_lambda);
    if (error != XML_NO_ERROR) m_lambda = -1.;
  }
}

//***********************************************************************

double Eos::computeEnthalpy(const double& density, const double& pressure) const
{
  return this->computeEnergy(density, pressure) + pressure / std::max(density, epsilonAlphaNull);
}

//***********************************************************************

double Eos::computeTotalEnthalpy(const double& density, const double& pressure, const double& velX, const double& velY, const double &velZ) const
{
  return this->computeEnthalpy(density, pressure) + 0.5*(velX*velX + velY*velY + velZ*velZ);
}

//***********************************************************************

double Eos::computeTotalEnthalpy(const double& density, const double& pressure, const Coord& velocity) const
{
  return this->computeEnthalpy(density, pressure) + 0.5*velocity.squaredNorm();
}

//***********************************************************************

void Eos::assignEpsilonForAlphaNull(bool alphaNull) const
{
  if (alphaNull) {
    epsilonAlphaNull = 1.e-15;
  }
  else {
    epsilonAlphaNull = 0.;
  }
}

//***********************************************************************