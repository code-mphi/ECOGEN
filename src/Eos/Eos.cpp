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

//! \file      Eos.cpp
//! \author    F. Petitpas, K. Schmidmayer, E. Daniel
//! \version   1.0
//! \date      May 14 2018

#include <iostream>
#include <string>
#include "Eos.h"

using namespace std;
using namespace tinyxml2;

double epsilon;

//***********************************************************************

Eos::Eos(){}

//***********************************************************************

Eos::Eos(int &number) :
m_number(number), m_mu(-1.), m_lambda(-1.)
{
  number++;
}

//***********************************************************************

Eos::~Eos(){}

//***********************************************************************

void Eos::readPhysicalParameter(XMLNode *element, string fileName)
{
  XMLError error;

  XMLElement *sousElement(element->FirstChildElement("physicalParameters"));
  if (sousElement != NULL) {
    //Recuperation des donnees
    error = sousElement->QueryDoubleAttribute("mu", &m_mu);
    if (error != XML_NO_ERROR) m_mu = -1.;
    error = sousElement->QueryDoubleAttribute("lambda", &m_lambda);
    if (error != XML_NO_ERROR) m_lambda = -1.;
  }
}

//***********************************************************************

void Eos::display() const
{
    cout << "Fluid : " << m_name << endl;
}

//***********************************************************************

string Eos::getName() const
{
  return m_name;
}

//***********************************************************************

int Eos::getNumber() const
{
  return m_number;
}

//***********************************************************************

double Eos::computeTotalEnthalpy(const double &density, const double &pressure, const double &velocity) const
{
  return this->computeEnergy(density, pressure) + pressure / density + 0.5*velocity*velocity;
}

//***********************************************************************

double Eos::getMu() const { return m_mu; }

//***********************************************************************

double Eos::getLambda() const { return m_lambda; }

//***********************************************************************

void Eos::assignEpsilonForAlphaNull(bool alphaNull, string fileName) const
{
  if (alphaNull) {
    epsilon = 1.e-15;
  }
  else {
    epsilon = 0.;
  }
}

//***********************************************************************