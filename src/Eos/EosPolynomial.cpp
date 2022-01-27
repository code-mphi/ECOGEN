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

#include "EosPolynomial.h"

//***********************************************************************

EosPolynomial::EosPolynomial(std::vector<std::string>& nameParameterEos, int& number) :
    Eos(number)
{
  nameParameterEos.push_back("a");
  nameParameterEos.push_back("b");
  nameParameterEos.push_back("c");
}

//***********************************************************************

EosPolynomial::~EosPolynomial(){}

//***********************************************************************

void EosPolynomial::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
  m_name = name;
  assert(parametersEos.size() == 3);
  m_a = parametersEos[0];
  m_b = parametersEos[1];
  m_c = parametersEos[2];
}

//***********************************************************************

double EosPolynomial::computePressure(const double& density, const double& /*temperature*/) const
{
  double epsilon = 1./density - 1.;
  return - (4.*m_a*(epsilon*epsilon*epsilon) + 3.*m_b*(epsilon*epsilon) + 2.*m_c*epsilon);
}

//***********************************************************************

double EosPolynomial::dedrho(const double& density, const double& temperature) const
{
  return this->computePressure(density, temperature) / (density*density);
}

//***********************************************************************

double EosPolynomial::dedrhoSecond(const double& density, const double& temperature) const
{
  double epsilon = 1./density - 1.;
  double dpdrho = (12.*m_a*(epsilon*epsilon) + 6.*m_b*epsilon + 2.*m_c) / (density*density);
  return dpdrho / (density*density) - 2. * this->dedrho(density, temperature) / density;
}

//***********************************************************************