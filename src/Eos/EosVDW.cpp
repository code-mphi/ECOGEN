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

#include "EosVDW.h"

//***********************************************************************

EosVDW::EosVDW(std::vector<std::string>& nameParameterEos, int& number) :
    Eos(number)
{
  nameParameterEos.push_back("Tc");
  nameParameterEos.push_back("Pc");
  nameParameterEos.push_back("r");
  nameParameterEos.push_back("gamma");
  nameParameterEos.push_back("p0");
  nameParameterEos.push_back("rho0");
}

//***********************************************************************

EosVDW::~EosVDW(){}

//***********************************************************************

void EosVDW::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
  m_name   = name;
  assert(parametersEos.size() == 6);
  double Tc = parametersEos[0]; // Critical temperature of the fluid 
  double Pc = parametersEos[1]; // Critical pressure of the fluid
  m_r = parametersEos[2];       // Ideal-gas constant
  m_a = 27. * m_r * m_r * Tc * Tc / (64. * Pc);
  m_b = m_r * Tc / (8. * Pc);
  m_gamma = parametersEos[3];  //Polytropic constant
  m_p0 = parametersEos[4];     //Reference pressure
  m_rho0 = parametersEos[5];   //Reference density
}

//***********************************************************************

double EosVDW::computeEnergy(const double& density, const double& temperature) const
{
  return -m_r * temperature * log(1./density - m_b) - density * m_a;
}

//***********************************************************************

double EosVDW::computePressure(const double& density, const double& /*temperature*/) const
{
  return (m_p0 + m_a * m_rho0*m_rho0) * std::pow((1. / m_rho0 - m_b) / (1. / density - m_b), m_gamma)
          - m_a * density*density;
}

//***********************************************************************

double EosVDW::computeTemperature(const double& density, const double& energy) const
{
  return (energy + density * m_a) / (-m_r * log(1./density - m_b));
}

//***********************************************************************

double EosVDW::dedrho(const double& density, const double& temperature) const
{
  return this->computePressure(density, temperature) / (density*density);
}

//***********************************************************************

double EosVDW::dedrhoSecond(const double& density, const double& temperature) const
{
  double dpdrho = (this->computePressure(density, temperature) + m_a * density*density) / (1. / density - m_b)
                    * m_gamma / (density*density) - 2. * m_a * density;
  return dpdrho / (density*density) - 2. * this->dedrho(density, temperature) / density;
}

//***********************************************************************

void EosVDW::verifyAndCorrectDensityMax(double& density) const
{
  if (density > 1. / m_b - 1.e-10) {
    density = 1. / m_b - 1.e-10;
  }
}

//***********************************************************************