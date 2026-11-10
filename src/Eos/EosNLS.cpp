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

#include "EosNLS.h"

//***********************************************************************

EosNLS::EosNLS(std::vector<std::string>& nameParameterEos, int& number) : Eos(number) { nameParameterEos.push_back("no_params"); }

//***********************************************************************

EosNLS::~EosNLS() {}

//***********************************************************************

void EosNLS::assignParametersEos(std::string name, std::vector<double> parametersEos)
{
  assert(parametersEos.size() == 1);

  m_name   = name;
  m_params = parametersEos[0];
}

//***********************************************************************

double EosNLS::computeEnergy(const double& density, const double& /*pressure*/) const { return 0.5 * density; }

//***********************************************************************

double EosNLS::computePressure(const double& density, const double& /*energy*/) const { return 0.5 * density * density; }

//***********************************************************************

double EosNLS::computePressureFromDensityAndEntropy(const double& density, const double& /*entropy*/) const { return 0.5 * density * density; }

//***********************************************************************

double EosNLS::computeTemperature(const double& /*density*/, const double& /*pressure*/) const
{
  return 1.; //This value does not matter
}

//***********************************************************************

double EosNLS::computeEntropyFromDensityAndEnergy(const double& /*density*/, const double& /*energy*/) const
{
  return 1.; //This value does not matter
}

//***********************************************************************

double EosNLS::computeEntropyFromDensityAndPressure(const double& /*density*/, const double& /*pressure*/) const
{
  return 1.; //This value does not matter
}

//***********************************************************************

double EosNLS::computeSquareSoundSpeed(const double& density, const double& /*energy*/) const { return density; }

//***********************************************************************

double EosNLS::kappa(const double& density) const { return 0.25 / density; }

//****************************************************************************

double EosNLS::kappaPrime(const double& density) const { return -0.25 / (density * density); }

//****************************************************************************

double EosNLS::kappaSecond(const double& density) const { return 0.5 / (density * density * density); }

//***********************************************************************
